// enamine_patch_features.cpp
// Resumable/multithreaded Enamine 64.9B -> 90B feature patcher.
// Identity is literal column-1 SMILES/CXSMILES, NOT ligand name.
// Existing feature batches are never modified; novel molecules are appended as
// new batch_XXXXXX directories using the existing Python/RDKit worker.
//
// Build:
//   g++ -O3 -std=c++17 -pthread enamine_patch_features.cpp -lbz2 -o enamine_patch_features
//
// Example:
// ./enamine_patch_features \
//   --old-root /old64/M/H32 \
//   --new-root /new90/M/H32 \
//   --feature-root /features/M/H32 \
//   --work-dir /scratch/$USER/patch_M_H32 \
//   --pipeline-script /path/enamine_recipe_feature_pipeline_with_query.py \
//   --threads 32 --feature-workers 32 --buckets 512 --batch-size 3000000 --resume

#include <bzlib.h>
#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <optional>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>
#include <unistd.h>

namespace fs = std::filesystem;

struct Args {
    fs::path old_root, new_root, feature_root, work_dir, pipeline_script;
    std::string python = "python";
    size_t threads = 1, feature_workers = 1, buckets = 512, buffer_bytes = 4096;
    uint64_t batch_size = 3000000;
    int linker_max_heavy_atoms = 4, minhash_perm = 64, lsh_rows_per_band = 8, feature_hash_seed = 0;
    bool resume=false, keep_new_duplicates=false, cleanup_on_complete=false;
    bool map_include_ligand_name=false, map_include_smiles=false;
    bool write_recipes=false, write_readable_features=false, write_readable_lsh=false;
    bool write_summary_csvs=false, write_feature_token_map=false, write_lsh_memberships_bin=false;
};

[[noreturn]] static void usage(const char* p,int rc=1){
    std::cerr
      <<"Usage: "<<p<<" --old-root DIR --new-root DIR --feature-root DIR --work-dir DIR\n"
      <<"  --pipeline-script FILE [--python python] [--threads N] [--feature-workers N]\n"
      <<"  [--buckets 512] [--batch-size 3000000] [--resume] [--cleanup-on-complete]\n"
      <<"  [--keep-new-duplicates]\n"
      <<"Feature options (must match old run):\n"
      <<"  --linker-max-heavy-atoms N --minhash-perm N --lsh-rows-per-band N\n"
      <<"  --feature-hash-seed N --map-include-ligand-name --map-include-smiles\n"
      <<"  --write-recipes --write-readable-features --write-readable-lsh\n"
      <<"  --write-summary-csvs --write-feature-token-map --write-lsh-memberships-bin\n";
    std::exit(rc);
}
static Args parse_args(int argc,char** argv){
    Args a;
    unsigned hc=std::thread::hardware_concurrency();
    a.threads=hc?hc:1; a.feature_workers=a.threads;
    for(int i=1;i<argc;++i){
        std::string x=argv[i];
        auto need=[&](const char* f)->std::string{ if(i+1>=argc) throw std::runtime_error(std::string("Missing ")+f); return argv[++i]; };
        if(x=="--old-root") a.old_root=need("--old-root");
        else if(x=="--new-root") a.new_root=need("--new-root");
        else if(x=="--feature-root") a.feature_root=need("--feature-root");
        else if(x=="--work-dir") a.work_dir=need("--work-dir");
        else if(x=="--pipeline-script") a.pipeline_script=need("--pipeline-script");
        else if(x=="--python") a.python=need("--python");
        else if(x=="--threads") a.threads=std::stoull(need("--threads"));
        else if(x=="--feature-workers") a.feature_workers=std::stoull(need("--feature-workers"));
        else if(x=="--buckets") a.buckets=std::stoull(need("--buckets"));
        else if(x=="--buffer-kb-per-bucket") a.buffer_bytes=std::stoull(need("--buffer-kb-per-bucket"))*1024ULL;
        else if(x=="--batch-size") a.batch_size=std::stoull(need("--batch-size"));
        else if(x=="--linker-max-heavy-atoms") a.linker_max_heavy_atoms=std::stoi(need("--linker-max-heavy-atoms"));
        else if(x=="--minhash-perm") a.minhash_perm=std::stoi(need("--minhash-perm"));
        else if(x=="--lsh-rows-per-band") a.lsh_rows_per_band=std::stoi(need("--lsh-rows-per-band"));
        else if(x=="--feature-hash-seed") a.feature_hash_seed=std::stoi(need("--feature-hash-seed"));
        else if(x=="--resume") a.resume=true;
        else if(x=="--keep-new-duplicates") a.keep_new_duplicates=true;
        else if(x=="--cleanup-on-complete") a.cleanup_on_complete=true;
        else if(x=="--map-include-ligand-name") a.map_include_ligand_name=true;
        else if(x=="--map-include-smiles") a.map_include_smiles=true;
        else if(x=="--write-recipes") a.write_recipes=true;
        else if(x=="--write-readable-features") a.write_readable_features=true;
        else if(x=="--write-readable-lsh") a.write_readable_lsh=true;
        else if(x=="--write-summary-csvs") a.write_summary_csvs=true;
        else if(x=="--write-feature-token-map") a.write_feature_token_map=true;
        else if(x=="--write-lsh-memberships-bin") a.write_lsh_memberships_bin=true;
        else if(x=="-h"||x=="--help") usage(argv[0],0);
        else throw std::runtime_error("Unknown option: "+x);
    }
    if(a.old_root.empty()||a.new_root.empty()||a.feature_root.empty()||a.work_dir.empty()||a.pipeline_script.empty()) usage(argv[0]);
    if(!a.threads||!a.feature_workers||!a.buckets||!a.batch_size) throw std::runtime_error("threads/workers/buckets/batch-size must be >0");
    if(a.minhash_perm%a.lsh_rows_per_band) throw std::runtime_error("minhash-perm must be divisible by lsh-rows-per-band");
    return a;
}

static uint64_t fnv1a64(const std::string& s){
    uint64_t h=14695981039346656037ULL; for(unsigned char c:s){h^=c;h*=1099511628211ULL;} return h;
}
static std::string shell_quote(const std::string& s){
    std::string o="'"; for(char c:s){ if(c=='\'') o+="'\\''"; else o+=c; } return o+"'";
}
static void atomic_text(const fs::path& p,const std::string& s){
    fs::create_directories(p.parent_path()); fs::path t=p; t+=".tmp."+std::to_string(getpid());
    { std::ofstream f(t,std::ios::binary|std::ios::trunc); if(!f) throw std::runtime_error("Cannot create "+t.string()); f<<s; }
    fs::rename(t,p);
}
static std::vector<fs::path> bz2_files(const fs::path& root){
    std::vector<fs::path> v;
    for(auto const& e:fs::recursive_directory_iterator(root)){
        if(!e.is_regular_file()) continue; auto n=e.path().filename().string();
        if(n.size()>=4&&n.substr(n.size()-4)==".bz2") v.push_back(fs::absolute(e.path()));
    }
    std::sort(v.begin(),v.end());
    if(v.empty()) throw std::runtime_error("No .bz2 files under "+root.string());
    return v;
}
static bool parse_smiles(const std::string& line,std::string& smiles){
    if(line.empty()) return false;
    size_t t=line.find('\t');
    if(t!=std::string::npos) smiles=line.substr(0,t);
    else{
        size_t s=line.find_first_not_of(" \r\n"); if(s==std::string::npos) return false;
        size_t e=line.find_first_of(" \t\r\n",s); if(e==std::string::npos) return false;
        smiles=line.substr(s,e-s);
    }
    if(smiles=="smiles"||smiles=="SMILES"||smiles=="cxsmiles"||smiles=="CXSMILES") return false;
    return !smiles.empty();
}

class BzReader{
    FILE* fp_=nullptr; BZFILE* bz_=nullptr; bool eof_=false; std::string pending_; fs::path p_;
public:
    explicit BzReader(const fs::path& p):p_(p){
        fp_=std::fopen(p.string().c_str(),"rb"); if(!fp_) throw std::runtime_error("Cannot open "+p.string());
        int e=BZ_OK; bz_=BZ2_bzReadOpen(&e,fp_,0,0,nullptr,0); if(e!=BZ_OK) throw std::runtime_error("BZ2 open failed "+p.string());
    }
    ~BzReader(){ if(bz_){int e=BZ_OK;BZ2_bzReadClose(&e,bz_);} if(fp_)std::fclose(fp_); }
    bool getline(std::string& out){
        out.clear();
        for(;;){
            auto n=pending_.find('\n');
            if(n!=std::string::npos){out=pending_.substr(0,n);pending_.erase(0,n+1);if(!out.empty()&&out.back()=='\r')out.pop_back();return true;}
            if(eof_){if(pending_.empty())return false;out.swap(pending_);if(!out.empty()&&out.back()=='\r')out.pop_back();return true;}
            char b[1<<20]; int e=BZ_OK; int got=BZ2_bzRead(&e,bz_,b,sizeof(b)); if(got>0)pending_.append(b,(size_t)got);
            if(e==BZ_STREAM_END)eof_=true; else if(e!=BZ_OK)throw std::runtime_error("BZ2 read failed "+p_.string());
        }
    }
};
class BzWriter{
    FILE* fp_=nullptr; BZFILE* bz_=nullptr; fs::path p_;
public:
    explicit BzWriter(const fs::path& p):p_(p){
        fs::create_directories(p.parent_path()); fp_=std::fopen(p.string().c_str(),"wb"); if(!fp_)throw std::runtime_error("Cannot create "+p.string());
        int e=BZ_OK; bz_=BZ2_bzWriteOpen(&e,fp_,9,0,30); if(e!=BZ_OK)throw std::runtime_error("BZ2 write open failed");
    }
    ~BzWriter(){close();}
    void line(const std::string& s){std::string t=s+"\n";int e=BZ_OK;BZ2_bzWrite(&e,bz_,const_cast<char*>(t.data()),(int)t.size());if(e!=BZ_OK)throw std::runtime_error("BZ2 write failed");}
    void close(){if(bz_){int e=BZ_OK;BZ2_bzWriteClose(&e,bz_,0,nullptr,nullptr);bz_=nullptr;}if(fp_){std::fclose(fp_);fp_=nullptr;}}
};

static fs::path src_dir(const fs::path& r,size_t i){std::ostringstream o;o<<"source_"<<std::setw(6)<<std::setfill('0')<<i;return r/o.str();}
static fs::path src_done(const fs::path& r,size_t i){std::ostringstream o;o<<"source_"<<std::setw(6)<<std::setfill('0')<<i<<".done";return r/o.str();}
static fs::path bucket_part(const fs::path& r,size_t b){std::ostringstream o;o<<"bucket_"<<std::setw(6)<<std::setfill('0')<<b<<".dat";return r/o.str();}

class SourceBucketWriter{
    fs::path r_; std::vector<std::string> buf_; size_t limit_;
public:
    SourceBucketWriter(const fs::path& r,size_t n,size_t l):r_(r),buf_(n),limit_(l){fs::create_directories(r_);}
    ~SourceBucketWriter(){try{flush_all();}catch(...){ }}
    void add(size_t b,const std::string& s){auto& x=buf_.at(b);x+=s;x+='\n';if(x.size()>=limit_)flush(b);}
    void flush(size_t b){auto& x=buf_.at(b);if(x.empty())return;std::ofstream f(bucket_part(r_,b),std::ios::binary|std::ios::app);if(!f)throw std::runtime_error("Cannot append bucket");f.write(x.data(),(std::streamsize)x.size());x.clear();}
    void flush_all(){for(size_t b=0;b<buf_.size();++b)flush(b);}
};

template<class F> static void parallel_for(size_t n,size_t nt,F fn){
    std::atomic<size_t> next{0}; std::atomic<bool> fail{false}; std::exception_ptr ep; std::mutex em;
    std::vector<std::thread> pool;
    for(size_t t=0;t<std::min(nt,std::max<size_t>(1,n));++t) pool.emplace_back([&]{
        while(!fail){
            size_t i=next.fetch_add(1); if(i>=n)break;
            try{fn(i);}catch(...){fail=true;std::lock_guard<std::mutex>lk(em);if(!ep)ep=std::current_exception();}
        }
    });
    for(auto& t:pool)t.join(); if(ep)std::rethrow_exception(ep);
}

enum class Kind{Old,New};
static uint64_t partition_source(const Args& a,const fs::path& phase,const fs::path& in,size_t idx,Kind kind){
    auto done=src_done(phase,idx);
    if(fs::exists(done)){std::ifstream f(done);uint64_t n=0;f>>n;return n;}
    auto out=src_dir(phase,idx);std::error_code ec;fs::remove_all(out,ec);fs::create_directories(out);
    SourceBucketWriter w(out,a.buckets,a.buffer_bytes);BzReader r(in);std::string line,sm;uint64_t n=0;
    while(r.getline(line)){if(!parse_smiles(line,sm))continue;size_t b=fnv1a64(sm)%a.buckets;w.add(b,kind==Kind::Old?sm:(sm+"\t"+line));++n;}
    w.flush_all();atomic_text(done,std::to_string(n)+"\n");return n;
}
static void partition_phase(const Args& a,const std::vector<fs::path>& files,const fs::path& phase,Kind kind,const char* label){
    fs::create_directories(phase);std::atomic<size_t> fin{0};std::mutex lm;
    parallel_for(files.size(),a.threads,[&](size_t i){
        auto n=partition_source(a,phase,files[i],i,kind);auto k=++fin;std::lock_guard<std::mutex>lk(lm);
        std::cerr<<"["<<label<<"] "<<k<<"/"<<files.size()<<" rows="<<n<<" "<<files[i]<<"\n";
    }); atomic_text(phase/"PHASE.done","complete\n");
}

static fs::path novel_bz2(const fs::path& r,size_t b){std::ostringstream o;o<<"novel_bucket_"<<std::setw(6)<<std::setfill('0')<<b<<".cxsmiles.bz2";return r/o.str();}
static fs::path cmp_done(const fs::path& r,size_t b){std::ostringstream o;o<<"bucket_"<<std::setw(6)<<std::setfill('0')<<b<<".done";return r/o.str();}
struct Cnt{uint64_t new_rows=0,present=0,novel=0,dup=0;};
static Cnt read_cnt(const fs::path& p){Cnt c;std::ifstream f(p);f>>c.new_rows>>c.present>>c.novel>>c.dup;return c;}

static Cnt compare_bucket(const Args& a,const fs::path& oldp,const fs::path& newp,size_t oldn,size_t newn,const fs::path& cmp,size_t b){
    auto done=cmp_done(cmp,b); if(fs::exists(done))return read_cnt(done);
    fs::create_directories(cmp);std::unordered_set<std::string> old;old.max_load_factor(.8f);
    for(size_t i=0;i<oldn;++i){auto p=bucket_part(src_dir(oldp,i),b);if(!fs::exists(p))continue;std::ifstream f(p);std::string s;while(std::getline(f,s)){if(!s.empty()&&s.back()=='\r')s.pop_back();if(!s.empty())old.insert(std::move(s));}}
    std::optional<std::unordered_set<std::string>> seen;if(!a.keep_new_duplicates){seen.emplace();seen->max_load_factor(.8f);}
    auto outp=novel_bz2(cmp,b);auto tmp=outp;tmp+=".partial."+std::to_string(getpid())+"."+std::to_string(b);std::error_code ec;fs::remove(tmp,ec);
    Cnt c;{BzWriter out(tmp);
      for(size_t i=0;i<newn;++i){auto p=bucket_part(src_dir(newp,i),b);if(!fs::exists(p))continue;std::ifstream f(p);std::string rec;
        while(std::getline(f,rec)){if(!rec.empty()&&rec.back()=='\r')rec.pop_back();auto t=rec.find('\t');if(t==std::string::npos)continue;std::string key=rec.substr(0,t);++c.new_rows;
          if(old.find(key)!=old.end()){++c.present;continue;}
          if(seen&&!seen->insert(key).second){++c.dup;continue;}
          out.line(rec.substr(t+1));++c.novel;
        }}
      out.close();
    }
    fs::rename(tmp,outp);std::ostringstream m;m<<c.new_rows<<"\t"<<c.present<<"\t"<<c.novel<<"\t"<<c.dup<<"\n";atomic_text(done,m.str());return c;
}
static void compare_phase(const Args& a,const fs::path& oldp,const fs::path& newp,size_t oldn,size_t newn,const fs::path& cmp){
    fs::create_directories(cmp);std::atomic<size_t> fin{0};std::mutex lm;
    parallel_for(a.buckets,a.threads,[&](size_t b){auto c=compare_bucket(a,oldp,newp,oldn,newn,cmp,b);auto k=++fin;std::lock_guard<std::mutex>lk(lm);
      std::cerr<<"[compare] "<<k<<"/"<<a.buckets<<" bucket="<<b<<" present="<<c.present<<" novel="<<c.novel<<" new_dup="<<c.dup<<"\n";
    });atomic_text(cmp/"PHASE.done","complete\n");
}

static uint64_t max_batch(const fs::path& root){
    uint64_t m=0;std::regex re("^batch_([0-9]+)$");for(auto const& e:fs::directory_iterator(root)){if(!e.is_directory())continue;std::smatch x;auto n=e.path().filename().string();if(std::regex_match(n,x,re))m=std::max<uint64_t>(m,static_cast<uint64_t>(std::stoull(x[1].str())));}return m;
}
static fs::path stats_path(const fs::path& root,uint64_t id){
    std::ostringstream o;o<<"batch_"<<std::setw(6)<<std::setfill('0')<<id;auto d=root/o.str();return d/(o.str()+".stats.json");
}
static bool valid_stats(const fs::path& p){
    std::error_code ec;if(!fs::is_regular_file(p,ec)||fs::file_size(p,ec)<2)return false;std::ifstream f(p);char c=0,first=0,last=0;while(f.get(c))if(!std::isspace((unsigned char)c)){if(!first)first=c;last=c;}return first=='{'&&last=='}';
}

static fs::path build_manifest(const Args& a,const fs::path& cmp){
    auto man=a.work_dir/"patch_manifest.tsv",meta=a.work_dir/"patch_manifest.meta";if(fs::exists(man)&&fs::exists(meta))return man;
    uint64_t base=max_batch(a.feature_root),task=base+1,first=task,cur=0,total=0;
    auto tmp=man;tmp+=".tmp."+std::to_string(getpid());std::ofstream o(tmp);o<<"task_id\tinput_path\tstart_row\tend_row\tnum_ligands\n";
    for(size_t b=0;b<a.buckets;++b){auto c=read_cnt(cmp_done(cmp,b));uint64_t s=0;
      while(s<c.novel){if(cur>=a.batch_size){++task;cur=0;}uint64_t take=std::min<uint64_t>(a.batch_size-cur,c.novel-s);
        o<<task<<"\t"<<fs::absolute(novel_bz2(cmp,b)).string()<<"\t"<<s<<"\t"<<s+take<<"\t"<<take<<"\n";cur+=take;s+=take;total+=take;
      }}
    o.close();if(total==0){first=task=0;}fs::rename(tmp,man);std::ostringstream m;m<<"base_existing_batch="<<base<<"\nfirst_patch_task="<<first<<"\nlast_patch_task="<<task<<"\ntotal_novel="<<total<<"\n";atomic_text(meta,m.str());
    std::cerr<<"[manifest] total novel="<<total<<(total?(" patch IDs "+std::to_string(first)+"-"+std::to_string(task)):"")<<"\n";return man;
}
static std::vector<uint64_t> manifest_ids(const fs::path& p){
    std::ifstream f(p);std::string line;std::getline(f,line);std::vector<uint64_t> ids;uint64_t last=0;
    while(std::getline(f,line)){auto t=line.find('\t');if(t==std::string::npos)continue;uint64_t id=std::stoull(line.substr(0,t));if(id!=last)ids.push_back(id);last=id;}return ids;
}
static std::string worker_cmd(const Args& a,const fs::path& man,uint64_t id){
    std::ostringstream c;c<<shell_quote(a.python)<<" "<<shell_quote(fs::absolute(a.pipeline_script).string())<<" worker"
      <<" --manifest "<<shell_quote(fs::absolute(man).string())<<" --task-id "<<id<<" --output-root "<<shell_quote(fs::absolute(a.feature_root).string())
      <<" --linker-max-heavy-atoms "<<a.linker_max_heavy_atoms<<" --minhash-perm "<<a.minhash_perm<<" --lsh-rows-per-band "<<a.lsh_rows_per_band<<" --feature-hash-seed "<<a.feature_hash_seed;
    if(a.map_include_ligand_name)c<<" --map-include-ligand-name";if(a.map_include_smiles)c<<" --map-include-smiles";if(a.write_recipes)c<<" --write-recipes";
    if(a.write_readable_features)c<<" --write-readable-features";if(a.write_readable_lsh)c<<" --write-readable-lsh";if(a.write_summary_csvs)c<<" --write-summary-csvs";
    if(a.write_feature_token_map)c<<" --write-feature-token-map";if(a.write_lsh_memberships_bin)c<<" --write-lsh-memberships-bin";return c.str();
}
static void feature_phase(const Args& a,const fs::path& man){
    auto ids=manifest_ids(man);std::vector<uint64_t> rem;for(auto id:ids)if(!valid_stats(stats_path(a.feature_root,id)))rem.push_back(id);
    std::cerr<<"[features] total="<<ids.size()<<" remaining="<<rem.size()<<" workers="<<a.feature_workers<<"\n";std::atomic<size_t> fin{0};std::mutex lm;
    parallel_for(rem.size(),a.feature_workers,[&](size_t i){auto id=rem[i];{std::lock_guard<std::mutex>lk(lm);std::cerr<<"[features] start "<<id<<"\n";}
      int rc=std::system(worker_cmd(a,man,id).c_str());if(rc!=0)throw std::runtime_error("Feature worker failed task "+std::to_string(id));
      if(!valid_stats(stats_path(a.feature_root,id)))throw std::runtime_error("Missing/invalid stats for task "+std::to_string(id));
      auto k=++fin;std::lock_guard<std::mutex>lk(lm);std::cerr<<"[features] completed "<<k<<"/"<<rem.size()<<" task="<<id<<"\n";
    });
    for(auto id:ids)if(!valid_stats(stats_path(a.feature_root,id)))throw std::runtime_error("Incomplete feature task "+std::to_string(id));
    atomic_text(a.work_dir/"features.done","complete\n");
}

static std::string signature(const Args& a,const std::vector<fs::path>& oldf,const std::vector<fs::path>& newf){
    std::ostringstream o;o<<"format=v2\nold="<<fs::absolute(a.old_root)<<"\nnew="<<fs::absolute(a.new_root)<<"\nfeatures="<<fs::absolute(a.feature_root)
      <<"\npipeline="<<fs::absolute(a.pipeline_script)<<"\nbuckets="<<a.buckets<<"\nbatch="<<a.batch_size<<"\nidentity=raw_smiles_col1\nkeep_dup="<<a.keep_new_duplicates
      <<"\nlinker="<<a.linker_max_heavy_atoms<<"\nminhash="<<a.minhash_perm<<"\nrowsband="<<a.lsh_rows_per_band<<"\nseed="<<a.feature_hash_seed
      <<"\nmap_name="<<a.map_include_ligand_name<<"\nmap_smiles="<<a.map_include_smiles<<"\nrecipes="<<a.write_recipes<<"\nread_feat="<<a.write_readable_features
      <<"\nread_lsh="<<a.write_readable_lsh<<"\nsummaries="<<a.write_summary_csvs<<"\ntoken_map="<<a.write_feature_token_map<<"\nlsh_mem="<<a.write_lsh_memberships_bin<<"\n[old]\n";
    for(auto&p:oldf)o<<p<<"\n";o<<"[new]\n";for(auto&p:newf)o<<p<<"\n";return o.str();
}
static void init_state(const Args& a,const std::vector<fs::path>& oldf,const std::vector<fs::path>& newf){
    fs::create_directories(a.work_dir);auto p=a.work_dir/"patch_config.txt";auto s=signature(a,oldf,newf);
    if(fs::exists(p)){std::ifstream f(p);std::stringstream b;b<<f.rdbuf();if(b.str()!=s)throw std::runtime_error("Work-dir checkpoint settings/files differ; use a fresh work-dir.");if(!a.resume)throw std::runtime_error("Checkpoint exists; use --resume.");}
    else atomic_text(p,s);
}

int main(int argc,char** argv){
    try{
        Args a=parse_args(argc,argv);fs::create_directories(a.feature_root);auto oldf=bz2_files(a.old_root),newf=bz2_files(a.new_root);init_state(a,oldf,newf);
        auto oldp=a.work_dir/"old_partitions",newp=a.work_dir/"new_partitions",cmp=a.work_dir/"novel_buckets";
        std::cerr<<"Patch run: raw SMILES identity; old files="<<oldf.size()<<" new files="<<newf.size()<<" buckets="<<a.buckets<<" threads="<<a.threads<<" feature_workers="<<a.feature_workers<<"\n";
        if(!fs::exists(oldp/"PHASE.done"))partition_phase(a,oldf,oldp,Kind::Old,"old");else std::cerr<<"[old] checkpoint complete\n";
        if(!fs::exists(newp/"PHASE.done"))partition_phase(a,newf,newp,Kind::New,"new");else std::cerr<<"[new] checkpoint complete\n";
        if(!fs::exists(cmp/"PHASE.done"))compare_phase(a,oldp,newp,oldf.size(),newf.size(),cmp);else std::cerr<<"[compare] checkpoint complete\n";
        auto man=build_manifest(a,cmp);feature_phase(a,man);
        if(a.cleanup_on_complete){std::error_code ec;fs::remove_all(oldp,ec);fs::remove_all(newp,ec);fs::remove_all(cmp,ec);atomic_text(a.work_dir/"cleanup.done","complete\n");}
        std::cerr<<"\nPATCH COMPLETE\nexisting batches untouched\nmanifest="<<man<<"\nfeatures="<<a.feature_root<<"\n";return 0;
    }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 1;}
}
