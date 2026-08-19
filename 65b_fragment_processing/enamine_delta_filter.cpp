// enamine_delta_filter.cpp
// Compare an old Enamine shard to the corresponding new shard by ligand ID
// (column 2), emitting only novel records. Uses disk hash partitions so the
// old shard does not need to fit in RAM.
//
// Build:
// g++ -O3 -std=c++17 enamine_delta_filter.cpp -lbz2 -o enamine_delta_filter
//
// Example:
// ./enamine_delta_filter \
//   --old-root /old/64b/M/H32 \
//   --new-root /new/90b/M/H32 \
//   --work-dir /scratch/$USER/M_H32_delta \
//   --output /delta/M/H32/H32_novel.cxsmiles.bz2 \
//   --buckets 4096

#include <bzlib.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
#include <algorithm>

namespace fs = std::filesystem;

struct Args {
    fs::path old_root, new_root, work_dir, output;
    std::size_t buckets = 4096;
    std::size_t bucket_buffer_bytes = 16384;
    bool keep_work = false;
};

[[noreturn]] void usage(const char* p, int rc=1) {
    std::cerr <<
      "Usage: " << p << " --old-root DIR --new-root DIR --work-dir DIR "
      "--output FILE.bz2 [--buckets N] [--keep-work]\n";
    std::exit(rc);
}

Args parse_args(int argc, char** argv) {
    Args a;
    for (int i=1; i<argc; ++i) {
        std::string x=argv[i];
        auto val=[&](const char* f)->std::string {
            if (i+1>=argc) throw std::runtime_error(std::string("Missing value for ")+f);
            return argv[++i];
        };
        if (x=="--old-root") a.old_root=val("--old-root");
        else if (x=="--new-root") a.new_root=val("--new-root");
        else if (x=="--work-dir") a.work_dir=val("--work-dir");
        else if (x=="--output") a.output=val("--output");
        else if (x=="--buckets") a.buckets=std::stoull(val("--buckets"));
        else if (x=="--keep-work") a.keep_work=true;
        else if (x=="-h" || x=="--help") usage(argv[0],0);
        else throw std::runtime_error("Unknown option: "+x);
    }
    if (a.old_root.empty() || a.new_root.empty() || a.work_dir.empty() || a.output.empty())
        usage(argv[0]);
    return a;
}

uint64_t fnv1a64(const std::string& s) {
    uint64_t h=14695981039346656037ULL;
    for (unsigned char c: s) { h^=c; h*=1099511628211ULL; }
    return h;
}

// Same basic interpretation as the Python pipeline:
// column 1 = smiles/cxsmiles, column 2 = ligand_name.
bool ligand_key(const std::string& line, std::string& key) {
    if (line.empty()) return false;

    std::size_t p1=line.find('\t');
    if (p1!=std::string::npos) {
        std::size_t p2=line.find('\t',p1+1);
        std::string smiles=line.substr(0,p1);
        key=(p2==std::string::npos) ? line.substr(p1+1)
                                    : line.substr(p1+1,p2-p1-1);
        if (smiles=="smiles" || smiles=="SMILES" ||
            smiles=="cxsmiles" || smiles=="CXSMILES") return false;
        return !key.empty();
    }

    std::size_t s1=line.find_first_not_of(" \r\n");
    if (s1==std::string::npos) return false;
    std::size_t e1=line.find_first_of(" \t\r\n",s1);
    if (e1==std::string::npos) return false;
    std::string smiles=line.substr(s1,e1-s1);
    std::size_t s2=line.find_first_not_of(" \t\r\n",e1);
    if (s2==std::string::npos) return false;
    std::size_t e2=line.find_first_of(" \t\r\n",s2);
    key=line.substr(s2,e2==std::string::npos ? std::string::npos : e2-s2);
    if (smiles=="smiles" || smiles=="SMILES" ||
        smiles=="cxsmiles" || smiles=="CXSMILES") return false;
    return !key.empty();
}

class BzReader {
    FILE* fp_=nullptr;
    BZFILE* bz_=nullptr;
    fs::path path_;
    bool eof_=false;
    std::string pending_;
public:
    explicit BzReader(const fs::path& p):path_(p) {
        fp_=std::fopen(p.string().c_str(),"rb");
        if (!fp_) throw std::runtime_error("Cannot open "+p.string());
        int e=BZ_OK;
        bz_=BZ2_bzReadOpen(&e,fp_,0,0,nullptr,0);
        if (e!=BZ_OK) throw std::runtime_error("BZ2 open failed "+p.string());
    }
    ~BzReader() {
        if (bz_) { int e=BZ_OK; BZ2_bzReadClose(&e,bz_); }
        if (fp_) std::fclose(fp_);
    }
    bool getline(std::string& out) {
        out.clear();
        for (;;) {
            auto n=pending_.find('\n');
            if (n!=std::string::npos) {
                out=pending_.substr(0,n);
                pending_.erase(0,n+1);
                if (!out.empty() && out.back()=='\r') out.pop_back();
                return true;
            }
            if (eof_) {
                if (pending_.empty()) return false;
                out.swap(pending_);
                if (!out.empty() && out.back()=='\r') out.pop_back();
                return true;
            }
            char buf[1<<20];
            int e=BZ_OK;
            int got=BZ2_bzRead(&e,bz_,buf,sizeof(buf));
            if (got>0) pending_.append(buf,(std::size_t)got);
            if (e==BZ_STREAM_END) eof_=true;
            else if (e!=BZ_OK) throw std::runtime_error("BZ2 read failed "+path_.string());
        }
    }
};

class BzWriter {
    FILE* fp_=nullptr;
    BZFILE* bz_=nullptr;
    fs::path path_;
public:
    explicit BzWriter(const fs::path& p):path_(p) {
        if (!p.parent_path().empty()) fs::create_directories(p.parent_path());
        fp_=std::fopen(p.string().c_str(),"wb");
        if (!fp_) throw std::runtime_error("Cannot create "+p.string());
        int e=BZ_OK;
        bz_=BZ2_bzWriteOpen(&e,fp_,9,0,30);
        if (e!=BZ_OK) throw std::runtime_error("BZ2 output open failed");
    }
    ~BzWriter() { close(); }
    void line(const std::string& s) {
        std::string t=s; t.push_back('\n');
        int e=BZ_OK;
        BZ2_bzWrite(&e,bz_,const_cast<char*>(t.data()),(int)t.size());
        if (e!=BZ_OK) throw std::runtime_error("BZ2 write failed "+path_.string());
    }
    void close() {
        if (bz_) { int e=BZ_OK; BZ2_bzWriteClose(&e,bz_,0,nullptr,nullptr); bz_=nullptr; }
        if (fp_) { std::fclose(fp_); fp_=nullptr; }
    }
};

fs::path bucket_path(const fs::path& root, std::size_t b) {
    char name[64];
    std::snprintf(name,sizeof(name),"bucket_%06zu.dat",b);
    return root/name;
}

class BucketWriter {
    fs::path root_;
    std::vector<std::string> buf_;
    std::size_t limit_;
public:
    BucketWriter(const fs::path& r, std::size_t n, std::size_t lim)
      : root_(r),buf_(n),limit_(lim) { fs::create_directories(root_); }
    ~BucketWriter(){ try{ flush_all(); }catch(...){ } }
    void add(std::size_t b,const std::string& s) {
        auto& x=buf_.at(b);
        x+=s; x.push_back('\n');
        if (x.size()>=limit_) flush(b);
    }
    void flush(std::size_t b) {
        auto& x=buf_.at(b);
        if (x.empty()) return;
        std::ofstream f(bucket_path(root_,b),std::ios::binary|std::ios::app);
        if (!f) throw std::runtime_error("Cannot append partition");
        f.write(x.data(),(std::streamsize)x.size());
        x.clear();
    }
    void flush_all(){ for(std::size_t i=0;i<buf_.size();++i) flush(i); }
};

std::vector<fs::path> bzfiles(const fs::path& root) {
    std::vector<fs::path> v;
    for (auto const& e: fs::recursive_directory_iterator(root)) {
        if (!e.is_regular_file()) continue;
        std::string s=e.path().filename().string();
        if (s.size()>=4 && s.substr(s.size()-4)==".bz2") v.push_back(e.path());
    }
    std::sort(v.begin(),v.end());
    return v;
}

uint64_t partition_old(const Args& a,const fs::path& root) {
    BucketWriter w(root,a.buckets,a.bucket_buffer_bytes);
    uint64_t n=0;
    std::string line,key;
    auto files=bzfiles(a.old_root);
    for (std::size_t i=0;i<files.size();++i) {
        BzReader r(files[i]);
        while(r.getline(line)) {
            if (!ligand_key(line,key)) continue;
            w.add(fnv1a64(key)%a.buckets,key);
            if (++n%10000000ULL==0) std::cerr<<"[old] "<<n<<"\n";
        }
        std::cerr<<"[old] file "<<i+1<<"/"<<files.size()<<" "<<files[i]<<"\n";
    }
    w.flush_all();
    return n;
}

uint64_t partition_new(const Args& a,const fs::path& root) {
    BucketWriter w(root,a.buckets,a.bucket_buffer_bytes);
    uint64_t n=0;
    std::string line,key;
    auto files=bzfiles(a.new_root);
    for (std::size_t i=0;i<files.size();++i) {
        BzReader r(files[i]);
        while(r.getline(line)) {
            if (!ligand_key(line,key)) continue;
            w.add(fnv1a64(key)%a.buckets,key+"\t"+line);
            if (++n%10000000ULL==0) std::cerr<<"[new] "<<n<<"\n";
        }
        std::cerr<<"[new] file "<<i+1<<"/"<<files.size()<<" "<<files[i]<<"\n";
    }
    w.flush_all();
    return n;
}

struct Counts { uint64_t compared=0,present=0,novel=0; };

Counts compare(const Args& a,const fs::path& oldp,const fs::path& newp) {
    BzWriter out(a.output);
    std::unordered_set<std::string> old_ids;
    old_ids.max_load_factor(0.8f);
    Counts c;

    for(std::size_t b=0;b<a.buckets;++b) {
        old_ids.clear();

        fs::path op=bucket_path(oldp,b), np=bucket_path(newp,b);

        if (fs::exists(op)) {
            std::ifstream f(op);
            std::string id;
            while(std::getline(f,id)) {
                if (!id.empty() && id.back()=='\r') id.pop_back();
                if (!id.empty()) old_ids.insert(id);
            }
        }

        if (fs::exists(np)) {
            std::ifstream f(np);
            std::string rec;
            while(std::getline(f,rec)) {
                if (!rec.empty() && rec.back()=='\r') rec.pop_back();
                auto t=rec.find('\t');
                if (t==std::string::npos) continue;
                std::string id=rec.substr(0,t);
                std::string original=rec.substr(t+1);
                ++c.compared;
                if (old_ids.find(id)!=old_ids.end()) ++c.present;
                else { out.line(original); ++c.novel; }
            }
        }

        if (!a.keep_work) {
            std::error_code ec;
            fs::remove(op,ec); fs::remove(np,ec);
        }

        if ((b+1)%32==0 || b+1==a.buckets)
            std::cerr<<"[compare] "<<b+1<<"/"<<a.buckets
                     <<" present="<<c.present<<" novel="<<c.novel<<"\n";
    }
    out.close();
    return c;
}

int main(int argc,char** argv) {
    try {
        Args a=parse_args(argc,argv);
        fs::create_directories(a.work_dir);
        fs::path oldp=a.work_dir/"old", newp=a.work_dir/"new";

        auto nonempty=[](const fs::path& p){
            return fs::exists(p) &&
                   fs::directory_iterator(p)!=fs::directory_iterator();
        };
        if (nonempty(oldp)||nonempty(newp))
            throw std::runtime_error("Temporary partition directory is not empty; use a fresh --work-dir.");

        std::cerr<<"old="<<a.old_root<<"\nnew="<<a.new_root
                 <<"\nwork="<<a.work_dir<<"\nout="<<a.output
                 <<"\nbuckets="<<a.buckets<<"\n";

        uint64_t nold=partition_old(a,oldp);
        uint64_t nnew=partition_new(a,newp);
        Counts c=compare(a,oldp,newp);

        if (!a.keep_work) {
            std::error_code ec;
            fs::remove_all(oldp,ec);
            fs::remove_all(newp,ec);
        }

        std::cerr<<"\nDONE\n"
                 <<"old records indexed: "<<nold<<"\n"
                 <<"new records scanned: "<<nnew<<"\n"
                 <<"already present: "<<c.present<<"\n"
                 <<"novel records: "<<c.novel<<"\n"
                 <<"output: "<<a.output<<"\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<"ERROR: "<<e.what()<<"\n";
        return 1;
    }
}
