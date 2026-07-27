#!/usr/bin/env python3
import argparse, csv, errno, io, os, re, shutil, subprocess, sys, tarfile, tempfile, time
from pathlib import Path
from typing import Iterable, Optional, Sequence, Set, Tuple, List

LIBRARY_ROOT_DEFAULT = Path('/pi/summer.thyme-umw/enamine-REAL-2.6billion')
ARCHIVE_TEMPLATE = 'condensed_params_and_db_{}.tar.gz'
EXPECTED_SUBCHUNKS = tuple(range(10))
SIMILARITY_REPORT_TEMPLATE = 'chunk_{:05d}_similarity_report.csv'
NAME_COLUMN_CANDIDATES = ('name','ligand_name','lig_name','ligand','title','id','identifier')

def parse_args():
    p=argparse.ArgumentParser(description='Filter Enamine parameter archives using blacklist and similarity-report ligand names.')
    p.add_argument('chunk', type=int)
    p.add_argument('similarity_root', type=Path)
    p.add_argument('output_root', type=Path)
    p.add_argument('--library-root', type=Path, default=LIBRARY_ROOT_DEFAULT)
    p.add_argument('--blacklist-name-column', default=None)
    p.add_argument('--similarity-name-column', default=None)
    p.add_argument('--similarity-status-column', default='status')
    p.add_argument('--removed-status', default='removed')
    p.add_argument('--allow-missing-archives', action='store_true')
    p.add_argument('--overwrite', action='store_true')
    p.add_argument('--keep-exclude-manifest', action='store_true')
    return p.parse_args()

def log(msg):
    print('[{}] {}'.format(time.strftime('%Y-%m-%d %H:%M:%S'), msg), flush=True)

def chunk_paths(root: Path, chunk: int) -> Tuple[int, Path]:
    s=chunk//100
    return s, root/str(s)/'{:05d}'.format(chunk)

def detect_column(fieldnames: Sequence[str], requested: Optional[str]) -> str:
    norm={n.strip().lower():n for n in fieldnames if n}
    if requested:
        k=requested.strip().lower()
        if k not in norm:
            raise ValueError('Requested column {!r} not found. Available: {}'.format(requested, ', '.join(fieldnames)))
        return norm[k]
    for c in NAME_COLUMN_CANDIDATES:
        if c in norm:
            return norm[c]
    raise ValueError('Could not auto-detect ligand-name column. Available: {}'.format(', '.join(fieldnames)))

def read_csv_names(handle, requested_name_column, status_column=None, required_status=None) -> Set[str]:
    r=csv.DictReader(handle)
    if not r.fieldnames: raise ValueError('CSV has no header row.')
    name_col=detect_column(r.fieldnames, requested_name_column)
    status_col=None
    if status_column is not None:
        norm={n.strip().lower():n for n in r.fieldnames if n}
        k=status_column.strip().lower()
        if k not in norm: raise ValueError('Status column {!r} not found.'.format(status_column))
        status_col=norm[k]
    names=set()
    for row in r:
        if status_col is not None and required_status is not None:
            if (row.get(status_col) or '').strip().lower()!=required_status.strip().lower():
                continue
        name=(row.get(name_col) or '').strip()
        if name:
            if '\n' in name or '\r' in name: raise ValueError('Ligand name contains newline: {!r}'.format(name))
            names.add(name)
    return names

def read_blacklist(path, requested_column):
    with path.open('r',encoding='utf-8-sig',newline='') as h:
        return read_csv_names(h,requested_column)

def locate_similarity_report(chunk_dir, chunk):
    base=SIMILARITY_REPORT_TEMPLATE.format(chunk)
    for c in (chunk_dir/(base+'.tar.gz'), chunk_dir/base):
        if c.is_file(): return c
    raise FileNotFoundError('Similarity report not found under {}'.format(chunk_dir))

def read_similarity_names(path, requested_name_column, status_column, removed_status):
    if path.name.endswith('.tar.gz'):
        with tarfile.open(str(path),'r:gz') as a:
            members=[m for m in a.getmembers() if m.isfile() and m.name.lower().endswith('.csv')]
            if len(members)!=1: raise ValueError('Expected exactly one CSV in {}, found {}.'.format(path,len(members)))
            ex=a.extractfile(members[0])
            if ex is None: raise ValueError('Could not read CSV from {}'.format(path))
            with ex, io.TextIOWrapper(ex,encoding='utf-8-sig',newline='') as h:
                return read_csv_names(h,requested_name_column,status_column,removed_status)
    with path.open('r',encoding='utf-8-sig',newline='') as h:
        return read_csv_names(h,requested_name_column,status_column,removed_status)

def glob_escape_literal(value):
    out=[]
    for ch in value:
        out.append('[*]' if ch=='*' else '[?]' if ch=='?' else '[[]' if ch=='[' else ch)
    return ''.join(out)

def write_exclude_manifest(path, names: Iterable[str]):
    with path.open('w',encoding='utf-8',newline='\n') as h:
        h.write('db.db\n')
        h.write('*lig_name_lis*\n')
        for n in sorted(set(names)):
            h.write('*{}*\n'.format(glob_escape_literal(n)))

def list_members(path):
    r=subprocess.run(['tar','-tzf',str(path)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    return [x for x in r.stdout.splitlines() if x]

def find_single_conf_root(members: Sequence[str]) -> str:
    roots=set(); pat=re.compile(r'^(?:\./)?condensed_params_and_db_[^/]+/single_conf_params/?$')
    for m in members:
        if pat.match(m.rstrip('/')): roots.add(m.rstrip('/'))
    if not roots:
        pat2=re.compile(r'^((?:\./)?condensed_params_and_db_[^/]+/single_conf_params)/')
        for m in members:
            z=pat2.match(m)
            if z: roots.add(z.group(1))
    if len(roots)!=1: raise ValueError('Expected one single_conf_params root, found {}: {}'.format(len(roots), ', '.join(sorted(roots))))
    return next(iter(roots))

def run_checked(cmd):
    log('Running: {}'.format(' '.join(str(x) for x in cmd)))
    subprocess.run(list(cmd),check=True)

def remove_tree_best_effort(path: Path):
    if not path.exists(): return
    def onerror(function, failed_path, exc_info): return
    try: shutil.rmtree(str(path),onerror=onerror)
    except Exception as e: log('WARNING: cleanup issue ignored: {}'.format(e))
    if path.exists(): log('WARNING: temporary directory remains for later cleanup: {}'.format(path))

def validate_output_archive(path: Path):
    r=subprocess.run(['tar','-tzf',str(path)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    if r.returncode!=0: raise RuntimeError('Archive validation failed for {}: {}'.format(path,r.stderr))
    members=[x for x in r.stdout.splitlines() if x]
    forbidden=[m for m in members if m.endswith('/db.db') or m=='db.db' or 'lig_name_lis' in m]
    outside=[m for m in members if '/single_conf_params/' not in m and not m.rstrip('/').endswith('/single_conf_params')]
    if forbidden: raise RuntimeError('Forbidden members remain: {}'.format(', '.join(forbidden[:10])))
    if outside: raise RuntimeError('Members outside single_conf_params remain: {}'.format(', '.join(outside[:10])))

def main():
    a=parse_args()
    if a.chunk<0 or a.chunk>53083:
        print('ERROR: chunk must be 0-53083.',file=sys.stderr); return 2
    superchunk, source_chunk=chunk_paths(a.library_root.resolve(),a.chunk)
    _, similarity_chunk=chunk_paths(a.similarity_root.resolve(),a.chunk)
    _, output_chunk=chunk_paths(a.output_root.resolve(),a.chunk)
    if not source_chunk.is_dir(): print('ERROR: source chunk missing: {}'.format(source_chunk),file=sys.stderr); return 2
    output_chunk.mkdir(parents=True,exist_ok=True)
    try:
        sim_report=locate_similarity_report(similarity_chunk,a.chunk)
        blacklist=read_blacklist(source_chunk/'blacklist_file.csv',a.blacklist_name_column)
        similar=read_similarity_names(sim_report,a.similarity_name_column,a.similarity_status_column,a.removed_status)
    except Exception as e:
        print('ERROR: {}'.format(e),file=sys.stderr); return 1
    excluded=blacklist|similar
    log('Chunk {:05d}: {:,} blacklist + {:,} similarity = {:,} unique exclusions'.format(a.chunk,len(blacklist),len(similar),len(excluded)))
    manifest=output_chunk/'chunk_{:05d}_tar_excludes.txt'.format(a.chunk)
    write_exclude_manifest(manifest,excluded)
    archives=[]; missing=[]
    for i in EXPECTED_SUBCHUNKS:
        p=source_chunk/ARCHIVE_TEMPLATE.format(i)
        (archives if p.is_file() else missing).append((i,p) if p.is_file() else p)
    if missing and not a.allow_missing_archives:
        print('ERROR: missing archives:\n  {}'.format('\n  '.join(str(x) for x in missing)),file=sys.stderr); return 1
    temp_dir=Path(tempfile.mkdtemp(prefix='.chunk_{:05d}_tar_filter_'.format(a.chunk),dir=str(output_chunk)))
    try:
        for pos,(subchunk,src) in enumerate(archives,1):
            out=output_chunk/src.name; partial=output_chunk/(src.name+'.partial')
            if out.exists() and not a.overwrite:
                log('Exists; skipping ({}/{}): {}'.format(pos,len(archives),out)); continue
            if partial.exists(): partial.unlink()
            exdir=temp_dir/'subchunk_{}'.format(subchunk); exdir.mkdir(parents=True,exist_ok=True)
            root=find_single_conf_root(list_members(src))
            run_checked(['tar','-xzf',str(src),'-C',str(exdir),'--wildcards','--no-anchored','--exclude-from',str(manifest),root])
            extracted=exdir/root
            if not extracted.is_dir(): raise RuntimeError('Expected extracted directory missing: {}'.format(extracted))
            run_checked(['tar','-czf',str(partial),'-C',str(exdir),root])
            validate_output_archive(partial)
            os.replace(str(partial),str(out))
            remove_tree_best_effort(exdir)
            log('Wrote {}'.format(out))
        marker=output_chunk/'chunk_{:05d}_archives.complete'.format(a.chunk)
        with marker.open('w',encoding='utf-8') as h:
            h.write('chunk={:05d}\n'.format(a.chunk)); h.write('superchunk={}\n'.format(superchunk)); h.write('unique_excluded_names={}\n'.format(len(excluded))); h.write('archives_written={}\n'.format(len(archives)))
    except Exception as e:
        print('ERROR: {}'.format(e),file=sys.stderr); return 1
    finally:
        remove_tree_best_effort(temp_dir)
        if not a.keep_exclude_manifest:
            try: manifest.unlink()
            except OSError: pass
    log('Completed chunk {:05d}.'.format(a.chunk)); return 0

if __name__=='__main__': sys.exit(main())
