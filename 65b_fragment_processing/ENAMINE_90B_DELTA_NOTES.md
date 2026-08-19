# 90B delta-filter workflow

The current Python feature pipeline interprets column 1 as SMILES/CXSMILES and
column 2 as the ligand name/identifier. The C++ filter follows the same format.

Recommended flow:

1. Compare a corresponding old/new shard, e.g. old `M/H32` vs new `M/H32`.
2. Emit only exact ligand IDs absent from the 64.9B shard.
3. Feed that `.cxsmiles.bz2` file into the existing Python/RDKit feature pipeline.
4. Keep the same RDKit version and feature/hash settings so the new feature
   hashes remain compatible with the existing 64.9B feature library.

Critical validation before scaling:
- Verify Enamine ligand IDs are stable between releases.
- Verify the H-sharding is stable enough that old `M/H32` compounds retained in
  the 90B release are still found in new `M/H32`.

If H-sharding changed, comparing only H32-to-H32 would create false "novel"
calls. In that case, build a broader old-library membership index instead.

Build:
```bash
g++ -O3 -std=c++17 enamine_delta_filter.cpp -lbz2 -o enamine_delta_filter
```

Example:
```bash
./enamine_delta_filter \
  --old-root /work/umassmed/thymelab_umassmed/compressed_unmodified_64b_library/M/H32 \
  --new-root /PATH/TO/90b_library/M/H32 \
  --work-dir /scratch/$USER/enamine_90b_delta/M_H32 \
  --output /work/umassmed/thymelab_umassmed/90b_delta/M/H32/H32_novel.cxsmiles.bz2 \
  --buckets 4096
```

Then point the existing node runner at the directory containing
`H32_novel.cxsmiles.bz2`.
