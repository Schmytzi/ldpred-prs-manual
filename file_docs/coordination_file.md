# Documentation of Coordination File Format

LDpred generates a coordination file at the end of the `coord` step.
It's a Hierarchical Data Format Version 5 (HDFv5) file.

This is a binary format.
You can read HDFv5 files in python using the package `h5py`.

```
/
├── cord_data [sic!]
│   ├── chrom_1 [table of SNP infos]
│   ├── chrom_2
│   └── ...
└── sum_stats
    ├── chrom_1 [table of GWAS statistics]
    ├── chrom_2
    └── ...
```