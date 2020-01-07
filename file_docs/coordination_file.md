# Documentation of Coordination File Format

## General Information
LDpred generates a coordination file at the end of the `coord` step.
It's a Hierarchical Data Format Version 5 (HDFv5) file.

## Reading the File
This is a binary format.
You can read HDFv5 files in python using the package `h5py`:

```python
import h5py

# File acts as a dict with additional methods
file = h5py.File('filename', 'r')

# Do stuff...

# It's still a file object! Close it when you're done!
file.close()
```

## Data Structure
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