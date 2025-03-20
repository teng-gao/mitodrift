# MitoDrift

# Installation
For now:
```R
devtools::install_local('.')
```

# Usage
## Input format
1. `mut_dat` includes the below columns, where `a` is the alternate allele count and `d` is the total depth.

| cell      | variant       | a    | d    |
|-----------|---------------|------|------|
| PD45534aj | MT_10448_T_C  | 0    | 2348 |
| PD45534aj | MT_11787_T_C  | 1462 | 2348 |
| PD45534aj | MT_11790_T_C  | 0    | 2338 |
| PD45534aj | MT_12242_A_G  | 0    | 2277 |
| PD45534aj | MT_12630_G_A  | 0    | 2382 |
| PD45534aj | MT_15791_A_G  | 0    | 2459 |

Note that depth information should be provided even if there is no alternate allele detected at the variant position (a=0), so that there is a row for each cell x variant combination.

## ML tree search
First, run `inst/bin/build_tree.r` to perform ML tree search.

```
Usage: build_tree.r [options]


Options:
	-m CHARACTER, --mut_dat=CHARACTER
		Mutation data file

	-i INTEGER, --max_iter=INTEGER
		Maximum number of iteration

	-o CHARACTER, --outfile=CHARACTER
		Output file .RDS

	-p INTEGER, --ncores=INTEGER
		Number of cores to use

	-e DOUBLE, --eps=DOUBLE
		mutation rate

	-s DOUBLE, --seq_err=DOUBLE
		sequencing error rate

	-n INTEGER, --n_pop=INTEGER
		Population size

	-g INTEGER, --n_gen=INTEGER
		Number of generations

	-f CHARACTER, --freq_dat=CHARACTER
		Mutation frequency table

	-r LOGICAL, --resume=LOGICAL
		Whether to resume from output file

	-h, --help
		Show this help message and exit
```
Example command:
```
Rscript $home/mitodrift/mitodrift/inst/bin/build_tree.r \
    -m $infile \
    -o $outfile \
    -p $ncores -i 1000 -e 0.005 -n 600 -g 100
```

## MCMC
Once ML tree search finishes, use `inst/bin/run_mcmc.r` to perform MCMC in order to estimate clade confidence.

```
Usage: run_mcmc.r [options]


Options:
	-t CHARACTER, --tree_file=CHARACTER
		Tree list file from mitodrift ML module

	-i INTEGER, --max_iter=INTEGER
		Maximum number of iteration

	-c INTEGER, --nchains=INTEGER
		Maximum number of iteration

	-o CHARACTER, --outfile=CHARACTER
		Output file .RDS

	-p INTEGER, --ncores=INTEGER
		Number of cores to use

	-h, --help
		Show this help message and exit
```
Example command:
```
Rscript $home/mitodrift/mitodrift/inst/bin/run_mcmc.r \
    -t $infile \
    -o $outfile \
    -i 1000 \
    -c 50 \
    -p $ncores
```
