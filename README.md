# MitoDrift

# Installation
For now:
```R
devtools::install_local('{repo_dir}')
```

# Usage
```
Usage: run_mitodrift_cpp.r [options]

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

	-h, --help
		Show this help message and exit
```


