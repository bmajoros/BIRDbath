# BIRDbath
Bayesian model for estimating variant effects from structured, multi-pool STARR-seq experiments

## Installation
To use the BIRDbath model, first install CmdStan from the STAN website:
https://mc-stan.org/tools/index.html


Clone this repo via the command GitHub clone:

```
git clone --recursive \
  https://github.com/bmajoros/BIRDbath.git
```

Compile the STAN model using the following command:

```
make all
```

## Usage
To run BIRDbath, run the wrapper script on your data (wrapper.sh).


The data should be imported from a LISP S-Expression format for each variant. Each variant has information regarding the pools where it is present, the size of the pools, the expected alternate allele frequency in the pool, and the counts of the alleles in the DNA and RNA libraries:
