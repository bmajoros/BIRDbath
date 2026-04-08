# BIRDbath
Bayesian model for estimating variant effects from structured, multi-pool STARR-seq experiments

## Installation
To use the BIRDbath model, first install CmdStan from the STAN website:
https://mc-stan.org/tools/index.html

Clone the supporting python repository via the GitHub clone:

```
git clone --recursive \
  https://github.com/bmajoros/python.git
```

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
To run BIRDbath, run the wrapper script on your data (`wrapper.sh`), modifying the variables which indicate file and path names:
- The directory with the BIRDbath input files, where the BIRDbath output files will be generated
- The directory with the BIRDbath model
- The directory with the GitHub repository: https://github.com/bmajoros/python
- The directory where STAN is installed
- The name of the input file
- The zero-based index of the first variant
- The zero-based index of the last variant

The directory where BIRDbath runs should include a subdirectory called "outputs" and a subdirectory called "temp".

The data should be imported from a LISP S-Expression format for each variant. Each variant has information regarding the pools where it is present, the size of the pools, the expected alternate allele frequency in the pool, and the counts of the alleles in the DNA and RNA libraries:

```
(variant (id chr1:729089:G:C)
        (pool 1
                (freq 0.5)
                (size 5)
                (DNA (ref 48) (alt 44))
                (RNA (ref 82) (alt 36))
        )
        (pool 5
                (freq 0.4)
                (size 5)
                (DNA (ref 20) (alt 14))
                (RNA (ref 11) (alt 2))
        )
        (pool 10
                (freq 0.3)
                (size 5)
                (DNA (ref 18) (alt 10))
                (RNA (ref 21) (alt 3))
        )
        (pool 20
                (freq 0.4)
                (size 5)
                (DNA (ref 17) (alt 13))
                (RNA (ref 34) (alt 13))
        )
)
```

After running BIRDbath, a tab-separated file is generated with 5 fields:
- VariantID: identifier of the variant tested.
- Effect size (theta): posterior median of effect size under the model. 
- LowerCI: lower bound of 95% credible interval. 
- UpperCI: upper bound of 95% credible interval. 
- P<sub>reg<sub>: posterior probability that the variant has a regulatory effect. 

For the above variant data, the output line might look like this:

```
chr1:729089:G:C 0.482   0.351   0.673   0.998
```
