# UTR-Seq

UTR-Seq is a Matlab package designed for identifying short sequence elements
that are involved in controling the stability of RNA transcripts from data 
collected in Massively Parallel Reporter Assays (MPRA). 

You can find the full method description in:

[A Massively Parallel Reporter Assay of 3' UTR Sequences Identifies In Vivo
Rules for mRNA Degradation](https://www.ncbi.nlm.nih.gov/pubmed/29727622).
Rabani M, Pieper L, Chew GL, Schier AF.
Mol Cell. 2018 May 3;70(3):565.

## Installation

These instructions will get you a copy of the UTR-Seq analysis software 
to use on your with Matlab local machine. 

### Prerequisites

UTR-Seq was written and tested for Matlab R2018b. Matlab can be obtained and
installed from [Mathworks](https://www.mathworks.com/products/matlab.html).

### Installing

The package source code is downloaded from GitHub. Unzip the package into a
directory of your choise (this will be the installation directory).

```
unzip UTR-Seq-master.zip
```

Import the intstallation directory into matlab, by running (in Matlab):

```
addpath <installation path>
```

Now you can run all the commands in the installation directory (see example
below).

## Example data

The downloaded package includes an example of reporter sequences and RNA-Seq
counts that can be used for running and testing:

```
example/reporter_counts.txt
example/reporter_seq.txt
example/bg_seq.txt
example/test_seq.txt
```

### Fitting model on training set

In matlab, browse to the installation directory. Then run:

```
cd <installation directory>;
utrseq_fit_all(3, 'example/reporter_seq.txt', 'example/reporter_counts.txt',
[0 1 2 3 4 5 6 7 8 10 12], 'example/reporter_');
```

This will run all 3 steps of optimization, and place the results in the example
directory.

### Running on test set

After optimizing the model on the training set, you can test its predictions
on a new set of sequences (test set).
To do that, in matlab, browse to the same installation directory. Then run:

```
cd <installation directory>;
utrseq_predict_all('example/test_seq.txt', 'example/bg_seq.txt', ...
'example/repoters_kmer_regression/run_linear.out.mat', ...
'example/repoters_fit_peaks/peaks.mat', 'example/test_');
```

This will run all 2 steps of predictions, and place the results in the example
directory.


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

