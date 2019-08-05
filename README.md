# UTR-Seq

UTR-Seq is a Matlab package for identifying short sequence elements
that regulate the stability of RNA transcripts.
It is designed to work on data that was collected by Massively Parallel
Reporter Assays (MPRA) in order to define the elements. The resulting
model can then apply to any input sequence of choice.

You can find the full method description in:

[A Massively Parallel Reporter Assay of 3' UTR Sequences Identifies In Vivo
Rules for mRNA Degradation](https://www.ncbi.nlm.nih.gov/pubmed/29727622).
Rabani M, Pieper L, Chew GL, Schier AF.
Mol Cell. 2018 May 3;70(3):565.

## Installation

These instructions will get you a copy of the UTR-Seq analysis software 
to use with Matlab on your local machine. 

### Prerequisites

We developed and tested UTR-Seq with Matlab R2018b. Matlab can be obtained and
installed from [Mathworks](https://www.mathworks.com/products/matlab.html).

### Installing

Determine in which directory you would like to install the package. 
This will be the installation directory.

Download the package source code from GitHub. Copy the zip file into the
installation directory, and unzip the package.

```
unzip UTR-Seq-master.zip
```

Import the intstallation directory into matlab, by running (in Matlab):

```
addpath <path to installation directory>;
```

Now you can run in matlab all the source code that is saved in the 
installation directory (see example below).

## Example data

The package that you downloaded includes an example for using with the
UTR-Seq package.

There are 4 files in the example directory: a set of reporter sequences, RNA-Seq
counts for those reporters, a set of test sequences and a test of background
sequences:

```
example/reporter_counts.txt
example/reporter_seq.txt
example/test_seq.txt
example/bg_seq.txt
```

In the example, we will apply two analysis steps to this data:
first, we will fit a new model by using the reporter data,
and second, we will use this model to predict sequence
elements (motifs) in the set of test sequences.


### Step 1: fitting a model on the training set

In matlab, browse to the installation directory. Then run:

```
cd <installation directory>;
utrseq_fit_all(3, 'example/reporter_seq.txt', 'example/reporter_counts.txt',
[0 1 2 3 4 5 6 7 8 10 12], 'example/reporter_');
```

This will run all 3 steps of optimization, and will place the results
in the example directory.

Note that some of the optimization steps can take a long time
(up to several hours) to complete.

### Step 2: predicting motifs in the test set

After optimizing the model on the training set, you can test its predictions
on a new set of sequences (the test set).
To do that, in matlab, browse to the same installation directory. Then run:

```
cd <installation directory>;
utrseq_predict_all('example/test_seq.txt', 'example/bg_seq.txt', ...
'example/repoters_kmer_regression/run_linear.out.mat', ...
'example/repoters_fit_peaks/peaks.mat', 'example/test_');
```

This will run 2 steps of predictions, and place the results in the example
directory.


## License

This project is licensed under the MIT License - see source files for details.

