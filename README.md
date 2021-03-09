1) Prerequisites

Install the Julia programming language, https://julialang.org/

Add a few packages: Optim, Statistics, Distributions, EmpiricalCDFs, Random, DelimitedFiles, SpecialFunctions, Random. 
To add a Julia package, start Julia in the command line, then type "]" to enter package management mode, then "add PackageName", finally press Backspace to exit package management and type "using PackageName" to precompile the package.


2) Input data

There are 3 types of input files: one for the mutant sites, one for the read depths of all sites, and those containing the tree branch lengths.

The mutated sites file contains 5 columns, separated by tabs:
1st column: read depth of site.
2nd-4th columns: readcounts of mutant types (like A, C, T for a site with wild type G). If less than 3 mutant types are observed, the other types should read "0". The order of the mutant types is not important, even consistent ordering for different lines is not necessary.
5th column: number of sites corresponding to the given mutant & total readcounts.

The readdepth file should contain 2 columns, separated by tabs: 
1st column: read depth value
2nd column: number of sites
The sum of the second column should be equal to the number of all sequenced sites.

Tree files:
Files containing the trees should be named "tree_1.txt", "tree_2.txt", ... However, the filename scheme can be easily edited to taste in the read_trees() function, line "treefile = ..."
Inside the files, there should be 3 columns, separated by tabs:
1st column: number of leaves in the branch
2nd column: number of all leaves in the tree
3rd column: length of the branch
Branches having the same number of leaves can be represented by a single line, the third column containing the sum of their branch lengths. The first column should still show the number of leaves in ONE branch.


3) Using the program

Running the code without arguments "julia estimate_parameters.jl" displays the possible arguments.
There are 3 mandatory parameters: the name of the file containing the mutant site readcounts, the name of the file containing the read depths of all sites, and the number of trees to use.
The software can calculate the loglikelihood value without any parameters optimization, i.e., the error rate and clonal peak size fixed to the input values ("-o off). Alternatively, either the error rate or the clonal peak size can be optimized for maximal loglikelihood ("-o error" and "-o clonalpeak"). In these cases, the input value of the corresponding parameter is neglected. Minimum and maximum values of the optimized parameter can be given by the -mine, -maxe, -minc, -maxc arguments.


4) Output

2 files are outputted. One contains the input parameters, the resulted mutation rate and loglikelihood values. If the error rate or the clonal peak size was optimized, the optimized value is shown, otherwise, the input value.
The other file contains the loglikelihood and mutation rate values of each individual tree.
Mutation rates are for the whole genome, to get sitewise values, they should be divided by the number of sequenced sites.
