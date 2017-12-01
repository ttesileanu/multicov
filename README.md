MultiCOV
========

This is a package for performing covariance analysis on alignments of biological sequences. In particular, it allows one to run various versions of the statistical coupling analysis (SCA) framework, as well as several different approaches to direct coupling analysis (DCA). The package can work with arbitrary alphabets (e.g., DNA, RNA, or protein) or with mixed alignments containing sequences from several alphabets.

This is also the software used to generate the analyses and figures used in the paper

Teşileanu, T., Colwell, L. J., & Leibler, S. (2015). Protein Sectors: Statistical Coupling Analysis versus Conservation. PLOS Computational Biology, 11(2), e1004091. doi:10.1371/journal.pcbi.1004091

The relevant files for that paper can be found in the `papers/sca_conservation` folder.

Installation
------------

Simply download and unzip the files in a folder of your choice. Most of the code is pure Matlab so it should work without any further trouble. It is recommended to run the `setup_paths.m` script before using any of the `multiCOV` functions. This makes sure to add all of the relevant folders to the paths so that all the functions can be properly accessed.

Some of the features of the package require compilation of C++ extensions. Estimating sequence weights, for example, is greatly accelerated by use of the C++ extensions, though it does work even without them. The functions related to Monte Carlo simulations rely on the C++ code and will not work if the extensions are not compiled.

The code for the extensions is in the `cpp` folder. On Unix-like machines (Mac or Linux), go to each of the folders and either run `make` (if you have it installed), or run the `bash` scripts called `mexify.sh`. Note that on Mac OS, you will need the Xcode command-line tools installed. On Windows, you may be able to compile the extensions using Cygwin (https://www.cygwin.com).

Usage
-----

All the functions in the `multiCOV` package are extensively documented. After running the `setup_paths.m` script, you can use Matlab's `help` command to get an explanation for what a function does and how it's used. You can, for example, type `help loadfasta` to learn about the function used to load FASTA alignments. You can also type `help alignment` to see short descriptions of all the functions in the `alignment` folder.

`multiCOV` contains a large number of functions, which can make it hard to navigate the package. A good way to get started is to look at the scripts that are included with the package. A short tutorial, `example.m`, on using `multiCOV` can be found in the `doc` folder. The scripts in `papers/sca_conservation`, that were used for the protein sectors paper, are also a good way to get an introduction to the code. You can also look at the scripts in the `experiments` folder.

Trouble?
--------

If you have any trouble or you think you've found a bug, please file a bug report at https://bitbucket.org/ttesileanu/multicov/issues.
