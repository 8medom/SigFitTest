# SigFitTest
SigFitTest is a package to create synthetic mutational catalogs and to evaluate tools for fitting mutational signatures.

## How to use this package
This package is written in Python. Code snippets to run individual fitting tools are primarily written in R.

There are three directories:
* `code`: the actual scripts
* `input`: reference signature catalogs
* `cosmic tissue data`: signature activity in tissue data analyzed by the COSMIC website

To run the package, you need to: (1) download it, (2) go to the directory `code`, (3) run `main.py`.

In `main.py', there are three main functions that you can use:
* `generate_catalogs(...)`: to only generate mutational catalogs based on single base substitutions (SBS) where the mutations are classified in 96 different contexts.
* `fit_and_evaluate(...)`: to generate mutational catalogs, use them as input for a signature fitting tool, and evaluate the results produced by this tool. All COSMICv3 signatures are used as a reference for fitting.
* `two_step_fit_and_evaluate()`: uses results of previously run `fit_and_evaluate()` to prune the list of reference signatures (only signatures that are active in sufficiently many samples are kept).

To run a fitting tool, the tool needs to be installed first following the instructions that come with each tool.
