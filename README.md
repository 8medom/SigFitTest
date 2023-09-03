# SigFitTest
SigFitTest is a package to create synthetic mutational catalogs and to evaluate tools for fitting mutational signatures.


## How to use this package
This package is written in Python. Code snippets to run individual fitting tools are primarily written in R.

There are three directories:
* `code`: the actual scripts
* `input`: reference signature catalogs
* `cosmic tissue data`: signature activity in tissue data analyzed by the COSMIC website

To run the package, you need to: (1) download it, (2) go to the directory `code`, (3) run `main.py`.

In `main.py', there are three main functions that you can use (comment or uncomment the corresponding lines at the end of this script according to which functions you want to run):
* `generate_synthetic_catalogs(...)`: to only generate mutational catalogs based on single base substitutions (SBS) where the mutations are classified in 96 different contexts. These synthetic catalogs can be later used for your own evaluation of signature fitting (or signature extraction).
* `fit_cosmic3_and_evaluate(...)`: to generate mutational catalogs, use them as input for a signature fitting tool, and evaluate the results produced by this tool. All COSMICv3 signatures are used as a reference for fitting.
* `prune_reference_fit_evaluate()`: uses results of `fit_and_evaluate()` run previously to prune the list of reference signatures (only signatures that are active in sufficiently many samples are kept), runs a fitting tool and evaluates the results.

Important points:
* To run a fitting tool, the tool needs to be installed first following the instructions that come with each tool.
* Several important variables are set in the file `MS_config.py` in the directory `code`:
  * `tool` (tool that will be run and evaluated),
  * `N_samples` (how many synthetic samples there are in each cohort; default: 100),
  * `num_realizaions` (how many synthetic cohorts will be generated; default: 50),
  * `timeout_time` (fitting tools are stopped if they do not finish until this time [in seconds]; you may need to increase it for slow tools such as mmsig, for example; default: 1800),
  * `num_muts_list_short` (how many mutations there will be in the synthetic catalogs; this is a list, catalog generation & fitting & evaluation are run separately for each number of mutations provided; default: [100, 2000, 50000]),
  * `input_signatures` (the reference signatures that are used to generate the catalogs; default: COSMICv3)
