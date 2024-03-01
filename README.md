# SigFitTest
SigFitTest is a package to create realistic synthetic mutational catalogs and to evaluate tools for fitting mutational signatures.


## How to use this package
This package is written in Python. Code snippets to run individual fitting tools are primarily written in R.

There are four directories:
* `code`: the actual scripts,
* `input`: reference signature catalogs,
* `cosmic tissue data`: signature activity (using the COSMICv3 catalog) in tissue data provided by the COSMIC website,
* `real mutational catalogs`: SBS96 mutational catalogs of four samples from the PCAWG project and the average signature weights estimated for them by sigLASSO, SigProfilerAssignment, and MuSiCal.

To run the package, you need to: (1) download and unpack it, (2) go to the directory `code`, (3) run `main.py`.

In `main.py`, there are five main functions that you can use (comment or uncomment the corresponding lines at the end of this script according to which functions you want to run):
* `generate_synthetic_catalogs(...)`: Generates mutational catalogs based on single base substitutions (SBS) where the mutations are classified in 96 different contexts. Signature weights in the synthetic data are obtained from the real tissue results in the directory `cosmic tissue data`. These synthetic catalogs can be later used for your own evaluation of signature fitting (or signature extraction). The output data (synthetic catalogs and true signature weights) are stored in the directory `generated_data`.
* `fit_with_cosmic3_synthetic_simple(...)`: Generates simple mutational catalogs with pre-defined signature contributions (same for every sample), uses them as input for a signature fitting tool, and evaluates the results produced by this tool. All COSMICv3 signatures are used as a reference for fitting.
* `fit_with_cosmic3_synthetic(...)`: Generates realistic mutational catalogs, uses them as input for a signature fitting tool, and evaluates the results produced by this tool. All COSMICv3 signatures are used as a reference for fitting. Sample output of this function for two tools (SigsPack and SigProfilerAssignment) for two cancer types (Head-SCC and ColoRect-AdenoCA) for two cohorts and three different mutation counts (100, 2000, 50000) is included in this repository (files reference_results-WGS-set6-SigsPack.dat and reference_results-WGS-set6-SPA.dat). When parameter `out_of_reference_weights` is provided, empirical signature weights are scaled down appropriately and the remaining mutations are driven by out-of-reference signatures (signatures from COSMICv3.3.1 that are absent in COSMICv3 that is used for fitting).
* `prune_reference_and_fit_synthetic()`: Uses results of `fit_with_cosmic3_synthetic()` run previously to prune the list of reference signatures (only signatures that are active in sufficiently many samples are kept), runs a fitting tool and evaluates the results. Note that the zip files with the previously computed results (e.g., `signature_results-WGS-SigsPack-set12-Head-SCC.zip`) must be located in the directory `code`.
* `fit_with_cosmic3_subsampled_real_catalogs()`: Generates mutational catalogs by subsampling from the real mutational catalogs stored in the directory `real mutational catalogs`. In this case, the fitting results are evaluated by comparing them with the ground truth result stored in the same directory (this result is obtained by analyzing the complete real catalogs using three different fitting methods and averaging the obtained results).


### Practical points:
* To run a fitting tool, the tool needs to be installed first following the instructions that come with that tool.
* To find the available cancer types, see the directory `cosmic tissue data`.
* Several important variables are set in the file `MS_config.py` in the directory `code`:
  * `tool` (tool that will be run and evaluated),
  * `N_samples` (how many synthetic samples there are in each cohort; default: 100),
  * `num_realizaions` (how many synthetic cohorts will be generated; default: 5),
  * `timeout_time` (fitting tools are stopped if they do not finish until this time [in seconds]; you may need to increase it for slow tools such as mmsig, for example; default: 1800),
  * `num_muts_list_short` (how many mutations there will be in the synthetic catalogs; this is a list, catalog generation & fitting & evaluation are run separately for each number of mutations provided; default: [100, 2000, 50000]),
  * `input_signatures` (the reference signatures that are used to generate the catalogs; default: WGS COSMICv3.3)
* Variable `code_name` is used to easily recognize the results obtained by the provided functions under various settings.
* The content of `stdout` and `stderr` is saved in the files `stdout-tool_name.txt` and `stderr-tool_name.txt`. The content of these files can help you with troubleshooting.
