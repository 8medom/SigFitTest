# SigFitTest
SigFitTest is a Python package to create realistic synthetic mutational catalogs and to evaluate tools for fitting mutational signatures.


## Requirements
To run SigFitTest, you need a standard computer and a working Python installation with the following packages: numpy, scipy, pandas (standard libraries for scientific computing), os (operating system utilities), time, subprocess, threading, shlex (to execute external scripts), and lzma (for compression). To use a fitting tool, the tool needs to be installed first following the instructions that come with that tool. For tools written in R, a working R installation is also necessary. Depending on the chosen tool, the generation and analysis of a cohort with 100 samples can take from a few seconds to more than an hour (slowest methods: mmsig, signature.tools.lib).


## How to use this package
This package is written in Python. The code snippets for running individual fitting tools are written in R and Python.

There are four directories:
* `code`: scripts to run,
* `input`: reference signature catalogs,
* `cosmic tissue data`: signature activity (using the COSMICv3 catalog) in the tissue data provided by the COSMIC website,
* `real mutational catalogs`: SBS96 mutational catalogs of four samples from the PCAWG project and the average signature weights estimated for them by sigLASSO, SigProfilerAssignment, and MuSiCal.

To run the package, you need to: (1) download and unpack it (or clone it: `git clone https://github.com/8medom/SigFitTest`), (2) go to the directory `code`, (3) run `main.py`.

In `main.py`, there are seven main functions that you can use (comment or uncomment the corresponding lines in this script according to which functions you want to run):
* `generate_synthetic_catalogs(...)`: Generates mutational catalogs based on single base substitutions (SBS) where the mutations are classified in 96 different contexts. Signature weights in the synthetic data are obtained from the real tissue results in the directory `cosmic tissue data`. These synthetic catalogs can be later used for your evaluation of signature fitting (or signature extraction). The output data (synthetic catalogs and true signature weights) are stored in the directory `generated_data`.
* `fit_with_cosmic3_synthetic_simple(...)`: Generates simple mutational catalogs with pre-defined signature contributions (same for every sample), uses them as input for a signature fitting tool, and evaluates the results produced by this tool. All COSMICv3 signatures are used as a reference for fitting. The produced results are shown in Figure 1c in the manuscript referenced below.
* `fit_with_cosmic3_synthetic(...)`: Generates realistic mutational catalogs with signature activities driven by empirical results from various cancer types, uses them as input for a signature fitting tool and evaluates the results produced by this tool. All COSMICv3 signatures are used as a reference for fitting. Sample output of this function for two tools (SigsPack and SigProfilerAssignment) for two cancer types (Head-SCC and ColoRect-AdenoCA) for two cohorts and three different mutation counts (100, 2000, 50000) is included in this repository (files reference_results-WGS-SET3-MutationalPatterns.dat, reference_results-WGS-SET3-SPA.dat, and reference_results-WGS-SET3-SigsPack.dat). When parameter `out_of_reference_weights` is provided, empirical signature weights are scaled down appropriately and the remaining mutations are driven by out-of-reference signatures (12 signatures from COSMICv3.3.1 that are absent in COSMICv3 that are used for fitting). The produced results are shown in Figures 2 and 4d,e.
* `prune_reference_and_fit_synthetic()`: Uses results of `fit_with_cosmic3_synthetic()` run previously to prune the list of reference signatures (only signatures that are active in sufficiently many samples are kept), runs a fitting tool, and evaluates the results. Note that the zip files with the previously computed results (e.g., `signature_results-WGS-SPA-SET3-Head-SCC.zip`) must be located in the directory `code`. The produced results are shown in Supplementary Figure 21.
* `fit_with_cosmic3_synthetic_compare_groups()`: Generates realistic mutational catalogs with signature activities driven by empirical results from various cancer types where systematic differences in the activity of a chosen signature are introduced between odd and even samples, uses the catalogs as input for a signature fitting tool (all COSMICv3 signatures are used as a reference for fitting), and applies Wilcoxon rank-sum test to compare true as well as estimated signature weights between the two groups of samples. When parameter `out_of_reference_weights` is provided, empirical signature weights are scaled down appropriately and the remaining mutations are driven by out-of-reference signatures (12 signatures from COSMICv3.3.1 that are absent in COSMICv3 that are used for fitting), which additionally complicates the search for significant differences between the two groups of samples. The produced results are shown in Figure 3b.
* `fit_with_cosmic3_subsampled_real_catalogs()`: Generates mutational catalogs by subsampling from real mutational catalogs stored in the directory `real mutational catalogs` and runs a fitting tool on them. The estimated signature activities are evaluated by comparing them with the ground truth result stored in the same directory (the ground truth is obtained by analyzing the complete real catalogs using three different fitting tools and averaging the contributions of all signatures that are identified by at least two tools). The produced results are shown in Figure 4a,b.
* `fit_external()`: Fits a provided mutational catalog (path set by input_catalog) and evaluates the estimated signature activities by comparing them with the provided ground truth activities (path set by catalog_GT). Variable reference_signatures can contain a subset of COSMICv3.3.1 signatures that will be used as a reference instead of COSMICv3 (default).


### Practical points
* To use a fitting tool, the tool needs to be installed first following the instructions that come with that tool.
* To find the cancer types for which empirical signature activities are available, see the directory `cosmic tissue data`.
* Several important variables are set in the file `MS_config.py` in the directory `code`:
  * `tool` (tool that will be run and evaluated),
  * `N_samples` (how many synthetic samples there are in each cohort; default: 100),
  * `num_realizaions` (how many synthetic cohorts will be generated; default: 5),
  * `timeout_time` (fitting tools are stopped if they do not finish until this time (in seconds); you may need to increase it for slow tools such as mmsig, for example; default: 1800),
  * `num_muts_list_short` (how many mutations there will be in the synthetic catalogs; this is a list, catalog generation & fitting & evaluation are run separately for each number of mutations provided; default: [100, 2000, 50000]),
  * `input_signatures` (the reference signatures that are used to generate the catalogs; default: WGS COSMICv3.3.1; note that the empirical signature activities obtained from the COSMIC catalog contain only COSMICv3 signatures; the remaining signatures that are in COSMICv3.3.1 thus have zero activity unless the out_of_reference_weights variable is set for some of the functions described above).
* The variable `code_name` is used to easily recognize the results obtained by the provided functions under various settings.
* The content of `stdout` and `stderr` is saved in the files `stdout-tool_name.txt` and `stderr-tool_name.txt`. The content of these files can sometimes help you with troubleshooting.
* There are now scripts for 12 fitting tools: deconstructSigs, mmsig, MuSiCal, MutationalPatterns, sigfit, sigLASSO, sigminer (in three variants: NNLS, QP, and SA), signature.tools.lib (referred to as signature_tools), SigProfilerAssignment (referred to as SPA), SigProfilerSingleSample (referred to as SPSS), SigsPack, and YAPSA.


### Evaluation metrics
Output files `results-*.dat` present the summary results for each synthetic cohort. These files contain the following columns:
* `cancer_type`: The cancer type whose empirical signature weights in the COSMIC catalog were used to generate signature weights in synthetic samples,
* `weights`: Realization counter,
* `samples`: Number of samples in each cohort,
* `muts`: Number of mutations in each sample,
* `TAE`: Total Absolute Error (referred to as the fitting error in the article) averaged over all samples; this quantifies the difference between the true and estimated signature weights,
* `TAE_std`: Standard deviation of the total absolute error (all other *_std columns are standard deviations of respective metrics),
* `RMSE`: Root mean square error averaged over all samples,
* `wT`: Total relative weight assigned to any signature, averaged over all samples,
* `n_eff`: The effective number of signatures (inverse Herfindahl index) to which weights are assigned, averaged over all samples,
* `TAE_TP`: Total Absolute Error computed only for true positives (i.e., the estimated active signatures that are active in a given sample),
* `wT_FP`: Total relative weight assigned to false positives (i.e., the estimated active signatures that are *not* active in a given sample),
* `n_FP`: Number of false positives per sample,
* `wT_FN`: Total relative weight of false negatives (i.e., the estimated inactive signatures that are active in a given sample),
* `n_FN`: Number of false negatives per sample,
* `P`: Precision (i.e., the fraction of estimated active signatures that are true positives),
* `R`: Recall/Sensitivity (i.e., the fraction of signatures active in a sample that are estimated to be active),
* `S`: Specificity (i.e., the fraction of signatures inactive in a sample that are estimated to be inactive),
* `F1`: F1 score, the harmonic mean of precision and recall,
* `MCC`: Matthews correlation coefficient, a correlation between the true and estimated signature activities,
* `Pearson`: Pearson correlation between the true and estimated signature weights, including only the signatures for which either true or estimated weight is positive, averaged over all samples,
* `r_S*`: Pearson correlation between the true and estimated signature weights for all samples and one chosen signature (the list of five chosen signatures for each cancer type is in `MS_config.py` in the dictionary `top_sigs`).


## Citation
M. Medo, Charlotte Ng, M. Medová, A comprehensive comparison of tools for fitting mutational signatures, Nature Communications (in press).


## Contact Information
Please address any queries to Matúš Medo at matus.medo@unibe.ch.
