# Analytical code for studying fossil hyracoid teeth
 
## For Locus Example Publication

Scripts associated with manuscript: Vitek, N., and P. Princehouse. 2024. Evaluating the utility of linear measurements to identify isolated tooth loci of extinct Hyracoidea. Acta Palaeontologica Polonica 69.

### Organization

Start at `hyracoid_base.R`. This is the entry point that loads librariers, sets file paths, etc.

The next steps (which is automatically sourced in the base script) are in `hyracoid_merge_data.R`. The first parts of this script source the formatting scripts and get the data in a state to be analyzed. The data merge script calls from among the following five scripts:  `literature_format_linear_data.R`, `procavia_format_linear_data_lowers.R`, `fayum_format_linear_data_lower.R`, `procavia_format_linear_data_uppers.R`, and `fayum_format_linear_data_uppers.R`.

The next part (`hyracoid_analyses.R`) produces descriptive statistics, then conduct analyses, including calling scripts that are specific to the upper or lower tooth arcade (because slightly different sets of variables are evaluated in each). The analytical code calls either of two sets: 

* `hyracoid_analyses_lowers.R` & `hyracoid_analyses_LDA_lowers.R`
* `hyracoid_analyses_uppers.R` & `hyracoid_analyses_LDA_uppers.R` 

as well as `hyracoid_ratio_physig.R` which itself calls `retention_borges_delta.R`

## For Species Description and Phylogenetic Analysis Publication

Scripts associated with manuscript: 
Each of the following R scripts can be run on its own. Each calls the data file `Topernawi_hyracoid_published_tooth_sizes.csv`, reposited in an associated Dryad data repository [https://doi.org/10.5061/dryad.2jm63xsw2](https://doi.org/10.5061/dryad.2jm63xsw2).

* `hyracoid_mass_estimate.R`: This script estimates body masses of hyracoid species from Topernawi based on the lengths of their second molars using previously published equations.
* `topernawi_dimensions.R`: This script formats and plots length and width information, and calculates some ratios used for species differentiation.
* `topernawi_association_context.R`: This scripts formats and plots additional metrics that can be used to evaluate associations between isolated specimens, such as coefficients of variation for collected specimens, proportions of upper vs lower molars, and proportions of teeth from different loci in the same arcade (ex: m1 vs m3 size). It calls additional data files associated with a separate publication (Vitek and Princehouse, 2024), reposited here: [https://doi.org/10.5061/dryad.n2z34tn40](ttps://doi.org/10.5061/dryad.n2z34tn40)
 

### Organization

Each of the R scripts can be run on their own. The following call the data file `Topernawi_hyracoid_published_tooth_sizes.csv`
* `hyracoid_mass_estimate.R`
* `topernawi_dimensions.R`
* `topernawi_association_context.R`

The remaining R script can be run on its own, calling phylogenetic analysis output files also reposited in the same Dryad data repository. 
* `plot_tree-JVP.R`
