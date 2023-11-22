# hyracoid-locus-example
 Scripts associated with manuscript
 
 ## Organization

Start at `hyracoid_base.R`. This is the entry point that loads librariers, sets file paths, etc.

The next steps (which is automatically sourced in the base script) are in `hyracoid_merge_data.R`. The first parts of this script source the formatting scripts and get the data in a state to be analyzed. The data merge script calls from among the following five scripts:  `literature_format_linear_data.R`, `procavia_format_linear_data_lowers.R`, `fayum_format_linear_data_lower.R`, `procavia_format_linear_data_uppers.R`, and `fayum_format_linear_data_uppers.R`.

The next part (`hyracoid_analyses.R`) produces descriptive statistics, then conduct analyses, including calling scripts that are specific to the upper or lower tooth arcade (because slightly different sets of variables are evaluated in each). The analytical code calls either of two sets: 

	- `hyracoid_analyses_lowers.R` & `hyracoid_analyses_LDA_lowers.R` 
	- `hyracoid_analyses_uppers.R` & `hyracoid_analyses_LDA_uppers.R` 

