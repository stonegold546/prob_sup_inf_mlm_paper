# Prob-Sup for clustered data

The **ps_functions_new.R** file is always required as it contains the functions for different estimation/inference techniques.

#### Examples

If wanting to try out examples, download ToolsK_archive.csv, share.dta and examples.R

The examples.R file include ordinal example 4.2 in Zou (2021), as a check that Zou's approach was correctly implemented - this is unreported in our manuscript. In Zou's manuscript, Zou used z critical values for this example. However, note that the function we built for Zou's approach uses t critical values rather than z. The function includes directions on how to use z critical values in place of t.

#### Simulation studies

If wanting to run the simulation studies, download:

* sim_new_0_null.R for study 1
* sim_new_1_beta.R for study 2
* sim_new_2_bin.R for study 3
