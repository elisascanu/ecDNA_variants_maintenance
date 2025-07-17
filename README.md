# ecDNA_variants_maintenance
Code for the manuscript "The maintenance and de-mixing of extrachromosomal1 DNA variants in single cells".

## Main function for simulating ecDNA proliferation in a growing population of cells, considering selection and switching between variants

The main function is called "generate" and it is included in most of the MATLAB scripts. The explanation step by step is present in `mixmaintainance.m`.

## The script `distribution_variants.m`

This script generates Figure 2 and it tracks the distribution of copies among the subpopulations.

## The script `switching_times.m`

This script generates Figure 3 and it tracks the fraction of different subpopulations when the switching among variants happens at different stages of the population. This script includes the case of neutral selection.

## The script `mixmaintainance.m`

This script generates Figure 4 and it tracks the fraction of different subpopulations under different scenarios of identical selection and switching.

## The script `probabilitytable.m`

This script generates a table of probabilities for the state of both mother cells and daughter cells. Precisely, it tracks the probability that a mother cell picked for division is either pure, mixed or ecDNA-free, and the probability of the relative outcome of this division, meaning the probability that a daughter is either mixed, pure or ecDNA-free. This table can be easily used to produce Figure 5 and the Supplementary Figures S5-S8.

## The script `shannon_singlecell.m`

This script generates Figure 6, thus a scatterplot for the Shannon index at single-cell level for different values of selection and switching probability.
