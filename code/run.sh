#!/usr/bin/env bash
set -ex

# This is the master script for the capsule. When you click "Reproducible Run", the code in this file will execute.

Rscript fig1_HHcohort_analysis.R
# Figure 1C
# Figure 1D
# Table 2
# Table S1

Rscript fig2_tcga_immune.R
# Figure 2A, B, C
# Figure S4

Rscript fig3a_HHcorhot_ciliated.R
# Figure 3A

Rscript fig3_caps_analysis.R
# Figure 3C, D, S6

Rscript figs5_repeat_association.R
# Figure S5

Rscript figs2_technical.R
# Figure S2

Rscript figs3_brca.R
# Figure S3; Tables S3,4
