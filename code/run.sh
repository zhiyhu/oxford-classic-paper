#!/usr/bin/env bash
set -ex

# This is the master script for the capsule. When you click "Reproducible Run", the code in this file will execute.


Rscript fig1_HHcohort_analysis.R
# Figure 1C
# Figure 1D
# Table 2
# Table S1


Rscript fig2d_repeat_association.R
# Figure 2D


Rscript fig3a_HHcorhot_ciliated.R
# Figure 3A

Rscript fig2_tcga_immune.R
# Figure 2A, B, C
# Figure S2

Rscript fig3_caps_analysis.R
# Figure 3C, D