#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate ./r_env
Rscript generate_review_plots.R
