# BB Bottleneck Estimator
## Description
This is the implementation for [transmission bottleneck estimation based on beta-binomial sampling](https://www.biorxiv.org/content/10.1101/101790v1).  

This code is adapted from:
https://github.com/weissmanlab/BB_bottleneck/blob/master/Bottleneck_size_estimation_exact.r
For details of how the code works and examples, see the ReadMe of this repo.
I vectorised and parallelised the functions, so it runs about 20x faster on 8 cores than the original code.

The input data for this project are a metadata xlsx file with data on each sequence, and
NGS diversity data output by VSensus (https://github.com/rjorton/VSensus).


## Requirements
- R 3.6.2+
- rmutil
- parallel
- openxlsx
