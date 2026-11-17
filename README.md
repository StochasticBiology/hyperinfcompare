# hyperinfcompare
Compare inferences from hypercubic EvAM models

<img width="747" height="261" alt="image" src="https://github.com/user-attachments/assets/515f1e8c-fd1e-4b59-90e4-d8895015adca" />

This R library provides functionality for a statistical comparison between the outputs of hypercubic evolutionary accumulation models (EvAM), based on the size and robustness of inferred differences in evolutionary dynamics.

The associated preprint is here https://arxiv.org/abs/2608.02781

Installation
----
Install with `remotes::install_github("StochasticBiology/hyperinfcompare", build_vignettes = TRUE)`. You'll also need the `hyperinf` package from here https://github.com/StochasticBiology/hyperinf , which wraps a lot of EvAM functionality. `hyperinf` also contains the functionality from this repo; it is included here for independent review.

Contents
----
The core function here is `compare_orderings`, which takes two fitted EvAM models and tries to identify differences in relative and/or absolute orderings of acquisition dynamics. There is also functionality to visualise these differences (`plot_hyperinf_compare_orderings`) and other associated statistics.

The code and data for the case studies in the manuscript, with sources, is in `inst/`. 

A vignette (`browseVignettes("hyperinfcompare")`) illustrates the shape of an analysis with this method.
