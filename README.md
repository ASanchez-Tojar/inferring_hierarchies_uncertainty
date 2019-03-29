# Inferring reliable dominance hierarchies: randomized Elo-rating method

This repository contains the R scripts use in the following study:

[Alfredo Sánchez-Tójar, Julia Schroeder, Damien R. Farine. 2018. *A practical guide for inferring reliable dominance hierarchies and estimating their uncertainty*. Journal of Animal Ecology, 87(3):594-608. DOI: 10.1111/1365-2656.12776](https://doi.org/10.1111/1365-2656.12776).

More information and materials (e.g. data) available at the [OSF project](http://doi.org/10.17605/OSF.IO/9GYEK). For any further information, please contact: [Alfredo Sánchez-Tójar](https://scholar.google.co.uk/citations?hl=en&user=Sh-Rjq8AAAAJ&view_op=list_works&sortby=pubdate), email: alfredo.tojar@gmail.com

## Scripts:

`Figure1_Example_of_an_interaction_dataset_for_diagram.R`: script to simulate some dummy dominance interactions to be used in **Figure 1** of [our publication](https://doi.org/10.1111/1365-2656.12776).

`Figure2_Exploring_parameter_space.R`: script to explore some of the parameter space of our approach (i.e. hierarchy scenarios simulated) - **Figure 2** of [our publication](https://doi.org/10.1111/1365-2656.12776).

`Figure3_Elo-rating_and_steepness.R`: script to show how the latent hierarhcy affects the performance of the original (non-randomized) Elo-rating - **Figure 4** of [our publication](https://doi.org/10.1111/1365-2656.12776).

`Figure4_Method_comparison_and_sampling_effort.R`: script to compare the performance of different methods for inferring dominance hierarchies and explorating the influence of sampling effort on the performance - **Figure 5** of [our publication](https://doi.org/10.1111/1365-2656.12776).

`Figure5_Randomized_Elo-rating_repeatability.R`: script to show how Elo-rating repeatability changes depending on the steepness of the hierarchy - **Figure 6** of [our publication](https://doi.org/10.1111/1365-2656.12776).

`Figure6_Halve_comparison.R`: script to show how the agreement between thetwo halves of the data changes depending on the steepness of the hierarchy - **Figure 7** of [our publication](https://doi.org/10.1111/1365-2656.12776).

`Supplementary_Information-Very_steep_hierarchies.R`: script to generate the supplementary figures of [our publication](https://doi.org/10.1111/1365-2656.12776).

`worked_examples_markdown.Rmd`: Rmarkdown script providing an easy step-by-step guide on how to infer dominance hierarchies, and their uncertainty, steepness and transitivity from dyadic interactions.


### Notes:

Notice that some of the figure numbers referred to in the name of the scripts do not match the figures numbers of the final publication.

11th March 2018: a new version of the R package '[aniDom](https://cran.r-project.org/web/packages/aniDom/index.html)' (v.0.1.3) has been uploaded to the CRAN.
