# Removing the influence of a group variable in high-dimensional predictive modelling
This repository is associated with the article: [Aliverti, Lum, Johndrow and Dunson (2018). *Removing the influence of a group variable in high-dimensional predictive modelling*](https://arxiv.org/abs/1810.08255).

### Abstract
In many application areas, predictive models are used to support or make important decisions. There is increasing awareness that these models may contain spurious or otherwise undesirable correlations. Such correlations may arise from a variety of sources, including batch effects,  systematic measurement errors, or sampling bias.
Without explicit adjustment, machine learning algorithms trained using these data can produce poor out-of-sample predictions which propagate these undesirable correlations.
We propose a method to pre-process the training data, producing an adjusted dataset that is statistically independent of the nuisance variables with minimum information loss.
We develop a conceptually simple approach for creating an adjusted dataset in high-dimensional settings based on a constrained form of matrix decomposition. 
The resulting dataset can then be used in any predictive algorithm with the guarantee that predictions will be statistically independent of the group variable.
We develop a scalable algorithm for implementing the method, along with theory support in the form of independence guarantees and optimality.
The method is illustrated on some simulation examples and applied to two case studies: removing machine-specific correlations from brain scan data, and removing race and ethnicity information from a dataset used to predict recidivism. That the motivation for removing undesirable correlations is quite different in the two applications illustrates the broad applicability of our approach.



The R package can be installed as follows.

```R
devtools::install_github("emanuelealiverti/SOG")
```

Simulation studies included in Section 3 of the paper can be reproduced with the scripts contained in the [folder](./SIMULATIONS) `SIMULATIONS`
