# Removing the influence of a group variable in high-dimensional predictive modelling


Predictive modelling relies on the assumption that observations used for training are representative of the data that will be encountered in future samples. In a variety of applications, this assumption is severely violated, since observational training data are often collected under sampling processes which are systematically biased with respect to group membership. Without explicit adjustment, machine learning algorithms can produce predictions that have poor generalization error with performance that varies widely by group. We propose a method to pre-process the training data, producing an adjusted dataset that is independent of the group variable with minimum information loss. We develop a conceptually simple approach for creating such a set of features in high dimensional settings based on a constrained form of principal components analysis. The resulting dataset can then be used in any predictive algorithm with the guarantee that predictions will be independent of the group variable. We develop a scalable algorithm for implementing the method, along with theory support in the form of independence guarantees and optimality. The method is illustrated on some simulation examples and applied to two real examples: removing machine-specific correlations from brain scan data, and removing race and ethnicity information from a dataset used to predict recidivism.

This repository is associated with the article: [Aliverti, Lum, Johndrow and Dunson (2018). *Removing the influence of a group variable in high-dimensional predictive modelling*](https://arxiv.org/abs/1810.08255).

The R package can be installed as follows.

```R
devtools::install_github("emanuelealiverti/SOG")
```
