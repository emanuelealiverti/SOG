# Removing the influence of a group variable in high-dimensional predictive modelling

Content:

```
├── README.md
├── SIMULATION_plot_rep.R
├── SIMULATION_rank.R
└── SIMULATION_y_rep.R
```
Files of this folder reproduce the simulation studies from Section 3 of the paper. 
Each result is averaged over 50 random splits intro training and test test.

- `SIMULATION_plot_rep.R` illustrate the empirical correlation between out of sample predictions and the group variable $Z$
- `SIMULATION_y_rep.R` illustrate the performance in predicting the response variable $Y$ in out-of-sample predictions
- `SIMULATION_rank.R` illustrate the performance in reconstructing the matrix $X$ for increasing values of the approximation rank 

