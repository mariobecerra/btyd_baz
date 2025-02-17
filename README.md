# btyd_baz

Fork from the [BTYDplus](https://github.com/mplatzer/BTYDplus/) R package. The BTYDplus package provides advanced statistical methods to describe and predict customer's purchase behavior. It uses historic transaction records to fit a probabilistic model, which then allows to compute quantities of managerial interest on a cohort- as well as on a customer level (Customer Lifetime Value, Customer Equity, P(alive), etc.).

The BTYDplus package complements the [BTYD](https://cran.r-project.org/package=BTYD) package by providing several additional buy-till-you-die models, that have been published in the marketing literature, but whose implementation are complex and non-trivial. These models are: NBD, MBG/NBD, BG/CNBD-k, MBG/CNBD-k, Pareto/NBD (HB), Pareto/NBD (Abe) and Pareto/GGG.

This package was forked to improve the memory usage and parallel processing of certain functions, especially when working with millions of customers.

## Installation

```
# install.packages("devtools")
devtools::install_github("mariobecerra/btyd_baz", dependencies=TRUE)
library(btyd_baz)
```

## Getting Started

```
demo("cdnow")        # Demonstration of fitting various models to the CDNow dataset
demo("mbg-cnbd-k")   # Demonstration of MBG/CNBD-k model with grocery dataset
demo("pareto-abe")   # Demonstration of Abe's Pareto/NBD variant with CDNow dataset
demo("pareto-ggg")   # Demonstration of Pareto/NBD (HB) & Pareto/GGG model with grocery dataset
```

## Implemented Models

These R source files extend the functionality of the BTYD package by providing functions for parameter estimation and scoring for NBD, MBG/NBD, BG/CNBD-k, MBG/CNBD-k, Pareto/NBD (HB), Pareto/NBD (Abe) and Pareto/GGG.

* **NBD** Ehrenberg, Andrew SC. "The pattern of consumer purchases." Applied Statistics (1959): 26-41.
* **MBG/NBD** Batislam, E.P., M. Denizel, A. Filiztekin. 2007. Empirical validation and comparison of models for customer base analysis. International Journal of Research in Marketing 24(3) 201–209.
* **(M)BG/CNBD-k** Platzer, Michael, and Thomas Reutterer (submitted)
* **Pareto/NBD (HB)** Ma, Shao-Hui, and Jin-Lan Liu. "The MCMC approach for solving the Pareto/NBD model and possible extensions." Natural Computation, 2007. ICNC 2007. Third International Conference on. Vol. 2. IEEE, 2007.
* **Pareto/NBD (Abe)** Abe, Makoto. "Counting your customers one by one: A hierarchical Bayes extension to the Pareto/NBD model." Marketing Science 28.3 (2009): 541-553.
* **Pareto/GGG** Platzer, Michael, and Thomas Reutterer. "Ticking Away the Moments: Timing Regularity Helps to Better Predict Customer Activity." Marketing Science (2016).

