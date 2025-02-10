# Lifted Atomic Norm Recovery

Signal recovery methods from incomplete linear measurements often make assumptions about the structure of the underlying signal. For instance, it is often assumed that the underlying signal is a sparse combination of elementary signals called *atoms*. If we let $\mathcal{A} := \{\boldsymbol{a}_1, \ldots, \boldsymbol{a}_n\}$ be an *atomic set* in $\mathbb{R}^d$ then a measure of complexity of a signal $\boldsymbol{x}_0 \in \mathbb{R}^d$ is the value of its *atomic norm*
$$
    \rho_{\mathcal{A}}(\boldsymbol{x}_0) = \inf\left\{\sum_{i=1}^{n} c_i:\,\, \boldsymbol{x}_0 = \sum_{i=1}^{n} c_i\boldsymbol{a}_i,\,\,c_i\geq 0\right\}.
$$
If the incomplete linear measurements are given by $\boldsymbol{y}_0 = \boldsymbol{\Phi}\boldsymbol{x}_0$ for some $m\times d$ matrix $\boldsymbol{\Phi}$, where typically $m \ll d$, then *atomic norm recovery* consists on solving
$$
    \text{minimize}_{\boldsymbol{x}\in\mathbb{R}^d} \quad \rho_{\mathcal{A}}{\boldsymbol{x}} \quad \text{subject to} \quad \boldsymbol{\Phi}\boldsymbol{x} = \boldsymbol{y}_0.
$$

However, there is an *implicit* geometric constraint, namely, not all atoms in $\mathcal{A}$ are used to represent the optimal solution $\boldsymbol{x}^{\star}$ to the above problem. In fact, only the *exposed atoms* will be selected to represent the optimal solution, while the *masked atoms* will never be selected. In fact, they may as well never have been included in the signal model. 

This repository contains the implementation of an alternative approach, called *lifted atomic norm recovery*, that enables the selection of masked atoms to represent the recovered signal.

Version 1.0

Date: 10.02.2025

## References

The following manuscripts describes the main motivation behind the method.

> P. Rademacher and C. A. Sing Long, "*Lift and unmask: Sparse models by lifting the atomic norm*." February 10, 2025.

The code reproducing the results of their numerical experiments can be found on ``experiments``.

### Citation

You can use the following code to cite our work.

```bibtex
@misc{rademacher_lift_2025,
	title     = {Lift and unmask: Sparse models by lifting the atomic norm},
	author    = {Rademacher, Pablo and Sing Long, Carlos A.},
	date      = {2025-02-10}}
```

## Repository

### Dependencies

The dependencies for the Python implementation are:
* ``numpy``
* ``scipy``
* ``cvxpy``
* ``matplotlib``

### Structure

``liftedAtomicNormRecovery\``
* ``solver\``
    * ``cvxpy.py``: implementation of lifted atomic norm recovery using ``cvxpy``

``experiments\``
  * ``conference\``: experiments in "*Lift and unmask: Sparse models by lifting the atomic norm*"
    * ``FIGS\``: folder in which the figures will be saved
    * ``DAT\``: folder in which the precomputed data is saved
    * ``01_Phase_Diagrams.ipynb``: computation of phase diagrams for several recovery strategies
 