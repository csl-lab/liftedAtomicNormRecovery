# Lifted Atomic Norm Recovery

Atomic norm minimization is a popular approach to perform signal recovery from incomplete linear measurements when the underlying signal is a sparse combination of elementary building blocks called *atoms*. However, there is an *implicit* geometric constraint in this strategy, namely, not all atoms are used to represent the recovered signal. In fact, only the *exposed atoms* will be selected to represent the optimal solution, while the *masked atoms* will never be selected. These masked atoms may as well never have been included in the signal model. 

This repository contains the implementation of an alternative approach to atomic norm minimization, called *lifted atomic norm recovery*, that enables the selection of masked atoms to represent the recovered signal.

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
 