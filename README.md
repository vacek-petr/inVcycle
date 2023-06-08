# inVcycle
Code for reproducing results in P. Vacek, E. C. Carson and K. M. Soodhalter, The effect of approximate coarsest-level solves on the convergence of multigrid V-cycle methods.
MATLAB R2023a is used.

* **experiments/variousRelativeResidualTolerances.m**: generates data included in Figure 1 of the paper.
* **experiments/BLRVariousLowRankThresholdParameters.m**: generates data included in Figure 2 of the paper.
* **experiments/stoppingCriteriaComparison.m**: generates data included in Figure 3 of the paper.
* **experiments/cgStoppingStrategy.m**: generates data included in Figures 4 and 5 of the paper.
* **experiments/BLRTestEstimateAccuracy.m**: generates data included in Figure 6 of the paper.


The BLR LU decompositions of system matrices, required in **experiments/BLRVariousLowRankThresholdParameters.m** and  **experiments/BLRTestEstimateAccuracy.m** are not included in the repository. They can by generated by the script **computeBLRLU.m**.

The script **computeBLRLU.m** and funciton **functions/blr.m** contain modified parts of code taken from the [BLRstability repository](https://gitlab.com/theo.andreas.mary/BLRstability).

## Matrices
The system matrices were generated by the FE software [FEniCS](image.png) (version 2019.1.0), see **matrixGeneration/assembleMatricesFenics.py**
