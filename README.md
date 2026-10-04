# ReVMD

This repository contains the implementation code for Reassembly variational mode decomposition.

## Usage

`x` is the input signal.

```matlab
Alpha =1000000;
Tau = 0;
K=5
Tol = 1e-7;
StopCriterion = 4;
[ReVMDModes, ~, ~] = ReVMD(x, Alpha,Tau, K, Tol,StopCriterion);
```
