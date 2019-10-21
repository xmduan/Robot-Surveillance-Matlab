Robot-Surveillance-Matlab
======
This is a free open source matlab toolbox for helping with calculating and optimizing some stochastic parameters on markov chain. It is basicly used in robot application when one or a group of robots randomly move on a graph to perform a surveillance task. Details of the algorithms and mathematic forms can be found in [professor Francesco Bullo's publications](http://motion.me.ucsb.edu/papers/index.html).
# Installation
* To install the toolbox into your matlab, please download [Robot Surveillance.mltbx](https://github.com/SJTUHan/Robot-Surveillance-Matlab/blob/master/Robot%20Surveillance.mltbx) into your workspace of matlab, then double click the file. Now all the function in the toolbox can be used.
* The other way to use the toolbox is download and decompress the software package, then add all files into the workspace. 
# Usage
## Notice
Functions in the package can be separated into two types: evaluation and optimization. For convex optimization problems, [CVX](http://cvxr.com/cvx/) is used, which means that CVX is necessary when using this kind of functions. For non-convex problems, [fmincon](https://www.mathworks.com/help/optim/ug/fmincon.html) is used. It should be noted that users can change the option of fmincon in the functions by themselves, no additional interfaces for changing them are offered in this version.
## Example
Detailed user instructions can be found in documentations of the functions. See below for an example. 

```
A=[1 1 0;
   1 0 1;
   0 1 1];
W=[1 2 3;
   4 5 6;
   7 8 9];
tau=10;
[F,K]=MC_OP(P,W,tau,'HittingTimeOp');
```
