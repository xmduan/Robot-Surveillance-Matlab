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
## Comparison with Julia
Julia code can be seen here [Robot-Surveillance-Julia](https://github.com/SJTUHan/Robot-Surveillance-julia), test code can be found in the packages. We test star graph, line graph, ring graph, complete graph and random graph with n=8 points, respectively. It should be noted that in the test file, every graph is randomly set(except for complete graph). Thus, in order to compare the efficiency, particular graphs are used instead of random ones. Graph structures are shown as follow.

Code efficiency are shown in the below tablet.

| Graph Type | Calculated Parameter |Computing Time (Matlab) | Optimization Time(Matlab)|  Computing Time (Julia)|Optimization Time(Julia) |
|:-:|:-:|:-:|:-:| :-:|:-:|
|Star Graph|Hitting Time|4.8650e-04s|0.0566s|0.0033s|0.215s|
|Ring Graph|Hitting Time|3.9410e-04s|0.0184s|0.0039s|0.203s|
|Complete Graph|Hitting Time|6.5570e-04s|0.0251s|0.0041s|17.26s|
|Random Graph|Hitting Time|6.6380e-04s|0.0247s|0.0039s|0.270s|
|Line Graph|Hitting Time|7.3720e-04s|0.0190s|0.0038s|0.275s|
|Star Graph|Entropy Rate|2.8690e-04s|0.5397s|9.99e-04s|23.08s|
|Ring Graph|Entropy Rate|3.5170e-04s|1.1642s|9.99e-04s|0.939s|
|Complete Graph|Entropy Rate|2.1030e-04s|0.0574s|1.00e-03s|0.239s|
|Random Graph|Entropy Rate|3.4160e-04s|0.1854s|9.99e-04s|1.934s|
|Line Graph|Entropy Rate|2.6370e-04s|0.1758s|9.99e-04s|1.823s|
|Star Graph|Kemeny|7.8910e-04s| 0.9552s|0.0016s|40.52s|
|Ring Graph|Kemeny|4.1860e-04s|0.7187s|0.0013s|39.24s|
|Complete Graph|Kemeny|4.3780e-04s|0.3816s|0.0017s|27.39s|
|Random Graph|Kemeny|3.8760e-04s|0.5134s|0.0016s|27.23s|
|Line Graph|Kemeny|3.8330e-04s|0.5225s|0.0017s|26.99s|
|Star Graph|Return Time Entropy|0.0513s|5.1133s|50%&50%|1 ms|
|Ring Graph|Return Time Entropy|0.1183s|30.2850s|50%&50%|1 ms|
|Complete Graph|Return Time Entropy|0.0097s|16.0799s|50%&50%|1 ms|
|Random Graph|Return Time Entropy|0.0147s|20.0534s|50%&50%|1 ms|
|Line Graph|Return Time Entropy|0.0188s|21.4639s|50%&50%|1 ms|
|Star Graph|Mixing Time|3.8840e-04s|3.4760s|0.0018s|0.0899 ms|
|Ring Graph|Mixing Time|2.1500e-05s|0.7130s|0.0017s|0.1909 ms|
|Complete Graph|Mixing Time|2.1500e-05s|0.5978s|0.0016s|0.0749 ms|
|Random Graph|Mixing Time|2.2000e-05s|0.5635s|0.0019s|0.2400 ms|
|Line Graph|Mixing Time|2.2500e-05s|0.5645s|0.0016s|0.1570 ms|
