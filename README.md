<link href="https://gist.githubusercontent.com/tuzz/3331384/raw/94f2380c2b798fab2139fd0a8f478c4f2d642e3b/github.css" rel="stylesheet"></link>

## Overview
This is a project that tries to solve a variant of Jermann Quadrini (2012 AER) without debt and adjustment. I try to solve it in two ways: Adrian's method and projection method. My guess is that projection method should be much faster than Adrian's method. It should solve the model within 5 seconds if given a good initial guess.

## Folder and Files
+ /cppcode.cpp: source code that **has** to be compiled by GCC4.8+
+ /cudacode.cu: source code that **has** to be compiled by NVCC 5.5+
+ /Model/: contains LyX and PDF that describe the model.
+ /MATLAB/: contains some codes written in MATLAB
+ /Dynare/: contains some codes written in Dynare

## Goal
+ Both methods yields the same solution. Can also use Value Function Iteration to check?
+ Implement linearization solution as a initial guess.
+ Implement Newton's method as generic as possible. Maybe use a `class` to contain function and derivative, if not too slow.
+ Make everything more portable/systematic:
	+ Use a class to contain model parameters and solution specific parameters.
	+ Put helper functions into appropriate header.
	+ Figure out a way to separate model "things" from solution "things". Model things like steady states, eureka should be put in a different header file. But the difficulty lies in the unpredictable (even at compile time) number/dimension of objects. For example, different models have different number of endogenous state variables and shocks, then how should we deal with grid point creation and initial guess? Maybe a multidimensional matrix/array can deal with this, but this would be henious to read and Juan wouldn't like it.

## To-do List
+ Use parameter structure to replace `#define` parameters.
+ A generic way to compute steady states
+ Writes the eureka function, assuming n = 0.3 always
+ Shamelessly clone the safe Newton's method in *Numerical Recipes*. Hopefully the function object thing isn't too slow. If it is then I function pointer it is. I just need to wrap another layer of function around the one that's depend on endogenuous variables.
+ Write objective and derivative functions in both cases, and finishes eureka
