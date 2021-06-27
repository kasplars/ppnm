This set of files has been created by Kasper Larsen, student nummber: 201806672.

PROBLEM STATEMENT:

- 6 - Adaptive integration of complex-valued functins

Implement an adaptive integrator which calculates the integral of a complex-valued function f(z) of a complex variable z 
along a straight line between two points in the complex plane. 

-----

This folder contains the following files (besides this file):
	auxfile.c		This file contains the functions which together make the adaptive integrator for complex
				functions. The Clenshaw-Curtis variable transformation is also implemented. Furthermore,
				the integrator estimates and returns the integration error. The problem statement is also given in here.
	auxfile.h		This header file holds the function declarations.
	main.c			This is the main file which produces the output in out.txt. Here, the adaptive integrator is tested
				for a specific function given a set of parameters.
	out.txt			This file contains the output.
	Makefile		The makefile produces the out.txt and error.txt files upon writing 'make' in the terminal. Moreover,
				writing 'make clean' removes the object and executable code and the errors.txt file.
	EXAMPROJECTS.md		This file styles the folder in github. Do not mind it.
	
Elaboration on some of the code:
	As mentioned above, I made the integrator estimate and return the integration error. This was estimated using eq. (46) in the notes,
	using the complex norm, |z|=sqrt(x*x+y*y) with z = x+iy. I could have taken the real and imaginary parts separately and found the
	error on each of these, but for simplicity, I did not.
	
	In main.c, I tested the integrator on a function which has a branch point. The csqrt function chooses the negative real axis as a branch
	cut, and my results are consistent with this. 
	
	I did not implement infinte limits as the problem was to integrate along a line between two points in the complex plane. As far as I recall,
	the coordinates of the points must be well-defined. If one coordinate is infinite, then it is not well-defined. 
	
Comment on the result: 
	My calculated integral values are consistent with the true integral value within the given precision.
	
Self-evalutaion:
	I believe that I have completed the task and done a little more than what was asked in the problem statement, i.e., making the integrator
	return an estimate of the error and having included the Clenshaw-Curtis variable transformation.
	
	Therefore, I would give myself 9 or 10 out of 10.

	   

