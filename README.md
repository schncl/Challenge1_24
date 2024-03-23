# Challenge1_24

This program read all the parameters for performing optimization from a getpot file, so we need to pass the file via -f filename
In the GetPot file you can also choose which optimization strategy you want between gradient descent (with or without momentum)
and nesterov. Furthermore you can also choose the update strategy for the learning rate. The options in this case are Armijo
Inverse or Exponential Decay or constant. 
The functions are hardcoded as lambda functions because i was not able to figure out muparserx in time
The choice of the learning rate and initial point are rather important, usually the bigger the modulus of x and y 
the more likely is the algorithm to not converge, and also the learning rate has to be tuned via trial and error.
You can also decide whether to use the hard-coded version of the gradient or the one computed via finite differences.

The program prints the minimum point, as well as the value of the function in said point, and write in a file
for each iteration the coordinate of the point at iteration k and the computed value of the function.

I created different test, and i also added a simple python notebook to visualize the convergence history over a contour plot of f, 
since i was not able to do that with gnuplot...
