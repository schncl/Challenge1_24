# Challenge1_24

This program read all the parameters for performing optimization from a getpot file, so we need to pass the file via -f filename
In the GetPot file you can also choose which optimization strategy you want between gradient descent (with or without momentum)
and nesterov. Furthermore you can also choose the update strategy for the learning rate. The options in this case are Armijo
Inverse or Exponential Decay or constant. 
The functions are hardcoded as lambda functions because i was not able to figure out muparserx in time
The choice of the learning rate and initial point are rather important, usually the bigger the modulus of x and y 
the more likely is the algorithm to not converge, and also the learning rate has to be tuned via trial and error

I have created various test and i will add also a way to plot the results either a python notebook or a script
