# COMSOLGeneticSolver
Genetic solver that uses COMSOL simulations to determine the fitness of simulation parameters.

My dissertation work had to do with acoustic localization of the source of a thermal shock inside a metal cylinder due to an electric arc. This code attempts to determine the exact location of the microphones and the material properties and component coupling spring constants that best reproduce the acoustic signal data obtained during experiments using a genetic solver. Fitness is determined by running a COMSOL simulation, generating plot data, and comparing it with experimental data.

The base genetic algorithm was tested against at least the Sphere function and the Rosenbrock function. I may have done others, but I can't seem to find the test code that I used in graduate school. At any rate, the algorithm handled the Sphere function just fine. It was inconsistent in finding the best solution for the Rosenbrock function. Considering Rosenbrock is notoriously difficult to optimize for and I coded the genetic algorithm from scratch, I considered this code a noble effort.

In the end the solver didn't yield a positive result, but it did convince me that my simulation model needed more tweaking. It was also a fun, worthwhile project that taught me a potentially useful skill.
