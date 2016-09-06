# Presolve
Presolve Routines for Optimization Problems

Presolver(s::AbtractMathProgBaseSolver) returns a conic solver.


This is an implementation of LP Presolving based on ideas from Andersen & Andersen paper - http://www.turma-aguia.com/davi/doc/Andersen.pdf
Presolving removes redundancies from the original problem given by the user and constructs a smaller equivalent problem
which is then fed to the solver.
The input to the presolver! function is the data fed into the linprog() function of MathProgBase API.
It returns the presolved problem along with a "Presolve" stack which contains information on each redundancy that was removed.
This was the objective of a 2016 Google Summer of Code project under Julia-Opt.
A blog describing the progress and issues can be found here - https://ramcha24.github.io/gsoc
Author : Ramchandran (https://github.com/ramcha24)
Mentor : Madeleiene Udell (https://github.com/madeleineudell/)
All queries and comments are welcome
