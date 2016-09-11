# Presolve
## Presolve Routines for Optimization Problems
This was part of a 2016 Google Summer of Code project under Julia-Opt.

A blog describing the progress and issues can be found here - [GSoC Blog](https://ramcha24.github.io/gsoc)

Author : [Ramchandran](https://github.com/ramcha24)

GSoC Mentor : [Madeleiene Udell](https://github.com/madeleineudell/)

Currently the package contains an implementation of LP Presolving based on ideas from [Andersen & Andersen paper](http://www.turma-aguia.com/davi/doc/Andersen.pdf) and the GLPK solver.

Presolving removes redundancies from the original problem given by the user and constructs a smaller equivalent problem which is then fed to the solver.
