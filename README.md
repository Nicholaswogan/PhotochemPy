# PhotochemPy
PhotochemPy will be a photochemical model of terrestrial atmospheres. Given inputs, like the stellar UV flux, the atmospheric temperature structure, etc., this code will find the **steady-state** chemical composition of an atmosphere.

I'm working on it right now. It isn't quite done, but is usable. The current version of the code can reproduce older calculations.

## New things

1. **Parallel Processing:** This version of the code can be run in Parallel. On my computer a calculation take about 60% of the time it did in serial. This speedup will be better for some problems and worse for others.

<p align="center">
<img src="Parallel_speed.jpg " width="450">
</p>


2. **Integration with Python:** The code interfaces with Python really well. To try this out run `./compile_parallel` (or `./compile` for serial processing). This will use the tool `f2py` to generate a python module of the fortran code. To test the code run `python test_integrate.py` to integrate to an atmosphere at photochemical equilibrium. Also check out the file `PhotochemPy.py`

3. **Fortran 90 features:** The code is more modern syntax. Variables are dynamically allocated (the code only needs to be compiled once for all problems). No variables are implicit.

## Ideas that I will work on

1. **The way hazes grow and shrink needs to change.** The radius of haze particles on a given timestep depends on the radius of particles from the previous timestep, but not by an ordinary-differential equation. This is a problem because it can cause numerical instability. I need to work something else out.

2. **Using optimization to solve the problem.** This photochemical model seeks a state of photochemical equilibrium. It so by starting with some atmospheric composition, then evolving it forward in time, governed by a large system of ordinary differential equations. After evolving for billions of years, we assume a steady-state is reached. *But there is another way to approach the problem!* The problem can be framed as an optimization problem, where the goal is to find the atmospheric composition that minimizes change in the atmosphere (thus finding the steady-state). It is possible this approach will converge to the steady-state much faster. I'll probably end up using the IPOPT algorithm, or Newthon's method.

3. **Estimating the solution with a 0-D model.** To find a atmosphere in steady-state, you need to start with an intial atmosphere. You will find the steady state faster if the initial atmosphere is close to the solution. If the intial condtions are far from the solution, then the solution may never be found. This is often a problem! One possible way to deal with it though is to estimate the intial atmosphere with simple 0-D photochemical model.

4. **Time-accurate ODE integration**. The code uses the backward euler method to evolve foward in time. This is fine for finding the steady-state, but it isn't ideal for tracking the atmosphere accurately with time. To deal with this I'd like to evolve the atmosphere with the CVODE integration method.

## History
PhotochemPy is an updated version of the `Atmos` photochemical model, originally developed by Jim Kasting and Kevin Zahnle and further developed by many of their students and colleges. The code was originally written in Fortran 77. This version of the code is written in Fortran 90.
