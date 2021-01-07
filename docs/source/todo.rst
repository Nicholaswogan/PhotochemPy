Improving PhotochemPy
=====================

There are a number of things that could be done to improve PhotochemPy. Here, I'll list a few of them. Any help is very appreciated. If you have any questions please email me at wogan@uw.edu.


- **Generalizing Xsections.f90:** Right now each species has its own subroutine to read in photolysis cross sections in Xsections.f90. This is pretty dumb. I think it is possible to make a single generalized subroutine which reads in cross sections data for any species.

- **Higher wavelength resolution:** PhotochemPy uses 118 wavelength bins. This low wavelength bin resolution requires strange "hacks" to get reasonable photolysis rates for species like NO (see Photo.f90). I believe these strange hacks can be ignored if the bin resolution is increased. This would require changed a fair number of files.

- **New method for determining particle size:** The radius of haze particles on a given timestep depends on the radius of particles from the previous timestep, but not by an ordinary-differential equation (in Sedmnt.f90). This is a problem because it can cause numerical instability. Something else needs to be worked out. This would require changing Sedmnt.f90 and probably integrate.f90.

- **Using optimization to solve the problem:** This photochemical model seeks a state of photochemical equilibrium. It does this by starting with some atmospheric composition, then evolving it forward in time, governed by a large system of ordinary differential equations. After evolving for billions of years, we assume a steady-state is reached. *But there is another way to approach the problem!* The problem can be framed as an optimization problem, where the goal is to find the atmospheric composition that minimizes change in the atmosphere (thus finding the steady-state). It is possible this approach will converge to the steady-state much faster. For example, check out `this paper <https://www.sciencedirect.com/science/article/pii/S2405896319321135>`_.

- **Estimating the solution with a 0-D model:** To find a atmosphere in steady-state, you need to start with an initial atmosphere. You will find the steady state faster if the initial atmosphere is close to the solution. If the initial condtions are far from the solution, then the solution may never be found. This is often a problem! One possible way to deal with it is to estimate the initial atmosphere with simple 0-D photochemical model.

- **Improving the setup.py for the Spike option:** Right now installing PhotochemPy with Spike is kind-of cumbersome and requires many manual steps. It could probably be improved by someone who is good at crafting setup.py files.
