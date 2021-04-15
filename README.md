# N-Link Pendulum

[![DOI](https://zenodo.org/badge/140887199.svg)](https://zenodo.org/badge/latestdoi/140887199)

This set of MATLAB scripts allows one to numerically simulate a two-dimensional N-link pendulum in time. The equations of motion are generated using one of:

* Newton-Euler approach (F = ma) with minimal coordinates (polar coordinates),
* Newton-Euler approach (F = ma) with maximal coordinates (Cartesian coordinates and differential kinematic constraints),
* Lagrange equations approach with minimal coordinates (polar coordinates).

p48_ROOT.m is the main scripts, which calls the other scripts/functions. Parameters are set here, and include

* number of links,
* timespan,
* option to run p48_simple_system,
* option to have nice initial conditions:
  - if true, links start out in-line
  - if false, links start out in a jumbled-fashion (chaotic),
* option to apply torques at the hinges,
* option to include torsional springs at the hinges,
* option to include linear viscous tosional friction at the hinges,
* plot fontsize,
* ODE45 tolerance.

Choosing to run the function p48_simple_system.m will assign simple values to these constants, equal amongst all links/hinges:

* gravitational constant,
* link length,
* distance from upper hinge to link's center of mass,
* link mass,
* link moment of inertia about its center of mass,
* inital angle for nice initial conditions,
* hinge torque,
* torsional spring constant,
* linear viscous friction constant.

The user may add to the appropriate area in p48_ROOT.m if this is not desired. 

The function p48_solve.m derives the symbolic equations of motion according to the approach specified. It then numerically integrates the equations of motion using ODE45, calling the function p48_RHS.m in the process. While solving, the progress is shown using odeprog.m, and is able to be aborted using odeabort.m; both scripts were created by Tim Franklin and retrieved from the MATLAB Central at https://www.mathworks.com/matlabcentral/fileexchange/9904-ode-progress-bar-and-interrupt. 

Instead of running the solver, the user may choose to use a pre-solved set of data in the next section of p48_ROOT.m, from the folder Saved Workspaces. There are quite a few cases included there.

After obtaining a solution to the equations of motion in time, the function p48_plot.m will plot

* the angle of each link in time, wrapped to [-pi,pi],
* the trajectory in (x,y) space of the ends of each link,
* the trajectory in (x,y) space of the end of the last link.

The function p48_energy.m plots the change in kinetic, potential, and total energy as a function of time. If friction and applied torques are omitted, the system is conservative and the total change in energy should always be zero. If large numbers of links are included, finite numerical precision tends to cause the total energy to drift and increase over time. If the maximal coordinates method is used, this also happens since the kinematic constraints are differential. 

The function p48_animate.m animates the system in real time. In the same section, the user can specify whether to include a timer and hinges in the animation, as well as a timescale to slow down or speed up the animation with respect to real time.

The function p48_gif.m then creates an .gif of the animation. This function calls gif.m created by Chad Greene, obtained from MATLAB central at https://www.mathworks.com/matlabcentral/fileexchange/63239-gif.



## History

This was my final project for the course "MAE 5730: Intermediate Dynamics and Vibrations" at Cornell University in the fall of 2017. On the set of problems provided for the semester, this was number 48. My final report is included here, which discusses some of the mathematical framework necessary to derive the equations of motion. I also included a problem from our final exam (3 links), as well as a bonus part (N-links) I added afterward. These two problems explain why for certain initial conditions, link masses, and link moments of inertia, the links remain in-line, and the set of pendulum links behave as one large pendulum. The stability of this case is not analyzed, but may be in the future. 
