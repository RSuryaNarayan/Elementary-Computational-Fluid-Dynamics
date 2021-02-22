# Elementary-Computational-Fluid-Dynamics
Hello there! If you are learning CFD from ground zero this is a great place for you to start. This is a repository of MATLAB codes that contains some popular CFD problems one usually encounters when learning CFD for the first time. This repository includes beginner-friendly and readable codes (with adequate comments) for some elementary CFD problems. As one would normally encounter in any CFD course, I have included the finite difference method (FDM) for the Poisson stencil, marched on to diffusion and convection equations with Finite Volume Method(FVM) and finally the Navier stokes itself!(with something called SIMPLE that you can find out more below). Since this is absolutely elementary, I have considered all problems at steady state. Though, it might be a great idea to include a Crank-Nicholson Scheme for hyperbolic and parabolic PDEs. It shall soon be added, stay tuned! For details on which code to refer for a specific problem, read on further below:
# Finite Difference Method-flow through a square pipe
This is a very simple finite difference formulation of the poisson's equation. You can refer to the `pipeflow.m` file for this implementation. Ideally you should get a contour something like this:\
![plot!](https://github.com/RSuryaNarayan/Elementary-Computational-Fluid-Dynamics/blob/master/Pipeflow.PNG)
# FVM-2D-Diffusion
Using a finite volume method on a 2-D structured mesh, the diffusion equation is solved in the `diffusionEquation.m` file. The results look something like so:\
![plot!](https://github.com/RSuryaNarayan/Elementary-Computational-Fluid-Dynamics/blob/master/Diffusion_contours.PNG)
# FVM-2D-Convection
The same except we use a convection stencil. Find it in `convection_equation.m`. Results should like so:\ 
![plot!](https://github.com/RSuryaNarayan/Elementary-Computational-Fluid-Dynamics/blob/master/convection_contours.PNG)
# Navier Stokes - Driven Cavity
A FVM staggered mesh based Navier Stokes solution is implemented in `NavierStokes_Driven_cavity.m`\
![plot!](https://github.com/RSuryaNarayan/Elementary-Computational-Fluid-Dynamics/blob/master/U_velocity_contours.PNG)
