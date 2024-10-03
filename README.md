# Multiparticle-Collision-Dynamics

## Technique

The MPCD technique models the fluid using particles,
whose positions and velocities are treated as continuous variables.
The system is divided up into cells that have no restriction on the number of particles,
each of the cells is part of a regular lattice.
The dynamics is split into two parts: Particle streaming and multiparticle collision dynamics.
Particle streaming is treated exactly for each particle in the system,
while the particle collisions are approximated on a cell level.
The multiparticle collision dynamics conserves mass,
momentum and energy and leads to the correct hydrodynamical equations. (Malevanets and
Kapral 1999)

The system we want to study consists of a solvent into which polymers are introduced.
The solvent is modelled using $N$ particles with mass $m$, continuous position **r**
and velocity **v**.
One timestep $\Delta t$ shall correspond to having calculated all the new particle positions
and velocities in the streaming and collision steps, respectively.
For each of the $N$ particles, the streaming and collision steps are applied,
and this pattern is repeated until the simulation is terminated.

## Results

![image](https://user-images.githubusercontent.com/10268570/196809846-a06b96d3-cdfd-4c03-a34b-3f08841c6c16.png)

![image](https://user-images.githubusercontent.com/10268570/196809816-200a7af2-1bca-4418-85d0-7fb927f93509.png)

Study of the polymers is more complicated and can be found in the thesis.

Thesis: [Polymer Behavior in narrow-channels blocked by circular obstacles.pdf](https://github.com/SaphCode/MPCD/files/9824520/Polymer.Behavior.in.narrow-channels.blocked.by.circular.obstacles.pdf)
