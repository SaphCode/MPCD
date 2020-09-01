



\newcommand{\matr}[1]\textbf{#1}
\newcommand{\vect}[1]{\vec{#1}}
\usepackage[table,xcdraw]{xcolor}



# To my dear supervisors

Dear Professor Likos!
Dear Dr. Bianchi!

Professor Likos, I hope you're enjoying your vacation.

You haven't heard from me in a while because I was on vaction too, but I have started my work on the thesis. 
I will review what we discussed in person and online with Dr. Bianchi, and I hope you are as interested in the
research question as I am. Dr. Bianchi, if you have suggestions for adapting dimensions/parameters of the situation (fluid flow through 2d cross section of pipe with cylindrical obstacles) or even tweak the situation, I would love to hear your suggestions too.
I hope I depicted the situation somewhat well below. I thought it would be interesting to not only try obstacles that are
regularly spaced, all the same size, but try a more heterogeneous situation. But that is still a little ways down the road,
so I'll have to see about that.

Progress-wise I'd say I'm still at the beginning. All the work here has been done in about 1 week. I'm working on the collision step now. It took a lot of time to get familiar with C++ again, and to find out how best to test my results. I think I have found a good way to do this, called unit testing, which basically means you test every unit of your program separately and make sure they work under various circumstances, to make debugging easier later down the road and also to add more credibility to my results.

Please do not worry if the writing style is sometimes informal, and at other times formal. I meant it to be that way, as some sections are unfinished/jumbled down thoughts/still expanding/experimenting. Also, the citing style is not the final version either, unfortunately, I cannot change it at the moment.

I'd be happy to hear your guidance and if everything is looking good, or if you suggest I focus more on other areas. If you have a lot to tell me I'd be happy to meet over Zoom/Teams or in person.

Where possible, I've not cited Wikipedia, or added other sources. Is it still okay if I cite Wikipedia, or is that completely undesirable?

All the best,

Chris

# MPCD simulation of polymers in solution

This notebook will serve as the documentation of our efforts and results for my Bachelor's Thesis. The goal of this thesis will be the study of short/long-chained polymers in a liquid thats flowing around obstacles. The liquid will be simulated with MPCD, or "Multi Particle Collision Dynamics".

I chose the language C++ for its familiar object-oriented nature and its proven execution-time. Output of the simulation will mostly be analysed in Python, specifically with Jupyter Notebook for a blend of beautiful visualizations and convenience. For simplification, the situation will be studied in 2D. The situation we specifically discussed is shown below.

![Situation](Assets/MPCD_Situation.png)

### Open questions & Description of situation

For these questions, I need to expand the program to even consider them. If I get close and still have no idea how to go about them, I will contact Maximilian Liebetreu (I think he was there when we first met, I'm sorry I can't remember exactly). Are there other people I can contact for the programming/implementation details?

The fluid flows from the -x direction through the 2D cross-section of some vessel. In this vessel there are cylindrically shaped objects, through which the fluid cannot pass. (->How will I simulate the obstacles, are they solid non-deformable, are they soft & elastic, ..?) (-> Methods: create a lot of particles and immobilize them?, check for collision between particles and obstacles on a per particle basis? (this will probably be too expensive computationally))

What happens when fluid goes out of bounds? (->Should particles that go too far in the +x direction be destroyed and new ones created from the -x direction, simulating continuous flow?)

Finally, polymers will be introduced into the fluid and their behavior studied. We might imagine this "vessel" as a blood vessel and the obstacles as cells moving around in the bloodstream, for example. (-> For the polymers, you told me I should wait until the MPCD is really working)

# Introduction

Multiparticle collision dynamics (MPCD), also known as Stochastic Rotation Dynamics (SRD)[@winkl2009] is a technique originally introduced to study the dynamics of complex fluids such as polymers in solution. Besides MPCD, there exist other mesoscopic models that have been constructed for this purpose, such as Langevin, Direct Simulation Monte Carlo and lattice Boltzmann methods.[@malev1999] We only concern ourselves with the application of MPCD, it follows that any comparison between methods are out of the scope of this thesis.

The MPCD technique models the fluid using particles, their positions and velocities are treated as continuous variables. The system is divided up into cells that have no restriction on the number of particles, each of the cells is part of a regular lattice. The dynamics is split into two parts: Particle streaming and multiparticle collision dynamics. Particle streaming is treated exactly for each particle in the system, while the collision step is approximated on a cell level. The multiparticle collision dynamics conserves mass, momentum and energy and leads to the correct hydrodynamical equations.[@malev1999] The streaming and collision step are described in more detail in (TODO).

# The MPCD algorithm

The system we are modelling consists of $N$ particles with mass $m$, continuous position $\vec{r_{i}}$ and velocity $\vec{v_{i}}$, where $i \in \{1, 2, \dots, N\}$. One timestep $\Delta t$ shall correspond to having calculated all the new particle positions and velocities in the streaming and collision steps, respectively. For each of the N particles, the streaming and collision steps are applied, and this pattern is repeated until the wanted number of timesteps have elapsed.

## The streaming step

The streaming step is very straightforward. The particle positions are simply updated according to

\begin{equation}
\vect{r_{i}} \rightarrow \vect{r_{i}} + \Delta t \cdot \vect{v_{i}}\textrm{,}
\end{equation}

where $\Delta t$ is a small time interval.[@winkl2009][@malev1999]

## The collision step

The collision step is somewhat more complicated. It involves the mean velocity of all particles in a particular cell, $\vect{V_c}$, the velocity of the particle $i$ $\vec{v_i}$ and a rotation matrix $\matr{R}(\alpha)$. The vector $\vect{v_i}$ is rotated relative to the mean velocity $\vect{V_c}$ of all particles in cell $c$, cell $c$ being the cell which particle $i$ belongs to. It is shown in [@malev1999] that the rule,

\begin{equation}
\vect{v_i} \rightarrow \vect{V_c} + \matr{R}(\alpha) [\vect{v_i} - \vect{V_c}] \textrm{,}
\end{equation}

conserves mass, momentum and energy under the molecular chaos assumption[@malev1999? cite differently][@winkl2009, molecular chaos, p.7]. The rotation matrix $\matr{R}(\alpha)$ is a simple 2d rotation matrix

\begin{equation}
R(\alpha) = 
\left[ \begin{array}{rr}
cos(\alpha) & -sin(\alpha) \\
sin(\alpha) & cos(\alpha) \\
\end{array}\right],
\end{equation}

where $\alpha$ is sampled randomly on a per-cell basis. Furthermore, for each particle in the cell $\alpha$ flips its sign with probability $\frac{1}{2}$.[@winkl2009, (p.6)]
The mean velocity of a cell is defined as

\begin{equation}
\vect{V_c} = \frac{1}{N_c} \sum_{i=1}^{N_c} \vect{v_i} \textrm{,}
\end{equation}

where $N_c$ is the number of particles in cell c.[@malev1999]

The original MPCD algorithm was not Galilean invariant. The problem lay in the "molecular chaos" assumption, which means that particles involved in a collision have no memory of earlier encounters when colliding. This assumption is problematic when the mean free path $\lambda = \Delta t \sqrt{\frac{k_{B}T}{m}}$ is small compared to the cell size $a$, since the same particles collide with each other repeatedly and thus build up correlations. When $\lambda \gg a$ Ihle and Kroll have shown that the molecular chaos assumption holds and the simulated results deviate from experimental ones only negligibly.[@ihlekroll2001, p.2][@winkl2009]

The solution to this problem is to shift all particles by the same random vector $s$ before the collision step. The components of $s$ are sampled randomly from a uniform distribution in the interval $[\frac{-a}{2}, \frac{a}{2}]$. After the collision, the particles are shifted back by the same amount.[@ihlekroll2001]

# Implementing the MPCD Algorithm

## (Pseudo) Random Number Generation

One of the pillars of this thesis is the generation of random rotation angles for the rotation matrix needed in the collision step. This proved to be somewhat difficult. First, the standard algorithm of the C++ standard library was tried, but it didn't qualify because it performed poorly in comparison to the second and third algorithms tried, which are called "Mersenne Twister" and "xoshiro256++", respectively.[@wiki:mersennetwister][@cppreference:prng][@unimi:xoshiro]

The Mersenne Twister was implemented using the C++ standard library. The xoshiro256++ was implemented using Sebastian Vigna's code with some additions.[@unimi:xoshiro]

To compare algorithms, and also to make sure that the implementation of the xoshiro256++ is right, a $\chi^2$ test for discrete observations was used. The generated angles in the interval $[0, 2\pi)$ were split into $k+1$ buckets, where $k$ is the number of degrees of freedom of the $\chi^2$ distribution. The test error

\begin{equation}
T = \sum_{b=1}^{k+1}{\frac{(N_o - E[N_b])^2}{E[N_b]}},
\end{equation}

where $E[N_b] = \frac{N}{b}, b \in \{1, 2, \dots , k+1\}$ is the expected bucket size, is compared to $\chi^2_{1-\alpha, k}$, where $\alpha$ is the signifigance level. The null hypothesis

$$
H_0: \textrm{The angles are distributed uniformly in the interval } [0, 2 \pi)
$$

is tested against the alternative hypothesis

$$
H_1: \textrm{The angles are not distributed uniformly in the interval } [0, 2 \pi) \textrm{.}
$$

If the test should have significance level $\alpha$, $H_0$ is rejected if $T \ge \chi^2_{1-\alpha, k}$.[@fruehwirthstat][@wiki:chisquaredtest][@wiki:goodnessoffit]

The results of the $\chi^2$ test are summarised in [TODO: Table, and table formatting].






\begin{table}[]
\begin{tabular}{llll}
\rowcolor[HTML]{4472C4} 
{\color[HTML]{FFF} \textbf{k}} & {\color[HTML]{FFF} \textbf{chi\textasciicircum{}2   probability}} & {\color[HTML]{FFF} \textbf{observed   MT}} & {\color[HTML]{FFF} \textbf{observed XS256++}} \\
\rowcolor[HTML]{D9E1F2} 
1                              & 5.991                                                             & 0.292                                      & 0.068                                         \\
2                              & 7.815                                                             & 1.626                                      & 3.682                                         \\
\rowcolor[HTML]{D9E1F2} 
3                              & 9.488                                                             & 3.735                                      & 2.124                                         \\
4                              & 11.07                                                             & 4.255                                      & 4.525                                         \\
\rowcolor[HTML]{D9E1F2} 
5                              & 12.592                                                            & 2.86                                       & 6.345                                         \\
6                              & 14.067                                                            & 4.071                                      & 6.377                                         \\
\rowcolor[HTML]{D9E1F2} 
7                              & 15.507                                                            & 8.662                                      & 11.731                                        \\
8                              & 16.919                                                            & 14.426                                     & 8.693                                         \\
\rowcolor[HTML]{D9E1F2} 
9                              & 18.307                                                            & 11.732                                     & 6.932                                         \\
10                             & 19.675                                                            & 9.857                                      & 5.939                                         \\
\rowcolor[HTML]{D9E1F2} 
11                             & 21.026                                                            & 10.438                                     & 11.73                                         \\
12                             & 22.362                                                            & 10.985                                     & 15.543                                        \\
\rowcolor[HTML]{D9E1F2} 
13                             & 23.685                                                            & 22.603                                     & 12.022                                        \\
14                             & 24.996                                                            & 13.988                                     & 12.923                                        \\
\rowcolor[HTML]{D9E1F2} 
15                             & 26.296                                                            & 15.526                                     & 20.71                                         \\
16                             & 27.587                                                            & 16.288                                     & 11.501                                        \\
\rowcolor[HTML]{D9E1F2} 
17                             & 28.869                                                            & \cellcolor[HTML]{F8CBAD}32.26              & 12.478                                        \\
18                             & 30.144                                                            & \cellcolor[HTML]{F8CBAD}31.787             & 13.882                                       
\end{tabular}
\end{table}



As we can see, both generators pass the $\chi^2$ test and we do not have to reject our null hypothesis $H_0$.

Visually, we can examine the generated buckets of both random generators in [TODO] the following plot. 



![A histogram of different bucket sizes generated by MT and xoshiro256++](Assets/angle_buckets.png)

The Mersenne Twister has been known to fail certain statistical tests since its inception, by virtue of its mathematical characteristics. There exist other algorithms that are designed to be faster and that do not fail any known statistical tests.[@vigna2019] Ultimately, the xoshiro256++, developed by Sebastian Vigna and David Blackman, was used. It is a variant of the xorshift algorithm, which extends the bit-shift and xor methods by bitrotation, making it still very fast, and more "random" than the xorshift.[@wiki:xorshift][@unimi:xoshiro]

Note that testing a (pseudo) random number generator is usually much more involved than this, but since this has already been done extensively by other authors, we are satisfied with the $\chi^2$ test, simply to test the implementation of the xoshiro256++, since it plays an important part.[@wiki:prng][@vigna2019]

## Particle Streaming

The particle positions were drawn from a uniform real distribution in the interval $[0, 5)$ for the $x$-coordinate, and $[0, 1)$ for the $y$-coordinate. The velocities were initialized to be in the interval $[-1\% \cdot 5, 1\% \cdot 5)$ for the $v_x$ component, and $[-1\% \cdot 1, 1\% \cdot 1)$ for the $v_y$ component. The results can be seen in figure [TODO] below. From the positions in the first row, the velocities in the second row, particle streaming is applied for 1 and 10 timesteps, according to equation (TODO).



![Particle streaming without collision with MT and with xoshiro](Assets/particle_streaming.png)

We see the x and y coordinates are randomly initialized according to the shape of the container. Looking closely, one can see that our particles look very much like noise. The absolute value of the velocity components are initialized to at most 1% of their respective dimensions. After one timestep, some of the particles on the outer ranges have moved out of bounds, and after ten timesteps, the particles have thinned out considerably along the edges.

## The collision step

### Grid

For the implementation of the collision step, a regular lattice, is needed. This grid has lattice constant $a$, where $a$ is a function of the desired average number of particles per cell, which is typically initialised to between three and 20[@winkl2009], although it can be as high as 70[@ihlekroll2001].

The lattice constant $a$ is calculated as follows. The average number of particles per cell, $\bar{N_c}$, is decided upon. If the grid has fixed dimensions, the number of cells is $n = \frac{N}{\bar{N_c}}$. This assumption is valid when only looking at flow through the region of interest and simulating continuous, neverending flow. Practically, this means ignoring particles that go too far out of this region of interest and creating particles that flow into it. If we assign a width $w$ and height $h$ to our region of interest, its area is $A = wh$. But since $n a^{2} = A$, this means that $a = \sqrt{\frac{wh}{n}}$.













## Velocity updating

The velocity of a particle updates according to
\begin{equation}
\vect{v_{i}} \rightarrow \vect{V_c} + \matr{R}(\alpha)[\vect{v_i} - \vect{V_c}],
\end{equation}

where $\vect{v_{i}}$ is the velocity of the $i-th$ particle, $\vect{V_c}$ is the mean velocity of all particles belonging to cell $c$, specified by $i$'s position. The matrix

$$
R(\alpha) = 
\left[ \begin{array}{rr}
cos(\alpha) & -sin(\alpha) \\
sin(\alpha) & cos(\alpha) \\
\end{array}\right]
$$

is a simple 2d-rotation matrix. The angle $\alpha$ is uniformly sampled from the interval $[0, 2\pi)$ on a per-cell basis.[@malev1999]

<!--

CONSIDER THIS. SHOULD BECOME MAXWELL WHEN IMPLEMENTING COLLISION STEP TOO
After having implemented the streaming step and that other one: after time driftting, looks like this: (+velocity distribution = maxwell?)

### some equations to copy

$$
r(t + \Delta t) = r(t) + \Delta t \cdot v(t)
$$

## Converting to Word doc (others possible too, f.ex. .tex)



## Equation Numbering jupyter extension
conda install -c conda-forge jupyter_contrib_nbextensions

jupyter contrib nbextension install --user

jupyter nbextension enable equation-numbering/main

### Turn equation numbering on/off



### Renumber equations



-->
