```latex
%%latex
\newcommand{\matr}[1]\textbf{#1}
\newcommand{\vect}[1]{\vec{#1}}
```


\newcommand{\matr}[1]\textbf{#1}
\newcommand{\vect}[1]{\vec{#1}}




```python
import numpy as np
import matplotlib.pyplot as plt
from math import floor
import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import os.path

```


```python
new_analysis = True
with_obstacles = False
```


```python
timesteps = []
I = []
J = []
U = []
V = []
pivot = []
```


```python

path = "Data/"
csv = ".csv"
av = 10
constants = 'constants_av' + str(av) + csv

constants = pd.read_csv(path + constants)
num_timesteps = int(constants['timesteps'])
```


```python
obstacles = None
if (with_obstacles):
    obstacles_path = path + 'constants_obstacles' + csv
    obstacles_csv = pd.read_csv(obstacles_path)
    obstacles = obstacles_csv[['x', 'y']]
    obstacle_radius = float(obstacles_csv['r'][0])
```


```python
shown_x = float(constants['width'])
shown_y = float(constants['height'])

saved = './Saved'
saved_particles = '/particles.pkl'
file_particles = saved + saved_particles

def load_particles(path):
    columns = []
    x_columns = ['x{}'.format(it) for it in range(0, num_timesteps)]
    vx_columns = ['vx{}'.format(it) for it in range(0, num_timesteps)]
    y_columns = ['y{}'.format(it) for it in range(0, num_timesteps)]
    vy_columns = ['vy{}'.format(it) for it in range(0, num_timesteps)]
    columns.extend(x_columns)
    columns.extend(y_columns)
    columns.extend(vx_columns)
    columns.extend(vy_columns)
    particles = pd.DataFrame(columns = columns)

    if (os.path.isfile(file_particles) and not new_analysis):
        print('Found saved particles_x and particles_y files!')
        particles = pd.read_pickle(file_particles)
        print('Loaded particles files.')
    else:
        # Loading particles
        print('Loading particles ..')
        filenames_particles = glob.glob('{}*.csv'.format(path))
        it = 0
        for file in filenames_particles:
            df = pd.read_csv(file)
            particles[['x{}'.format(it), 'y{}'.format(it), 'vx{}'.format(it), 'vy{}'.format(it)]] = df[['x', 'y', 'vx', 'vy']]
            it += 1
            if (it % 10 == 0):
                print('--loaded {}'.format(it))
        particles.to_pickle(file_particles)
        it = 0
        print('Particles loaded and saved!\n')
        # Particles loaded
    return particles


particles_path = '{parent}particles_av{av}'.format(parent = path, av=av)    
particles = load_particles(path = particles_path)
```

    Loading particles ..
    --loaded 10
    --loaded 20
    --loaded 30
    --loaded 40
    --loaded 50
    --loaded 60
    --loaded 70
    --loaded 80
    --loaded 90
    Particles loaded and saved!
    
    


```python
cell_dim = float(constants['cell_dim'])
width = float(constants['width'])
height = float(constants['height'])
shown_cols = floor(width / cell_dim)
shown_rows = floor(height / cell_dim)

saved = './Saved'
saved_cells = '/cells.pkl'
file_cells = saved + saved_cells

columns = ['i', 'j', 'meanX', 'meanY', 'num']

def load_cells(path, columns, shown_rows, shown_cols):
    
    index = ['t']
    if (os.path.isfile(file_cells)):
        cells = pd.read_pickle(file_cells)
        #print(cells)
        #cells.set_index(index, inplace=True)
    else:
        # Loading cells
        print('Loading cells')
        
        cells_timesteps = []
        
        it = 0
        filenames_cells = glob.glob('{}*.csv'.format(path))
        for file in filenames_cells:
            df = pd.read_csv(file)
            cells_timesteps.append(df)
            #df[[i,j]] = (df[[i,j]] + 1/2) * cell_dim
            #print(df.head())
            it += 1
            if (it % 10 == 0):
                print('--loaded {}'.format(it))
        #cells.to_pickle(file_cells)
        print('Cells loaded and saved!\n')
        # Cells loaded
        return cells_timesteps
    
def prepare_cells(cells_timesteps, columns, shown_rows, shown_cols):
    # Preparing cell values
    print('Preparing cell values ..')

    array_i = np.arange(0, shown_rows)
    array_j = np.arange(0, shown_cols)
    I,J = np.meshgrid(array_j, array_i)

    U = []
    V = []

    i = columns[0]
    j = columns[1]
    vx = columns[2]
    vy = columns[3]
    num = columns[4]

    pivots = []
    for df in cells_timesteps:
        # only the rows and cols above 0
        # and below shown_rows, shown_cols
        # this is to 
        # --1. no vaccuum around simulated region
        # --2. I,J are fixed size
        U_inner = []
        V_inner = []
        for it in array_i: # TODO: check this code something seems foul (row, cols, but only using rows)
            temp = df.loc[df[i] == it]
            u = np.array(temp[vx])
            U_inner.append(u)
            v = np.array(temp[vy])
            V_inner.append(v)
        U.append(np.array(U_inner))#, dtype = object))
        V.append(np.array(V_inner))#, dtype = object))

        pivot = df.pivot(index = i, columns = j, values = num)
        pivots.append(pivot)

    print('Cell preparation complete!')
    return I,J,U,V, pivots
    # Cell preparation complete


cells_path = path + 'cells_av{}'.format(av)
cells_timesteps = load_cells(cells_path, columns, shown_rows, shown_cols)
I,J,U,V,pivots = prepare_cells(cells_timesteps, columns, shown_rows, shown_cols)
```

    Loading cells
    --loaded 10
    --loaded 20
    --loaded 30
    --loaded 40
    --loaded 50
    --loaded 60
    --loaded 70
    --loaded 80
    --loaded 90
    --loaded 100
    --loaded 110
    --loaded 120
    --loaded 130
    --loaded 140
    Cells loaded and saved!
    
    Preparing cell values ..
    Cell preparation complete!
    


```python
# Plotting
#with (sns.plotting_context(sns.set())):
x_0_region = 0
x_max_region = int(constants['width'])
y_0_region = 0
y_max_region = int(constants['height'])

streamplot_density = [0.5, 1]

print('Plotting data ..')
fig = plt.figure(figsize=(20,20))
ax = [fig.add_subplot(2,2,i+1) for i in range(4)]

for a in ax:
    a.set_xticklabels([])
    a.set_yticklabels([])
    #a.set_aspect('equal')
    
fig.subplots_adjust(wspace=0, hspace=0)

color = np.sqrt(U[0]**2 + V[0]**2)
point_size = 0.3

#circles = []
if with_obstacles:
    for index, o in obstacles.iterrows():
        circle1 = plt.Circle((o['x'], o['y']), obstacle_radius, color = 'gray')
        circle2 = plt.Circle((o['x'], o['y']), obstacle_radius, color = 'gray')
        #circle3 = plt.Circle((o['x'], o['y']), obstacle_radius, color = 'black')
        #circles.append(circle)
        ax[1].add_artist(circle1)
        ax[2].add_artist(circle2)
        #ax[3].add_artist(circle3)
        #ax[0].xaxis.set_ticks([])
#ax[0].yaxis.set_ticks([])
ax[0].quiver(I, J, U[0], V[0], color)
ax[0].set(xlim=(x_0_region-1,x_max_region), ylim=(y_0_region-1,y_max_region))

#ax[1].xaxis.set_ticks([])
#ax[1].yaxis.set_ticks([])
ax[1].streamplot(I, J, U[0], V[0], color=color, density = streamplot_density) # grid
ax[1].set(xlim=(x_0_region-1,x_max_region), ylim=(y_0_region-1,y_max_region))

#ax[2].xaxis.set_ticks([])
#ax[2].yaxis.set_ticks([])
ax[2].plot(particles['x0'], particles['y0'], "o", markersize = point_size)
if with_obstacles:
    sns.scatterplot(ax = ax[2], x = obstacles['x'], y = obstacles['y'], s = obstacle_radius)
ax[2].set(xlim=(x_0_region, x_max_region), ylim=(y_0_region,y_max_region))

#ax[3].xaxis.set_ticks([])
#ax[3].yaxis.set_ticks([])
#img = ax[3].imshow(pivot, cmap='hot')
#fig.colorbar(img, ax=ax[3], fraction=0.046, pad=0.005)
sns.heatmap(pivots[0], ax=ax[3], xticklabels = False, yticklabels = False, cbar_kws={"fraction": 0.046, "pad": 0.01})
ax[3].set(xlim=(x_0_region-1,x_max_region), ylim=(y_0_region,y_max_region))
#ax[3].xticks('')
#ax[3].yticks('')
ax[3].set_ylabel('')
ax[3].set_xlabel('')
#ax[1,1].imshow(pivot, cmap='hot')

plt.savefig("Assets/initial_region.png")
#plt.close()
print('Data plotted and saved!')
```

    Plotting data ..
    Data plotted and saved!
    


    
![png](thesis_files/thesis_8_1.png)
    



```python
print('Plotting data ..')
fig = plt.figure(figsize=(20,20))
ax = [fig.add_subplot(2,2,i+1) for i in range(4)]

timesteps = int(constants['timesteps']) - 1


for a in ax:
    a.set_xticklabels([])
    a.set_yticklabels([])
    #a.set_aspect('equal')
    
fig.subplots_adjust(wspace=0, hspace=0)

color = np.sqrt(U[-1]**2 + V[-1]**2)
point_size = 0.3

if with_obstacles:
    for index, o in obstacles.iterrows():
        circle1 = plt.Circle((o['x'], o['y']), obstacle_radius, color = 'orange')
        circle2 = plt.Circle((o['x'], o['y']), obstacle_radius, color = 'orange')
        #circle3 = plt.Circle((o['x'], o['y']), obstacle_radius, color = 'black')
        #circles.append(circle)
        ax[1].add_artist(circle1)
        ax[2].add_artist(circle2)
        #ax[3].add_artist(circle3)
        #ax[0].xaxis.set_ticks([])

ax[0].xaxis.set_ticks([])
ax[0].yaxis.set_ticks([])
ax[0].quiver(I, J, U[-1], V[-1], color)
ax[0].set(xlim=(x_0_region-1,x_max_region), ylim=(y_0_region-1,y_max_region))

ax[1].xaxis.set_ticks([])
ax[1].yaxis.set_ticks([])
ax[1].streamplot(I, J, U[-1], V[-1], color=color, density=streamplot_density) # grid
ax[1].set(xlim=(x_0_region-1,x_max_region), ylim=(y_0_region-1,y_max_region))

ax[2].xaxis.set_ticks([])
ax[2].yaxis.set_ticks([])
ax[2].plot(particles['x{}'.format(timesteps)], particles['y{}'.format(timesteps)], "o", markersize = point_size)
ax[2].set(xlim=(x_0_region,x_max_region), ylim=(y_0_region,y_max_region))

ax[3].xaxis.set_ticks([])
ax[3].yaxis.set_ticks([])
#img = ax[3].imshow(pivot, cmap='hot')
#fig.colorbar(img, ax=ax[3], fraction=0.046, pad=0.005)
sns.heatmap(pivots[-1], ax=ax[3], xticklabels = False, yticklabels = False, cbar_kws={"fraction": 0.046, "pad": 0.005})
ax[3].set(xlim=(x_0_region,x_max_region), ylim=(y_0_region,y_max_region))
#ax[3].xticks('')
#ax[3].yticks('')
ax[3].set_ylabel('')
ax[3].set_xlabel('')
#ax[1,1].imshow(pivot, cmap='hot')

plt.savefig("Assets/stationary_region.png")
#plt.close()
print('Data plotted and saved!')
```

    Plotting data ..
    Data plotted and saved!
    


    
![png](thesis_files/thesis_9_1.png)
    


# Todo

1. C:/Users/chris/OneDrive/Documents/Studium/Studiengänge/Bachelor%20Physik/Bachelorarbeit/Cell-level%20canonical%20sampling%20by%20velocity%20scaling%20for%20multiparticlecollision%20dynamics%20simulations.pdf formula (2), $t+\Delta t$ The velocity scaling might not be right because of this mistake?
2. Add virtual particles where the obstacles are. Done
3. Clean code. It's hardly readable. Half done
4. Maybe remove the shared pointers, theyre stupid actually.
5. For \_obstacles, you can actually do it. (make unique f.ex? or just use two vectors for wall and obstacles) They only appear in pipe and simulation, but particles need to stay in simulation, and particles need obstacles to check if they are in bounds initially.
6. For interactors too. (also pipe)
8. Search for TODO
9. check velocity for before and after thermostat/ maybe adjust temperature?
10. at() of map now throws. this should be in grid or smth. fix it.

# Introduction

The aim of this thesis is a study of polymer characterisitcs and behavior in solvent that flows around circular shapes. The study is done coupling two techniques, namely Multiparticle Collision Dynamics (MPCD) and Molecular Dynamics (MD). Multiparticle collision dynamics (MPCD), also known as Stochastic Rotation Dynamics (SRD)[@winkl2009] is a technique originally introduced to study the dynamics of complex fluids. Besides MPCD, there exist other mesoscopic models that have been constructed for this purpose, such as Langevin, Direct Simulation Monte Carlo and lattice Boltzmann methods.[@malev1999] We only concern ourselves with the application of MPCD, it follows that any comparison between methods are out of the scope of this thesis. To model a polymer in this solution, we employed a velocity Verlet-algorithm from the field of Molecular Dynamics and couple it to the MPCD by letting it participate in the collision step. The gyration tensor and the center of mass behavior is recorded and studied for different locations in the flow field. Since this is a work that aims to introduce the reader to the techniques, the simulation is simplified to 2 dimensions.

#### maybe remove

The MPCD technique models the fluid using particles, their positions and velocities are treated as continuous variables. The system is divided up into cells that have no restriction on the number of particles, each of the cells is part of a regular lattice. The dynamics is split into two parts: Particle streaming and multiparticle collision dynamics. Particle streaming is treated exactly for each particle in the system, while the collision step is approximated on a cell level. The multiparticle collision dynamics conserves mass, momentum and energy and leads to the correct hydrodynamical equations.[@malev1999] The streaming and collision step are described in more detail in (TODO: section numbering).

# Description

The solvent flows through a narrow slit-like channel, which can be seen in \ref{fig:situation}. The solvent is a Newtonian liquid which flows according to a pressure gradient in the $\vect{x}$ direction. The collision of MPCD particles with the walls and obstacles was modelled with no-slip boundary conditions. Without the presence of the obstacles, this leads to a parabolic flow profile.

After the simulation has run $60000$ MPCD timesteps, a long-chained polymer of $n=50$ beads is inserted into the solvent and the first $1000$ MPCD timesteps are used for the polymer to get into a polymer like shape to not introduce artificial shapes into the gyration tensor analysis, because the polymer starts out on a line with a random angle. The polymer interacts with the solvent by participating in the collision step, this is described in more detail in section \ref{sec:MolecDyn}. The polymer interacts with the walls and obstacles using a Lennard-Jones potential, spring-like behavior is modelled using a FENE potential. More information can be found in section \ref{sec:molecDyn}.



# The MPCD algorithm

The system we are modelling consists of $N$ particles with mass $m$, continuous position $\vec{r_{i}}$ and velocity $\vec{v_{i}}$, where $i \in \{1, 2, \dots, N\}$. One timestep $\Delta t$ shall correspond to having calculated all the new particle positions and velocities in the streaming and collision steps, respectively. For each of the $N$ particles, the streaming and collision steps are applied, and this pattern is repeated until the simulation is terminated.

## The streaming step

The streaming step is very straightforward. The particle positions are simply updated according to

\begin{equation}
\vect{r_{i}} \rightarrow \vect{r_{i}} + \Delta t \cdot \vect{v_{i}}\textrm{,}
\end{equation}

where $\Delta t$ is a small time interval.[@winkl2009][@malev1999]

## The collision step

The collision step is somewhat more complicated. It involves the mean velocity of all particles in a particular cell, $\vect{V_c}$, the velocity of the particle $i$, $\vec{v_i}$, and a rotation matrix $\matr{R}(\alpha)$. The vector $\vect{v_i}$ is rotated relative to the mean velocity $\vect{V_c}$ of all particles in cell $c$, cell $c$ being the cell which particle $i$ belongs to. It is shown in [@malev1999] that the rule,

\begin{equation}
\vect{v_i} \rightarrow \vect{V_c} + \matr{R}(\alpha) [\vect{v_i} - \vect{V_c}] \textrm{,}
\end{equation}

conserves mass, momentum and energy under the molecular chaos assumption[@malev1999][@winkl2009, molecular chaos, p.7]. The rotation matrix $\matr{R}(\alpha)$ is a simple 2d rotation matrix

\begin{equation}
R(\alpha) = 
\left[ \begin{array}{rr}
cos(\alpha) & -sin(\alpha) \\
sin(\alpha) & cos(\alpha) \\
\end{array}\right],
\end{equation}

where $\alpha$ is sampled on a per-cell basis from the interval $[0, 2\pi)$ with a uniform distribution. Furthermore, for each particle in the cell $\alpha$ flips its sign with probability $\frac{1}{2}$.[@winkl2009, (p.6)]
The mean velocity of a cell is defined as

\begin{equation}
\vect{V_c} = \frac{1}{N_c} \sum_{i=1}^{N_c} \vect{v_i} \textrm{,}
\end{equation}

where $N_c$ is the number of particles in cell c.[@malev1999]

The original MPCD algorithm was not Galilean invariant. The problem lay in the "molecular chaos" assumption, which means that particles involved in a collision have no memory of earlier encounters when colliding. This assumption is problematic when the mean free path 

\begin{equation}
\lambda = \Delta t \sqrt{\frac{k_{B}T}{m}}
\end{equation}

is small compared to the cell size $a$, since the same particles collide with each other repeatedly and thus build up correlations. When $\lambda \gg a$ Ihle and Kroll have shown that the molecular chaos assumption holds and the simulated results deviate from experimental ones only negligibly.[@ihlekroll2001, p.2][@winkl2009]

The solution to this problem is to shift all particles by the same random vector $s$ before the collision step. The components of $s$ are sampled randomly from a uniform distribution in the interval $[-\frac{a}{2}, \frac{a}{2}]$. After the collision, the particles are shifted back by the same amount.[@ihlekroll2001]

## Additions to the original MPCD Algorithm

### Grid shift

After particle streaming, the grid is shifted. The components of this shift, $\vect s$ are sampled from a uniform random distribution in the interval $[-\frac{a}{2},\frac{a}{2}]$. This step is necessary to restore Galilean invariance, which is violated when the molecular chaos assumption does not hold. This happens when simulating cold fluids or when using very small timesteps.[@ihlekroll2001]
The shift is undone at the end of the collision step, after the velocity of all particles has been updated.

If particles go out of the simulation bounds due to the grid shift, they reappear at the other end as described in section (TODO: section numbering). Because the wall will not align with the grid anymore, virtual particles as described in section (TODO: section numbering) have to be introduced.

### Interaction of fluid particles with obstacles

The fluid flows through a pipe that will be setup somewhere between 100 to 400 width, and 20 - 50 height, as can be seen in figure (TODO: figure number). In SI units, we might imagine .. (TODO: expand this section). The pipe has two parallel walls, the fluid-wall interaction is modeled using stick (or no-slip) boundary conditions. (TODO: figure). When a particle hits the wall, it goes back the same way it came there, which means that the sign of the velocity vector is flipped. Stick boundary conditions are shown in figure. (TODO: figure numbering) The fluid interacts with obstacles, which are modeled exactly like the wall, with a small complication from the geometry. To simplify this problem, the approximate collision process found in [@nikoubashman2013] is used.

![Container of MPCD fluid, with obstacles](Release/Assets/MPCD_Pipe.png)

![Reflection of particle from wall, no slip boundary conditions](Release/Assets/Wall_noslip_boundary_conditions_reflection.png)

Because of the shifting of the grid before the collision step, the cells next to the walls might be partially blocked by the wall. This partial blocking by the wall causes the cell to have, on average, less particles than it would have, had the grid not been shifted. The change in average particles distorts the collision step of particles near the wall. For more complex geometries than the one used in this thesis, the wall wont even be parallel to the grid lines, which makes this a problem even without a grid shift. To compensate this, it has been shown in (TODO: find a source) that the following process undoes the distortion.

As is shown in figure (TODO: figure numbering), imagine a cell being blocked a little bit by the wall. Let's assume that the average number of particles per cell is 4. In the first cell, we count 3 particles. What is now done, is to introduce a "virtual" particle, which we might imagine as being behind the wall, though the position of it does not matter. Jumping to the second cell, we introduce 2 "virtual particles". In general, $\bar N_c - N_c$ particles are introduced, where $\bar N_c$ is the average number of particles per cell, chosen to be an integer, and $N_c$ is the actual number of particles found in cell $c$. Their velocities are sampled from two independent normal distributions with mean $\mu = 0$ and variation $\sigma^2 = \frac{k_{B}T}{m}$, where $k_{B}$ is the boltzmann constant, $T$ is the temperature of the fluid, and $m$ is the mass of one particle.

![Undoing the distortion caused by the grid shift](Release/Assets/Wall_stick_boundary_conditions_virtual_particles.png)

### The ballistics step

The ballistics step might be called a substep of particle streaming. After the particles are moved, their position is checked for interaction with the wall or obstacles. If particles overshoot the bounds set by either the wall or obstacles, their position is set to the collision point, their velocity is reversed, subsequently they are moved for the rest of the distance they would have travelled, as explained in section. (TODO: section numbering) This means we are assuming elastic collision between the particles, the walls and obstacles.

### Modelling poseuille flow

Several methods exist to model poseuille flow.[@nikoubashman2013] Here, a gravitational approach is chosen. This means an external force acts on the unit volume of the fluid, which is given by,

\begin{equation}
\vect{F} = \rho_S g \hat x,
\end{equation}

where $\rho_S$ is the mass density of the solvent, $g$ is the acceleration constant and $\hat x$ is the unit vector in the $x$-direction. This force is incorporated into the simulation by updating the position and velocity of a particle according to the solution for constant force Newton's equations.[@nikoubashman2013]

Because now a force acts on the particles, as seen in figure (TODO figure#), additional energy is coming into the system which needs to be controlled. The tool to do this in MPCD is called a thermostat.

The constant force in combination with no-slip boundary conditions applied by the bounce-back rule leads to a parabolic velocity profile

\begin{equation}
v_x(y) = \frac{4 v_{max} (D - y) y}{D^2},
\end{equation}

where $D$ is the distance between walls. The maximum velocity of this profile is given by

\begin{equation}
v_{max} = \frac{m N_c g D^2}{8 \eta},[@huang2010]
\end{equation}

where $\eta$ refers to the viscosity of the fluid.
The profile obtained from the simulation is plotted against the expected profile in figure (TODO: figure #)

### Thermostat

To counteract the heating up of the fluid by the constant force on every particle, a thermostat is needed. As described in Winkler[@winkl2009], the cell-level thermostat was used. For this, we need to sample a kinetic energy from the gamma distribution

\begin{equation}
P(E_k) = \frac{q}{E_k\Gamma(f/2)}\left(\frac{E_k}{k_BT}\right)^{f/2} e^{-\frac{E_k}{k_B T}},
\end{equation}

where $f = 2 (N_c - 1)$ are the degrees of freedom of the system. The velocities (TODO: which velocities? I scaled the $\Delta v_i$ and the $v_i$, and both gave the same effect) are then scaled by the factor 

\begin{equation}
\kappa = \left(\frac{2E_k}{\sum_{i=1}^{N_c}m\Delta v_i^2}\right)^{1/2},
\end{equation}

where $\Delta v_i = v_i - V_c$. The effect of the thermostat can be seen in figure (TODO: figure#), in which the theoretical distribution[@winkl2009]

\begin{equation}
P(\Delta v) = \sum_{N_c = 2}^{\infty} e^{-\bar{N_c}} \frac{\bar{N_c}^{N_c}}{N_c!} \frac{P(\Delta v, N_c)} { \left( 1 - (\bar{N_c} + 1) \cdot e^{-\bar{N_c}} \right)},
\end{equation}
where 

\begin{equation}
P(\Delta v, N_c) = \left(\frac{m}{2 \pi k_B T (1 - 1/N_c)}\right)^{3/2} exp\left(-\frac{m}{2 k_B T (1 - 1/N_c)} \Delta v^2 \right),
\end{equation}

is plotted against the distribution of $|\Delta v|$ in a collision cell.


### Collision with circular obstacles

The collision rules are detailed in figure (TODO: figure#).

![Collision with circular obstacle](Assets/MPCD_obstacle_collision_rule.png)

From the particle\'s position before the streaming step \vect{p_0}, to the particle\'s position after moving, \vect{p_1}, a line with equation

\begin{equation}
(y - y_0) = k (x - x_0),
\end{equation}

where $k = \frac{y_1 - y_0}{x_1 - x_0}$ is drawn. This line intersects with a circle with equation

\begin{equation}
(x - x_c)^2 + (y - y_c)^2 = r^2,
\end{equation}

where $(x_c, y_c)$ is the center position of the circle. By solving the equation, we get 2 collision points, the rule for choosing the right one is then to take the one that's on the direct path from the new position to the old. We calculate the quantity by which the particle has "overshot" the collision point, and correct its position to $\vect{p_{correct}} = \vect{p_1} - 2 \vect{o},$ where $\vect{o}$ is the amount which the particle has gone too far.

# Molecular Dynamics

## Verlet algorithm
To simulate the physics of a polymer dissolved in our solvent, we employ Molecular Dynamics. In particular, the velocity verlet-algorithm 

\begin{equation}
\vect{r}(t + \Delta t) = \vect{r}(t) + \Delta t \vect{v}(t) + \frac{1}{2} \Delta t^{2} \vect{a}(t)
\vect{v}(t + \Delta t) = \vect{v}(t) + \frac{1}{2} \Delta t (\vect{a}(t) + \vect{a}(t + \Delta t)
\end{equation}

was used. [@allenTildesley] 

## Polymers

### Interaction

For the simulation, we used long-chained polymers with a chain length of $n_c = 50$, the interaction of the monomers is modelled by a truncated-shifted Lennard Jones potential[@weiss2017]

\begin{equation}
V_{LJ} = \sum_{i=1}^{n_c} \sum_{j=i+1}^{n_c} 4 \epsilon \left[\left(\frac{\sigma}{r_{i,j}}\right)^{12} - \left(\frac{\sigma}{r_{i,j}}\right)^6 + \frac{1}{4}\right] \Theta(2^{1/6} \sigma - r_{i,j})
\end{equation}

where $\epsilon$ controls the interaction strength, $r_{i,j}$ is the distance between two monomers $i$ and $j$, $\sigma$ is the mean of the diameters of the two interaction partners $i$ and $j$ and $\Theta$ is the Heaviside function. This same function was also used to model the monomer-obstacle interaction.

The attraction between neighboring beads of the polymer was modelled using a FENE potential [@weiss2017]

\begin{equation}
V_{FENE} = \sum_{i=1}^{n_c-1} \frac{k}{2} R_0^2 ln\left(1-\left(\frac{r_{i,i+1}}{R_0}\right)^2\right)
\end{equation},

where $k$ denotes the stiffness of the bond and $R_0$ is the maximum bond extension.

The interaction of the polymers with the walls is controlled again by a truncated-shifted Lennard Jones potential[@weiss2017]

\begin{equation}
V_{W} = \sum_{i=1}^{n_c} \epsilon_{W} \left[ \frac{2}{15} \left( \frac{\sigma_{W}}{y_i} \right)^9 - \left( \frac{\sigma_W}{y_i} \right)^3 + \frac{\sqrt{10}}{3} \right] \Theta\left(\left(\frac{2}{5}\right)^{1/6} \sigma_W - y_i\right)
\end{equation},

where $y_i$ is the vertical distance of monomer $i$ to the wall, $\sigma_W = \sigma$ from before and $\epsilon_W$ controls the interaction strength with the wall.

The coupling of the MPCD algorithm to the MD algorithm is achieved by letting the monomers participate in the MPCD collision step.[@weiss2017]



```python
def viscosity(v_max, m = 1, av_N_c = 10, g = 0.01, D = 20):
    return (m * av_N_c * g * D**2) / (8 * v_max)

def parabolic_flow(y, v_max = 4, D = 20):
    return (4 * v_max * (D - y) * y)/(D**2)

```


```python
y = np.linspace(0, 20, 100)
plt.plot(y, parabolic_flow(y))
plt.title('Parabolic flow')

rowsums = [row.sum()/len(row) for row in U[-1]]
plt.plot(range(len(rowsums)), rowsums, "o-")
plt.xlabel('x-Velocity')
plt.ylabel('Cell number')
plt.title('Velocity profile, averaged over all rows, for the last timestep ({})'.format(timesteps))
plt.savefig("Assets/velocity_profile.png")
#plt.close()
```


    
![png](thesis_files/thesis_14_0.png)
    



```python
import math

def theoretical_dist(Delta_v, N_c, av_N_c = 10, m = 1, k_B = 1, T = 1):
    A = (m/(2 * math.pi * k_B * T * (1 - 1/N_c)))**(3/2)
    B = np.exp(-(m/(2 * k_B * T * (1 - 1/N_c))) * Delta_v**2)
    return A * B

def poisson_avg(Delta_v, max_sum = 1000, av_N_c = 10, m = 1, k_B = 1, T = 1):
    e = math.exp(-av_N_c)
    total = 0
    for N_c in range(2, max_sum):
        P = (av_N_c**N_c)/(math.factorial(N_c))
        N = theoretical_dist(Delta_v, N_c, m, k_B, T)
        total += P * (N)
    D = 1 - (av_N_c + 1) * e
    return e * total/D

def plot_theoretical():
    Delta_v = np.linspace(0, 4, 1000)
    plt.plot(Delta_v, poisson_avg(Delta_v))
    plt.show()
    
plot_theoretical()
```


    
![png](thesis_files/thesis_15_0.png)
    


![Initial State of region with barebones MPCD implementation](Assets/initial_region.png)

![Ending (maybe stationary) state of region with barebones MPCD implementation](Assets/stationary_region.png)

![Ending velocity profile](Release/Assets/velocity_profile.png)

# Histograms velocity

Function for animation of histogram


```python
from matplotlib import animation

def update_hist(num, data):
    if (num % 10 == 0):
        print('Frame: {}'.format(num))
    plt.cla()
    plt.hist(data[num], bins = 20)
```

## vx


```python
x_velocities = [particles['vx{}'.format(i)] for i in range(timesteps)]

n_frames = num_timesteps - 1

fig = plt.figure()
x_velocities_hist = plt.hist(x_velocities[0])

anim = animation.FuncAnimation(fig, update_hist, n_frames, fargs=(x_velocities,))
anim.save('./Assets/x_vel_hist.mp4', fps = 10)
plt.show()
```

    Frame: 0
    Frame: 0
    Frame: 10
    Frame: 20
    Frame: 30
    Frame: 40
    Frame: 50
    Frame: 60
    Frame: 70
    Frame: 80
    Frame: 90
    

    C:\Users\chris\anaconda3\envs\datascience\lib\site-packages\matplotlib\axes\_axes.py:6620: RuntimeWarning: All-NaN slice encountered
      xmin = min(xmin, np.nanmin(xi))
    C:\Users\chris\anaconda3\envs\datascience\lib\site-packages\matplotlib\axes\_axes.py:6621: RuntimeWarning: All-NaN slice encountered
      xmax = max(xmax, np.nanmax(xi))
    


    ---------------------------------------------------------------------------

    ValueError                                Traceback (most recent call last)

    <ipython-input-15-841a2dc48ddb> in <module>
          7 
          8 anim = animation.FuncAnimation(fig, update_hist, n_frames, fargs=(x_velocities,))
    ----> 9 anim.save('./Assets/x_vel_hist.mp4', fps = 10)
         10 plt.show()
    

    ~\anaconda3\envs\datascience\lib\site-packages\matplotlib\animation.py in save(self, filename, writer, fps, dpi, codec, bitrate, extra_args, metadata, extra_anim, savefig_kwargs, progress_callback)
       1139                 for anim, d in zip(all_anim, data):
       1140                     # TODO: See if turning off blit is really necessary
    -> 1141                     anim._draw_next_frame(d, blit=False)
       1142                     if progress_callback is not None:
       1143                         progress_callback(frame_number, total_frames)
    

    ~\anaconda3\envs\datascience\lib\site-packages\matplotlib\animation.py in _draw_next_frame(self, framedata, blit)
       1174         # post- draw, as well as the drawing of the frame itself.
       1175         self._pre_draw(framedata, blit)
    -> 1176         self._draw_frame(framedata)
       1177         self._post_draw(framedata, blit)
       1178 
    

    ~\anaconda3\envs\datascience\lib\site-packages\matplotlib\animation.py in _draw_frame(self, framedata)
       1724         # Call the func with framedata and args. If blitting is desired,
       1725         # func needs to return a sequence of any artists that were modified.
    -> 1726         self._drawn_artists = self._func(framedata, *self._args)
       1727         if self._blit:
       1728             if self._drawn_artists is None:
    

    <ipython-input-14-29b224c71586> in update_hist(num, data)
          5         print('Frame: {}'.format(num))
          6     plt.cla()
    ----> 7     plt.hist(data[num], bins = 20)
    

    ~\anaconda3\envs\datascience\lib\site-packages\matplotlib\pyplot.py in hist(x, bins, range, density, weights, cumulative, bottom, histtype, align, orientation, rwidth, log, color, label, stacked, data, **kwargs)
       2667         orientation='vertical', rwidth=None, log=False, color=None,
       2668         label=None, stacked=False, *, data=None, **kwargs):
    -> 2669     return gca().hist(
       2670         x, bins=bins, range=range, density=density, weights=weights,
       2671         cumulative=cumulative, bottom=bottom, histtype=histtype,
    

    ~\anaconda3\envs\datascience\lib\site-packages\matplotlib\__init__.py in inner(ax, data, *args, **kwargs)
       1436     def inner(ax, *args, data=None, **kwargs):
       1437         if data is None:
    -> 1438             return func(ax, *map(sanitize_sequence, args), **kwargs)
       1439 
       1440         bound = new_sig.bind(ax, *args, **kwargs)
    

    ~\anaconda3\envs\datascience\lib\site-packages\matplotlib\axes\_axes.py in hist(self, x, bins, range, density, weights, cumulative, bottom, histtype, align, orientation, rwidth, log, color, label, stacked, **kwargs)
       6646             # this will automatically overwrite bins,
       6647             # so that each histogram uses the same bins
    -> 6648             m, bins = np.histogram(x[i], bins, weights=w[i], **hist_kwargs)
       6649             tops.append(m)
       6650         tops = np.array(tops, float)  # causes problems later if it's an int
    

    <__array_function__ internals> in histogram(*args, **kwargs)
    

    ~\anaconda3\envs\datascience\lib\site-packages\numpy\lib\histograms.py in histogram(a, bins, range, normed, weights, density)
        790     a, weights = _ravel_and_check_weights(a, weights)
        791 
    --> 792     bin_edges, uniform_bins = _get_bin_edges(a, bins, range, weights)
        793 
        794     # Histogram is an integer or a float array depending on the weights.
    

    ~\anaconda3\envs\datascience\lib\site-packages\numpy\lib\histograms.py in _get_bin_edges(a, bins, range, weights)
        424             raise ValueError('`bins` must be positive, when an integer')
        425 
    --> 426         first_edge, last_edge = _get_outer_edges(a, range)
        427 
        428     elif np.ndim(bins) == 1:
    

    ~\anaconda3\envs\datascience\lib\site-packages\numpy\lib\histograms.py in _get_outer_edges(a, range)
        321         first_edge, last_edge = a.min(), a.max()
        322         if not (np.isfinite(first_edge) and np.isfinite(last_edge)):
    --> 323             raise ValueError(
        324                 "autodetected range of [{}, {}] is not finite".format(first_edge, last_edge))
        325 
    

    ValueError: autodetected range of [nan, nan] is not finite



    
![png](thesis_files/thesis_23_3.png)
    


## vy


```python
last_vels = particles['vy{}'.format(0)]
sns.distplot(last_vels)
plt.xlabel('vy, timestep 0')
plt.ylabel('density')
plt.title('distribution of y velocity at t=0')
```


```python
last_vels = particles['vy{}'.format(timesteps)]
sns.distplot(last_vels)
plt.xlabel('vy, last timestep')
plt.ylabel('density')
plt.title('distribution of y velocity at t=last')
```

# Temperature of solvent


```python
v_2_t = [particles[['vx{}'.format(t), 'vy{}'.format(t)]]**2 for t in range(round(len(particles.columns)/4 - 1))]
#print(v_2_t)
avg_v_2_t = [v_2_particles.sum().sum()/len(v_2_particles) for v_2_particles in v_2_t]# first col sum, then row sum (sum vxs, then sum vys)
print('Average sqared velocities after first timestep:\n{}'.format(v_2_t[0].sum()/len(v_2_t[0])))
print('Average v_x2 after last timestep:\n{}'.format(v_2_t[-1].sum()/len(v_2_t[-1])))

particle_mass = float(constants['particle_mass'])

temp_t = np.array(avg_v_2_t)
plt.plot(temp_t)
#plot = sns.lineplot(data=temp_t)
plt.xlabel('Timestep')
plt.ylabel('average of v2')
plt.title('Temperature is increasing')
plt.xlim(0,1000)
plt.savefig('Assets/temperature.png')
```

    Average sqared velocities after first timestep:
    vx0    2.133656
    vy0    1.120110
    dtype: float64
    Average v_x2 after last timestep:
    vx4998    0.0
    vy4998    0.0
    dtype: float64
    


    
![png](thesis_files/thesis_28_1.png)
    


# Animations


```python
from matplotlib import animation
print('Animating Quiver ...\n')

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(x_0_region, x_max_region), ylim=(y_0_region, y_max_region))

quiv = ax.quiver(I, J, 10**-12*U[0], 10**-12*V[0], scale = 1) #10**-5*U10**-5*V

# initialization function: plot the background of each frame
def init():
    #quiv.set_data([], [], [], [])
    return quiv,

# animation function.  This is called sequentially
def animate(it, quiv, I, J):
    #color = np.sqrt(U[it]**2 + V[it]**2)
    quiv.set_UVC(10**-12*U[it], 10**-12*V[it]) #*10**-5*
    if (it % 10 == 0):
        print('--Created {} frame.\n'.format(it))
    #quiv.set_color(color)
    return quiv,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, fargs = (quiv, I, J),#init_func=init,
                               frames=timesteps, blit=False)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('./Assets/quiver_animation.mp4', fps=5) #extra_args=['-vcodec', 'libx264'])
print('Animated and saved!')

plt.close()
```


```python
from matplotlib import animation
import seaborn as sns

print('Animating Heatmap ...\n')

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(x_0_region, x_max_region), ylim=(y_0_region, y_max_region))
sns.heatmap(pivots[0], square = True, xticklabels = False, yticklabels = False, cbar = False)#,cbar_kws={"fraction": 0.046, "pad": 0.005})

# initialization function: plot the background of each frame
def init():
    plt.clf()
    ax = sns.heatmap(pivots[0], square = True, xticklabels = False, yticklabels = False, cbar = False)#,cbar_kws={"fraction": 0.046, "pad": 0.005})
    #return heatmap,

# animation function.  This is called sequentially
def animate(it):
    plt.clf()
    ax = sns.heatmap(pivots[it], square = True, xticklabels = False, yticklabels = False, cbar = False)#,cbar_kws={"fraction": 0.046, "pad": 0.005})
    ax.set_xlim((x_0_region, x_max_region))
    ax.set_ylim((y_0_region, y_max_region))
    if (it % 10 == 0):
        print('--Created {} frame.\n'.format(it))
    #return heatmap

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=timesteps)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('./Assets/heatmap_animation.mp4', fps=5) #extra_args=['-vcodec', 'libx264'])
print('Animated and saved!')

plt.close()
```


```python
from matplotlib import animation

print('Animating particles ...')
# First set up the figure, the axis, and the plot element we want to animate
point_size = 0.1
fig = plt.figure()
ax = plt.axes(xlim=(x_0_region, x_max_region), ylim=(0, y_max_region))
scatter, = ax.plot(particles['x0'], particles['y0'], "o", markersize = point_size)

x = 'x'
y = 'y'

# initialization function: plot the background of each frame
def init():
    #quiv.set_data([], [], [], [])
    return scatter,

# animation function.  This is called sequentially
def animate(it):
    scatter.set_xdata(particles[x + str(it)])
    scatter.set_ydata(particles[y + str(it)])
    if (it % 10 == 0):
        print('--Created {} frame.\n'.format(it))
    return scatter,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, #init_func=init,
                               frames=timesteps, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('./Assets/scatterplot_animation.mp4', fps=5, extra_args=['-vcodec', 'libx264'])
print('Animated and saved!')

plt.close()
```


```python
from matplotlib import animation

print('Animating Streamplot ...\n')

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(x_0_region, x_max_region), ylim=(y_0_region, y_max_region))
color = np.sqrt(U[0]**2 + V[0]**2)
stream = ax.streamplot(I, J, U[0], V[0], color=color, density = streamplot_density)

# initialization function: plot the background of each frame
def init():
    #quiv.set_data([], [], [], [])
    return stream

# animation function.  This is called sequentially
def animate(it):
    ax.collections = [] # clear lines streamplot
    ax.patches = [] # clear arrowheads streamplot
    color = np.sqrt(U[it]**2 + V[it]**2)
    stream = ax.streamplot(I, J, U[it], V[it], color=color, density = streamplot_density)
    if (it % 10 == 0):
        print('--Created {} frame.\n'.format(it))
    return stream

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, #init_func=init,
                               frames=timesteps, blit=False)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('./Assets/streamplot_animation.mp4', fps=5) #extra_args=['-vcodec', 'libx264'])
print('Animated and saved!')

plt.close()
```

#### Conservation of number of particles

The number of particles is (just for convenience) plotted and it stays constant.


```python
import matplotlib.pyplot as plt
nums = []
for pivot in pivots:
    nums.append(sum(sum(lis) for lis in pivot.values))

plt.plot(nums)
plt.title('Variation in number of particles with timesteps')
plt.ylabel('Number of Particles')
#plt.ylim((79990, 80010))
plt.xlabel('Timestep')
plt.savefig('./Assets/number_particles.png')
plt.close()
```

![Constant number of particles](Assets/number_particles.png)

#### Conservation of momentum

Because the angle of rotation in the collision step is chosen randomly, and the streaming step does not change momentum, in a large system momentum should be conserved. The velocities of particles were taken and added up. The base momentum is the initial momentum, the error (or variation) from this base is plotted below.


```python
momentum_x = []
momentum_y = []
for it in range(0, timesteps-1):
    vx = sum(particles['vx{}'.format(it)])
    vy = sum(particles['vy{}'.format(it)])
    momentum_x.append(vx)
    momentum_y.append(vy)

base_momentum_x = momentum_x[0]
base_momentum_y = momentum_y[0]
error_x = [base_momentum_x - m for m in momentum_x]
error_y = [base_momentum_y - m for m in momentum_y]
plt.plot(error_x)
plt.title('Variation in x-Velocity, initially: {}'.format(round(base_momentum_x, 2)))
plt.xlabel('Timestep')
plt.ylabel('x-Velocity')
plt.savefig('./Assets/x_velocity_variation.png')
plt.show()
#plt.close()
plt.plot(error_y)
plt.title('Variation in y-Velocity, initally: {}'.format(round(base_momentum_y, 2)))
plt.xlabel('Timestep')
plt.ylabel('y-Velocity')
plt.savefig('./Assets/y_velocity_variation.png')
plt.show()
#plt.close()
plt.plot(error_x, error_y, "o", markersize = 3)
plt.title('Variation in velocity')
plt.xlabel('Variation in x-velocity')
plt.ylabel('Variation in y-velocity')
plt.savefig('./Assets/velocity_variation.png')
plt.show()
#plt.close()
# -----------------------------------


#plt.close()
```

![Variation in x-velocity throughout simulation](Assets/x_velocity_variation.png)

![Variation in y-velocity throughout simulation](Assets/y_velocity_variation.png)

![Variation in velocity throughout simulation](Assets/velocity_variation.png)

Error $\sim 10^{-5}$

#### Energy

The collision step of the MPCD algorithm conserves energy _locally_, which is to say on a cell level. [@winkl2009] The energy should also be conserved globally, since no force is acting upon the particles _yet_, the streaming and collision steps conserve energy, and the particle number remains the same.

To inspect this, the energy of every particle is added up. The base energy is the initial energy, the error (or variation from this base) is calculated and plotted below.


```python
'''
def square(lis):
    for e in lis:
        yield e**2

xvels_squared = []
for cell_xvels in U:
    xvels_squared.append(sum(sum(square(lis)) for lis in cell_xvels))
#print(xvels)
yvels_squared = []
for cell_yvels in V:
    yvels_squared.append(sum(sum(square(lis)) for lis in cell_yvels))
assert(len(xvels) == len(yvels))
it = 0
energy_cell_level = []
mass = 2.988e-26
for xvel_squared in xvels_squared:
    yvel_squared = yvels_squared[it]
    energy_cell_level.append((xvel_squared + yvel_squared)) # might add mass here
    
plt.plot(energy_cell_level)
plt.title('Variation in energy, cell method')
plt.xlabel('Timestep')
plt.ylabel('Energy')
plt.savefig('./Assets/constant_energy_cellcalc.png')
plt.close()
'''
```

Error $\sim 10^{-8}$

## Converting to Word doc (others possible too, f.ex. .tex)


```python
# just do it manually, it works on anaconda env datascience

import subprocess
#automatic document conversion to markdown and then to word
#first convert the ipython notebook paper.ipynb to markdown
subprocess.run("jupyter nbconvert --to markdown thesis.ipynb --output-dir='./Generated'") #--output-dir='./Generated'
#next remove code
path = "./Generated/thesis.md"
with open(path, "r") as f:
    lines = f.readlines()
    idx = []
    idx_files = []
    for i, line in enumerate(lines):
        if (line.startswith("```")):
            idx.append(i)
        if ("thesis_files" in line):
            c = line.find("thesis_files")
            lines[i] = line[0:c] + "" + line[c:]

idx = sorted(idx, reverse=True) # reverse order so not deleting lines and then missing others
for current, previous in zip(idx[::2], idx[1::2]):
    print("Deleting {p}:{c}".format(p=previous, c=current+1))
    print('\n'.join(lines[previous:current+1]))
    del lines[previous:current+1]
    
with open(path, "w") as f:
    #f.write("\\newcommand{\matr}[1]\\textbf{#1}")
    #f.write("\\newcommand{\\vect}[1]{\\vec{#1}}")
    for line in lines:
        f.write("%s" % line)
#next convert markdown to ms word
conversion_tex = "pandoc -s ./Generated/thesis.md -o ./Generated/thesis.tex --filter pandoc-citeproc --bibliography=\"list.bib\" --csl=\"apa.csl\""
subprocess.run(conversion_tex)
conversion_pdf = "pandoc -s ./Generated/thesis.md -o ./Generated/thesis.pdf --filter pandoc-citeproc --bibliography=\"list.bib\" --csl=\"apa.csl\""
subprocess.run(conversion_pdf)
# LATEX TO DOCX pandoc -s math.tex -o example30.docx
```

    Deleting 2686:2689
    ```python
    
    
    
    ```
    
    Deleting 2674:2682
    ```javascript
    
    %%javascript
    
    MathJax.Hub.Queue(
    
      ["resetEquationNumbers", MathJax.InputJax.TeX],
    
      ["PreProcess", MathJax.Hub],
    
      ["Reprocess", MathJax.Hub]
    
    );
    
    ```
    
    Deleting 2664:2670
    ```javascript
    
    %%javascript
    
    MathJax.Hub.Config({
    
        TeX: { equationNumbers: { autoNumber: "AMS" } }
    
    });
    
    ```
    
    Deleting 1178:1216
    ```python
    
    # just do it manually, it works on anaconda env datascience
    
    
    
    import subprocess
    
    #automatic document conversion to markdown and then to word
    
    #first convert the ipython notebook paper.ipynb to markdown
    
    subprocess.run("jupyter nbconvert --to markdown thesis.ipynb --output-dir='./Generated'") #--output-dir='./Generated'
    
    #next remove code
    
    path = "./Generated/thesis.md"
    
    with open(path, "r") as f:
    
        lines = f.readlines()
    
        idx = []
    
        idx_files = []
    
        for i, line in enumerate(lines):
    
            if (line.startswith("```")):
    
                idx.append(i)
    
            if ("thesis_files" in line):
    
                c = line.find("thesis_files")
    
                lines[i] = line[0:c] + "" + line[c:]
    
    
    
    idx = sorted(idx, reverse=True) # reverse order so not deleting lines and then missing others
    
    for current, previous in zip(idx[::2], idx[1::2]):
    
        print("Deleting {p}:{c}".format(p=previous, c=current+1))
    
        print('\n'.join(lines[previous:current+1]))
    
        del lines[previous:current+1]
    
        
    
    with open(path, "w") as f:
    
        #f.write("\\newcommand{\matr}[1]\\textbf{#1}")
    
        #f.write("\\newcommand{\\vect}[1]{\\vec{#1}}")
    
        for line in lines:
    
            f.write("%s" % line)
    
    #next convert markdown to ms word
    
    conversion_tex = "pandoc -s ./Generated/thesis.md -o ./Generated/thesis.tex --filter pandoc-citeproc --bibliography=\"list.bib\" --csl=\"apa.csl\""
    
    subprocess.run(conversion_tex)
    
    conversion_pdf = "pandoc -s ./Generated/thesis.md -o ./Generated/thesis.pdf --filter pandoc-citeproc --bibliography=\"list.bib\" --csl=\"apa.csl\""
    
    subprocess.run(conversion_pdf)
    
    # LATEX TO DOCX pandoc -s math.tex -o example30.docx
    
    ```
    
    Deleting 1143:1172
    ```python
    
    '''
    
    def square(lis):
    
        for e in lis:
    
            yield e**2
    
    
    
    xvels_squared = []
    
    for cell_xvels in U:
    
        xvels_squared.append(sum(sum(square(lis)) for lis in cell_xvels))
    
    #print(xvels)
    
    yvels_squared = []
    
    for cell_yvels in V:
    
        yvels_squared.append(sum(sum(square(lis)) for lis in cell_yvels))
    
    assert(len(xvels) == len(yvels))
    
    it = 0
    
    energy_cell_level = []
    
    mass = 2.988e-26
    
    for xvel_squared in xvels_squared:
    
        yvel_squared = yvels_squared[it]
    
        energy_cell_level.append((xvel_squared + yvel_squared)) # might add mass here
    
        
    
    plt.plot(energy_cell_level)
    
    plt.title('Variation in energy, cell method')
    
    plt.xlabel('Timestep')
    
    plt.ylabel('Energy')
    
    plt.savefig('./Assets/constant_energy_cellcalc.png')
    
    plt.close()
    
    '''
    
    ```
    
    Deleting 1088:1127
    ```python
    
    momentum_x = []
    
    momentum_y = []
    
    for it in range(0, timesteps-1):
    
        vx = sum(particles['vx{}'.format(it)])
    
        vy = sum(particles['vy{}'.format(it)])
    
        momentum_x.append(vx)
    
        momentum_y.append(vy)
    
    
    
    base_momentum_x = momentum_x[0]
    
    base_momentum_y = momentum_y[0]
    
    error_x = [base_momentum_x - m for m in momentum_x]
    
    error_y = [base_momentum_y - m for m in momentum_y]
    
    plt.plot(error_x)
    
    plt.title('Variation in x-Velocity, initially: {}'.format(round(base_momentum_x, 2)))
    
    plt.xlabel('Timestep')
    
    plt.ylabel('x-Velocity')
    
    plt.savefig('./Assets/x_velocity_variation.png')
    
    plt.show()
    
    #plt.close()
    
    plt.plot(error_y)
    
    plt.title('Variation in y-Velocity, initally: {}'.format(round(base_momentum_y, 2)))
    
    plt.xlabel('Timestep')
    
    plt.ylabel('y-Velocity')
    
    plt.savefig('./Assets/y_velocity_variation.png')
    
    plt.show()
    
    #plt.close()
    
    plt.plot(error_x, error_y, "o", markersize = 3)
    
    plt.title('Variation in velocity')
    
    plt.xlabel('Variation in x-velocity')
    
    plt.ylabel('Variation in y-velocity')
    
    plt.savefig('./Assets/velocity_variation.png')
    
    plt.show()
    
    #plt.close()
    
    # -----------------------------------
    
    
    
    
    
    #plt.close()
    
    ```
    
    Deleting 1066:1080
    ```python
    
    import matplotlib.pyplot as plt
    
    nums = []
    
    for pivot in pivots:
    
        nums.append(sum(sum(lis) for lis in pivot.values))
    
    
    
    plt.plot(nums)
    
    plt.title('Variation in number of particles with timesteps')
    
    plt.ylabel('Number of Particles')
    
    #plt.ylim((79990, 80010))
    
    plt.xlabel('Timestep')
    
    plt.savefig('./Assets/number_particles.png')
    
    plt.close()
    
    ```
    
    Deleting 1020:1060
    ```python
    
    from matplotlib import animation
    
    
    
    print('Animating Streamplot ...\n')
    
    
    
    # First set up the figure, the axis, and the plot element we want to animate
    
    fig = plt.figure()
    
    ax = plt.axes(xlim=(x_0_region, x_max_region), ylim=(y_0_region, y_max_region))
    
    color = np.sqrt(U[0]**2 + V[0]**2)
    
    stream = ax.streamplot(I, J, U[0], V[0], color=color, density = streamplot_density)
    
    
    
    # initialization function: plot the background of each frame
    
    def init():
    
        #quiv.set_data([], [], [], [])
    
        return stream
    
    
    
    # animation function.  This is called sequentially
    
    def animate(it):
    
        ax.collections = [] # clear lines streamplot
    
        ax.patches = [] # clear arrowheads streamplot
    
        color = np.sqrt(U[it]**2 + V[it]**2)
    
        stream = ax.streamplot(I, J, U[it], V[it], color=color, density = streamplot_density)
    
        if (it % 10 == 0):
    
            print('--Created {} frame.\n'.format(it))
    
        return stream
    
    
    
    # call the animator.  blit=True means only re-draw the parts that have changed.
    
    anim = animation.FuncAnimation(fig, animate, #init_func=init,
    
                                   frames=timesteps, blit=False)
    
    
    
    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    
    # installed.  The extra_args ensure that the x264 codec is used, so that
    
    # the video can be embedded in html5.  You may need to adjust this for
    
    # your system: for more information, see
    
    # http://matplotlib.sourceforge.net/api/animation_api.html
    
    anim.save('./Assets/streamplot_animation.mp4', fps=5) #extra_args=['-vcodec', 'libx264'])
    
    print('Animated and saved!')
    
    
    
    plt.close()
    
    ```
    
    Deleting 978:1018
    ```python
    
    from matplotlib import animation
    
    
    
    print('Animating particles ...')
    
    # First set up the figure, the axis, and the plot element we want to animate
    
    point_size = 0.1
    
    fig = plt.figure()
    
    ax = plt.axes(xlim=(x_0_region, x_max_region), ylim=(0, y_max_region))
    
    scatter, = ax.plot(particles['x0'], particles['y0'], "o", markersize = point_size)
    
    
    
    x = 'x'
    
    y = 'y'
    
    
    
    # initialization function: plot the background of each frame
    
    def init():
    
        #quiv.set_data([], [], [], [])
    
        return scatter,
    
    
    
    # animation function.  This is called sequentially
    
    def animate(it):
    
        scatter.set_xdata(particles[x + str(it)])
    
        scatter.set_ydata(particles[y + str(it)])
    
        if (it % 10 == 0):
    
            print('--Created {} frame.\n'.format(it))
    
        return scatter,
    
    
    
    # call the animator.  blit=True means only re-draw the parts that have changed.
    
    anim = animation.FuncAnimation(fig, animate, #init_func=init,
    
                                   frames=timesteps, blit=True)
    
    
    
    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    
    # installed.  The extra_args ensure that the x264 codec is used, so that
    
    # the video can be embedded in html5.  You may need to adjust this for
    
    # your system: for more information, see
    
    # http://matplotlib.sourceforge.net/api/animation_api.html
    
    anim.save('./Assets/scatterplot_animation.mp4', fps=5, extra_args=['-vcodec', 'libx264'])
    
    print('Animated and saved!')
    
    
    
    plt.close()
    
    ```
    
    Deleting 935:976
    ```python
    
    from matplotlib import animation
    
    import seaborn as sns
    
    
    
    print('Animating Heatmap ...\n')
    
    
    
    # First set up the figure, the axis, and the plot element we want to animate
    
    fig = plt.figure()
    
    ax = plt.axes(xlim=(x_0_region, x_max_region), ylim=(y_0_region, y_max_region))
    
    sns.heatmap(pivots[0], square = True, xticklabels = False, yticklabels = False, cbar = False)#,cbar_kws={"fraction": 0.046, "pad": 0.005})
    
    
    
    # initialization function: plot the background of each frame
    
    def init():
    
        plt.clf()
    
        ax = sns.heatmap(pivots[0], square = True, xticklabels = False, yticklabels = False, cbar = False)#,cbar_kws={"fraction": 0.046, "pad": 0.005})
    
        #return heatmap,
    
    
    
    # animation function.  This is called sequentially
    
    def animate(it):
    
        plt.clf()
    
        ax = sns.heatmap(pivots[it], square = True, xticklabels = False, yticklabels = False, cbar = False)#,cbar_kws={"fraction": 0.046, "pad": 0.005})
    
        ax.set_xlim((x_0_region, x_max_region))
    
        ax.set_ylim((y_0_region, y_max_region))
    
        if (it % 10 == 0):
    
            print('--Created {} frame.\n'.format(it))
    
        #return heatmap
    
    
    
    # call the animator.  blit=True means only re-draw the parts that have changed.
    
    anim = animation.FuncAnimation(fig, animate, init_func=init,
    
                                   frames=timesteps)
    
    
    
    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    
    # installed.  The extra_args ensure that the x264 codec is used, so that
    
    # the video can be embedded in html5.  You may need to adjust this for
    
    # your system: for more information, see
    
    # http://matplotlib.sourceforge.net/api/animation_api.html
    
    anim.save('./Assets/heatmap_animation.mp4', fps=5) #extra_args=['-vcodec', 'libx264'])
    
    print('Animated and saved!')
    
    
    
    plt.close()
    
    ```
    
    Deleting 895:933
    ```python
    
    from matplotlib import animation
    
    print('Animating Quiver ...\n')
    
    
    
    # First set up the figure, the axis, and the plot element we want to animate
    
    fig = plt.figure()
    
    ax = plt.axes(xlim=(x_0_region, x_max_region), ylim=(y_0_region, y_max_region))
    
    
    
    quiv = ax.quiver(I, J, 10**-12*U[0], 10**-12*V[0], scale = 1) #10**-5*U10**-5*V
    
    
    
    # initialization function: plot the background of each frame
    
    def init():
    
        #quiv.set_data([], [], [], [])
    
        return quiv,
    
    
    
    # animation function.  This is called sequentially
    
    def animate(it, quiv, I, J):
    
        #color = np.sqrt(U[it]**2 + V[it]**2)
    
        quiv.set_UVC(10**-12*U[it], 10**-12*V[it]) #*10**-5*
    
        if (it % 10 == 0):
    
            print('--Created {} frame.\n'.format(it))
    
        #quiv.set_color(color)
    
        return quiv,
    
    
    
    # call the animator.  blit=True means only re-draw the parts that have changed.
    
    anim = animation.FuncAnimation(fig, animate, fargs = (quiv, I, J),#init_func=init,
    
                                   frames=timesteps, blit=False)
    
    
    
    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    
    # installed.  The extra_args ensure that the x264 codec is used, so that
    
    # the video can be embedded in html5.  You may need to adjust this for
    
    # your system: for more information, see
    
    # http://matplotlib.sourceforge.net/api/animation_api.html
    
    anim.save('./Assets/quiver_animation.mp4', fps=5) #extra_args=['-vcodec', 'libx264'])
    
    print('Animated and saved!')
    
    
    
    plt.close()
    
    ```
    
    Deleting 857:875
    ```python
    
    v_2_t = [particles[['vx{}'.format(t), 'vy{}'.format(t)]]**2 for t in range(round(len(particles.columns)/4 - 1))]
    
    #print(v_2_t)
    
    avg_v_2_t = [v_2_particles.sum().sum()/len(v_2_particles) for v_2_particles in v_2_t]# first col sum, then row sum (sum vxs, then sum vys)
    
    print('Average sqared velocities after first timestep:\n{}'.format(v_2_t[0].sum()/len(v_2_t[0])))
    
    print('Average v_x2 after last timestep:\n{}'.format(v_2_t[-1].sum()/len(v_2_t[-1])))
    
    
    
    particle_mass = float(constants['particle_mass'])
    
    
    
    temp_t = np.array(avg_v_2_t)
    
    plt.plot(temp_t)
    
    #plot = sns.lineplot(data=temp_t)
    
    plt.xlabel('Timestep')
    
    plt.ylabel('average of v2')
    
    plt.title('Temperature is increasing')
    
    plt.xlim(0,1000)
    
    plt.savefig('Assets/temperature.png')
    
    ```
    
    Deleting 846:853
    ```python
    
    last_vels = particles['vy{}'.format(timesteps)]
    
    sns.distplot(last_vels)
    
    plt.xlabel('vy, last timestep')
    
    plt.ylabel('density')
    
    plt.title('distribution of y velocity at t=last')
    
    ```
    
    Deleting 837:844
    ```python
    
    last_vels = particles['vy{}'.format(0)]
    
    sns.distplot(last_vels)
    
    plt.xlabel('vy, timestep 0')
    
    plt.ylabel('density')
    
    plt.title('distribution of y velocity at t=0')
    
    ```
    
    Deleting 700:712
    ```python
    
    x_velocities = [particles['vx{}'.format(i)] for i in range(timesteps)]
    
    
    
    n_frames = num_timesteps - 1
    
    
    
    fig = plt.figure()
    
    x_velocities_hist = plt.hist(x_velocities[0])
    
    
    
    anim = animation.FuncAnimation(fig, update_hist, n_frames, fargs=(x_velocities,))
    
    anim.save('./Assets/x_vel_hist.mp4', fps = 10)
    
    plt.show()
    
    ```
    
    Deleting 687:696
    ```python
    
    from matplotlib import animation
    
    
    
    def update_hist(num, data):
    
        if (num % 10 == 0):
    
            print('Frame: {}'.format(num))
    
        plt.cla()
    
        plt.hist(data[num], bins = 20)
    
    ```
    
    Deleting 644:669
    ```python
    
    import math
    
    
    
    def theoretical_dist(Delta_v, N_c, av_N_c = 10, m = 1, k_B = 1, T = 1):
    
        A = (m/(2 * math.pi * k_B * T * (1 - 1/N_c)))**(3/2)
    
        B = np.exp(-(m/(2 * k_B * T * (1 - 1/N_c))) * Delta_v**2)
    
        return A * B
    
    
    
    def poisson_avg(Delta_v, max_sum = 1000, av_N_c = 10, m = 1, k_B = 1, T = 1):
    
        e = math.exp(-av_N_c)
    
        total = 0
    
        for N_c in range(2, max_sum):
    
            P = (av_N_c**N_c)/(math.factorial(N_c))
    
            N = theoretical_dist(Delta_v, N_c, m, k_B, T)
    
            total += P * (N)
    
        D = 1 - (av_N_c + 1) * e
    
        return e * total/D
    
    
    
    def plot_theoretical():
    
        Delta_v = np.linspace(0, 4, 1000)
    
        plt.plot(Delta_v, poisson_avg(Delta_v))
    
        plt.show()
    
        
    
    plot_theoretical()
    
    ```
    
    Deleting 623:636
    ```python
    
    y = np.linspace(0, 20, 100)
    
    plt.plot(y, parabolic_flow(y))
    
    plt.title('Parabolic flow')
    
    
    
    rowsums = [row.sum()/len(row) for row in U[-1]]
    
    plt.plot(range(len(rowsums)), rowsums, "o-")
    
    plt.xlabel('x-Velocity')
    
    plt.ylabel('Cell number')
    
    plt.title('Velocity profile, averaged over all rows, for the last timestep ({})'.format(timesteps))
    
    plt.savefig("Assets/velocity_profile.png")
    
    #plt.close()
    
    ```
    
    Deleting 613:621
    ```python
    
    def viscosity(v_max, m = 1, av_N_c = 10, g = 0.01, D = 20):
    
        return (m * av_N_c * g * D**2) / (8 * v_max)
    
    
    
    def parabolic_flow(y, v_max = 4, D = 20):
    
        return (4 * v_max * (D - y) * y)/(D**2)
    
    
    
    ```
    
    Deleting 317:377
    ```python
    
    print('Plotting data ..')
    
    fig = plt.figure(figsize=(20,20))
    
    ax = [fig.add_subplot(2,2,i+1) for i in range(4)]
    
    
    
    timesteps = int(constants['timesteps']) - 1
    
    
    
    
    
    for a in ax:
    
        a.set_xticklabels([])
    
        a.set_yticklabels([])
    
        #a.set_aspect('equal')
    
        
    
    fig.subplots_adjust(wspace=0, hspace=0)
    
    
    
    color = np.sqrt(U[-1]**2 + V[-1]**2)
    
    point_size = 0.3
    
    
    
    if with_obstacles:
    
        for index, o in obstacles.iterrows():
    
            circle1 = plt.Circle((o['x'], o['y']), obstacle_radius, color = 'orange')
    
            circle2 = plt.Circle((o['x'], o['y']), obstacle_radius, color = 'orange')
    
            #circle3 = plt.Circle((o['x'], o['y']), obstacle_radius, color = 'black')
    
            #circles.append(circle)
    
            ax[1].add_artist(circle1)
    
            ax[2].add_artist(circle2)
    
            #ax[3].add_artist(circle3)
    
            #ax[0].xaxis.set_ticks([])
    
    
    
    ax[0].xaxis.set_ticks([])
    
    ax[0].yaxis.set_ticks([])
    
    ax[0].quiver(I, J, U[-1], V[-1], color)
    
    ax[0].set(xlim=(x_0_region-1,x_max_region), ylim=(y_0_region-1,y_max_region))
    
    
    
    ax[1].xaxis.set_ticks([])
    
    ax[1].yaxis.set_ticks([])
    
    ax[1].streamplot(I, J, U[-1], V[-1], color=color, density=streamplot_density) # grid
    
    ax[1].set(xlim=(x_0_region-1,x_max_region), ylim=(y_0_region-1,y_max_region))
    
    
    
    ax[2].xaxis.set_ticks([])
    
    ax[2].yaxis.set_ticks([])
    
    ax[2].plot(particles['x{}'.format(timesteps)], particles['y{}'.format(timesteps)], "o", markersize = point_size)
    
    ax[2].set(xlim=(x_0_region,x_max_region), ylim=(y_0_region,y_max_region))
    
    
    
    ax[3].xaxis.set_ticks([])
    
    ax[3].yaxis.set_ticks([])
    
    #img = ax[3].imshow(pivot, cmap='hot')
    
    #fig.colorbar(img, ax=ax[3], fraction=0.046, pad=0.005)
    
    sns.heatmap(pivots[-1], ax=ax[3], xticklabels = False, yticklabels = False, cbar_kws={"fraction": 0.046, "pad": 0.005})
    
    ax[3].set(xlim=(x_0_region,x_max_region), ylim=(y_0_region,y_max_region))
    
    #ax[3].xticks('')
    
    #ax[3].yticks('')
    
    ax[3].set_ylabel('')
    
    ax[3].set_xlabel('')
    
    #ax[1,1].imshow(pivot, cmap='hot')
    
    
    
    plt.savefig("Assets/stationary_region.png")
    
    #plt.close()
    
    print('Data plotted and saved!')
    
    ```
    
    Deleting 238:305
    ```python
    
    # Plotting
    
    #with (sns.plotting_context(sns.set())):
    
    x_0_region = 0
    
    x_max_region = int(constants['width'])
    
    y_0_region = 0
    
    y_max_region = int(constants['height'])
    
    
    
    streamplot_density = [0.5, 1]
    
    
    
    print('Plotting data ..')
    
    fig = plt.figure(figsize=(20,20))
    
    ax = [fig.add_subplot(2,2,i+1) for i in range(4)]
    
    
    
    for a in ax:
    
        a.set_xticklabels([])
    
        a.set_yticklabels([])
    
        #a.set_aspect('equal')
    
        
    
    fig.subplots_adjust(wspace=0, hspace=0)
    
    
    
    color = np.sqrt(U[0]**2 + V[0]**2)
    
    point_size = 0.3
    
    
    
    #circles = []
    
    if with_obstacles:
    
        for index, o in obstacles.iterrows():
    
            circle1 = plt.Circle((o['x'], o['y']), obstacle_radius, color = 'gray')
    
            circle2 = plt.Circle((o['x'], o['y']), obstacle_radius, color = 'gray')
    
            #circle3 = plt.Circle((o['x'], o['y']), obstacle_radius, color = 'black')
    
            #circles.append(circle)
    
            ax[1].add_artist(circle1)
    
            ax[2].add_artist(circle2)
    
            #ax[3].add_artist(circle3)
    
            #ax[0].xaxis.set_ticks([])
    
    #ax[0].yaxis.set_ticks([])
    
    ax[0].quiver(I, J, U[0], V[0], color)
    
    ax[0].set(xlim=(x_0_region-1,x_max_region), ylim=(y_0_region-1,y_max_region))
    
    
    
    #ax[1].xaxis.set_ticks([])
    
    #ax[1].yaxis.set_ticks([])
    
    ax[1].streamplot(I, J, U[0], V[0], color=color, density = streamplot_density) # grid
    
    ax[1].set(xlim=(x_0_region-1,x_max_region), ylim=(y_0_region-1,y_max_region))
    
    
    
    #ax[2].xaxis.set_ticks([])
    
    #ax[2].yaxis.set_ticks([])
    
    ax[2].plot(particles['x0'], particles['y0'], "o", markersize = point_size)
    
    if with_obstacles:
    
        sns.scatterplot(ax = ax[2], x = obstacles['x'], y = obstacles['y'], s = obstacle_radius)
    
    ax[2].set(xlim=(x_0_region, x_max_region), ylim=(y_0_region,y_max_region))
    
    
    
    #ax[3].xaxis.set_ticks([])
    
    #ax[3].yaxis.set_ticks([])
    
    #img = ax[3].imshow(pivot, cmap='hot')
    
    #fig.colorbar(img, ax=ax[3], fraction=0.046, pad=0.005)
    
    sns.heatmap(pivots[0], ax=ax[3], xticklabels = False, yticklabels = False, cbar_kws={"fraction": 0.046, "pad": 0.01})
    
    ax[3].set(xlim=(x_0_region-1,x_max_region), ylim=(y_0_region,y_max_region))
    
    #ax[3].xticks('')
    
    #ax[3].yticks('')
    
    ax[3].set_ylabel('')
    
    ax[3].set_xlabel('')
    
    #ax[1,1].imshow(pivot, cmap='hot')
    
    
    
    plt.savefig("Assets/initial_region.png")
    
    #plt.close()
    
    print('Data plotted and saved!')
    
    ```
    
    Deleting 127:215
    ```python
    
    cell_dim = float(constants['cell_dim'])
    
    width = float(constants['width'])
    
    height = float(constants['height'])
    
    shown_cols = floor(width / cell_dim)
    
    shown_rows = floor(height / cell_dim)
    
    
    
    saved = './Saved'
    
    saved_cells = '/cells.pkl'
    
    file_cells = saved + saved_cells
    
    
    
    columns = ['i', 'j', 'meanX', 'meanY', 'num']
    
    
    
    def load_cells(path, columns, shown_rows, shown_cols):
    
        
    
        index = ['t']
    
        if (os.path.isfile(file_cells)):
    
            cells = pd.read_pickle(file_cells)
    
            #print(cells)
    
            #cells.set_index(index, inplace=True)
    
        else:
    
            # Loading cells
    
            print('Loading cells')
    
            
    
            cells_timesteps = []
    
            
    
            it = 0
    
            filenames_cells = glob.glob('{}*.csv'.format(path))
    
            for file in filenames_cells:
    
                df = pd.read_csv(file)
    
                cells_timesteps.append(df)
    
                #df[[i,j]] = (df[[i,j]] + 1/2) * cell_dim
    
                #print(df.head())
    
                it += 1
    
                if (it % 10 == 0):
    
                    print('--loaded {}'.format(it))
    
            #cells.to_pickle(file_cells)
    
            print('Cells loaded and saved!\n')
    
            # Cells loaded
    
            return cells_timesteps
    
        
    
    def prepare_cells(cells_timesteps, columns, shown_rows, shown_cols):
    
        # Preparing cell values
    
        print('Preparing cell values ..')
    
    
    
        array_i = np.arange(0, shown_rows)
    
        array_j = np.arange(0, shown_cols)
    
        I,J = np.meshgrid(array_j, array_i)
    
    
    
        U = []
    
        V = []
    
    
    
        i = columns[0]
    
        j = columns[1]
    
        vx = columns[2]
    
        vy = columns[3]
    
        num = columns[4]
    
    
    
        pivots = []
    
        for df in cells_timesteps:
    
            # only the rows and cols above 0
    
            # and below shown_rows, shown_cols
    
            # this is to 
    
            # --1. no vaccuum around simulated region
    
            # --2. I,J are fixed size
    
            U_inner = []
    
            V_inner = []
    
            for it in array_i: # TODO: check this code something seems foul (row, cols, but only using rows)
    
                temp = df.loc[df[i] == it]
    
                u = np.array(temp[vx])
    
                U_inner.append(u)
    
                v = np.array(temp[vy])
    
                V_inner.append(v)
    
            U.append(np.array(U_inner))#, dtype = object))
    
            V.append(np.array(V_inner))#, dtype = object))
    
    
    
            pivot = df.pivot(index = i, columns = j, values = num)
    
            pivots.append(pivot)
    
    
    
        print('Cell preparation complete!')
    
        return I,J,U,V, pivots
    
        # Cell preparation complete
    
    
    
    
    
    cells_path = path + 'cells_av{}'.format(av)
    
    cells_timesteps = load_cells(cells_path, columns, shown_rows, shown_cols)
    
    I,J,U,V,pivots = prepare_cells(cells_timesteps, columns, shown_rows, shown_cols)
    
    ```
    
    Deleting 66:111
    ```python
    
    shown_x = float(constants['width'])
    
    shown_y = float(constants['height'])
    
    
    
    saved = './Saved'
    
    saved_particles = '/particles.pkl'
    
    file_particles = saved + saved_particles
    
    
    
    def load_particles(path):
    
        columns = []
    
        x_columns = ['x{}'.format(it) for it in range(0, num_timesteps)]
    
        vx_columns = ['vx{}'.format(it) for it in range(0, num_timesteps)]
    
        y_columns = ['y{}'.format(it) for it in range(0, num_timesteps)]
    
        vy_columns = ['vy{}'.format(it) for it in range(0, num_timesteps)]
    
        columns.extend(x_columns)
    
        columns.extend(y_columns)
    
        columns.extend(vx_columns)
    
        columns.extend(vy_columns)
    
        particles = pd.DataFrame(columns = columns)
    
    
    
        if (os.path.isfile(file_particles) and not new_analysis):
    
            print('Found saved particles_x and particles_y files!')
    
            particles = pd.read_pickle(file_particles)
    
            print('Loaded particles files.')
    
        else:
    
            # Loading particles
    
            print('Loading particles ..')
    
            filenames_particles = glob.glob('{}*.csv'.format(path))
    
            it = 0
    
            for file in filenames_particles:
    
                df = pd.read_csv(file)
    
                particles[['x{}'.format(it), 'y{}'.format(it), 'vx{}'.format(it), 'vy{}'.format(it)]] = df[['x', 'y', 'vx', 'vy']]
    
                it += 1
    
                if (it % 10 == 0):
    
                    print('--loaded {}'.format(it))
    
            particles.to_pickle(file_particles)
    
            it = 0
    
            print('Particles loaded and saved!\n')
    
            # Particles loaded
    
        return particles
    
    
    
    
    
    particles_path = '{parent}particles_av{av}'.format(parent = path, av=av)    
    
    particles = load_particles(path = particles_path)
    
    ```
    
    Deleting 56:64
    ```python
    
    obstacles = None
    
    if (with_obstacles):
    
        obstacles_path = path + 'constants_obstacles' + csv
    
        obstacles_csv = pd.read_csv(obstacles_path)
    
        obstacles = obstacles_csv[['x', 'y']]
    
        obstacle_radius = float(obstacles_csv['r'][0])
    
    ```
    
    Deleting 44:54
    ```python
    
    
    
    path = "Data/"
    
    csv = ".csv"
    
    av = 10
    
    constants = 'constants_av' + str(av) + csv
    
    
    
    constants = pd.read_csv(path + constants)
    
    num_timesteps = int(constants['timesteps'])
    
    ```
    
    Deleting 34:42
    ```python
    
    timesteps = []
    
    I = []
    
    J = []
    
    U = []
    
    V = []
    
    pivot = []
    
    ```
    
    Deleting 28:32
    ```python
    
    new_analysis = True
    
    with_obstacles = False
    
    ```
    
    Deleting 13:26
    ```python
    
    import numpy as np
    
    import matplotlib.pyplot as plt
    
    from math import floor
    
    import pandas as pd
    
    import glob
    
    import numpy as np
    
    import matplotlib.pyplot as plt
    
    import matplotlib.gridspec as gridspec
    
    import seaborn as sns
    
    import os.path
    
    
    
    ```
    
    Deleting 0:5
    ```latex
    
    %%latex
    
    \newcommand{\matr}[1]\textbf{#1}
    
    \newcommand{\vect}[1]{\vec{#1}}
    
    ```
    
    




    CompletedProcess(args='pandoc -s ./Generated/thesis.md -o ./Generated/thesis.pdf --filter pandoc-citeproc --bibliography="list.bib" --csl="apa.csl"', returncode=83)



## Equation Numbering jupyter extension
conda install -c conda-forge jupyter_contrib_nbextensions

jupyter contrib nbextension install --user

jupyter nbextension enable equation-numbering/main

### Turn equation numbering on/off


```javascript
%%javascript
MathJax.Hub.Config({
    TeX: { equationNumbers: { autoNumber: "AMS" } }
});
```

### Renumber equations


```javascript
%%javascript
MathJax.Hub.Queue(
  ["resetEquationNumbers", MathJax.InputJax.TeX],
  ["PreProcess", MathJax.Hub],
  ["Reprocess", MathJax.Hub]
);
```

-->


```python

```
