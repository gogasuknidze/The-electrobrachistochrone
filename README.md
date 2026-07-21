# Numerical Method for Solving the Classical and Relativistic Electrobrachistochrone


## Introduction
The brachistochrone problem is one of the most well-known optimization problems in mathematics and physics. The original problem asks for the path that minimizes the travel time of a particle moving in a uniform gravitational field between two fixed points. In this work, I consider an electric charge moving in an electric field and ask the same question: what path minimizes the travel time? I developed a numerical method in Python that solves the problem for a classical particle in both uniform and non-uniform electric fields, and for a relativistic particle in the uniform field case.



## Classical Electrobrachistochrone

We consider a free-space region with an electric field. There is a particle with an initial velocity of zero, $v_0 = 0$, and charge $q = 1$ without loss of generality. Due to the influence of the field, it starts to move from the initial point $A$ to the final point $B$ along a wire. We aim to find the path that minimizes the travel time for different kinds of electric potentials.

### Numerical Implementation

We start by discretizing the plane into $N$ points: $x \in [0,N], \qquad y\in [0,N], \qquad (x_i,y_i)$
Once we have the points, we can define the electric potential as $N \times N$ matrix whose entries will be the values of the potential energy:          
$$ \phi_{ij} = \phi(x_i,y_j) $$

We will see that when we run the algorithm for optimization, knowing the values of the electric potential on the grid points won't be enough, and we will need the values inside the cells also, so for that reason I use the weighted average. The idea is simple because we know the field value and cell nodes, and a point with coordinates $(x_0,y_0)$ is somewhere inside the cell we define $ t_x =  
