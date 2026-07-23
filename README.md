# Numerical Method for Solving the Classical and Relativistic Electrobrachistochrone


## Introduction
The brachistochrone problem is one of the most well-known optimization problems in mathematics and physics. The original problem asks for the path that minimizes the travel time of a particle moving in a uniform gravitational field between two fixed points. In this work, I consider an electric charge moving in an electric field and ask the same question: what path minimizes the travel time? I developed a numerical method in Python that solves the problem for a classical particle in both uniform and non-uniform electric fields, and for a relativistic particle in the uniform field case.



## Classical Electrobrachistochrone

We consider a free-space region with an electric field. There is a particle with an initial velocity of zero, $v_0 = 0$, and charge $q = 1$ without loss of generality. Due to the influence of the field, it starts to move from the initial point $A$ to the final point $B$ along a wire. We aim to find the path that minimizes the travel time for different kinds of electric potentials.

### Numerical Implementation
We start by discretizing the plane into \(N\) points:

$$
x\in[0,N],\qquad
y\in[0,N],\qquad
(x_i,y_j).
$$

Once we have the grid points, we define the electric potential as an \(N \times N\) matrix whose entries are the values of the electric potential,

$$
\phi_{ij}=\phi(x_i,y_j).
$$

During the optimization, the trajectory does not necessarily pass through the grid points. Instead, for every point \((x,y)\) on the trajectory, the surrounding grid cell is identified, consisting of the four neighboring points

$$
(x_i,y_j),\;
(x_{i+1},y_j),\;
(x_i,y_{j+1}),\;
(x_{i+1},y_{j+1}).
$$

The relative position of the point inside the cell is then calculated as

$$
t_x=\frac{x-x_i}{x_{i+1}-x_i},
\qquad
t_y=\frac{y-y_j}{y_{j+1}-y_j}.
$$

Using these parameters, the electric potential at the point \((x,y)\) is obtained by bilinear interpolation of the potential values at the four surrounding grid points:

$$
\begin{aligned}
\phi(x,y)=
&(1-t_x)(1-t_y)\phi_{i,j}
+t_x(1-t_y)\phi_{i+1,j}\\
&+(1-t_x)t_y\phi_{i,j+1}
+t_xt_y\phi_{i+1,j+1}.
\end{aligned}
$$

We calculate the velocity of the particle from energy conservation:

$$
v(x,y)=\sqrt{\frac{2q}{m}\left(\phi_A-\phi(x,y)\right)}
$$

The length of each segment is calculated as

$$
\Delta S_k=\sqrt{(x_{k+1}-x_k)^2+(y_{k+1}-y_k)^2}.
$$

The time required to cover each segment is

$$
\Delta T_k=
\frac{\Delta S_k}
{v(x_k^{\mathrm{mid}},y_k^{\mathrm{mid}})},
\qquad
x_k^{\mathrm{mid}}=\frac{x_k+x_{k+1}}{2},
\qquad
y_k^{\mathrm{mid}}=\frac{y_k+y_{k+1}}{2}.
$$

The total travel time along the path is

$$
T\approx
\sum_{k=0}^{n-1}
\frac{\Delta S_k}
{v(x_k^{\mathrm{mid}},y_k^{\mathrm{mid}})}.
$$

## Optimization Algorithm

The type of algorithm we use to find the optimal path does not differentiate between physically allowed and physically forbidden paths. For this reason, we first need to introduce restrictions that guarantee physically realistic results.

Because of energy conservation, a particle cannot appear at points that correspond to higher potential energy than the initial point \(A\):

$$
q(\phi_A-\phi(x,y))\ge0.
$$

If the motion ends at point \(B\),

$$
q(\phi_A-\phi_B)\ge0.
$$

For the algorithm to run more efficiently, we fix the direction and step size of the particle along the \(x\)-axis. This limits the class of potentials that can be considered, but it does not reduce the generality of the method because the restriction is introduced only for computational reasons. If we want to calculate the optimal path in a potential where the particle moves from right to left, a simple modification of the code is sufficient.

$$
x_k=x_A+\frac{k}{n}(x_B-x_A)
$$

where \(n\) is the number of points satisfying the restriction.

Initially, we guess the solution, for example, a straight line,

$$
y_{\mathrm{line}}(s)=y_A+(y_B-y_A)s,
\qquad
0\le s\le1,
$$

and then deform it according to:

$$
y_{\mathrm{guess}}=
y_A+\left(y_B-y_A\right)s
-0.4|y_A-y_B|\sin(\pi s).
$$

The L-BFGS-B algorithm then calculates the travel time along each path and selects the one that yields the smallest value.

When I ran the algorithm, I encountered some physically unreasonable results because the paths were non-smooth and had a zig-zag shape. Therefore, I introduced an additional restriction that guarantees smoothness:

$$
S=\sum_{k=1}^{n-1}
\left(y_{k+1}-2y_k+y_{k-1}\right)^2.
$$

The final objective function to be minimized is

$$
J=T+\lambda S,
$$

where \(\lambda\) determines how strongly the smoothness requirement is enforced.



## Relativistic Electrobrachistochrone

For the relativistic case we only consider moving charged particle in unfirom electric field, pointing downwards: $vec{E} = (0,-E)$. So, in the code I don't define the electric potential, I just define the electrical field strength and construct the suitable paramets for optimization based on that. 
Relativistic kinetik energy: $ K = (\gamma - 1)mc^2, \qquad \gamma = frac{1}{\sqrt{1-frac{v^2}{c^2}}}$.
At every point transfered energy is: $ K=qE_0(y_A-y) $
based on this we can express the $\gamma$ factor becomes: $\gamma=1-frac{qE_0(y_A-y)}{mc^2}$
relativistic velocity can be expressed as: $v_{\mathrm{rel}}(y)= c\times\sqrt{1 - frac{1}{(1+frac{qE_0(y_A-y)}{mc^2})^2}}$
from these definitions in order to control the regime we use the factor: $ \eta = frac{qE_0\Delta y}{mc^2}$
if $\eta \ll 1 $ then we have the non-relativstic case, if $\eta \approx 1$ we have quasyrelativistic case, and for $\eta \gg 1$ we will have relativstic regime. 
Algorithm used for optimization is almost same as in calssical case, only changes made are the parameters what we defined differently. 
Physical intuition tells us that in relativstic case straight path should be straight line. Because when velocities are reltaivstic there is no need of initial decline because increment of velocity becomes negligable so particle should choose the correct path. There exists also articles where analytic derivation of the same situation shows that path is straight line in relativistic regime. 










