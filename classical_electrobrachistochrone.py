import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, root_scalar


# Define the coordinate grid and physical parameters
st = 0
en = 100
N = 200
m = 1
q = 1

x = np.linspace(st, en, N)
y = np.linspace(st, en, N)


def polynomial_potential(x, y,
                         ax=0, ay=0,
                         bx=0, by=0, bxy=0,
                         cx=0, cy=0, cx2y=0, cxy2=0):

    # Construct the potential on every point of the two-dimensional grid
    phi = np.zeros((N, N))

    for i in range(len(x)):
        for j in range(len(y)):
            phi[i, j] = (
                ax * x[i] + ay * y[j]
                + bx * x[i]**2 + by * y[j]**2 + bxy * x[i] * y[j]
                + cx * x[i]**3 + cy * y[j]**3
                + cx2y * x[i]**2 * y[j]
                + cxy2 * x[i] * y[j]**2
            )

    return phi


# Different examples of the potential field
linear_phi = polynomial_potential(x, y, ay=2)
quadratic_phi = polynomial_potential(x, y, bx=-2, by=0.01)
cubic_phi = polynomial_potential(x, y, ax=-100, ay=0.001)


def check_endpoints(A, B, phi):

    # Check if the particle can gain kinetic energy while moving from A to B
    if q * (phi[A[0], A[1]] - phi[B[0], B[1]]) >= 0:
        print("Allowed: motion between these two points can happen")
    else:
        print("Denied: motion can not happen between these two points")


def phi_interpolate(x0, y0, x, y, phi):

    # Find the grid cell which contains the point (x0, y0)
    i = np.searchsorted(x, x0) - 1
    j = np.searchsorted(y, y0) - 1

    # Keep the indices inside the grid
    i = max(0, min(i, N - 2))
    j = max(0, min(j, N - 2))

    x1, x2 = x[i], x[i + 1]
    y1, y2 = y[j], y[j + 1]

    # Potential values at the four corners of the grid cell
    Q11 = phi[i, j]
    Q21 = phi[i + 1, j]
    Q12 = phi[i, j + 1]
    Q22 = phi[i + 1, j + 1]

    tx = (x0 - x1) / (x2 - x1)
    ty = (y0 - y1) / (y2 - y1)

    # Bilinear interpolation gives the potential between the grid points
    return (
        Q11 * (1 - tx) * (1 - ty)
        + Q21 * tx * (1 - ty)
        + Q12 * (1 - tx) * ty
        + Q22 * tx * ty
    )


def optimize_smooth_path(x, y, phi, in_point, fin_point, n_points=80):

    # Convert the endpoint indices into their physical coordinates
    xA, yA = x[in_point[0]], y[in_point[1]]
    xB, yB = x[fin_point[0]], y[fin_point[1]]

    # The x coordinates stay fixed, while the optimizer changes y coordinates
    x_path = np.linspace(xA, xB, n_points)

    phi_A = phi[in_point[0], in_point[1]]

    s = np.linspace(0, 1, n_points)

    # Start from a straight line and bend it downward using a sine function
    y_guess = yA + (yB - yA) * s
    y_guess = y_guess - 0.4 * abs(yA - yB) * np.sin(np.pi * s)

    # The endpoints are fixed, so only the inner points are optimized
    y_inner_guess = y_guess[1:-1]

    # Controls how strongly sharp bending of the path is penalized
    smooth_lambda = 0.05

    def travel_time(y_inner):

        # Reconstruct the full path by adding the fixed endpoints
        y_path = np.concatenate(([yA], y_inner, [yB]))

        T = 0

        for k in range(n_points - 1):

            # Use the midpoint of every small path segment
            x_mid = 0.5 * (x_path[k] + x_path[k + 1])
            y_mid = 0.5 * (y_path[k] + y_path[k + 1])

            phi_mid = phi_interpolate(x_mid, y_mid, x, y, phi)

            # Energy conservation: kinetic energy comes from potential-energy loss
            kinetic = q * (phi_A - phi_mid)

            # Such a path is forbidden if the kinetic energy becomes non-positive
            if kinetic <= 0:
                return 1e12

            v_mid = np.sqrt(2 * kinetic / m)

            # Length of the current small path segment
            ds = np.sqrt(
                (x_path[k + 1] - x_path[k])**2
                + (y_path[k + 1] - y_path[k])**2
            )

            # Add the travelling time of the segment
            T += ds / v_mid

        smooth_penalty = 0

        # Penalize large second differences, which correspond to sharp bending
        for k in range(1, n_points - 1):
            smooth_penalty += (
                y_path[k + 1] - 2 * y_path[k] + y_path[k - 1]
            )**2

        # Minimize travelling time while keeping the curve reasonably smooth
        return T + smooth_lambda * smooth_penalty

    # Every inner point must stay inside the coordinate grid
    bounds = [(y[0], y[-1]) for _ in range(n_points - 2)]

    result = minimize(
        travel_time,
        y_inner_guess,
        method="L-BFGS-B",
        bounds=bounds
    )

    y_opt = np.concatenate(([yA], result.x, [yB]))
    path = list(zip(x_path, y_opt))

    return path, result.fun, result


def cycloid_between_points(xA, yA, xB, yB, n=300):

    # Horizontal and vertical distances between the endpoints
    dx = xB - xA
    dy = yA - yB

    if dy <= 0:
        raise ValueError(
            "For the standard brachistochrone, the final point must be lower than the initial point."
        )

    ratio = dx / dy

    # Equation which determines the final cycloid parameter theta_B
    def equation(theta):
        return (theta - np.sin(theta)) / (1 - np.cos(theta)) - ratio

    sol = root_scalar(
        equation,
        bracket=[1e-6, 2 * np.pi - 1e-6]
    )

    theta_B = sol.root

    # Determine the cycloid radius from the endpoint positions
    R = dy / (1 - np.cos(theta_B))

    theta = np.linspace(0, theta_B, n)

    # Parametric equations of the exact cycloid
    x_cyc = xA + R * (theta - np.sin(theta))
    y_cyc = yA - R * (1 - np.cos(theta))

    return x_cyc, y_cyc, R, theta_B


# Endpoint indices on the coordinate grid
A = (20, 160)
B = (160, 100)

check_endpoints(A, B, linear_phi)

# Find the numerical path with the minimum travelling time
path, T, result = optimize_smooth_path(x, y, linear_phi, A, B)

path_x = [p[0] for p in path]
path_y = [p[1] for p in path]

xA, yA = x[A[0]], y[A[1]]
xB, yB = x[B[0]], y[B[1]]

# Calculate the exact cycloid for comparison with the numerical result
cyc_x, cyc_y, R, theta_B = cycloid_between_points(
    xA, yA, xB, yB
)

X, Y = np.meshgrid(x, y, indexing="ij")

# Plot the potential field, numerical path and exact cycloid
plt.figure(figsize=(7, 7))

plt.contourf(X, Y, linear_phi, levels=40, alpha=0.5)
plt.colorbar(label="Potential φ(x,y)")

plt.plot(
    path_x,
    path_y,
    color="black",
    label="Numerical optimized path",
    linewidth=2
)

plt.plot(
    cyc_x,
    cyc_y,
    "--",
    label="Exact cycloid brachistochrone",
    linewidth=2
)

plt.scatter(
    [xA, xB],
    [yA, yB],
    label="Endpoints",
    s=70
)

plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
plt.axis("equal")
plt.legend()
plt.show()
