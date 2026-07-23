import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, brentq


st = 0
en = 100
N = 200
c = 1
m = 1
q = 1

E0 = 2  # Electric field strength

x = np.linspace(st, en, N)
y = np.linspace(st, en, N)


def uniform_electric_field(x, y, E0):
    # Uniform electric field directed downward
    Ex = np.zeros((N, N))
    Ey = np.zeros((N, N))

    for i in range(len(x)):
        for j in range(len(y)):
            Ex[i, j] = 0
            Ey[i, j] = -E0

    E_strength = np.sqrt(Ex**2 + Ey**2)

    return Ex, Ey, E_strength


Ex, Ey, E_strength = uniform_electric_field(x, y, E0)


def optimize_smooth_path_uniform_E_relativistic(
    x, y, E0, c, in_point, fin_point, n_points=80
):
    xA, yA = x[in_point[0]], y[in_point[1]]
    xB, yB = x[fin_point[0]], y[fin_point[1]]

    # The optimizer changes only the y coordinates of the path
    x_path = np.linspace(xA, xB, n_points)

    s = np.linspace(0, 1, n_points)

    # Initial path is a straight line bent downward by a sine function
    y_guess = yA + (yB - yA) * s
    y_guess = y_guess - 0.4 * abs(yA - yB) * np.sin(np.pi * s)

    y_inner_guess = y_guess[1:-1]

    smooth_lambda = 0.05

    def travel_time(y_inner):
        y_path = np.concatenate(([yA], y_inner, [yB]))

        T = 0

        for k in range(n_points - 1):
            y_mid = 0.5 * (y_path[k] + y_path[k + 1])

            # Energy gained by moving downward in the electric field
            energy_gain = q * (-E0) * (y_mid - yA)

            if energy_gain <= 0:
                return 1e12

            # Relativistic energy determines gamma and therefore the velocity
            gamma = 1 + energy_gain / (m * c**2)
            v_mid = c * np.sqrt(1 - 1 / gamma**2)

            ds = np.sqrt(
                (x_path[k + 1] - x_path[k])**2
                + (y_path[k + 1] - y_path[k])**2
            )

            T += ds / v_mid

        # Penalize sharp and non-smooth changes in the path
        smooth_penalty = 0

        for k in range(1, n_points - 1):
            smooth_penalty += (
                y_path[k + 1] - 2 * y_path[k] + y_path[k - 1]
            )**2

        return T + smooth_lambda * smooth_penalty

    bounds = [(y[0], y[-1]) for _ in range(n_points - 2)]

    # Find the path which minimizes the total travel time
    result = minimize(
        travel_time,
        y_inner_guess,
        method="L-BFGS-B",
        bounds=bounds
    )

    y_opt = np.concatenate(([yA], result.x, [yB]))
    path = list(zip(x_path, y_opt))

    return path, result.fun, result


def exact_cycloid_between_points(A, B, n=500):
    xA, yA = A
    xB, yB = B

    L = xB - xA
    H = yA - yB

    if L <= 0:
        raise ValueError("Final point must be to the right of initial point.")

    if H <= 0:
        raise ValueError("Final point must be below initial point.")

    # Find the cycloid parameter which connects the selected endpoints
    def equation(theta):
        return H / L - (1 - np.cos(theta)) / (theta - np.sin(theta))

    theta_final = brentq(
        equation,
        1e-6,
        2 * np.pi - 1e-6
    )

    a = L / (theta_final - np.sin(theta_final))

    theta = np.linspace(0, theta_final, n)

    # Parametric equations of the exact classical cycloid
    X = a * (theta - np.sin(theta))
    Y = a * (1 - np.cos(theta))

    x_cyc = xA + X
    y_cyc = yA - Y

    return x_cyc, y_cyc


in_point = (20, 160)
fin_point = (180, 80)


path, T, result = optimize_smooth_path_uniform_E_relativistic(
    x, y, E0, c, in_point, fin_point, n_points=80
)

path = np.array(path)

A = (x[in_point[0]], y[in_point[1]])
B = (x[fin_point[0]], y[fin_point[1]])

# Classical cycloid is calculated only for comparison
x_cyc, y_cyc = exact_cycloid_between_points(A, B)


plt.figure(figsize=(8, 6))

plt.plot(path[:, 0], path[:, 1], "o-", label="Relativistic numerical path")
plt.plot(x_cyc, y_cyc, "-", label="Exact classical cycloid")

plt.scatter(A[0], A[1], s=80, label="Initial point")
plt.scatter(B[0], B[1], s=80, label="Final point")

plt.xlabel("x")
plt.ylabel("y")
plt.title("Relativistic numerical path vs classical cycloid")
plt.grid()
plt.legend()
plt.axis("equal")
plt.show()

print("Relativistic travel time:", T)
print("Optimization success:", result.success)
