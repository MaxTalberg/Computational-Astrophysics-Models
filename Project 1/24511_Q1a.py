#
#
# Coursework A
# Question 1
# ____a_____
#
import numpy as np
import matplotlib.pyplot as plt

# constants
G = 6.6743e-11
M = 1.9884e+30
c = G*M

# conditions
a = 0; b = 2.5e9
N = 1e4
h = (b-a)/N
x = 5.2e12; y = 0
vx = 0; vy = 880

tpoints = np.arange(a, b, h)    # time points
y0 = np.array([x, y, vx, vy])   # ODE values

# ODE equations


def orbits(u):
    x, y, vx, vy = u
    r = np.hypot(x, y)
    f = -c/r**3
    return np.array([vx, vy, f*x, f*y])

# RK4 algorithm


def RK4step(f, u, h):
    k1 = h*f(u)
    k2 = h*f(u+0.5*k1)
    k3 = h*f(u+0.5*k2)
    k4 = h*f(u+k3)
    return u + (k1+2*k2+2*k3+k4)/6

# RK4 integration

def RK4integrate(f, y0, tspan):
    z = np.zeros([len(tpoints), len(y0)])
    z[0, :] = y0
    for k in range(1, len(tspan)):
        z[k, :] = RK4step(f, z[k-1], tspan[k]-tspan[k-1])
    return z

# Run

sol_RK4 = RK4integrate(orbits, y0, tpoints)  # solve RK4
x, y, vx, vy = sol_RK4.T    # transpose matrix

v = np.hypot(vx, vy)
r = np.hypot(x, y)

# Calculate aphelion

#print(max(x), min(x))

# plotting results

# FIG2 and FIG3 (b*5)
plt.scatter(x, y, c=v, s=1, cmap='turbo')
plt.colorbar(label='$v$ [m/s]', orientation='vertical')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('RK4: Halleys Comet')
plt.scatter(0, 0, marker="o", c = 'black')      # Sun
#plt.show()
