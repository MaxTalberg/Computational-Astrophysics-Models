#
#
# Coursework A
# Question 2
# ____b_____
#
import numpy as np
import matplotlib.pyplot as plt

# constants
G = 6.6743e-11
m1M = 1e-3
m2M = 4e-2
M = 1.9884e+30
AU = 1.496e+11
c = G*M

P1 = (4*(np.pi**2)*(2.52*1.54960e11)**3)/(c*(m1M+1))    #sqrd
P2 = (4*(np.pi**2)*(5.24*1.54960e11)**3)/(c*(m2M+1))    #sqrd

# P1 = 133026545.51
# P2 = 391322484.047

# conditions
a = 0; b = 391322484.047 * 5
N = 1e4
h = (b-a)/N
r1 = np.asarray([2.52*AU, 0])
r2 = np.asarray([5.24*AU, 0])
v1 = (c/(2.52*AU))**0.5
v2 = (c/(5.24*AU))**0.5
v1 = np.asarray([0, v1])
v2 = np.asarray([0, v2])

tpoints = np.arange(a, b, h)    # time points
y0 = np.array([r1, r2, v1, v2])   # ODE values

# ODE equations

def orbits(u):
    r1, r2, v1, v2 = u
    r1_norm = np.linalg.norm(r1)
    r2_norm = np.linalg.norm(r2)
    r21_norm = np.linalg.norm(r1-r2)
    f1 = (-c/r1_norm**3)*r1 + ((c*m2M)/r21_norm**3)*(r1-r2)
    f2 = (-c/r2_norm**3)*r2 - ((c*m1M)/r21_norm**3)*(r1-r2)
    return np.array([v1, v2, f1, f2])

# RK4 algorithm


def RK4step(f, u, h):
    k1 = h*f(u)
    k2 = h*f(u+0.5*k1)
    k3 = h*f(u+0.5*k2)
    k4 = h*f(u+k3)
    return u + (k1+2*k2+2*k3+k4)/6

# RK4 integration


def RK4integrate(f, y0, tspan):
    y = np.zeros([len(tpoints), len(y0), len(y0[0])])
    y[0, :] = y0
    for k in range(1, len(tspan)):
        y[k, :] = RK4step(f, y[k-1], tspan[k]-tspan[k-1])
    return y


# Run

sol_RK4 = RK4integrate(orbits, y0, tpoints)  # solve RK4
X, Y = sol_RK4.T    # transpose matrix


r1x, r2x, v1x, v2x = X
r1y, r2y, v1y, v2y = Y

# Analysis

v1 = np.hypot(v1x, v1y)
v2 = np.hypot(v2x, v2y)
r1 = np.hypot(r1x, r1y)
r2 = np.hypot(r2x, r2y)

# plotting results
# FIG 7 and 8 varying b
plt.scatter(r1x, r1y, c=v1, s=5, cmap='turbo')
plt.scatter(r2x, r2y, c=v2, s=5, cmap='turbo')
plt.colorbar(label='$v$ [m/s]', orientation='vertical')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('RK4: 3-body problem')
plt.scatter(0, 0, marker="o", c = 'black')
#plt.show()

# FIG9
plt.plot(tpoints, r1-2.52*AU)
plt.plot(tpoints, r2-5.24*AU)
plt.xlabel('Time [s]')
plt.ylabel('Displacment [m]')
plt.title('RK4: 3-body problem')
plt.legend(['planet 1', 'planet 2'])
#plt.show()
