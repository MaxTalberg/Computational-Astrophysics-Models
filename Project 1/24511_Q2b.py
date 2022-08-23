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
m1 = 5.683e26
m2 = 1.898e27
m3 = 1.9884e+30
AU = 1.496e+11
M_sun = 1.9884e+30
c = G*M_sun


# conditions
a = 0; b = 3e9
N = 1e4
h = (b-a)/N
M = np.asarray([m1, m2, m3])
r1 = np.asarray([1429.258e9, 0])
r2 = np.asarray([777.258e9, 0])
r3 = np.asarray([-742e6, 0])
v1 = np.asarray([0, 9.09e3])
v2 = np.asarray([0, 12.44e3])
v3 = np.asarray([0, -12])


tpoints = np.arange(a, b, h)    # time points
y0 = np.array([r1, r2, r3, v1, v2, v3])   # ODE values


# CoM formula

r_com = (m1*r1+m2*r2+m3*r3) / (m1+m2+m3)
v_com = (m1*v1+m2*v2*m3*v3) / (m1+m2+m3)

# ODE equations


def orbits(u):
    r1, r2, r3, v1, v2, v3 = u
    r12 = np.linalg.norm(r2-r1)
    r13 = np.linalg.norm(r3-r1)
    r23 = np.linalg.norm(r3-r2)
    f1 = (-G*m2/r12**3)*(r1-r2) + (-G*m3/r13**3)*(r1-r3)
    f2 = (-G*m1/r12**3)*(r2-r1) + (-G*m3/r23**3)*(r2-r3)
    f3 = (-G*m1/r13**3)*(r3-r1) + (-G*m2/r23**3)*(r3 - r2)
    return np.array([v1, v2, v3, f1, f2, f3])

# RK4 algorithm


def RK4step(f, u, h):
    k1 = h*f(u)
    k2 = h*f(u+0.5*k1)
    k3 = h*f(u+0.5*k2)
    k4 = h*f(u+k3)
    return u + (k1+2*k2+2*k3+k4)/6

# RK4 integration


def RK4integrate(f, y0, tspan):
    y = np.zeros((len(tspan), len(y0), len(y0[0])))
    y[0] = y0
    for k in range(1, len(tspan)):
        y[k] = RK4step(f, y[k-1], tspan[k]-tspan[k-1])
    return y

# Run

sol_RK4 = RK4integrate(orbits, y0, tpoints)  # solve RK4
X, Y = sol_RK4.T    # transpose matrix

r1X, r2X, r3X, v1X, v2X, v3X = np.asarray(X)
r1Y, r2Y, r3Y, v1Y, v2Y, v3Y = np.asarray(Y)

# Analysis

v1 = np.hypot(v1X, v1Y)
v2 = np.hypot(v2X, v2Y)
v3 = np.hypot(v3X, v3Y)
r1 = np.hypot(r1X, r1Y)
r2 = np.hypot(r2X, r2Y)
r3 = np.hypot(r3X, r3Y)
r12 = np.linalg.norm(r2-r1)
r13 = np.linalg.norm(r3-r1)
r23 = np.linalg.norm(r3-r2)
dis1 = r1-1429.258e9    # displacment
dis2 = r2-777.258e9
dis3 = r3+742e6
r_com = (m1*r1+m2*r2+m3*r3) / (m1+m2+m3)
v_com = (m1*v1+m2*v2+m3*v3) / (m1+m2+m3)

# Calculate aphelion and perihelion
#print(max(r2X), min(r2X))
#print(max(r1X), min(r1X))

# FIG 10
plt.scatter(r1X, r1Y, c=v1, s=1, cmap='turbo')
plt.scatter(r2X, r2Y, c=v2, s=1, cmap='turbo')
plt.scatter(r3X, r3Y, c=v2, s=1, cmap='turbo')
plt.colorbar(label='$v$ [m/s]', orientation='vertical')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('RK4: 3-body problem')
#plt.show()

# FIG 11 and 12
plt.plot(tpoints, dis1)
plt.plot(tpoints, dis2)
plt.plot(tpoints, dis3)
plt.scatter(tpoints, r_com, s=0.05, c= "black")
plt.xlabel('Time [s]')
plt.ylabel('Displacment [m]')
plt.title('RK4: 3-body problem')
plt.legend(['Saturn', 'Jupiter', 'Sun', 'CoM'])
#plt.show()
