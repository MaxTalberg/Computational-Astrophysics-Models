#
#
# Coursework A
# Question 1
# ____b_____
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
h = (b-a)/N     #2.5e5
x = 5.2e12; y = 0
vx = 0; vy = 880

tpoints = np.arange(a, b, h)    # time points
y0 = np.array([x, y, vx, vy])   # ODE values

tol = 1e-6

# ODE equations


def orbits(u):
    x, y, vx, vy = u
    r = np.hypot(x, y)
    f = -c/r**3
    return np.array([vx, vy, f*x, f*y])

# Runge-Kutta Cash-Karp algorithm


def RKCK(f, t, u, h, tol):
    err = 2 * tol
    while (err > tol):
        tn = t + h
        k1 = h*f(u)
        k2 = h*f(u+(1/5)*k1)
        k3 = h*f(u+(3/40)*k1+((9/40)*k2))
        k4 = h*f(u+(3/10)*k1-((9/10)*k2)+((6/5)*k3))
        k5 = h*f(u-(11/54)*k1+((5/2)*k2)-((70/27)*k3)+((35/27)*k4))
        k6 = h*f(u+(1631/55296)*k1+((175/512)*k2)+((575/13824)*k3)+((44275/110592)*k4)+((253/4096)*k5))
        or4 = ((37/378)*k1)+((250/621)*k3)+((125/594)*k4)+((512/1771)*k6)
        or5 = ((2825/27648)*k1)+((18575/48384)*k3)+((13525/55296)*k4)+((277/14336)*k5)+((1/4)*k6)
        err = 1e2*tol+max(abs(or4-or5))
        h = 0.8 * h * (tol * h / err) ** (1 / 4)
        un = u + or4
        return un, h, tn

# RK4 integration


def RKCKintegrate(f, y0, a, b, h, tol):
    Z = [y0]        # initial values
    H = [h]
    T = [a]
    z = y0
    t = a           # initial time
    while (t < b):
        h = min(h, b - t)
        z, h, t = RKCK(f, t, z, h, tol)
        Z.append(z)
        H.append(h)
        T.append(t)
    return np.array(Z), np.asarray(H), np.asarray(T)

# Run


sol_RK4 = RKCKintegrate(orbits, y0, a, b, h, tol)  # solve RK4
x, y, vx, vy = sol_RK4[0].T    # transpose matrix
H = sol_RK4[1]
T = sol_RK4[2]

# Analysis

Npoints = list(range(len(H)))
v = np.hypot(vx, vy)
r = np.hypot(x, y)
KE = 0.5*v**2
GPE = -(c/r)

# calculating aphelion
#print(max(x), min(x))

# plotting results

# FIG4
plt.scatter(x, y, c=v, s=1, cmap='turbo')
plt.colorbar(label='$v$ [m/s]', orientation='vertical')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('RK CK: Halleys Comet')
plt.scatter(0, 0, marker="o", c = 'black')
#plt.show()

# FIG5  b*5
plt.scatter(Npoints, H, c=v, s=1, cmap='turbo')
plt.colorbar(label='$v$ [m/s]', orientation='vertical')
plt.xlabel('Step')
plt.ylabel('Step size [seconds]')
plt.title('RK CK: Halleys Comet')
#plt.show()

# FIG6
plt.plot(T, KE, c = "red")
plt.plot(T, GPE, c = "blue")
plt.plot(T, KE + GPE, c = "black")
plt.xlabel('Time')
plt.ylabel('Energy')
plt.legend(['KE', 'GPE', 'Etot'])
plt.title('RK CK: Halleys Comet')
#plt.show()
