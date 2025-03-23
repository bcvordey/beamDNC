import numpy as np
import math
import matplotlib.pyplot as plt

# Parameters
L = 1   # Length of the domain
T = 10   # Total time
M = 100  # Number of spatial points
N = 10000  # Number of time steps
c = 2  # Wave speed

# Other special parameters
beta = 0.1
h = 1
k = 0.1

#Velocity and RHS
v0 = 1
f = 10
t0 = -h/v0

#Tolerance
Tol = 1e-12



# Spatial and time steps
dx = L / M
dt = (T-t0) / N

# phi
phi = 1/(1 + k*dx + beta*dx/dt)

# Initialize solution arrays
u = np.zeros((N+1, M+1))


# Initial condition
u[0,:] = 0  
u[1,:] = v0 * dt


for n in range(1, N):
    for i in range(1, M):
        u[n+1,i] = 2*(1 - c**2 * (dt/dx**2)**2) * u[n,i] + \
                    c**2 * (dt/dx**2)**2 * (u[n,i+2] - 4*u[n,i+1]  - 4*u[n,i-1] + u[n,i-2]) - \
                    u[n-1,i] + \
                    (dt**2)*f
                    
    #left End
    u[n+1,0] = u[n+1,1]

    #Contact conditions at right end
    if abs(u[n,M] - h) < Tol or u[n,M] > h: 
        u[n+1,M] = phi * (u[n+1,M-1] + k*h*dx + beta*(dx/dt)*u[n,M])
 
    else:
        u[n+1,M] = u[n+1, M - 1]


# Plot results
plt.plot(u[:, 0] - L, np.linspace(t0, T, N+1), color = 'red', label = 'u(-l,t)') #left boundary
plt.plot(u[:, -1], np.linspace(t0, T, N+1), color = 'blue', label = 'u(0,t)') #right boundary
plt.axvline(x = h, color = 'g', ls = "--", label = 'h') # obstacle at h 
plt.xlabel('x')
plt.ylabel('u(-l, t) and u(0,t)')
plt.title("FDM")
plt.legend()
plt.show()