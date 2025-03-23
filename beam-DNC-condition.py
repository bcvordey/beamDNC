####### SOLVING THE BEAM EQUATION WITH DNC CONDITION ON THE BOUNDARY   ######

## Solve the EULER-BERNOULLI beam equation u_tt - alpha u_xxxx = f ##
# With Dirichlet boundary condition u(0,t) = phi(t) and u_x(0,t) = 0
# Neumann boundary condition on x = 1 and u(1,t) = -k(u(1,t) - y )^\delta - beta * u_t(1,t)



#  Known terms #
kappa = 0.1  #kappa
obstaclePosition = 0.1    #distance between the obstacle and the beam
beamLength = 1    #Length of the beam


# Calculating alpha #
youngModulus = 0.01 # denoted by E
crossSectionalInertia = 1 #denoted by I
linearDensity = 1 #denoted by rho 
crossSectionalArea = 1 # cross sectional area denoted by A
alpha = (youngModulus * crossSectionalInertia) / (linearDensity * crossSectionalArea)

#Time Intervals
T = 10   # Total time  


#force  
force = 1 #right hand side
beta = 0.1


delta = 1 ###


#Tolerance
Tol = 1e-5

###############################
N = 10 #number of time intervals

u0 = 0  #initial displacement of the beam
v0 = 1  #initial velocity of the beam

t0 = -obstaclePosition/v0
dt = (T - t0)/N
u_1 = v0 * dt

#number of elements
numberOfElements = 5

####Packages########
from netgen.meshing import *
from netgen.meshing import Element0D, Element1D, Element2D, MeshPoint, FaceDescriptor
from netgen.csg import Pnt


## Generate a 1D mesh ##
m = Mesh(dim = 1)
pnums = []
for i in range(0,numberOfElements+1):
    pnums.append(m.Add (MeshPoint( Pnt(i/numberOfElements, 0, 0))))

# add segments
for i in range(0,numberOfElements):
    m.Add(Element1D([pnums[i],pnums[(i+1)]], index=1))

m.SetMaterial(1,'material')

# Setting boundary elements and BC
m.Add (Element0D(pnums[0], index=1))
m.Add (Element0D(pnums[numberOfElements], index=2))

# set boundary condition names
m.SetBCName(0,'left')
m.SetBCName(1,'right')


## Ngsolve Packages ##
import ngsolve
from ngsolve import *
from hdg_method import hyb_Con_IntPenalty

## Ngsolve Meshing ##
ngsmesh = ngsolve.Mesh(m)

## Change the order of the FEM Method ##
order = 3

V1 = H1(ngsmesh, order=order, dirichlet="left")
V2 = NormalFacetFESpace(ngsmesh, order=order-1, dirichlet="left")
fes = V1 * V2

gfu = GridFunction(fes)
gfu0 = GridFunction(fes)
gfu1= GridFunction(fes)
gfu0.Set(u0)
gfu1.Set(u_1)

U_valsAtFreeEnd = []
U_valsAtClampedEnd = []

tn = 0

for i in range(N+1):
    gfu.Set(0)   #setting all displacement values to zero before computing new ones 

    #setting up DG to deal with grad^2u * grad^2v team

    u, uhat = fes.TnT() # making u and uhat symbolic objects  
    v, vhat = fes.TnT() # making u and vhat symbolic objects  
    
    # when beam is in Contact with obstacle
    if abs(gfu1.vec[numberOfElements] - obstaclePosition) < Tol or gfu1.vec[numberOfElements] > obstaclePosition:

        a = BilinearForm(fes)
    
        a += (
            1/dt**2*u*v*dx  
            + alpha**2*(hyb_Con_IntPenalty(u=u, uhat=uhat, v=v, vhat=vhat, order=order, dx=dx)) 
            - kappa*u*v*ds(definedon='right') 
            - (beta)/(dt)*u*v*ds(definedon='right') 
            ).Assemble()

        f = LinearForm(
            force*v*dx 
            + 1/(dt**2)*(2*gfu1 - gfu0)*v*dx 
            - kappa*obstaclePosition*v*ds(definedon='right')  
            - (beta)/(dt)*gfu1*v*ds(definedon='right')
            ).Assemble()
    # when beam is NOT in Contact with obstacle
    else:
        a = BilinearForm(fes)
        a += (
            1/dt**2*u*v*dx  
            + alpha**2*(hyb_Con_IntPenalty(u=u, uhat=uhat, v=v, vhat=vhat, order=order, dx=dx))
            ).Assemble()
        
        f = LinearForm( 
            force*v*dx 
            + 1/dt**2*(2*gfu1 - gfu0)*v*dx 
            ).Assemble()
             
    gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * f.vec 


## Resetting u values ##
    gfu0.vec.data = gfu1.vec.data
    gfu1.vec.data = gfu.vec.data
    
    tn += dt
    
    U_valsAtClampedEnd.append(gfu.vec[0]) 
    U_valsAtFreeEnd.append(gfu.vec[numberOfElements])
    

## Plotting ##
from numpy import *
import math
import matplotlib.pyplot as plt
t_vals = linspace(t0,T,N+1)



plt.plot(U_valsAtClampedEnd, t_vals, color = 'red', label = 'Clamped End of the Beam') #left boundary
plt.plot(U_valsAtClampedEnd, t_vals, color = 'blue', label = 'Free End of the Beam') #right boundary
plt.axvline(x = obstaclePosition, color = 'g', ls = "--", label = 'Obstacle') # obstacle Position 
plt.xlabel('x')
plt.ylabel('Left and Right ends of Rod')
plt.title("Finite Element Method ")
plt.legend()
plt.show()



