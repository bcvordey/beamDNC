from ngsolve import *
from ngsolve.webgui import Draw
mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))

order = 3

V = HDivDiv(mesh, order=order-1)
Q = H1(mesh, order=order, dirichlet="left|bottom|top|right")
X = V*Q

print ("ndof-V:", V.ndof, ", ndof-Q:", Q.ndof)

sigma, w = X.TrialFunction()
tau, v = X.TestFunction()

n = specialcf.normal(2)

def tang(u): return u-(u*n)*n

a = BilinearForm(X, symmetric=True)
a += (InnerProduct(sigma, tau) + div(sigma)*grad(v) \
      + div(tau)*grad(w) - 1e-10*w*v )*dx \
      + (-(sigma*n) * tang(grad(v)) - (tau*n)*tang(grad(w)))*dx(element_boundary=True)
a.Assemble()

f = LinearForm(-1 * v * dx).Assemble()

gfu = GridFunction(X)
gfu.vec.data = a.mat.Inverse(X.FreeDofs()) * f.vec
# Draw(gfu)
Draw (gfu.components[0], mesh, name="sigma")
Draw (gfu.components[1], mesh, name="disp")



