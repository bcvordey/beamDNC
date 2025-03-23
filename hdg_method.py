
from ngsolve import *
def hyb_Con_IntPenalty(u, uhat, v, vhat, order, dx):
    n = specialcf.normal(2)
    h = specialcf.mesh_size

    def jumpdn(v,vhat):
        return n*(grad(v)-vhat)
    def hesse(v):
        return v.Operator("hesse")
    def hessenn(v):
        return InnerProduct(n, hesse(v)*n)

    dS = dx(element_boundary=True)
    return InnerProduct (hesse(u), hesse(v)) * dx \
        - hessenn(u) * jumpdn(v,vhat) * dS \
        - hessenn(v) * jumpdn(u,uhat) * dS \
        + 3*order*order/h * jumpdn(u,uhat) * jumpdn(v,vhat) * dS
