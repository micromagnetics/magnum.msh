import magnummsh
from dolfin import *
from magnumfe import *

reader = magnummsh.GmshReader("stack.msh")
orig = reader.mesh()
File("origdomains.pvd") << MeshFunction("size_t", orig, 3, orig.domains())
File("origfacets.pvd") << MeshFunction("size_t", orig, 2, orig.domains())
#faces, size = reader._prepare_domains()
#reader.create_cuboid((100,100,5), (20,20,1))
reader.create_shell(2, 0.5, (10,10,10))

#e = reader.create_cuboid([1,1,1], [4,4,4])
#reader.print_domain_info()
m = reader.mesh()
File("mesh.pvd") << m
File("domains.pvd") << MeshFunction("size_t", m, 3, m.domains())

state = State(m, cell_domains={'magnetic':(1,2,3,4,5,6)})
sm = state['magnetic'].mesh
File("subdomains.pvd") << MeshFunction("size_t", sm, 3, sm.domains())
File("subfacets.pvd") << MeshFunction("size_t", sm, 2, sm.domains())
print orig.coordinates()[0:10]
print m.coordinates()[0:10]
