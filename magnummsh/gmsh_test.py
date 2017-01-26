import magnummsh
from dolfin import *

reader = magnummsh.GmshReader("input.msh")
#faces, size = reader._prepare_domains()
#reader.create_cuboid((100,100,5), (20,20,1))
reader.create_shell(2, 0.3, (10,10,10))

#e = reader.create_cuboid([1,1,1], [4,4,4])
#reader.print_domain_info()
m = reader.mesh()
File("mesh.pvd") << m
File("domains.pvd") << MeshFunction("size_t", m, 3, m.domains())
