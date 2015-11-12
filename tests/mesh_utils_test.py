import unittest
from dolfin import *
from magnummsh import *

class MeshUtilsTest(unittest.TestCase):
  def test_add(self):
    mt = MeshTool()

    mesh1 = UnitCubeMesh(1,2,2)
    mesh2 = UnitCubeMesh(1,1,1)
    mesh2.translate(Point(1.1,0.0,0.0))

    mt.add(mesh1)
    mt.add(mesh2)

    File("mesh.pvd") << mt.mesh()
    self.assertEqual(mt.mesh().num_cells(), 30)
    self.assertEqual(mt.mesh().num_vertices(), 26)

if __name__ == '__main__':
    unittest.main()
