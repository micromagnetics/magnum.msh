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

    self.assertEqual(mt.mesh().num_cells(), 30)
    self.assertEqual(mt.mesh().num_vertices(), 26)

  def test_add_domains(self):
    mesh1 = Mesh("mesh/mesh_with_domains.xml.gz")
    mesh1b = Mesh("mesh/mesh_with_domains.xml.gz")
    mesh1b.translate(Point(-2.1,1.1,0.0))
    mesh2 = UnitCubeMesh(1,1,1)
    mesh2.translate(Point(1.1,0.0,0.0))
    mesh3 = UnitCubeMesh(1,2,2)

    mt = MeshTool()
    mt.add(mesh1) # contains domains
    mt.add(mesh2, domain_id=4) # add domain 4

    self.assertEqual(mt.mesh().num_cells(), 6006)
    self.assertEqual(mt.mesh().num_vertices(), 1339)

    domains = MeshFunction("size_t", mt.mesh(), 3, mt.mesh().domains())
    self.assertEqual(len(np.unique(domains.array())), 4)

    with self.assertRaises(ValueError):
      mt.add(mesh3) # add domain without ID -> error

    with self.assertRaises(ValueError):
      mt.add(mesh1b, [5,6,7]) # use domain_id map of wrong type -> error
    with self.assertRaises(ValueError):
      mt.add(mesh1b, {1:5,2:6,4:7}) # use domain_id map with missing domain -> error
    mt.add(mesh1b, {1:5,2:6,3:7,4:8}) # use domain_id map (with additional domain)

    domains = MeshFunction("size_t", mt.mesh(), 3, mt.mesh().domains())
    self.assertEqual(mt.mesh().num_cells(), 12006)
    self.assertEqual(mt.mesh().num_vertices(), 2670)
    self.assertEqual(len(np.unique(domains.array())), 7)

    mt = MeshTool()
    mt.add(mesh2) # start with mesh without domain IDs
    with self.assertRaises(ValueError):
      mt.add(mesh3, domain_id=4) # add mesh with domain_id -> error

    mt.add(mesh3) # add mesh without id
    with self.assertRaises(ValueError): # add mesh with domain IDs -> error
      mt.add(mesh1b)

    self.assertEqual(mt.mesh().num_cells(), 30)
    self.assertEqual(mt.mesh().num_vertices(), 26)


if __name__ == '__main__':
    unittest.main()
