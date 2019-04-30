import unittest
from magnummsh import *
from magnumfe import *

set_fenics_log_level(40)
set_log_level(40)

class GmshUtilsTest(unittest.TestCase):
  def test_wrong_filename(self):
    with self.assertRaises(IOError):
      reader = GmshReader("mesh/nofile.msh")
      mesh = reader.mesh()

  def test_no_domains(self):
    reader = GmshReader("mesh/cube_nodomains.msh")
    mesh = reader.mesh()
    self.assertEqual(mesh.num_vertices(), 8)
    self.assertEqual(mesh.num_cells(), 6)
    self.assertEqual(set(mesh.domains().markers(2).values()), set({}))
    self.assertEqual(set(mesh.domains().markers(3).values()), set({}))
    state = State(mesh)
    self.assertAlmostEqual(state.volume(), 1.0)

  def test_with_domains(self):
    reader = GmshReader("mesh/multidomain.msh")
    mesh = reader.mesh()
    self.assertEqual(mesh.num_vertices(), 245)
    self.assertEqual(mesh.num_cells(), 921)
    self.assertEqual(set(mesh.domains().markers(2).values()), set({2,3,4}))
    self.assertEqual(set(mesh.domains().markers(3).values()), set({1,2}))

    mesh = Mesh("mesh.xml")
    state = State(mesh, cell_domains={'domain1':1, 'domain2':2}, facet_domains={'face1':2, 'face2':3, 'face3':4})
    self.assertAlmostEqual(state.volume(), 259.05740941566467)
    self.assertAlmostEqual(state['domain1'].volume(), 24.0)
    self.assertAlmostEqual(state['domain2'].volume(), 235.057409416)

#  def test_shell_simple(self):
#    reader = GmshReader("mesh/new.msh")
#    #reader = GmshReader("mesh/cube_nodomains.msh")
#    #reader = GmshReader("mesh/simple.msh")
#    #reader = GmshReader("mesh/sphere.msh")
#    reader.create_shell(2, 0.5, (3,3,3))
#    mesh = reader.mesh()
#    File("simple.pvd") << MeshFunction("size_t", mesh, 3, mesh.domains())
#    state = State(mesh, cell_domains={'magnetic':1, 'air':1000, 'shell':(1001,1002,1003)})
#    self.assertAlmostEqual(state['magnetic'].volume(), 1.0)
#    self.assertAlmostEqual(state['air'].volume(), 26.0)
#
##  @unittest.expectedFailure
#  def test_shell_simple_broken(self):
#    reader = GmshReader("mesh/cube_nodomains.msh")
#    reader.create_shell(2, 0.5, (3,3,3))
#    mesh = reader.mesh()
#    state = State(mesh, cell_domains={'magnetic':1, 'air':1000, 'shell':(1001,1002,1003)})
#    self.assertAlmostEqual(state['magnetic'].volume(), 1.0)
#    self.assertAlmostEqual(state['air'].volume(), 26.0)
#
#  def test_shell_domains(self):
#    reader = GmshReader("mesh/stack.msh")
#    reader.create_shell(2, 0.5, (10,10,10))
#    mesh = reader.mesh()
#    state = State(mesh, cell_domains={'magnetic':(1,2,3,4,5,6), 'air':1000, 'shell':(1001,1002,1003)})
#    self.assertAlmostEqual(state['magnetic'].volume(), 69160.89421653724)
#    self.assertAlmostEqual(state['air'].volume(), 26.0)

if __name__ == '__main__':
    unittest.main()
