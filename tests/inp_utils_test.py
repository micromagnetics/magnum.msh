import unittest
from dolfin import *
from magnummsh import *
from tempfile import mkstemp
from os import remove 

class InpUtilsTest(unittest.TestCase):
  def test_without_domain(self):
    (fd, filename) = mkstemp()
    try: 
      mesh = UnitCubeMesh(2,2,2)

      writer = INPWriter(mesh)
      writer.write_inp(filename)

      reader = INPReader(filename)
      mesh2 = reader.mesh()
      domains = reader.domains()

      self.assertEqual(mesh.num_cells(), mesh2.num_cells())
      self.assertEqual(mesh.num_vertices(), mesh.num_vertices())
    finally: 
      remove(filename)
    

if __name__ == '__main__':
    unittest.main()
