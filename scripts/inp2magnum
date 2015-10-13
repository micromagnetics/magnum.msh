#!/usr/bin/python

# Copyright (C) 2011-2015 Florian Bruckner
#
# This file is part of magnum.msh.
#
# magnum.msh is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# magnum.msh is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with magnum.msh. If not, see <http://www.gnu.org/licenses/>.
#
# Last modified by Florian Bruckner, 2015-10-13

import argparse
from dolfin import *
import numpy as np
from StringIO import StringIO

class INPConverter(object):
  def __init__(self):
    self._num_materials = None

  def read_inp(self, filename):
    # TODO: handle element data
    # TODO: fix problem with node_data datatypes
    with open(filename, "r") as f:
      lines = f.readlines()
      NNodes, NElements, NodeData, ElementData, TensorData = [int(x) for x in lines[0].split()]
      self.nodes = np.genfromtxt(StringIO(''.join(lines[1:NNodes+1])), dtype=(int, float, float, float))
      self.elements = np.genfromtxt(StringIO(''.join(lines[NNodes+1:NNodes+NElements+1])), dtype=(int, int, np.dtype((str, 10)), int, int, int, int))

      X = [int(x) for x in lines[NNodes+NElements+1].split()]
      NNodeDataComponents = X[0]
      NodeDataComponenetSizes = X[1:]
      NodeDataSize = sum(NodeDataComponenetSizes)
      NodeDataTags = [lines[i].split(',')[0] for i in range(NNodes+NElements+2,NNodes+NElements+NNodeDataComponents+2)]

      self.node_data = np.genfromtxt(StringIO(''.join(lines[NNodes+NElements+NNodeDataComponents+2:2*NNodes+NElements+NNodeDataComponents+2])), dtype=(float, float))
    self._need_update = True

  def mesh(self, include_domains=True):
    if self._need_update == True:
      self._need_update = False
      editor = MeshEditor()
      self._mesh = Mesh()
      editor.open(self._mesh, 3, 3)
      editor.init_vertices(len(self.nodes))
      editor.init_cells(len(self.elements))
      map(lambda x: editor.add_vertex(x[0]-1, np.array([x[1], x[2], x[3]])), self.nodes)
      map(lambda x: editor.add_cell(x[0]-1, np.array([x[3]-1, x[4]-1, x[5]-1, x[6]-1], dtype=np.uintp)), self.elements)
      editor.close()

      if include_domains == True:
        d = self._mesh.domains()
        d.init(3)
        map(lambda i: d.set_marker((i, self.elements[i][1]), 3), range(self.elements.shape[0]))
    return self._mesh

  def domains(self):
    mesh = self.mesh()
    mat = MeshFunction('size_t', mesh, 3)
    mat.array()[:] = map(lambda x: x[1], self.elements)
    self._num_materials = set(mat.array()[:])
    return mat

  def data(self, name='data'):
    mesh = self.mesh()
    V = FunctionSpace(mesh, "CG", 1)
    components = []
    for i in range(1, self.node_data.shape[1]):
      u = Function(V)
      u.vector()[vertex_to_dof_map(V)] = self.node_data[:, i].copy()
      components.append(u)

    if len(components) == 1:
      return u
    else:
      VV = VectorFunctionSpace(mesh, "CG", 1, self.node_data.shape[1]-1)
      res = Function(VV, name=name)
      assign(res, components)
      return res

  def print_stats(self):
    print "Found Mesh:"
    print "  Nodes:     %-10d" % len(self.nodes)
    print "  Elements:  %-10d" % len(self.elements)
    if self._num_materials is not None:
      print "  Materials: %-10d" % len(self._num_materials)
    print "  Node Data: %-10d" % (self.node_data.shape[1]-1)


if __name__ == "__main__":
  parser = argparse.ArgumentParser(prog='magnum-msh')
  parser.add_argument('-d', '--domains', help='Output file for domain information')
  parser.add_argument('-f', '--field', help='Output file for extracted data')
  parser.add_argument('source', help='Input INP file')
  parser.add_argument('target', help='Output mesh file(\'.pvd\', \'.xml\')')
  args = parser.parse_args()

  converter = INPConverter()
  converter.read_inp(args.source)
  File(args.target) << converter.mesh()
  if args.domains:
    File(args.domains) << converter.domains()
  if args.field:
    File(args.field) << converter.data()
  converter.print_stats()
