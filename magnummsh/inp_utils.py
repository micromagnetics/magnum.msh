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
# Last modified by Florian Bruckner, 2015-10-14

from dolfin import *
import numpy as np
from StringIO import StringIO

class INPReader(object):
  def __init__(self, filename):
    self._num_materials = None
    self._node_data = None

    # TODO: handle element data
    # TODO: fix problem with _node_data datatypes
    with open(filename, "r") as f:
      lines = f.readlines()
      NNodes, NElements, NodeData, ElementData, TensorData = [int(x) for x in lines[0].split()]
      self.nodes = np.genfromtxt(StringIO(''.join(lines[1:NNodes+1])), dtype=(int, float, float, float))
      self.elements = np.genfromtxt(StringIO(''.join(lines[NNodes+1:NNodes+NElements+1])), dtype=(int, int, np.dtype((str, 10)), int, int, int, int))

      if NodeData > 0:
        X = [int(x) for x in lines[NNodes+NElements+1].split()]
        NNodeDataComponents = X[0]
        NodeDataComponenetSizes = X[1:]
        NodeDataSize = sum(NodeDataComponenetSizes)
        NodeDataTags = [lines[i].split(',')[0] for i in range(NNodes+NElements+2,NNodes+NElements+NNodeDataComponents+2)]
        self._node_data = np.genfromtxt(StringIO(''.join(lines[NNodes+NElements+NNodeDataComponents+2:2*NNodes+NElements+NNodeDataComponents+2])), dtype=(float, float))
    self._need_update = True

  def mesh(self, include_domains=True):
    if self._need_update == True:
      self._need_update = False
      editor = MeshEditor()
      self._mesh = Mesh()
      editor.open(self._mesh, "tetrahedron", 3, 3)
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
    for i in range(1, self._node_data.shape[1]):
      u = Function(V)
      u.vector()[vertex_to_dof_map(V)] = self._node_data[:, i].copy()
      components.append(u)

    if len(components) == 1:
      return u
    else:
      VV = VectorFunctionSpace(mesh, "CG", 1, self._node_data.shape[1]-1)
      res = Function(VV, name=name)
      assign(res, components)
      return res

  def print_stats(self):
    print "Found Mesh:"
    print "  Nodes:     %-10d" % len(self.nodes)
    print "  Elements:  %-10d" % len(self.elements)
    if self._num_materials is not None:
      print "  Materials: %-10d" % len(self._num_materials)
    if self._node_data is not None:
      print "  Node Data: %-10d" % (self._node_data.shape[1]-1)


class INPWriter(object):
  def __init__(self, mesh, field=None):
    self._mesh = mesh
    self._domains = MeshFunction('size_t', self._mesh, 3, self._mesh.domains())
    if field is None:
      self._fields = None
    else:
      V = FunctionSpace(self._mesh, "CG", 1)
      VV = VectorFunctionSpace(self._mesh, "CG", 1)
      u, v, w = Function(V), Function(V), Function(V)
      assign([u,v,w], Function(VV, field))
      self._fields = (u, v, w)
      self._map = vertex_to_dof_map(V)

  def write_inp(self, filename):
    with open(filename, "w") as f:
      # NNodes, NElements, NodeData, ElementData, TensorData
      N = self._mesh.num_vertices()
      NE = self._mesh.num_cells()
      f.write("%d %d %d %d %d\n" % (N, NE, 0 if self._fields is None else 1, 0, 0))
      nodes = np.hstack((np.arange(1,N+1).reshape(N,1), self._mesh.coordinates()))
      np.savetxt(f, nodes, fmt="%8d %.15e %.15e %.15e")
      if self._mesh.domains().num_marked(3) != 0:
        self._mat = self._domains.array()
      else:
        self._mat = np.repeat(1, NE)
      elements = np.hstack((np.arange(1,NE+1).reshape(NE,1), self._mat.reshape(NE,1), self._mesh.cells()+1))
      np.savetxt(f, elements, fmt="%8d %4d tet %8d %8d %8d %8d")
      if self._fields is not None:
        f.write("3 1 1 1\n")
        f.write("dataX, none\n")
        f.write("dataY, none\n")
        f.write("dataZ, none\n")

        field = np.hstack((np.arange(1,N+1).reshape(N,1),
                           self._fields[0].vector()[self._map].reshape(N,1),
                           self._fields[1].vector()[self._map].reshape(N,1),
                           self._fields[2].vector()[self._map].reshape(N,1)))
        np.savetxt(f, field, fmt="%8d %.15e %.15e %.15e")

  def print_stats(self):
    print "Found Mesh:"
    print "  Nodes:     %-10d" % self._mesh.num_vertices()
    print "  Elements:  %-10d" % self._mesh.num_cells()
    print "  Materials: %-10d" % len(np.unique(self._mat))
    #print "  Node Data: %-10d" % (self._node_data.shape[1]-1)

