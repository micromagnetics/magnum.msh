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
# Last modified by Florian Bruckner, 2015-11-19

from dolfin import *
import numpy as np

class MeshTool(object):
  def __init__(self):
    """
    Class provided various operations with dolfin meshes.

    Examples:
      mesh1 = UnitCubeMesh(2,2,2)
      mesh2 = UnitCubeMesh(1,1,1)
      mesh3 = Mesh("mesh_with_domains.xml")

      mt = MeshTool()
      mt.add(mesh1, domain_id=1)
      mt.add(mesh2, domain_id=2)
      mt.add(mesh3, domain_id={1:3,2:4})
      File("mesh.pvd") << mt.mesh()
    """
    self._mesh = None

  def add(self, mesh, domain_id=None):
    """
    Add (non-overlapping) mesh

    *Arguments*
      mesh (:class:`Mesh`)
        Mesh that should be added
      domain_id (:class:`int` or :class:`dict`)
        domain ID if not included in mesh, or dictionary map if domains are included and should be re-mapped
    """
    editor = MeshEditor()
    newmesh = Mesh()
    editor.open(newmesh, 3, 3)
    if self._mesh is None:
      N = 0
      NE = 0
    else:
      N = self._mesh.num_vertices()
      NE = self._mesh.num_cells()

    editor.init_vertices(N+mesh.num_vertices())
    editor.init_cells(NE+mesh.num_cells())

    # add original mesh
    if self._mesh is not None:
      map(lambda n, x: editor.add_vertex(n, x), np.arange(N), self._mesh.coordinates())
      map(lambda n, c: editor.add_cell(n, np.array(c, dtype=np.uintp)), np.arange(NE), self._mesh.cells())

    # add given mesh (shift indices)
    map(lambda n, x: editor.add_vertex(N+n, x), np.arange(mesh.num_vertices()), mesh.coordinates())
    map(lambda n, c: editor.add_cell(NE+n, np.array(c+N, dtype=np.uintp)), np.arange(mesh.num_cells()), mesh.cells())
    editor.close()

    # handle domains
    d = newmesh.domains()
    d.init(3)

    if self._mesh is None or self._mesh.domains().num_marked(3) == 0:
      domains_orig = None
    else:
      domains_orig = MeshFunction('size_t', self._mesh, 3, self._mesh.domains()).array()
      map(lambda i: d.set_marker((i, domains_orig[i]), 3), range(len(domains_orig)))

    if mesh.domains().num_marked(3) == 0:
      if domain_id is None:
        if domains_orig is not None:
          raise ValueError("Specified Mesh contains no domains IDs, but original mesh does!")
        else:
          pass
      elif isinstance(domain_id, int):
        if domains_orig is None and self._mesh is not None:
          raise ValueError("Specified domain_id, but original mesh does not contains domain IDs!")
        else:
          map(lambda i: d.set_marker((NE+i, domain_id), 3), range(mesh.num_cells()))
      else:
        raise ValueError("Specified Mesh does not contain domains IDs. The specified domain_id needs to be 'Integer' or 'None'!")
    else:
      if self._mesh is not None and domains_orig is None:
        raise ValueError("Specified Mesh contains domains IDs, but original mesh does not!")
      else:
        domains = MeshFunction('size_t', mesh, 3, mesh.domains()).array()
        unique_domains = np.unique(domains)
        n_domains = len(unique_domains)
        if isinstance(domain_id, dict):
          if set(unique_domains) - set(domain_id.keys()) != set():
            raise ValueError("Specified domain_id map does contain all occuring domain IDs (%s != %s)!" % (set(unique_domains), set(domain_id.keys())))
          domain_map = domain_id
        elif domain_id is None:
          domain_map = range(n_domains+1)
        else:
          raise ValueError("Specified Mesh contains domains IDs, but specified domain_id-map has wrong type!")
        map(lambda i: d.set_marker((NE+i, domain_map[domains[i]]), 3), range(mesh.num_cells()))

    self._mesh = newmesh

  def mesh(self):
    return self._mesh
