#!/usr/bin/python

# Copyright (C) 2016-2017 Claas Abert, Florian Bruckner
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
# Last modified by Florian Bruckner, 2017-01-27

try:
  import gmsh
  _gmsh_found = True
except Exception as e:
  _gmsh_found = False
  _gmsh_exeption = e

import dolfin as df
import numpy as np
import os

class GmshReader(object):
  def __init__(self, filename):
    if not _gmsh_found: raise _gmsh_exeption

    # suppress loggin
    if not os.path.isfile(filename):
      raise IOError("No such file or directory: '%s'" % filename)
    gmsh.initialize()
    gmsh.open(filename)
    self._model = gmsh.model
    self._physical_cell_ids = None
    self._physical_face_ids = None

  def mesh(self):
    # prepare FEniCS mesh and mesh editor
    mesh = df.Mesh()
    mesh_editor = df.MeshEditor()
    mesh_editor.open(mesh, "tetrahedron", 3, 3);

    # prepare container for physical IDs
    self._physical_cell_ids = set()
    self._physical_face_ids = set()

    # get mesh vertices
    ids, nodes, _ = self._model.mesh.getNodes()
    ids = np.array(ids)-1
    nodes = np.array(nodes).reshape((-1,3))
    mesh_editor.init_vertices(len(ids))
    for n, xyz in zip(ids,nodes):
        mesh_editor.add_vertex(n, df.Point(xyz))

    # get mesh cells
    _, ids, cells = self._model.mesh.getElements(3)
    ids = np.arange(len(ids[0]))
    cells = np.array(cells,dtype=np.uintp).reshape((-1,4))-1
    mesh_editor.init_cells(len(ids))

    # read PhysicalGroup ids for cells
    i = 0
    cell_groups = self._model.getPhysicalGroups(3)
    if len(cell_groups) == 0:
      for n, ijkl in zip(ids,cells):
        mesh_editor.add_cell(n, ijkl)
    else:
      for g in cell_groups:
        entities = self._model.getEntitiesForPhysicalGroup(*g)
        for e in entities:
          _, ids, cells = self._model.mesh.getElements(3, e)
          cells = np.array(cells,dtype=np.uintp).reshape((-1,4))-1
          for ijkl in cells:
            mesh_editor.add_cell(i, ijkl)
            mesh.domains().set_marker((i, g[1]), 3)
            self._physical_cell_ids.add(g[1])
            i += 1

    # finish mesh creation
    mesh_editor.close()

    # set facet domains
    mesh.init(0, 2)
    facet_groups = self._model.getPhysicalGroups(2)
    for g in facet_groups:
      entities = self._model.getEntitiesForPhysicalGroup(*g)
      for e in entities:
        _, ids, cells = self._model.mesh.getElements(2, e)
        facets = np.array(cells,dtype=np.uintp).reshape((-1,3))-1
        for ijb in facets:
          i = reduce(lambda x,y: x.intersection(y), map(lambda i: set(mesh.topology()(0,2)(ijb[i])), range(3))).pop()
          mesh.domains().set_marker((i, g[1]), 2)
          self._physical_face_ids.add(g[1])
    return mesh

  def print_domain_info(self):
    # initialize physical IDs
    if self._physical_face_ids is None or self._physical_cell_ids is None: self.mesh()

    # domain name mappings
    print "========================================"
    print " DOMAIN NAME MAPPING"
    print " Type  |  ID  | Domain name"
    print "-------+------+-------------------------"

    for physical_id in self._physical_face_ids:
      print " Facet |% 5d | %s" % (physical_id, self._model.getPhysicalName(2, physical_id))
    for physical_id in self._physical_cell_ids:
      print " Cell  |% 5d | %s" % (physical_id, self._model.getPhysicalName(3, physical_id))

    print "=======+======+========================="

  def _create_cuboid_geo(self, size, n):
    dx, dy, dz = size
    lc = 1.0 # TODO what about this value?

    vertices = gp.GVertexVector()
    vertices.push_back(self._model.addVertex(-dx, -dy, -dz, lc))
    vertices.push_back(self._model.addVertex(-dx, -dy,  dz, lc))
    vertices.push_back(self._model.addVertex(-dx,  dy, -dz, lc))
    vertices.push_back(self._model.addVertex(-dx,  dy,  dz, lc))
    vertices.push_back(self._model.addVertex( dx, -dy, -dz, lc))
    vertices.push_back(self._model.addVertex( dx, -dy,  dz, lc))
    vertices.push_back(self._model.addVertex( dx,  dy, -dz, lc))
    vertices.push_back(self._model.addVertex( dx,  dy,  dz, lc))

    edges = gp.GEdgeVector()
    direction = [2, 1, 2, 1, 2, 1, 2, 1, 0, 0, 0, 0]
    for i in range(12):
      edge = self._model.addLine(vertices[_vertex_data[i][0]], vertices[_vertex_data[i][1]])
      edge.resetMeshAttributes()
      edge.setTransfinite(nbPointsTransfinite=n[direction[i]] + 1, typeTransfinite=1, coeffTransfinite=1)
      edges.push_back(edge)

    faces = gp.GFaceVector()
    for i in range(6):
      a = gp.GEdgeVector()
      for j in range(4):
        a.push_back(edges[_face_data[i][j]])
      aa = gp.GEdgeVectorVector()
      aa.append(a)
      face = self._model.addPlanarFace(aa)
      face.setTransfinite() #transfiniteArrangement = -1
      faces.push_back(face)

    return vertices, edges, faces

  def _prepare_domains(self):
    if self._model.noPhysicalGroups(): # handle simple mesh with just one volume and one surface
      # load entities of sample (edges, vertices, region)
      #self._model.createTopologyFromMesh()
      #self._model.removeDuplicateMeshVertices(1e-18);
      faces = self._model.bindingsGetFaces()
      if len(faces) == 0:
        raise RuntimeError("GmshModel: No faces found (maybe createTopologyFromMesh is needed)!")
      sample_faces = gp.GFaceVector()
      for f in faces:
        f.addPhysicalEntity(1)
        sample_faces.push_back(f)

      for region in self._model.bindingsGetRegions():
        region.addPhysicalEntity(1)
    else: # handle sophisticated meshes
      # automatically retrieve outer faces of sample for meshing of shell
      #self._model.createTopologyFromMesh();
      #self._model.removeDuplicateMeshVertices(1e-18);
      # TODO: handle holes
      all_regions = gp.GRegionVector()
      for region in self._model.bindingsGetRegions():
        all_regions.push_back(region)
      compound = gp.GRegionCompound(self._model, 0, all_regions)
      faces = compound.bindingsGetFaces()
      if len(faces) == 0:
        raise RuntimeError("GmshModel: No faces found (maybe createTopologyFromMesh is needed)!")

      sample_faces = gp.GFaceVector()
      for f in faces:
        #f.addPhysicalEntity(1)
        sample_faces.push_back(f)

    bounds = self._model.bounds()
    sample_size = np.array([ max(np.fabs(bounds.min().x()), np.fabs(bounds.max().x())),
                             max(np.fabs(bounds.min().y()), np.fabs(bounds.max().y())),
                             max(np.fabs(bounds.min().z()), np.fabs(bounds.max().z()))
                           ])
    return sample_faces, sample_size

  def create_shell(self, d, margin, n, shell_progression=1.0):
    # retrieve box size, add margin and create inner cuboid
    sample_faces, sample_size = self._prepare_domains()
    size = sample_size + margin
    inner_vertices, inner_edges, inner_faces = self._create_cuboid_geo(size, n)

    # retrieve shell width create outer cuboid
    shell_width = min(sample_size) + margin
    size = sample_size + shell_width
    outer_vertices, outer_edges, outer_faces = self._create_cuboid_geo(size, n)

    # setup connections between cuboids
    connection_edges = gp.GEdgeVector()
    for i in range(8):
      connection = self._model.addLine(inner_vertices[i], outer_vertices[i])
      connection.setTransfinite(nbPointsTransfinite=d+1, typeTransfinite=1, coeffTransfinite=shell_progression)
      connection_edges.push_back(connection)

    connection_faces = gp.GFaceVector()
    for i in range(12):
      a = gp.GEdgeVector()
      a.push_back(inner_edges[i])
      a.push_back(connection_edges[_vertex_data[i][0]])
      a.push_back(outer_edges[i])
      a.push_back(connection_edges[_vertex_data[i][1]])
      aa = gp.GEdgeVectorVector()
      aa.append(a)
      face = self._model.addPlanarFace(aa)
      face.setTransfinite()
      connection_faces.push_back(face)

    # add air region
    faces = gp.GFaceVectorVector()
    faces.push_back(inner_faces)
    faces.push_back(sample_faces)
    air_region = self._model.addVolume(faces)
    air_region.addPhysicalEntity(1000)

    # shell
    shell_regions = gp.GRegionVector()
    for i in range(6):
      a = gp.GFaceVector()
      a.push_back(inner_faces[i])
      a.push_back(outer_faces[i])
      for j in range(4):
        a.push_back(connection_faces[_face_data[i][j]]);
      aa = gp.GFaceVectorVector()
      aa.push_back(a)
      region = self._model.addVolume(aa)
      region.setTransfinite()
      shell_regions.push_back(region);

    for i in range(6):
      # i/2+2 leads to 1001, 1002, 1003 for x, y, z transformation area
      shell_regions[i].addPhysicalEntity(i/2+1001)

    self._model.mesh(3)
    self._model.save("shell.msh")
    self._model.writeGEO("shell.geo")

  def create_cuboid(self, size, n):
    sample_vertices, sample_edges, sample_faces = self._create_cuboid_geo(size, n)

    for f in sample_faces:
      f.addPhysicalEntity(1)
    faces = gp.GFaceVectorVector()
    faces.push_back(sample_faces)

    sample_region = self._model.addVolume(faces)
    sample_region.setTransfinite()
    sample_region.addPhysicalEntity(1)
    self._model.mesh(3)
    self._model.writeGEO("test.geo")
    self._model.save("test.msh")

_vertex_data = [
  [0,1],[1,3],[3,2],[2,0],
  [4,5],[5,7],[7,6],[6,4],
  [0,4],[1,5],[2,6],[3,7]
  ]

_face_data = [
  [0, 1, 2, 3],[4, 5, 6, 7],
  [0, 9, 4, 8],[10, 6,11,2],
  [8, 7,10, 3],[9, 5,11, 1]
  ]

