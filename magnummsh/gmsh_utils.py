import gmshpy as gp
import dolfin as df
import numpy as np

class GmshReader(object):
  def __init__(self, filename):
    # suppress loggin
    gp.Msg_SetVerbosity(0)

    # load GMSH model
    gp.GmshOpenProject(filename)
    self._model = gp.GModel.current()

    # init attributes
    self._physical_cell_ids = None
    self._physical_face_ids = None

  def mesh(self):
    # prepare FEniCS mesh and mesh editor
    mesh = df.Mesh()
    mesh_editor = df.MeshEditor()
    mesh_editor.open(mesh, "tetrahedron", 3, 3);
    mesh_editor.init_vertices(self._model.getNumMeshVertices())

    # prepare container for physical IDs
    self._physical_cell_ids = set()
    self._physical_face_ids = set()

    # get mesh vertices
    entities = gp.GEntityVector()
    self._model.getEntities(entities)

    for entity in entities:
      for vertex in entity.mesh_vertices:
        mesh_editor.add_vertex(vertex.getIndex() - 1, df.Point(vertex.x(), vertex.y(), vertex.z()))

    # get mesh cells
    mesh_editor.init_cells(
        reduce(lambda x,y: x+y, map(lambda r: r.tetrahedra.size(), self._model.bindingsGetRegions()))
    )

    i = 0
    for region in self._model.bindingsGetRegions():
      if len(region.physicals) == 0:
        region_id = None
      elif len(region.physicals) == 1:
        region_id = region.physicals[0]
      else:
        raise Exception("Tetrahedra should belong to only one physical domain.")
        
      for tet in region.tetrahedra:
        mesh_editor.add_cell(i, np.array([
          tet.getVertex(0).getIndex() - 1,
          tet.getVertex(1).getIndex() - 1,
          tet.getVertex(2).getIndex() - 1,
          tet.getVertex(3).getIndex() - 1
        ], dtype=np.uintp))

        if region_id is not None:
          mesh.domains().set_marker((i, region_id), 3)
          self._physical_cell_ids.add(region_id)

        i += 1

    # finish mesh creation
    mesh_editor.close()

    # set facet domains
    mesh.init(0, 2)
    for face in self._model.bindingsGetFaces():
      if len(face.physicals) == 0:
        continue
      elif len(face.physicals) == 1:
        face_id = face.physicals[0]
      else:
        raise Exception("Facets should belong to only one physical domain.")

      for tri in face.triangles:
        i = reduce(lambda x,y: x.intersection(y), map(lambda i: set(mesh.topology()(0,2)(tri.getVertex(i).getIndex() - 1)), range(3))).pop()
        mesh.domains().set_marker((i, face_id), 2)
        self._physical_face_ids.add(face_id)

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
