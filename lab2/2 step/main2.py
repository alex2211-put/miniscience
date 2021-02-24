import vtk
import numpy as np
import gmsh
import math
import os


# Класс расчётной сетки
class CalcMesh:

    # Конструктор сетки, полученной из stl-файла
    def __init__(self, nodes_coords, tetrs_points):
        # 3D-сетка из расчётных точек
        # Пройдём по узлам в модели gmsh и заберём из них координаты
        self.nodes = np.array([nodes_coords[0::3], nodes_coords[1::3], nodes_coords[2::3]])

        # Модельная скалярная величина распределена как-то вот так
        # np.power(self.nodes[0, :], 2) + np.power(self.nodes[1, :], 2)
        self.smth = (np.sin(np.power(self.nodes[0, :] + 1, 3) * 40) * 4000 + 5000) + (
                    np.cos(np.power(self.nodes[1, :] + 1, 3) * 50) * 3000 + 1000)
        w = np.array([1, 5, 2])

        self.velocity = np.cross(np.transpose(self.nodes), w).transpose()

        self.tetrs = np.array([tetrs_points[0::4], tetrs_points[1::4], tetrs_points[2::4], tetrs_points[3::4]])
        self.tetrs -= 1

    def move(self, tau, t):
        w = 60
        k = 25
        self.nodes += self.velocity*tau
        d = self.nodes.copy()
        self.smth = ((np.sin(np.power(self.nodes[0, :] + 1, 3) * 40) * 5000) + (
                    np.cos(np.power(self.nodes[1, :] + 1, 3) * 50) * 5000)) * np.cos(w * t +
                                                                                     k * np.power(
            np.dot(d[0], d[0].transpose()) + np.dot(d[1], d[1].transpose()) +
            np.dot(d[2], d[2].transpose()), 0.5))


    def snapshot(self, snap_number):
        unstructuredGrid = vtk.vtkUnstructuredGrid()
        points = vtk.vtkPoints()

        smth = vtk.vtkDoubleArray()
        smth.SetName("smth")

        vel = vtk.vtkDoubleArray()
        vel.SetNumberOfComponents(3)
        vel.SetName("vel")

        for i in range(0, len(self.nodes[0])):
            points.InsertNextPoint(self.nodes[0, i], self.nodes[1, i], self.nodes[2, i])
            smth.InsertNextValue(self.smth[i])
            vel.InsertNextTuple((self.velocity[0, i], self.velocity[1, i], self.velocity[2, i]))

        unstructuredGrid.SetPoints(points)

        unstructuredGrid.GetPointData().AddArray(smth)
        unstructuredGrid.GetPointData().AddArray(vel)

        for i in range(0, len(self.tetrs[0])):
            tetr = vtk.vtkTetra()
            for j in range(0, 4):
                tetr.GetPointIds().SetId(j, self.tetrs[j, i])
            unstructuredGrid.InsertNextCell(tetr.GetCellType(), tetr.GetPointIds())

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetInputDataObject(unstructuredGrid)
        writer.SetFileName("cake3d-step-" + str(snap_number) + ".vtu")
        writer.Write()


gmsh.initialize()

try:
    path = os.path.dirname(os.path.abspath(__file__))
    gmsh.merge(os.path.join(path, 'cake.stl'))
except:
    print("Could not load STL mesh: bye!")
    gmsh.finalize()
    exit(-1)

angle = 5
forceParametrizablePatches = False
includeBoundary = True
curveAngle = 180
gmsh.model.mesh.classifySurfaces(angle * math.pi / 180., includeBoundary, forceParametrizablePatches,
                                 curveAngle * math.pi / 180.)
gmsh.model.mesh.createGeometry()

# Зададим объём по считанной поверхности
s = gmsh.model.getEntities(2)
l = gmsh.model.geo.addSurfaceLoop([s[i][1] for i in range(len(s))])
gmsh.model.geo.addVolume([l])

gmsh.model.geo.synchronize()

# Зададим мелкость желаемой сетки
funny = False
f = gmsh.model.mesh.field.add("MathEval")
if funny:
    gmsh.model.mesh.field.setString(f, "F", "2*Sin((x+y)/5) + 3")
else:
    gmsh.model.mesh.field.setString(f, "F", "0.9")
gmsh.model.mesh.field.setAsBackgroundMesh(f)

# Построим сетку
gmsh.model.mesh.generate(3)

# Теперь извлечём из gmsh данные об узлах сетки
nodeTags, nodesCoord, parametricCoord = gmsh.model.mesh.getNodes()

# И данные об элементах сетки тоже извлечём, нам среди них нужны только тетраэдры, которыми залит объём
GMSH_TETR_CODE = 4
tetrsNodesTags = None
elementTypes, elementTags, elementNodeTags = gmsh.model.mesh.getElements()
for i in range(0, len(elementTypes)):
    if elementTypes[i] != GMSH_TETR_CODE:
        continue
    tetrsNodesTags = elementNodeTags[i]

if tetrsNodesTags is None:
    print("Can not find tetra data. Exiting.")
    gmsh.finalize()
    exit(-2)

print("The model has %d nodes and %d tetrs" % (len(nodeTags), len(tetrsNodesTags) / 4))

for i in range(0, len(nodeTags)):
    assert (i == nodeTags[i] - 1)
assert (len(tetrsNodesTags) % 4 == 0)

mesh = CalcMesh(nodesCoord, tetrsNodesTags)
mesh.snapshot(0)
tau = 0.0004
t = 0
for i in range(1, 150):
    t += tau
    mesh.move(tau, t)
    mesh.snapshot(i)
    print('.', end='')
gmsh.finalize()
