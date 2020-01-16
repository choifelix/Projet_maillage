import gmsh
import sys

gmsh.initialize(sys.argv)
filename = "square.msh"
# model = gmsh.model
# model = gmsh.merge(filename)

# list_point = gmsh.model.mesh.getNodes()

# #for i in range(len(list_point[0])):
# #	print(list_point[0][i]," : ", list_point[1][3*i],", ",list_point[1][3*i+1],", ",list_point[1][3*i+2])

# list_physical_entity = gmsh.model.getPhysicalGroups()
# for elem in list_physical_entity :
# 	list_entity = gmsh.model.getEntitiesForPhysicalGroup(elem[0],elem[1])
# 	for entity in list_entity :
# 		list_elements = gmsh.model.mesh.getElements(elem[0],entity)
# 		print(list_elements)



class Point :
	def __init__(self, index, x, y) :
		self.id = index
		self.x = x
		self.y = y

	def __str__(self):
		return "Point " + str(self.id) + " : (" + str(self.x) + " , " + str(self.y)+")"

class segment :
	def __init__(self, index, tag, Point) :
		self.id = index
		self.tag = tag
		self.Point = Point

	def __str__(self):
		return "Segment : " + str(self.id) + " : " + str(self.Point[0]) + " , " + str(self.Point[1])

	def area(self) :
		return sqrt( (Point[0].x-Point[1].x)**2 + (Point[0].y-Point[1].y)**2)

	def jac(self) :
		return sqrt( (Point[0].x-Point[1].x)**2 + (Point[0].y-Point[1].y)**2)


class triangle :
	def __init__(self, index, tag, Point) :
		self.id = index
		self.tag = tag
		self.Point = Point

	def __str__(self):
		return "Triangle : " + str(self.id) + " : " + str(self.Point[0]) + " , " + str(self.Point[1])


	def area(self) :
		return 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))

	def jac(self) :
		return 2 * area(self)


class Mesh :
	def __init__(self):
		self.points    = []
		self.segments  = []
		self.triangles = []

	def gmshToMesh(self, filename):
		model = gmsh.model
		model = gmsh.merge(filename)

		#------ remplissage points ---------
		list_point = gmsh.model.mesh.getNodes()

		for i in range(len(list_point[0])):
			print(list_point[0][i]," : ", list_point[1][3*i],", ",list_point[1][3*i+1],", ",list_point[1][3*i+2])
			self.points.append(Point(list_point[0][i],list_point[1][3*i],list_point[1][3*i+1]))

		#
		list_physical_entity = gmsh.model.getPhysicalGroups()
		for elem in list_physical_entity :
			list_entity = gmsh.model.getEntitiesForPhysicalGroup(elem[0],elem[1])
			for entity in list_entity :
				list_elements = gmsh.model.mesh.getElements(elem[0],entity)
				
				

point = Point(1,1,0)
print(point)





