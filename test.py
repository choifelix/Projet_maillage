import gmsh
import sys

gmsh.initialize(sys.argv)
filename = "square.msh"
class Point :
	def __init__(self, index, x, y) :
		self.id = index
		self.x = x
		self.y = y


	def __str__(self):
		return "Point " + str(self.id) + " : (" + str(self.x) + " , " + str(self.y)+")"

class Segment :
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


class Triangle :
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



	def getPoint_id(self, index) :

		#chercher 1 par un l id
		for i in self.points :
			if i.id == index :
				return i

		print("Not found")
		return 0

	def gmshToMesh(self, filename):
		model = gmsh.model
		model = gmsh.merge(filename)

		#------ remplissage points ---------
		list_point = gmsh.model.mesh.getNodes()

		for i in range(len(list_point[0])):
			#print(list_point[0][i]," : ", list_point[1][3*i],", ",list_point[1][3*i+1],", ",list_point[1][3*i+2])
			self.points.append(Point(list_point[0][i],list_point[1][3*i],list_point[1][3*i+1]))

		#
		list_physical_entity = gmsh.model.getPhysicalGroups()
		for elem in list_physical_entity :

			# elem = [dim,tag_physicalEntity]
			list_entity = gmsh.model.getEntitiesForPhysicalGroup(elem[0],elem[1])
			for entity in list_entity :
				# entity = [tag_entity]
				list_elements = gmsh.model.mesh.getElements(elem[0],entity)

				# list_elements = [[tag],[elements]]
				if(elem[0] == 1): #segment
		
					for i in range(len(list_elements[1][0])):
						#print(list_elements[1][0][i])
						#print(list_elements[2][0][i])

						self.segments.append(Segment(list_elements[1][0][i], entity, [ list_elements[2][0][2*i],list_elements[2][0][2*i+1] ] ))

				if(elem[0] == 2): #triangle, list_elements[1][0] -> indices des triangles
					# list_elements[2][0][3*i] -> indices des sommets composant le triangles -> 3 sommets a la suite dans la liste
					for i in range(len(list_elements[1][0])):
						#list_elements[2][0][3*i],list_elements[2][0][3*i+1],list_elements[2][0][3*i+2]]
						s1 = self.getPoint_id(list_elements[2][0][3*i])
						s2 = self.getPoint_id(list_elements[2][0][3*i+1])
						s3 = self.getPoint_id(list_elements[2][0][3*i+2])
						self.triangles.append(Triangle(list_elements[1][0][i], entity, [s1, s2, s3] ))

				# print(list_elements)

mesh = Mesh()
mesh.gmshToMesh(filename)


