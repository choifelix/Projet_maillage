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
		#return "Point " + str(self.id) + " : (" + str(self.x) + " , " + str(self.y)+")"
		return ""

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
	def __init__(self, index, tag, Point_list) :
		self.id = index
		self.tag = tag
		self.Point = Point_list

	def __str__(self):
		return "Triangle : " + str(self.id) + " : " + str(self.Point[0]) + " , " + str(self.Point[1])


	# def area(self) :
	# 	return 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))

	def gaussPoint(self,order=2):
		res = []
		if order == 1 : 
			#xi, eta, poids(omega)
			res.append([1./3., 1./3., 1./6.])
		if order == 2 :
			res.append([1./6., 1./6., 1./6.])
			res.append([4./6., 1./6., 1./6.])
			res.append([1./6., 4./6., 1./6.])

		if (order!=1 and order!=2 ) : 
			print("l ordre doit etre compris entre 1 et 2 -> ordre au dela non implementé!")
			return


		return res

	# def jac(self) :
	# 	return 2 * area(self)


class Mesh :
	def __init__(self):
		self.points    = []
		self.segments  = []
		self.triangles = []
		self.physical_tag = []



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
		size_seg = 0
		size_tri = 0

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

					# -- nb de points des segments - sans 'doublons'  --
					#for j in range(1,len(list_elements)-1) :
					#	size_seg = size_seg + len(list_elements[j][0])
					# -- -------------------------------------------- --
					
					for i in range(len(list_elements[1][0])):
						#print(list_elements[1][0][i])
						#print(list_elements[2][0][i])
						s1 = self.getPoint_id(list_elements[2][0][2*i])
						s2 = self.getPoint_id(list_elements[2][0][2*i+1])

						self.segments.append(Segment(list_elements[1][0][i], elem[1], [s1, s2] ))

				if(elem[0] == 2): #triangle, list_elements[1][0] -> indices des triangles
					# list_elements[2][0][3*i] -> indices des sommets composant le triangles -> 3 sommets a la suite dans la liste

					
					# -- nb de points des triangles - sans 'doublons' --
					#unique_list = [] 
					# traverse for all elements 
					#for x in list_elements[2][0]: 
					# check if exists in unique_list or not 
					#	if x not in unique_list: 
					#		unique_list.append(x)

					#size_tri = len(unique_list)
					# -- -------------------------------------------- --



					for i in range(len(list_elements[1][0])):
						#list_elements[2][0][3*i],list_elements[2][0][3*i+1],list_elements[2][0][3*i+2]]
						s1 = self.getPoint_id(list_elements[2][0][3*i])
						s2 = self.getPoint_id(list_elements[2][0][3*i+1])
						s3 = self.getPoint_id(list_elements[2][0][3*i+2])
						self.triangles.append(Triangle(list_elements[1][0][i], elem[1], [s1, s2, s3] )) #eleme[1] c'est le physical tag
			self.physical_tag.append((elem[0],elem[1]))
				# print(list_elements)




