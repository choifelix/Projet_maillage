#bibliothèques

import numpy as np
from math import *
import gmsh
import scipy.sparse as sp
from test import *
import matplotlib.pyplot as plt


###MASSE DE REF###
masse_ref = 1/24 * (np.ones((3,3)) + np.eye(3))



#####################FONCTION DE FORMES TRIANGLE DE REF##################



def phi(x,y,i) :
	#si on se situe au sommet 1
	if i == 1 :
		phi = 1-x-y 
	if i == 2 : 
		phi = x	
	if i == 3 : 
		phi = y
	if (i!=1 and i!=2 and i!=3) : 
		print("le dernier argument doit etre compris entre 1 et 3 -> numero du sommet!")
		return 

	return  phi

#res = phi(1,1,2)
#print(res)


#######################GRADIENT DE PHI########################



def gradient_phi(i) :
 	#triangle de reference P1
	#gradient vecteur
	grad_phi = []
 
	#sommet 1
	#grad = (-1, -1)
	if i == 1 : 
		grad_phi.append(-1)
		grad_phi.append(-1)  	
	if i == 2 : 
		grad_phi.append(1)
		grad_phi.append(0)
	if i == 3 : 
		grad_phi.append(0)
		grad_phi.append(1)
	if (i!=1 and i!=2 and i!=3) : 
		print("l argument doit etre compris entre 1 et 3 -> numero du sommet!")
		return
	grad_phi = np.array([grad_phi], dtype=np.float64)
	return grad_phi




#res = gradient_phi(1)
#print(res)
def area(x1, y1, x2, y2, x3, y3):
	return 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))


def Bp(x1, y1, x2, y2, x3, y3):
	b1 = y3-y1
	b2 = y1-y2
	b3 = x1-x3
	b4 = x2-x1
	detJ = 2 * area(x1, y1, x2, y2, x3, y3)
	return np.mat([[b1, b2], [b3, b4]], dtype=np.float64)







################MASSE ELEMENTAIRE##############


def matrice_masse_elem(x1, y1, x2, y2, x3, y3) : 
	#pour i, j = 1,2,3
	det = 2 * 0.5* abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
	M =  det * masse_ref
	return M

 
def matrice_masse_globale(mesh, physical_tag) : 
	row = []
	col = []
	val = []

	for triangle in mesh.triangles :
		if(triangle.tag == physical_tag):
			x1 = triangle.Point[0].x
			y1 = triangle.Point[0].y
			x2 = triangle.Point[1].x
			y2 = triangle.Point[1].y
			x3 = triangle.Point[2].x
			y3 = triangle.Point[2].y
			M_e = matrice_masse_elem(x1, y1, x2, y2, x3, y3)
			for i in range(3):
				I = triangle.Point[i].id
				for j in range(3):
					J = triangle.Point[j].id
					row.append(I)
					col.append(J)
					val.append(M_e[i][j])

	row = np.array(row)
	col = np.array(col)
	val = np.array(val)

	data = (val,(row,col))
	M = sp.coo_matrix(data, shape=(int(max(row)+1),int(max(row)+1)), dtype=np.float64 )

	return M




def matrice_rigidite_elem(x1, y1, x2, y2, x3, y3):
	K = np.zeros((3,3), dtype=np.float64)
	Tp = area(x1, y1, x2, y2, x3, y3)
	B_p = Bp(x1, y1, x2, y2, x3, y3)
	produit_B = B_p.transpose()*B_p

	for i in range(3):
		for j in range(3):
			temp = np.dot(gradient_phi(j+1), produit_B )
			K[i][j] = Tp * np.dot(temp, gradient_phi(i+1).T).item()

	return K

def matrice_rigidite_globale(mesh, physical_tag) : 
	row = []
	col = []
	val = []


	for triangle in mesh.triangles :
		if(triangle.tag == physical_tag):
			x1 = triangle.Point[0].x
			y1 = triangle.Point[0].y
			x2 = triangle.Point[1].x
			y2 = triangle.Point[1].y
			x3 = triangle.Point[2].x
			y3 = triangle.Point[2].y
			K_e = matrice_rigidite_elem(x1, y1, x2, y2, x3, y3)
			for i in range(3):
				I = triangle.Point[i].id
				#recup max id pour la taille de matrice
				for j in range(3):
					J = triangle.Point[j].id
					row.append(I)
					col.append(J)
					val.append(K_e[i][j])

	row = np.array(row)
	col = np.array(col)
	data = np.array(val)
	#data = (val,(row,col))

	K = sp.coo_matrix((data, (row,col)), shape=(int(max(row)+1),int(max(row)+1)), dtype=np.float64 )
	#K = sp.csr_matrix((data, (row,col)), shape=(int(max(row))+1,int(max(row))+1), dtype=np.float64 )


	return K



# def Sum(K,M):
# 	A = M
	

# 	for i_k in range(len(K[0])):
# 		for i in range(len(A[0])):
# 			if( K[1][0][i_k] == A[1][0][i] and K[1][1][i_k] == A[1][1][i]):
# 				A[0][i] += K[0][i_k]
# 			else:
# 				A[0].append(k[0][i_k])
# 				A[1][0].append(k[1][0][i_k])
# 				A[1][1].append(k[1][1][i_k])
# 	return A



def f(x,y) :
	return (1+2*(pi**2))*sin(pi*x)*sin(pi*y)

def calcul_X(m, triangle, xi, eta):
	x = 0
	y = 0
	for i in range(3):
		phi_ = phi(xi, eta, i+1)
		x += triangle.Point[i].x * phi_  #triangle.Point[i].x = S(i) en x
		y += triangle.Point[i].y * phi_  #triangle.Point[i].y = S(i) en y

	return x,y

def calcul_B(mesh, order, physical_tag, size):

	#nombre de point de gauss en fonctionde la precision
	if(order ==1 ):
		size_m = 1
	if(order == 2):
		size_m = 3

	B = np.zeros((size,1))


	for triangle in mesh.triangles :
		if(triangle.tag == physical_tag):

			gauss = triangle.gaussPoint(order)
			for i in range(3):
				I = triangle.Point[i].id
				I = int(I) - 1
				for m in range(size_m):
					xi    = gauss[m][0]
					eta   = gauss[m][1]
					poids = gauss[m][2]
					#print(poids)
					(x,y) = calcul_X(m, triangle, xi, eta)
					#print(f(x,y))
					B[I] += poids * f(x,y) * phi(xi,eta,i+1)
					#print(B[I])

	return B










#dim = 1 : segments bords
def Dirichlet(mesh, dim, physical_tag, g, triplets, B) :

	row,col,val = sp.find(triplets)
	row=row.tolist()
	col=col.tolist()
	val=val.tolist()
	# size = len(row) # matrice carre
	if(dim == 1) :	
		for segment in mesh.segments: #chaque segment
			if(segment.tag == physical_tag):
				for point in segment.Point: #chaque points du segment
					I = point.id

					#------ chercher les lignes d'indices I et les mettre a 0 -----
					for i in range(len(row)):
						if (row[i] == I):
							val[i] = 0.

					#------ ajout de 1 sur l'elem de diago I,I -> (I,I,1) -----
					# val.append(1)
					# row.append(I)
					# col.append(I)

					val.append(1.)
					row.append(I)
					col.append(I)



					#------ modif de B pour retouver u = B = g sur les noeuds du bords -----
					B[I] = g

	

	row = np.array(row)
	col = np.array(col)
	val = np.array(val)


	data = (val,(row,col))

	A = sp.coo_matrix(data, shape=(int(max(row)+1),int(max(row)+1)), dtype=np.float64 )

	return (A, B)




def solve(A,B):
	

	u = np.linalg.solve(A,B)
	return u











#######SCRIPT######

# ------ initialisation du maillage -----
gmsh.initialize(sys.argv)
filename = "square.msh"
mesh = Mesh()
mesh.gmshToMesh(filename)


order      = 2
dim_dirich = 1
g          = 0

# ------ construction du probleme  ------


# le choix du physical tag depend du .msh -> ici square.msh
for (dim, physical_tag) in mesh.physical_tag:
	if(dim == 2):
		physical_tag_Triangle = physical_tag
	if(dim == 1):
		physical_tag_segment = physical_tag

#print(physical_tag_Triangle)
#print(physical_tag_segment)

# -- calcul de A --
K = matrice_rigidite_globale(mesh, physical_tag_Triangle)
M = matrice_masse_globale(mesh, physical_tag_Triangle)

A = K+M
A_tmp = A.toarray()
#print(np.linalg.det(A_tmp))

size, sizecol = np.shape(A_tmp)

# print(np.linalg.det(A))


# -- calcul de B --
B = calcul_B(mesh, order, physical_tag_Triangle, size)


# -- Pose les conditions aux bords --
(A, B) =Dirichlet(mesh, dim_dirich, physical_tag_segment, g, A, B)
print(A)

# -- reshape A et B --
A = A.toarray()
print(A)
A = np.delete(A,0,0)
A = np.delete(A,0,1)

print(A)

B = np.array(B)
B = np.delete(B,0,0)
print(B)

# A_inv = np.linalg.det(A)
# print("inverse de A:")
# print(A_inv)

# -- resolution --
U = solve(A,B)
print("solution U :")
print(U)



# Visualisation
x= [pt.x for pt in mesh.points]
y= [pt.y for pt in mesh.points]
connectivity=[]
for tri in mesh.triangles:
  connectivity.append([ p.id for p in tri.Point]) 

plt.tricontourf(x, y, connectivity, U)
plt.colorbar()
plt.show()

### U de référence
# Uref = np.zeros((msh.Npts,))
# for pt in msh.points:
#   I = int(pt.id)
#   Uref[I] = g(pt.x, pt.y)
# plt.tricontourf(x, y, connectivity, Uref, 12)
# plt.colorbar()
# plt.show()

