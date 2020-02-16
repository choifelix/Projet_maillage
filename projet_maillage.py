#bibliothèques

import numpy as np
from math import *
import gmsh
import scipy.sparse as sp
import scipy.sparse.linalg as linalg
from test import *
import matplotlib.pyplot as plt
import matplotlib.tri as mtri


###MASSE DE REF###
masse_ref = 1/24 * (np.ones((3,3)) + np.eye(3))



#####################FONCTION DE FORMES TRIANGLE DE REF##################



def phi(x,y,i) :
	#si on se situe au sommet 1
	if i == 0 :
		phi = 1-x-y 
	if i == 1 : 
		phi = x	
	if i == 2 : 
		phi = y
	if (i!=0 and i!=1 and i!=2) : 
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
	if i == 0 : 
		grad_phi.append(-1)
		grad_phi.append(-1)  	
	if i == 1 : 
		grad_phi.append(1)
		grad_phi.append(0)
	if i == 2 : 
		grad_phi.append(0)
		grad_phi.append(1)
	if (i!=0 and i!=1 and i!=2) : 
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

 
def matrice_masse_globale(mesh, physical_tag, triplets) : 

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
					triplets.append( I,J,M_e[i][j] ) 


	return triplets




def matrice_rigidite_elem(x1, y1, x2, y2, x3, y3):
	K = np.zeros((3,3), dtype=np.float64)
	Tp = area(x1, y1, x2, y2, x3, y3)
	B_p = Bp(x1, y1, x2, y2, x3, y3)
	produit_B = np.dot(B_p.transpose(),B_p)
	print(produit_B)

	for i in range(3):
		for j in range(3):
			# temp = np.dot(gradient_phi(j).transpose(), produit_B )
			# K[i][j] = Tp * np.dot(temp, gradient_phi(i))
			# print(np.dot( np.dot(gradient_phi(j),produit_B),gradient_phi(i).transpose() ))

			K[i][j] = Tp * np.dot( np.dot(gradient_phi(j),produit_B),gradient_phi(i).transpose() )

	return K

def matrice_rigidite_globale(mesh, physical_tag, triplets) : 

	for triangle in mesh.triangles :
		if(triangle.tag == physical_tag):

			# -- coordonnées des points du triangle --
			x1 = triangle.Point[0].x
			y1 = triangle.Point[0].y
			x2 = triangle.Point[1].x
			y2 = triangle.Point[1].y
			x3 = triangle.Point[2].x
			y3 = triangle.Point[2].y

			# -- matrice elementaire --
			K_e = matrice_rigidite_elem(x1, y1, x2, y2, x3, y3)

			# -- remplissage matrice globale --
			for i in range(3):
				I = triangle.Point[i].id

				for j in range(3):
					J = triangle.Point[j].id
					triplets.append( I,J,K_e[i][j] )


	return triplets




def f(x,y) :
	return (1+2*(pi**2))*sin(pi*x)*sin(pi*y)

def calcul_X(m, triangle, xi, eta):
	x = 0
	y = 0
	for i in range(3):
		phi_ = phi(xi, eta, i)
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
			x1 = triangle.Point[0].x
			y1 = triangle.Point[0].y
			x2 = triangle.Point[1].x
			y2 = triangle.Point[1].y
			x3 = triangle.Point[2].x
			y3 = triangle.Point[2].y
			detJ = 2 * area(x1, y1, x2, y2, x3, y3)

			gauss = triangle.gaussPoint(order)

			for i in range(3):
				I = triangle.Point[i].id
				I = int(I) 
				for m in range(size_m):
					xi    = gauss[m][0]
					eta   = gauss[m][1]
					poids = gauss[m][2]
					#print(poids)
					(x,y) = calcul_X(m, triangle, xi, eta)
					#print(f(x,y))
					B[I] += poids * f(x,y) * phi(xi,eta,i)
					#print(B[I])
			B[I] = detJ * B[I]

	return B










#dim = 1 : segments bords
def Dirichlet(mesh, dim, physical_tag, g, triplets, B) :


	if(dim == 1) :	
		for segment in mesh.segments: #chaque segment
			if(segment.tag == physical_tag):
				for point in segment.Point: #chaque points du segment
					I = point.id

					#------ chercher les lignes d'indices I et les mettre a 0 -----
					for i in range(len(triplets.data[0])):
						if (triplets.data[1][0][i] == I):
							# val[i] = 0.
							triplets.data[0][i] = 0


					#------ ajout de 1 sur l'elem de diago I,I -> (I,I,1) -----
					triplets.append( I,I,1 )


					#------ modif de B pour retouver u = B = g sur les noeuds du bords -----
					B[I] = g


	return triplets, B




def solve(A,B):
	

	u = linalg.spsolve(A,B)
	return u


###g resolution U ref
def g_sol(x,y) :
	return np.sin(np.pi*x)*np.sin(np.pi*y)






def test():

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


	triplets = Triplet()
	triplets = matrice_rigidite_globale(mesh, physical_tag_Triangle, triplets)
	triplets = matrice_masse_globale(mesh, physical_tag_Triangle, triplets)
	A = sp.coo_matrix(triplets.data ).tocsr()
	size, sizecol = np.shape(A)
	print(A)



	# -- calcul de B --
	B = calcul_B(mesh, order, physical_tag_Triangle, size)


	# -- Pose les conditions aux bords --
	triplets, B =Dirichlet(mesh, dim_dirich, physical_tag_segment, g, triplets, B)
	# print(triplets)
	# -- reshape A et B --


	A = sp.coo_matrix(triplets.data ).tocsr()
	print(A)


	# B = np.array(B)


	# -- resolution --
	U = solve(A,B)


	# Visualisation
	x= [pt.x for pt in mesh.points]

	y= [pt.y for pt in mesh.points]

	#triangulation avec les coo des points
	triang = mtri.Triangulation(x, y)

	#id des triangles
	connectivity=[]
	for tri in mesh.triangles:
	  connectivity.append([ p.id for p in tri.Point]) 

	print(len(x))
	print(len(connectivity))

	# print(U)
	####U APPROX NOTRE SOLUTION
	plt.tricontourf(x,y,connectivity, U, 12)
	plt.colorbar()
	plt.show()


	############

	### U de REFERENCE
	Uref = np.zeros((len(mesh.points),))

	for pt in mesh.points:
	  I = int(pt.id) 
	  Uref[I] = g_sol(pt.x, pt.y)


	# # print(Uref)

	# plt.tricontourf(x,y,connectivity, Uref, 12)
	# plt.colorbar()
	# plt.show()


	return A,U,B

def affiche_matrice():
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


	K_triplets = Triplet()
	M_triplets = Triplet()
	K_triplets = matrice_rigidite_globale(mesh, physical_tag_Triangle, K_triplets)
	M_triplets = matrice_masse_globale(mesh, physical_tag_Triangle, M_triplets)
	K = sp.coo_matrix(K_triplets.data ).tocsr()
	M = sp.coo_matrix(M_triplets.data ).tocsr()



	print(K.toarray())
	print(M.toarray())

	return K,M


#######SCRIPT######


# ------ initialisation du maillage -----
gmsh.initialize(sys.argv)
filename = "square.msh"
mesh = Mesh()
mesh.gmshToMesh(filename)


K,M = affiche_matrice()


A,U,B = test()




