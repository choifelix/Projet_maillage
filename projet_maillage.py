#bibliothèques

import numpy as np
from math import *
import gmsh
import scipy as sp
from test import *


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

res = phi(1,1,2)
print(res)


#######################GRADIENT DE PHI########################



def gradient_phi(i) :

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
	return grad_phi




res = gradient_phi(1)
print(res)


################MASSE ELEMENTAIRE##############


def matrice_masse_elem(x1, y1, x2, y2, x3, y3) : 
	#pour i, j = 1,2,3
	det = 2 * 0.5* abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
	M =  det * masse_ref
	return M


def matrice_masse_globale(mesh) : 
	row = []
	col = []
	val = []

	for triangle in mesh.triangles :

		x1 = triangle.Point[0].x
		y1 = triangle.Point[0].y
		x2 = triangle.Point[1].x
		y2 = triangle.Point[1].y
		x3 = triangle.Point[2].x
		y3 = triangle.Point[2].y
		M_e = matrice_masse_elem(x1, y1, x2, y2, x3, y3)

		for i in range(3):
			I = triangle.Point[i]
			for j in range(3):
				J = triangle.Point[i]
				row.append(I)
				col.append(J)
				val.append(M_e[i][j])

	data = (val,(row,col))

	M = sp.coo_matrix(data)

	return M


#######SCRIPT######


gmsh.initialize(sys.argv)
filename = "square.msh"
mesh = Mesh()
mesh.gmshToMesh(filename)
res = matrice_masse_globale(mesh)

print(res.toarray())


#triangle p, 
#def matrice_rigidite(p) :




#return 
