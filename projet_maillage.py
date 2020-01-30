#bibliothèques

import numpy as np
from math import *
import gmsh
import scipy.sparse as sp
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
	M = sp.coo_matrix(data, shape=(len(mesh.triangles),len(mesh.triangles)), dtype=np.float64 )

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

def matrice_rigidite_globale(mesh) : 
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
		K_e = matrice_rigidite_elem(x1, y1, x2, y2, x3, y3)
		for i in range(3):
			I = triangle.Point[i].id
			for j in range(3):
				J = triangle.Point[j].id
				row.append(I)
				col.append(J)
				val.append(K_e[i][j])

	row = np.array(row)
	col = np.array(col)
	val = np.array(val)

	data = (val,(row,col))

	K = sp.coo_matrix(data, shape=(len(mesh.triangles),len(mesh.triangles)), dtype=np.float64 )

	return K



def f(x,y) :
	return (1+2*(pi**2))*sin(pi*x)*sin(pi*y)

def calcul_X(m, triangle, xi, eta):
	x = 0
	y = 0
	for i in range(3):
		phi_ = phi(xi, eta, i+1)
		x += triangle.Point[i].x * phi_
		y += triangle.Point[i].y * phi_

	return x,y

def calcul_B(mesh, order):

	#nombre de point de gauss en fonctionde la precision
	if(order ==1 ):
		size_m = 1
	if(order == 2):
		size_m = 3

	B = np.zeros((len(mesh.points)))
	print(len(mesh.points))
	print(B)

	for triangle in mesh.triangles :
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











# def quadrature(mesh) :






#######SCRIPT######


gmsh.initialize(sys.argv)
filename = "square.msh"
mesh = Mesh()
mesh.gmshToMesh(filename)
#res = matrice_rigidite_globale(mesh)
#print(res.toarray())

res = calcul_B(mesh,1)
print(res)


