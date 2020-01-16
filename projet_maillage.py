#bibliothèques

import numpy as np
from math import *
import gmsh



masse_ref = 1/24 * (np.ones((3,3)) + np.eye(3))
#print(masse_ref)



################GMSH -> RECUPERATION DES DONNEES #############

def gmshToMesh(filename) :

	points_list = []
	seg_list = []
	triangle_list = []

	#ouverture fichier gmsh
	#https://gitlab.onelab.info/gmsh/gmsh/blob/master/api/gmsh.py?fbclid=IwAR2os2AtzsAZa2Y-lOZSi0Jv2EXRSu4meq5OEAz4aBiYyv15iJ0S_YUVFHo
	ierr = c_int()
	lib.gmshOpen(
        	c_char_p(fileName.encode()),
        	byref(ierr))
	print("cool")
	print(ierr.value)
	if ierr.value != 0:
		raise ValueError(
            	"gmshOpen returned non-zero error code: ",
            	ierr.value)
	


	return points_list, seg_list, triangle_list

filename = "square.msh"
mesh = gmsh.merge(filename)
print("tkt\n")
gmshToMesh(filename)


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


#def matrice_masse_globale() : 




#triangle p, 
#def matrice_rigidite(p) :




#return 
