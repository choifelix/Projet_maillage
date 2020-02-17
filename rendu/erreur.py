import numpy as np
from math import *
import gmsh
from projet_maillage import *
import matplotlib.pyplot as plt

import gmsh
import sys

pas_x = []
erreur = []
diff = []

for pas in np.arange(0.1,1,0.1):


	gmsh.initialize(sys.argv)
	gmsh.option.setNumber("General.Terminal", 1)
	gmsh.option.setNumber("Mesh.CharacteristicLengthMin", pas);
	gmsh.option.setNumber("Mesh.CharacteristicLengthMax", pas);
	# Model
	model = gmsh.model
	model.add("Square")
	# Rectangle of (elementary) tag 1
	factory = model.occ
	factory.addRectangle(0,0,0, 1, 1, 1)
	# Sync
	factory.synchronize()
	# Physical groups
	gmsh.model.addPhysicalGroup(1, [1, 2, 3, 4], 1)
	# gmsh.model.addPhysicalGroup(1, [2,3,4], 2)
	gmsh.model.addPhysicalGroup(2, [1], 10)
	# Mesh (2D)
	model.mesh.generate(2)
	# ==============

	# ==============
	#Save mesh
	gmsh.write("square.msh")
	# Finalize GMSH
	gmsh.finalize()



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


	# -- calcul de A --
	K = matrice_rigidite_globale(mesh, physical_tag_Triangle)
	M = matrice_masse_globale(mesh, physical_tag_Triangle)

	A = K+M

	A_tmp = A.toarray()

	size, sizecol = np.shape(A_tmp)
	B = calcul_B(mesh, order, physical_tag_Triangle, size)


	# -- Pose les conditions aux bords --
	(A, B) =Dirichlet(mesh, dim_dirich, physical_tag_segment, g, A, B)

	# -- reshape A et B --
	A = A.toarray()

	B = np.array(B)

	# -- resolution --
	U = solve(A,B)

	# Visualisation
	x= [pt.x for pt in mesh.points]

	y= [pt.y for pt in mesh.points]


	U = U.flatten()

	############

	### U de REFERENCE
	Uref = np.zeros((len(mesh.points),))
	for pt in mesh.points:
	  I = int(pt.id) - 1
	  Uref[I] = g_sol(pt.x, pt.y)



	#####ERREUR RELATIVE

	for i in range(len(U)) :

		diff.append((U[i] - Uref[i])**2)


	## erreur relative
	err = sqrt(sum(diff)/sum(U**2))

	#erreur
	#err = sqrt(sum(diff))
	erreur.append(err)

	pas_x.append(pas)



###VISU DES ERREURS
plt.plot(pas_x, erreur)
plt.show()


###LOG LOG
plt.loglog(pas_x, erreur)
plt.show()

