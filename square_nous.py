# File square.py
import gmsh
import sys


#print('vous avez choisi un pas de maillage egal a ', str(sys.argv[1]))


gmsh.initialize(sys.argv)
gmsh.option.setNumber("General.Terminal", 1)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.1);
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.1);
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