# Projet_maillage
Approximation numérique par méthode des éléments finis


1. Dependances

Numpy 
Matplotlib
Scipy
Math
Gmsh
Sys

2. Utilisation

Pour lancer une visualisation, il faut avoir un maillage (fichier.msh) - nous avons pris ici l'exemple d'un carré (fichier square_nous.py).
Dans le fichier erreur.py, le pas du maillage varie entre 0.1 et 1 avec un pas valant 0.1 (permettant d'analyser l'erreur relative entre la valeur U calculée et la valeur de reference). Pour visualiser ces courbes d'erreur, veuillez lancer le programme erreur.py.
Pour ce qui est de la visualisation de la solution U trouvee, il faut lancer projet_maillage.py, programme appelant maillage.py (definissant les classes tq triangle, segment etc).

3. Details sur les fichiers
 
projet_maillage.py : fichier qui construit puis resout le probleme à partir d'un maillage gmsh existant. Ce fichier contient les fonctions de creations des differentes matrices et de la resolutions, puis plus loin du script correspondant à la mise en oeuvre de la resolution. C'est notamment ici que sont definies les variables comme les physical_tag ou encore la valeur de g(=0)

erreur.py : fait une boucle sur les differents pas (h), construit les maillages avec gmsh et resout le probleme sur ce maillage, l'erreur est ensuite calculée à partir de la solution excate. Les erreurs sont recuperées et visualisées.

4. Résultats

Le résultat que nous obtenons ne semble etre satisfaisant (reste nulle sur les bords - les conditions de Dirichlets sont satisfaites cependant). L'écart entre les valeur de U que nous calculons et les valeurs de référence est élevé. Néanmoins, l'erreur augmente de facon linéaire lorsque le pas augmente, résultat qui semble correct. 
Le calcul de A semble etre la cause de nos erreurs de calculs. En effet, les valeurs sont 






