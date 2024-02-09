# PROJET
 
Ce repo contient
- un dossier data contenant les instances
- un dossier results contenant nos resultats
- un dossier txtFiles contenant l'output du solveur de la résolution de certaines instances
- un dossier utils contenant :
     - les fonctions de parsing, de simplification des graphes et de calcul des partitions des sommets et des arêtes (variables d'agrégation)
     - un fichier utils_heuristic (déplacement vers un voisin, test sur les chemins, calcul des distances et poids robustes)
     - un fichier voisinages (recherche de voisinages améliorants)
 
 Des fichiers Julia pour chaque méthode
 - PLNE_compacte.jl (méthode statique)
 - branch_cut.jl (méthode Branch&Cut)
 - constrSol.jl qui construit une solution réalisable avec la variante de l'algorithme de Dirjska
 - heurDich.jl implémentant l'heuristique dichotomique
 - heurVois.jl qui implémente la recherche locale de voisinages
 - pipeline.jl pour tester les résultats
 - plans_coupants.jl (méthode Plans coupants)
 - PLNE_dual.jl (méthode Plans coupants)
 - des fichiers sandbox et des notebooks où nous avons fait des test

   Plusieurs méthodes ont des paramètres pour utiliser des variables d'agrégation, perturber l'instance par exemple.


 Un fichier Python analyse.py d'analyse des résultats. Un fichier sandbox.py est présent dans la branche exploitation-resultats pour les graphiques et comparaison des méthodes
