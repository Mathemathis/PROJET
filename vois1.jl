include("utils/parsing.jl")
include("utils/utils_heuristic.jl")
include("constrSol.jl")

"""Garde un nombre de noeud constants
chemin[i-1] -> chemin[i] = old_noeud -> chemin[i+1]
chemin[i-1] -> nv_noeud -> chemin[i+1]"""

function nvResPoids(i_ph_dec, i_to_i_ph_dec, sum_poids, nv_noeud, old_noeud, d2, ph, p, longueur)
    """test si le nouveau voisin est admissible pour les poids des sommets"""
    i_lim = min(div(d2, 2), longueur)
    i_old_noeud = i_to_i_ph_dec[old_noeud]
    nv_poids = sum_poids + p[nv_noeud] - p[old_noeud]
    if i_old_noeud <= i_lim # nv_noeud compte dans le poids
        nv_poids -= (2*ph[old_noeud]) # on enleve ce poids
        if ((i_lim == longueur) || (ph[nv_noeud] >= ph[i_ph_dec[(i_lim + 1)]])) # attention peut être placé juste à droite
            nv_poids += (2*ph[nv_noeud]) # ajout du noeud # a droite de i_lim reste pareil
        else # le nouveau noeud sature est a droite de i_lim
            nv_poids +=  2* ph[i_ph_dec[(i_lim + 1)]] # on sature a droite
            nv_poids -= ph[i_ph_dec[(i_lim + 1)]]* (d2 - 2 * i_lim)  # on enleve la variation du sommet à droite
            if (i_lim + 2 > longueur) || ph[nv_noeud] >= ph[i_ph_dec[(i_lim + 2)]] # deux cas ou nv_noeud devient le noeud partiellement sature
                nv_poids += (ph[nv_noeud]* (d2 - 2 * i_lim))
            else 
                nv_poids += (ph[i_ph_dec[(i_lim + 2)]]* (d2 - 2 * i_lim))
            end
        end 
    else  # ancien noeud pas enleve
        if ((i_lim > 0) && (ph[nv_noeud] > ph[i_ph_dec[i_lim]])) # saturation nv noeud
            nv_poids += (2*ph[nv_noeud]) # on sature le noeud
            nv_poids -= (2*ph[i_ph_dec[i_lim]]) # on enleve le noeud limite
            nv_poids +=(ph[i_ph_dec[i_lim]] * (d2 - 2 * i_lim)) # saturation partielle du noeud limite
            nv_poids -= (ph[i_ph_dec[(i_lim+1)]]* (d2 - 2 * i_lim)) # i_lim +1 existe forcement. Possible que i_lim =1 = i_old_noeud
        else 
            if ((i_old_noeud == i_lim+1) || (ph[nv_noeud] >  ph[i_ph_dec[(i_lim+1)]])) # il est placé juste à droite
                if (i_lim+2 <= longueur) && (ph[nv_noeud] <=  ph[i_ph_dec[(i_lim+2)]]) 
                    nv_poids += (ph[i_ph_dec[(i_lim+2)]]* (d2 - 2 * i_lim)) 
                else
                    nv_poids += (ph[nv_noeud] * (d2 - 2 * i_lim))
                end
                nv_poids -= (ph[i_ph_dec[(i_lim+1)]]* (d2 - 2 * i_lim))  # saturation partielle du noeud limite
            end
        end
    end
    return(nv_poids)
end

function RechLocEchange(chemin, i_ph_dec, i_to_i_ph_dec, sum_poids, d2, ph, p, d1, d, D, longueur, deltap, deltam, S, sum_arcs, s, t)
    i = 2
    i_lim = min(div(d2, 2), longueur)
    while (i <= longueur-1) 
        sommets_admissibles = intersect(deltap[chemin[i-1]], deltam[chemin[i+1]]) # il existe un chemin
        sommets_admissibles = filter(x -> x != chemin[i], sommets_admissibles) # enlever le chemin actuel
        println("sommets admissibles = ", sommets_admissibles, " pour le sommet ", chemin[i])
        for nv_noeud in collect(sommets_admissibles)
            nv_poids = nvResPoids(i_ph_dec, i_to_i_ph_dec, sum_poids, nv_noeud, chemin[i], d2, ph, p, longueur) # se fait en temps constant (youpi !)
            
            # verification
            nv_chemin = nvChemin(chemin, chemin[i], nv_noeud)
            nv_poids_lent = getInfoSommets(nv_chemin, p, ph, d2)
            if  nv_poids!= nv_poids_lent
                @warn("calcul poids incorrect")
               
            end
            # fin verification

            if nv_poids <= S
                #println("Solution admissible")
                if nvDist(chemin, nv_noeud, chemin[i], d1, d, D) < sum_arcs # temps lineaire en le nombre d'aretes
                    println("____________________________")
                    println("on a trouve une solution ameliorante")
                    println("nv_noeud admissible ? ", isAdmissible(chemin, nv_noeud, chemin[i], d2, p, ph, deltap, S, s, t))
                    println("nv distances ? ", nvDist(chemin, nv_noeud, chemin[i], d1, d, D))
                    println("____________________________")
                    return(true, chemin[i], nv_noeud)
                end
            end
        end
        i += 1
    end
    return(false, -1, -1)
end

function deplacementEch(chemin, i_ph_dec, i_to_i_ph_dec, sum_poids, sum_arcs, longueur, d1, d2, p, ph, d, D, old_noeud, nv_noeud)
    """on enleve old noeud du chemin et on met nv_noeud, retourne les informations"""
    # i_lim, longueur restent égaux
    sum_poids = nvResPoids(i_ph_dec, i_to_i_ph_dec, sum_poids, nv_noeud, old_noeud, d2, ph, p, longueur) # utiliser le temps constant

    chemin = nvChemin(chemin, old_noeud, nv_noeud)

    i_ph_dec =sort(chemin, lt = (x, y) -> ph[x] >= ph[y])
    i_to_i_ph_dec=Dict() # passer du numero du sommet a sa position dans i_ph_dec
    for i in collect(chemin)
        i_to_i_ph_dec[i]= findfirst(x -> x == i, i_ph_dec)
    end
    sum_arcs =  Dist(chemin, d1, d, D)
    return(chemin, i_ph_dec, i_to_i_ph_dec, sum_poids, sum_arcs)
end
