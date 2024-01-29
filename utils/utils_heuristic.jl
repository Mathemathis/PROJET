function nvChemin(chemin, old_noeud, nv_noeud)
    """Nouveau chemin ou on remplace old_noeud par un nouveau noeud"""
    nv_chemin = copy(chemin) # calcul du nouveau chemin
    for i in 1:length(chemin) 
        if nv_chemin[i] == old_noeud
            nv_chemin[i] = nv_noeud
        end
    end
    return(nv_chemin)
end

function cheminToAretes(chemin)
    """chemin : liste de sommets
    renvoie une liste d'aretes [(i,j)] correspondant au chemin"""
    aretes = []
    current_node = chemin[1]
    for i in 2:length(chemin) # calcul du nouveau chemin
        push!(aretes, (current_node, chemin[i]))
        current_node = chemin[i]
    end
    return aretes
end

function isChemin(chemin, deltap)
    """On peut tracer le chemin"""
    for i in 1:(length(chemin)-1)
        if !(in(chemin[i+1], deltap[chemin[i]]))
            return(false) # on ne peut pas tracer le chemin
        end
    end
    return(true)
end

function Dist(chemin, d1, d, D)
    "renvoie la distance d'un chemin"
    aretes = cheminToAretes(chemin)
    return(getInfoArcs(aretes, d, D, d1))
end


function nvDist(chemin, nv_noeud, old_noeud, d1, d, D)
    nv_chemin = nvChemin(chemin, old_noeud, nv_noeud)
    return(Dist(nv_chemin, d1, d, D))
end

function getInfoSommets(chemin, p, ph, d2)
    """Renvoie le poids robuste d'un chemin ainsi que i_res le sommet limite dans les sommets orientes par ordre decroissant des ph"""
    i_ph_dec =sort(chemin, lt = (x, y) -> ph[x] <= ph[y], rev = true)
    res=sum([p[i] for i in collect(i_ph_dec)]) # somme deterministe
    capa=0 # budget pour augmenter les delta^2
    i_res = 0 # dernier indice à être rempli totalement (indice dans i_poids_croissants) 
    for i in i_ph_dec
        if capa+2 <= d2
            capa+=2 # augmentation du budget
            res+=2*ph[i] # augmentation du poids total des sommets
            i_res += 1 # on rempli totatlement ce sommet
        else 
            res+=(d2-capa)*ph[i]
            capa = 0
            break
        end
    end
    return res, i_res # total des poids max, indice du dernier sommet rempli totalement
end


function getInfoArcs(aretes, d, D, d1)
    """Renvoie la distance d'une suite d'aretes"""
    i_aretes_d =sort(aretes, lt = (x, y) -> d[x] >= d[y])
    res=sum([d[i] for i in collect(i_aretes_d)]) # somme deterministe
    capa=0 # budget pour augmenter les delta^2
    
    for (i,j) in i_aretes_d
        if capa+D[i,j]<=d1
            capa+=D[i,j]
            res+=d[i,j]*D[i,j]
        else
            res+=d[i,j]*(d1-capa)
            capa = d1
            break
        end
    end
    return res
end


function isAdmissible(chemin, nv_noeud, old_noeud, d2, p, ph, deltap, S)
    """Renvoie true si la nouvelle solution est 
    -  bien un chemin
    -  de poids robuste admissible"""
    nv_chemin = nvChemin(chemin, old_noeud, nv_noeud) # calcul du nouveau chemin
    if isChemin(nv_chemin, deltap)
        res_poids, _ = getInfoSommets(chemin, p, ph, d2)
        return(res_poids <= S) # poids 
    else
        return(false)
    end
end

function initDelta(d, n)
    # init delta p delta m
    deltap=Dict()
    deltam=Dict()
    for i in 1:n
        deltap[i]=[]
        deltam[i]=[]
    end
    for (i,j) in keys(d)
        push!(deltap[i],j)
        push!(deltam[j],i)
    end
    return(deltap, deltam)
end
