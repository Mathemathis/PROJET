function cheminToAretes(chemin)
    """chemin : liste de sommets
    renvoie une liste d'aretes [(i,j)] correspondant au chemin"""
    aretes = []
    current_node = chemin[1]
    for (i, item) in enumerate(chemin[2:end]) # calcul du nouveau chemin # attention i commence à 1
        push!(aretes, (current_node, item))
        current_node = item
    end
    return aretes
end

function isChemin(chemin, deltap, s, t)
    """On peut tracer le chemin"""
    if chemin[1] != s
        @warn("Le chemin ne commence pas en s")
        return false
    else  
        if chemin[end] != t
            @warn("Le chemin ne finit pas en t")
            return false
        end
    end
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


function getInfoSommets(chemin, p, ph, d2)
    """Renvoie le poids robuste d'un chemin ainsi que i_res le sommet limite dans les sommets orientes par ordre decroissant des ph"""
    i_ph_dec =sort(chemin, lt = (x, y) -> ph[x] <= ph[y], rev = true)
    res=sum([p[i] for i in collect(i_ph_dec)]) # somme deterministe
    capa=0 # budget pour augmenter les delta^2
    for i in i_ph_dec
        if capa+2 <= d2
            capa+=2 # augmentation du budget
            res+=2*ph[i] # augmentation du poids total des sommets
        else 
            res+=(d2-capa)*ph[i]
            capa = 0
            break
        end
    end
    return res # total des poids max, indice du dernier sommet rempli totalement
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


function isAdmissible(chemin, nv_noeud, old_noeud, d2, p, ph, deltap, S, s, t)
    """Renvoie true si la nouvelle solution est 
    -  bien un chemin
    -  de poids robuste admissible"""
    nv_chemin = nvChemin(chemin, old_noeud, nv_noeud) # calcul du nouveau chemin
    if isChemin(chemin, deltap, s, t)
        res_poids = getInfoSommets(chemin, p, ph, d2)
        return(res_poids <= S) # poids 
    else
        return(false)
    end
end

function initDelta(d, n)
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

function initChemin(a, deltap, s, t)
    chemin = [s]
    current_node=s
    while current_node!=t
        for j in collect(deltap[current_node])
            if (a[j]>= 1 - 1e-5) && !(j in chemin)
                push!(chemin, j)
                current_node = j
                break
            end
        end
    end
    return chemin
end