using DataStructures


function deleteBoucles(chemin)
    #println("chemin = ", chemin)
    B = [(i, count(==(i), chemin)) for i in unique(chemin)]
    #println("B =", B)
    duplicates = filter(kv -> kv[2] > 1, B) # sommets comptés deux fois
    #println("duplicates = ", duplicates)
    if duplicates == []
        return chemin
    else
        println("chemin avec des doublons")
        doublon = duplicates[1][1]
        indices = findall(x -> x == doublon, chemin)
        #println("indices = ", indices)
        chemin = vcat(chemin[1:indices[1]], chemin[(indices[2]+1):end])
        #println("chemin = ", chemin)
        return deleteBoucles(chemin)
    end
end
                
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
    #println("aretes = ", aretes)
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

function isAdmissible(chemin, nv_noeud, old_noeud, d2, p, ph, deltap, S, s, t) # je ne me sers plus de cette fonction
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
    """Prend les valeurs de a du PLNE et retourne le chemin associe"""
    chemin = [s]
    current_node=s
    while current_node!=t
        for j in collect(deltap[current_node])
            if (a[j]>= 1 - 1e-5) && !(j in chemin) # attention à ne pas faire de boucles
                push!(chemin, j)
                current_node = j
                break
            end
        end
    end
    return chemin
end

# fonctions de deplacement associees aux voisinages

function nvCheminInf(chemin, noeud, nv_noeud)
    """Saute mouton / au lieu de chemin[i] -> chemin[i+1] -> chemin[i+2] -> chemin[i+3]
    on fait chemin[i] -> nv_noeud -> chemin[i+3]"""
    nv_chemin = [] # calcul du nouveau chemin
    i = 1
    while i <= length(chemin)
        push!(nv_chemin, chemin[i])
        if chemin[i] == noeud
            push!(nv_chemin, nv_noeud)
            i += 2
        end
        i += 1
    end
    return(nv_chemin)
end

function nvCheminSup(chemin, noeud, nv_noeud)
    """Saute mouton / au lieu de chemin[i] -> chemin[i+1] 
    on fait chemin[i] -> nv_noeud -> chemin[i+1]"""
    nv_chemin = [] # calcul du nouveau chemin
    i = 1
    while i <= length(chemin)
        push!(nv_chemin, chemin[i])
        if chemin[i] == noeud
            push!(nv_chemin, nv_noeud)
        end
        i += 1
    end
    return(nv_chemin)
end

function nvCheminEch(chemin, old_noeud, nv_noeud)
    #println("nouveau noeud = ", nv_noeud)
    #println("old noeud = ", old_noeud)
    """Nouveau chemin ou on remplace old_noeud par un nouveau noeud"""
    #println("chemin = ", chemin)
    nv_chemin = [] # calcul du nouveau chemin
    i = 1
    while i <= length(chemin)
        if chemin[i] == old_noeud
            push!(nv_chemin, nv_noeud)
        else
            push!(nv_chemin, chemin[i])
        end
        i += 1
    end
    #println("nv chemin = ", nv_chemin)
    return(nv_chemin)
end


function nvCheminInf2(chemin, noeud, k, l)
    """Saute mouton / au lieu de chemin[i] -> chemin[i+1] 
    on fait chemin[i] -> nv_noeud -> chemin[i+1]"""
    nv_chemin = [] # calcul du nouveau chemin
    i = 1
    while i <= length(chemin)
        push!(nv_chemin, chemin[i])
        if chemin[i] == noeud
            push!(nv_chemin, k)
            push!(nv_chemin, l)
            i += 3
        end
        i += 1
    end
    return(nv_chemin)
end

function nvCheminEch2(chemin, noeud, k, l)
    """Nouveau chemin ou on remplace old_noeud par un nouveau noeud"""
    nv_chemin = [] # calcul du nouveau chemin
    i = 1
    while i <= length(chemin) 
        push!(nv_chemin, chemin[i])
        if chemin[i] == noeud
            push!(nv_chemin, k)
            push!(nv_chemin, l)
            i += 2
        end
        i += 1
    end
    return(nv_chemin)
end

function nvCheminSup2(chemin, noeud, k, l)
    """Nouveau chemin ou on remplace old_noeud par un nouveau noeud"""
    nv_chemin = [] # calcul du nouveau chemin
    i = 1
    while i <=length(chemin) 
        push!(nv_chemin, chemin[i])
        if chemin[i] == noeud
            push!(nv_chemin, k)
            push!(nv_chemin, l)
            i += 1
        end
        i += 1
    end
    return(nv_chemin)
end