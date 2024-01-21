using JuMP
using CPLEX

function plans_coupants(n::Int64, s::Int64, t::Int64, S::Int64, d1::Int64, d2::Int64, p::Vector{Int64}, ph::Vector{Int64}, d::Dict{Any, Any}, D::Dict{Any, Any}, timelimit)
    """Résout la méthode des plans coupants (par callback)"""
    # initialisation des voisinages entrants et sortants pour chaque sommet
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

    order_ph=sortperm(ph, rev=true) # ph = p chapeau (incertitude pour le poids de chaque sommet) -> donne les indices dans l'ordre décroissant des ph
    Dh=[D[key] for key in collect(keys(D))] # stocke les distances dans un vecteur
    order_dh=sortperm(Dh, rev=true) # ordre décroissant des distances
    nD=size(order_dh)[1]
    Keys=collect(keys(D))

    m = Model(CPLEX.Optimizer)
    #set_silent(m)
    
    set_optimizer_attribute(m, "CPXPARAM_Preprocessing_Presolve", 0)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Strategy_FPHeur", -1)
    # set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)
    set_optimizer_attribute(m, "CPXPARAM_TimeLimit", timelimit)
    set_optimizer_attribute(m, "CPX_PARAM_THREADS", 1)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Display", 4)

    MOI.set(m, MOI.NumberOfThreads(), 1)


    # contraintes du PLNE maitre

    @variable(m, z >=0 )    
    @variable(m, x[i in 1:n, j in deltap[i]], Bin)
    #@variable(m, a[1:n]>=0)
    @variable(m, a[1:n], Bin)

    @constraint(m,  [i in 1:n; i!=s && i!=t], sum(x[i,j] for j in deltap[i]) - sum(x[j,i] for j in deltam[i])==0)
    @constraint(m,  sum(x[s,j] for j in deltap[s]) - sum(x[j,s] for j in deltam[s])==1)
    @constraint(m,  sum(x[t,j] for j in deltap[t]) - sum(x[j,t] for j in deltam[t])==-1)

    @constraint(m,  [i in 1:n; i!=t], sum(x[i,j] for j in deltap[i])==a[i])
    @constraint(m,  a[t]==1)

    @objective(m, Min, z)

    iter=0
    
    derniere_solution=[]
    x_last=Dict()

    function mon_super_callback(cb_data::CPLEX.CallbackContext, context_id::Clong)
        if isIntegerPoint(cb_data, context_id)

            iter+=1
    
            CPLEX.load_callback_variable_primal(cb_data, context_id) # on regardera les valeurs des variables dans la solution courante

            x_val=Dict(key => callback_value(cb_data, x[key]) for key in Keys) # on stocke x
            a_val = [callback_value(cb_data, a[i]) for i in 1:n] # on stocke a

            nodes_visited= [i for i in 1:n if a_val[i]>1e-5]

            z_etoile = callback_value(cb_data, z) # valeur objectif actuelle
              
            # poids du chemin
            """Idée : on ajoute une contrainte sur le poids des sommets en trouvant un augmentation des poids très élevée"""

            res=sum([a_val[i]*p[i] for i in 1:n])
            capa=0 # budget pour augmenter les delta^2
            I=1 # on regarde les indices dans l'ordre dans l'ordre des p^hat qui permettent d'augmenter le plus le poids du chemin
            delta_2_val = [0 for i in 1:n] # on stocke les delta^2

            while capa+2<=d2 && I<=n # tant que j'ai encore du budget pour augmenter le poids de mes sommets
                if a_val[order_ph[I]]>1e-5 # c'est bien un sommet que j'ai sélectionné 
                    capa+=2 # augmentation du budget
                    res+=2*ph[order_ph[I]] # augmentation du poids total des sommets
                    delta_2_val[order_ph[I]]=2 # stockage valeur delta^2
                end
                I+=1
            end

            while I<=n && a_val[order_ph[I]]<1e-5 # on passe les sommets restants non selectionnés
                I+=1
            end

            if I<=n # on augmente avec le budget qu'il nous reste
                res+=(d2-capa)*ph[order_ph[I]] 
                delta_2_val[order_ph[I]]=d2-capa
            end
        
            if res >= S + 1e-5 # si notre poids est plus grand que S -> ajout de la contrainte
                cstr = @build_constraint(sum(a[i]*(p[i]+delta_2_val[i]*ph[i]) for i in 1:n) <= S)
                MOI.submit(m, MOI.LazyConstraint(cb_data), cstr)
                # println("Ajout d'une contrainte de type 23")
            
            else
                # duree du chemin
                """Même idée sur les chemins : on augmente le poids d'un chemin en trouvant les arcs qui augmente le plus le coût"""

                res=sum([x_val[key]*d[key] for key in Keys])
                capa=0
                I=1
                delta_1_val = Dict(key => 0.0 for key in collect(keys(D))) # on initialise tous les delta^1 des arêtes à 0
                while capa+D[Keys[order_dh[I]]]<=d1 && I<nD # exactement même idée que pour les sommets
                    if x_val[Keys[order_dh[I]]]>1e-5
                        capa+=D[Keys[order_dh[I]]]
                        res+=d[Keys[order_dh[I]]]*D[Keys[order_dh[I]]]
                        delta_1_val[Keys[order_dh[I]]]=D[Keys[order_dh[I]]]
                    end
                    I+=1
                end
                while I<=nD && x_val[Keys[order_dh[I]]]<1e-5
                    I+=1
                end
                if I<=nD
                    res+=d[Keys[order_dh[I]]]*(d1-capa)
                    delta_1_val[Keys[order_dh[I]]]=d1-capa
                end

                if res >=  z_etoile + 1e-5 # si le coût du chemin est plus grand que la valeur objectif actuelle -> ajout de la contrainte
                    cstr2 = @build_constraint(sum(x[key]*d[key]*(1+delta_1_val[key]) for key in Keys) <= z)
                    MOI.submit(m, MOI.LazyConstraint(cb_data), cstr2)
                    # println("Ajout d'une contrainte de type 24 ")

                elseif derniere_solution!=nodes_visited # nodes visited : noeuds dans la solution que l'on regarde [ALH je ne comprends pas cette partie là]
                    # println("Ajout d'une contrainte de type 25 ")
                    # Contrainte à prouver: elle permet d'enlever de la symmétrie. 
                    # L'hypothèse est que si sum(x*d*(1+D))<= sum(y*d*(1+D)) alors v(x)<=v(y) (avec v(x) la vraie valeure de la solution)
                    # Une fois une solution y est trouvée, on va donc chercher des solutions vraiment plus petites (ça enlève des solutions de même coût)
                    # cstr3 = @build_constraint(sum(x[key]*d[key]*(1+D[key]) for key in Keys)+ sum(a_val[i]-a[i] for i in 1:n if a_val[i]>1e-5)/sum(a_val)<= sum([x_val[key]*d[key]*(1+D[key]) for key in Keys]))
                    
                    # Cette contrainte fonctionne aussi avec l'hypothèseque d[i,j]=d[i', j'] <=> D[i,j]=D[i', j'] ce qui est le cas dans nos données
                    cstr3 = @build_constraint(sum(x[key]*d[key] for key in Keys)+ sum(a_val[i]-a[i] for i in 1:n if a_val[i]>1e-5)/sum(a_val)<= sum([x_val[key]*d[key] for key in Keys]))
                    MOI.submit(m, MOI.LazyConstraint(cb_data), cstr3)
                    x_last=deepcopy(x_val)
                    derniere_solution=deepcopy(nodes_visited)
                end
            end
        end
    end

    logfile_name = "txtFiles/plans_coupants.txt"

    # Obtenir le chemin absolu du fichier de journal dans le répertoire actuel
    logfile_path = abspath(logfile_name)
    logfile = open(logfile_path, "w")
    redirect_stdout(logfile)

    # On précise que le modèle doit utiliser notre fonction de callback
    MOI.set(m, CPLEX.CallbackFunction(), mon_super_callback)
    optimize!(m)

    println(primal_status(m))
    println(termination_status(m))
    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL

    if feasibleSolutionFound
        # Récupération des valeurs d’une variable 
        # println("Value x : ", JuMP.value.(x))
        println("Valeur de l’objectif : ", JuMP.objective_value(m))
        println("Nombre d'itérations : ", iter)
        #return JuMP.value.(x), JuMP.value.(a)
    end
end

function isIntegerPoint(cb_data::CPLEX.CallbackContext, context_id::Clong)

    # context_id  == CPX_CALLBACKCONTEXT_CANDIDATE si le  callback est
    # appelé dans un des deux cas suivants :
    # cas 1 - une solution entière a été obtenue; ou
    # cas 2 - une relaxation non bornée a été obtenue
    if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
        return false
    end

    # Pour déterminer si on est dans le cas 1 ou 2, on essaie de récupérer la
    # solution entière courante
    ispoint_p = Ref{Cint}()
    ret = CPXcallbackcandidateispoint(cb_data, ispoint_p)

    # S'il n'y a pas de solution entière
    if ret != 0 || ispoint_p[] == 0
        return false
    else
        return true
    end
end
