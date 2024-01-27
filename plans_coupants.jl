using JuMP
using CPLEX

function calcul_d_D_k_set(d, D)
    d_D=[(d[key], D[key]) for key in collect(keys(D))]
    count_element=Dict(d_D .=> 0)

    for element in d_D
        count_element[element] += 1
    end
    d_D_set=sort([key for key in collect(keys(count_element))], lt = (x, y) -> (x[1] < y[1]) || (x[1] == y[1] && x[2] < y[2]))
    A_d_D_set=Dict()
    for d_D in d_D_set
        A_d_D_set[d_D]=[]
    end
    for key in collect(keys(D))
        push!(A_d_D_set[d[key], D[key]], key)
    end
    A_d_D_set
    return d_D_set , count_element, A_d_D_set
end

function plans_coupantsALG(n::Int64, s::Int64, t::Int64, S::Int64, d1::Int64, d2::Int64, p::Vector{Int64}, ph::Vector{Int64}, d::Dict{Any, Any}, D::Dict{Any, Any}, name_instance, var_y, var_eps, timelimit)
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
    Keys=collect(keys(D))
    d_D_set, K_set, A_d_D_set=calcul_d_D_k_set(d, D)

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
    @variable(m, a[1:n], Bin)
    if var_y=="var_y"
        @variable(m, y[d_D in d_D_set, k in 1:K_set[d_D]], Bin)
        @constraint(m,  [d_D in d_D_set], sum(k*y[d_D,k] for k in 1:K_set[d_D]) ==sum(x[i,j] for (i,j) in A_d_D_set[d_D]))
        @constraint(m,  [d_D in d_D_set], sum(y[d_D,k] for k in 1:K_set[d_D]) <=1)
    end

    if var_eps=="var_eps"
        @variable(m, eps[1:40, 1:30], Bin)
    end

    @constraint(m,  [i in 1:n; i!=s && i!=t], sum(x[i,j] for j in deltap[i]) - sum(x[j,i] for j in deltam[i])==0)
    @constraint(m,  sum(x[s,j] for j in deltap[s]) - sum(x[j,s] for j in deltam[s])==1)
    @constraint(m,  sum(x[t,j] for j in deltap[t]) - sum(x[j,t] for j in deltam[t])==-1)

    @constraint(m,  [i in 1:n; i!=t], sum(x[i,j] for j in deltap[i])==a[i])
    @constraint(m,  a[t]==1)
    @constraint(m, constraint_expr, sum(x[key]*d[key] for key in Keys) <= z)

    @objective(m, Min, z)

    iter=0
    iter24=0

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

            while I<=n # exactement même idée que pour les sommets
                if a_val[order_ph[I]]>1e-5
                    if capa+2<=d2
                        capa+=2 # augmentation du budget
                        res+=2*ph[order_ph[I]] # augmentation du poids total des sommets
                        delta_2_val[order_ph[I]]=2 # stockage valeur delta^2
                    else # on augmente avec le budget qu'il nous reste
                        res+=(d2-capa)*ph[order_ph[I]] 
                        delta_2_val[order_ph[I]]=d2-capa
                        I+=n #on sort de la boucle
                    end
                end
                I+=1
            end
        
            if res >= S + 1e-5 # si notre poids est plus grand que S -> ajout de la contrainte
                cstr = @build_constraint(sum(a[i]*(p[i]+delta_2_val[i]*ph[i]) for i in 1:n) <= S)
                MOI.submit(m, MOI.LazyConstraint(cb_data), cstr)
                println("Ajout d'une contrainte de type 23")
            
            else
                # duree du chemin
                """Même idée sur les chemins : on augmente le poids d'un chemin en trouvant les arcs qui augmente le plus le coût"""
                res=sum([x_val[key]*d[key] for key in Keys])
                capa=0
                delta_1_val = Dict(key => 0.0 for key in Keys) # on initialise tous les delta^1 des arêtes à 0
                for d_D in reverse(d_D_set)
                    if capa>d1-1e-5
                        break
                    end
                    for (i,j) in A_d_D_set[d_D]
                        if x_val[i,j]>1e-5
                            if capa+D[i,j]<=d1
                                capa+=D[i,j]
                                res+=d[i,j]*D[i,j]
                                delta_1_val[i,j]=D[i,j]
                            else
                                res+=d[i,j]*(d1-capa)
                                delta_1_val[i,j]=d1-capa
                                capa=d1
                                break
                            end
                        end
                    end
                end

                if res >=  z_etoile + 1e-5 # si le coût du chemin est plus grand que la valeur objectif actuelle -> ajout de la contrainte
                    iter24+=1
                    println("Ajout d'une contrainte de type 24 ", iter24)
                    cstr2 = @build_constraint(sum(x[key]*d[key]*(1+delta_1_val[key]) for key in Keys) <= z)
                    MOI.submit(m, MOI.LazyConstraint(cb_data), cstr2)
                    println(res, " ", z_etoile)
                    if var_y=="var_y"
                        y_val=Dict((d_D,k) => callback_value(cb_data, y[d_D, k]) for d_D in d_D_set for k in 1:K_set[d_D])
                        cstr3 = @build_constraint(sum(y[d_D, k] for d_D in d_D_set for k in 1:K_set[d_D] if y_val[d_D,k]>1e-5)<= sum(y_val[d_D, k] for d_D in d_D_set for k in 1:K_set[d_D] if y_val[d_D,k]>1e-5)-sum(x_val[i,j]-x[i,j] for i in 1:n for j in deltap[i] if x_val[i,j]>1e-5)/sum([x_val[i,j] for i in 1:n for j in deltap[i]]))
                        println(cstr3)
                        MOI.submit(m, MOI.LazyConstraint(cb_data), cstr3)
                        if var_eps=="var_eps"
                            K_set_star=Dict()
                            i=1
                            for (index, d_D) in enumerate(d_D_set)
                                for k in 1:K_set[d_D]
                                    if y_val[d_D,k]>1e-5
                                        K_set_star[i, index, d_D]=k
                                        i+=1
                                    end
                                end
                            end
                            Keys_K_set_star=collect(keys(K_set_star))
                            d_D_set_star=[key[3] for key in Keys_K_set_star]
                            for (i, index, d_D) in Keys_K_set_star
                                k_star=K_set_star[i, index, d_D]
                                Mi=res/(d_D[1])
                                cstr4 = @build_constraint(sum(k*y[d_D, k] for k in k_star+1:K_set[d_D])
                                                            + sum(K_set_star[key]*y[key[3], K_set_star[key]] for key in Keys_K_set_star if key[1]>=i) 
                                                            +1- sum(K_set_star[key] for key in Keys_K_set_star if key[1]>=i)
                                                            + sum(sum(k*y[d_Dprim,k] for k in 1:K_set[d_Dprim]) for d_Dprim in d_D_set if (!(d_Dprim in d_D_set_star) && d_Dprim[1]>=d_D[1] && d_Dprim[2] >= d_D[2])) 
                                                            <= Mi*eps[iter24, i])
                                MOI.submit(m, MOI.LazyConstraint(cb_data), cstr4)
                                # println(cstr4)
                            end
                            cstr5=@build_constraint(sum(eps[iter24, i] for (i, index, d_D) in Keys_K_set_star) <= sum(1 for (i, index, d_D) in Keys_K_set_star)-sum(x_val[i,j]-x[i,j] for i in 1:n for j in deltap[i] if x_val[i,j]>1e-5)/sum([x_val[i,j] for i in 1:n for j in deltap[i]]))
                            MOI.submit(m, MOI.LazyConstraint(cb_data), cstr5)
                        end
                    end
                end
            end
        end
    end

    if var_y=="var_y"
        if var_eps == "var_eps"
            logfile_name = "txtFiles/plans_coupantsALG/var_y/var_eps/$name_instance.txt"
        else
            logfile_name = "txtFiles/plans_coupantsALG/var_y/no_var_eps/$name_instance.txt"
        end
    else
        logfile_name = "txtFiles/plans_coupantsALG/no_var_y/$name_instance.txt"
    end

    # Obtenir le chemin absolu du fichier de journal dans le répertoire actuel
    logfile_path = abspath(logfile_name)
    logfile = open(logfile_path, "w")
    redirect_stdout(logfile)

    # JuMP.set_start_value.(eps, 1)

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
        
        return JuMP.value.(x), JuMP.value.(a)
    end
end

function plans_coupantsPLNE(n::Int64, s::Int64, t::Int64, S::Int64, d1::Int64, d2::Int64, p::Vector{Int64}, ph::Vector{Int64}, d::Dict{Any, Any}, D::Dict{Any, Any}, timelimit)
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

    m = Model(CPLEX.Optimizer)
    set_silent(m)
    
    set_optimizer_attribute(m, "CPXPARAM_Preprocessing_Presolve", 0)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Strategy_FPHeur", -1)
    set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)
    set_optimizer_attribute(m, "CPXPARAM_TimeLimit", timelimit)

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

    function mon_super_callback(cb_data::CPLEX.CallbackContext, context_id::Clong)
        if isIntegerPoint(cb_data, context_id)
    
            CPLEX.load_callback_variable_primal(cb_data, context_id)


            x_val = Matrix{Float64}(undef,n,n) # matrice

            for i in 1:n
                for j in deltap[i]
                    x_val[i,j] = callback_value(cb_data, x[i,j])
                end
            end

            a_val = [callback_value(cb_data, a[i]) for i in 1:n]

            z_etoile = callback_value(cb_data, z)
            println("current z* = ", z_etoile)
              
            # poids du chemin

            SP2 = Model(CPLEX.Optimizer)
            set_silent(SP2)

            @variable(SP2, delta_2[1:n]>=0)    

            @constraint(SP2,  [i in 1:n], delta_2[i]<=2)
            @constraint(SP2,  sum(delta_2[i] for i in 1:n)<=d2)
            
            set_optimizer_attribute(SP2, "CPXPARAM_Preprocessing_Presolve", 0)
            set_optimizer_attribute(SP2, "CPXPARAM_TimeLimit", timelimit)
            set_optimizer_attribute(SP2, "CPX_PARAM_THREADS", 1)
            set_optimizer_attribute(SP2, "CPXPARAM_MIP_Display", 4)

            @objective(SP2, Max, sum(a_val[i]*(p[i]+delta_2[i]*ph[i]) for i in 1:n))
            optimize!(SP2)
            delta_2_val = JuMP.value.(delta_2)

            if JuMP.objective_value(SP2) >= S + 1e-5
                cstr = @build_constraint(sum(a[i]*(p[i]+delta_2_val[i]*ph[i]) for i in 1:n) <= S)
                MOI.submit(m, MOI.LazyConstraint(cb_data), cstr)
                println("Ajout d'une contrainte de type 23")
            
            end

            # duree du chemin

            SP1 = Model(CPLEX.Optimizer)
            set_silent(SP1)
            @variable(SP1, delta_1[i in 1:n, j in deltap[i]]>=0)    

            @constraint(SP1,  [i in 1:n, j in deltap[i]], delta_1[i,j]<=D[i,j])

            @constraint(SP1,  sum(delta_1[i,j] for i in 1:n for j in deltap[i])<=d1)
            
            set_optimizer_attribute(SP1, "CPXPARAM_Preprocessing_Presolve", 0)
            set_optimizer_attribute(SP1, "CPXPARAM_TimeLimit", timelimit)
            set_optimizer_attribute(SP1, "CPX_PARAM_THREADS", 1)
            set_optimizer_attribute(SP1, "CPXPARAM_MIP_Display", 4)

            @objective(SP1, Max, sum(x_val[i,j]*d[i,j]*(1+delta_1[i,j]) for i in 1:n for j in deltap[i]))
            optimize!(SP1)
            delta_1_val = JuMP.value.(delta_1)

            if JuMP.objective_value(SP1) >=  z_etoile + 1e-5
                println("solution SP1 =  ", JuMP.objective_value(SP1))
                println([(key, delta_1_val[key], D[key], d[key]) for key in collect(keys(D)) if delta_1_val[key]>0.001])
                cstr2 = @build_constraint(sum(x[i,j]*d[i,j]*(1+delta_1_val[i,j]) for i in 1:n for j in deltap[i]) <= z)
                MOI.submit(m, MOI.LazyConstraint(cb_data), cstr2)
                println("Ajout d'une contrainte de type 24")
            end
        end
    end

    logfile_name = "txtFiles/plans_coupantsPLNE.txt"

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
        return JuMP.value.(x), JuMP.value.(a)
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
