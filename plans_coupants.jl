using JuMP
using CPLEX
include("./utils/parsing.jl")
include("PLNE_dual_poids.jl")
include("./constrSol.jl")

function plans_coupantsALG(n::Int64, s::Int64, t::Int64, S::Int64, d1::Int64, d2::Int64, p::Vector{Int64}, ph::Vector{Int64}, d, D, name_instance, symmetry, with_initial_values, timelimit)
    """Résout la méthode des plans coupants (par callback)"""
    if with_initial_values=="with initial values"
        Keys=collect(keys(D))
        d_D_set, K_set, A_d_D_set=calcul_d_D_k_set(p, d, D, S)
        p_set, ph_set, Kp_set, Kph_set, A_p_set, A_ph_set=calcul_p_ph_k_set(p, ph, S)
        # initial_values=plne_dual_poids(n, s, t, S, d1, d2, p, ph, d, D)
        initial_values=get_init_sol(name_instance)

        borne_sup=sum([initial_values[1][key]*d[key] for key in Keys])
        capa=0
        for d_D in reverse(d_D_set)
            if capa>d1-1e-5
                break
            end
            for (i,j) in A_d_D_set[d_D]
                if initial_values[1][i,j]>1e-5
                    if capa+D[i,j]<=d1
                        capa+=D[i,j]
                        borne_sup+=d[i,j]*D[i,j]
                    else
                        borne_sup+=d[i,j]*(d1-capa)
                        capa=d1
                        break
                    end
                end
            end
        end
        println(borne_sup)
        n, s, t, S, d1, d2, p, ph, d, D, initial_values = simplify_instance(name_instance, borne_sup, initial_values)
    end
    # initialisation des voisinages entrants et sortants pour chaque sommet
    deltap, deltam=initDelta(d, n)
    println(n)

    Keys=collect(keys(D))
    d_D_set, K_set, A_d_D_set=calcul_d_D_k_set(p, d, D, S)
    p_set, ph_set, Kp_set, Kph_set, A_p_set, A_ph_set=calcul_p_ph_k_set(p, ph, S)

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
    if symmetry=="no_symmetry"
        @variable(m, y[d_D in d_D_set, k in 1:K_set[d_D]], Bin)
        @constraint(m,  [d_D in d_D_set], sum(k*y[d_D,k] for k in 1:K_set[d_D]) ==sum(x[i,j] for (i,j) in A_d_D_set[d_D]))
        @constraint(m,  [d_D in d_D_set], sum(y[d_D,k] for k in 1:K_set[d_D]) <=1)

        @constraint(m, sum(d_D[1]*k*y[d_D,k] for d_D in d_D_set for k in 1:K_set[d_D]) <= z)
        
        @variable(m, wp[p_i in p_set, k in 1:Kp_set[p_i]], Bin)
        @constraint(m,  [p_i in p_set], sum(k*wp[p_i,k] for k in 1:Kp_set[p_i]) ==sum(a[i] for i in A_p_set[p_i]))
        @constraint(m,  [p_i in p_set], sum(wp[p_i,k] for k in 1:Kp_set[p_i]) <=1)

        @constraint(m, sum(p_i*k*wp[p_i,k] for p_i in p_set for k in 1:Kp_set[p_i]) <= S)

        @variable(m, wph[ph_i in ph_set, k in 1:Kph_set[ph_i]], Bin)
        @constraint(m,  [ph_i in ph_set], sum(k*wph[ph_i,k] for k in 1:Kph_set[ph_i]) ==sum(a[i] for i in A_ph_set[ph_i]))
        @constraint(m,  [ph_i in ph_set], sum(wph[ph_i,k] for k in 1:Kph_set[ph_i]) <=1)
    end

    @constraint(m,  [i in 1:n; i!=s && i!=t], sum(x[i,j] for j in deltap[i]) - sum(x[j,i] for j in deltam[i])==0)
    @constraint(m,  sum(x[s,j] for j in deltap[s]) - sum(x[j,s] for j in deltam[s])==1)
    @constraint(m,  sum(x[t,j] for j in deltap[t]) - sum(x[j,t] for j in deltam[t])==-1)

    @constraint(m,  [i in 1:n; i!=t], sum(x[i,j] for j in deltap[i])==a[i])
    @constraint(m,  a[t]==1)
    @constraint(m, constraint_expr, sum(x[key]*d[key] for key in Keys) <= z)
    @constraint(m, sum(a[i]*p[i] for i in 1:n) <= S)

    @objective(m, Min, z)

    iter=0
    iter24=0
    best_integer=-1

    function mon_super_callback(cb_data::CPLEX.CallbackContext, context_id::Clong)
        if isIntegerPoint(cb_data, context_id)

            iter+=1
            flag_23=false
            flag_24=false
    
            CPLEX.load_callback_variable_primal(cb_data, context_id) # on regardera les valeurs des variables dans la solution courante

            x_val=Dict(key => callback_value(cb_data, x[key]) for key in Keys) # on stocke x
            a_val = [callback_value(cb_data, a[i]) for i in 1:n] # on stocke a
            println([i for i in 1:n if a_val[i]>1e-5])

            z_etoile = callback_value(cb_data, z) # valeur objectif actuelle
              
            # poids du chemin
            """Idée : on ajoute une contrainte sur le poids des sommets en trouvant un augmentation des poids très élevée"""
            res=sum([a_val[i]*p[i] for i in 1:n])
            capa=0 # budget pour augmenter les delta^2
            delta_2_val = [0 for i in 1:n] # on stocke les delta^2
            delta_2_val_ph = Dict(ph_i => 0.0 for ph_i in ph_set)
            for ph_i in reverse(ph_set)
                if capa>d2-1e-5
                    break
                end
                for i in A_ph_set[ph_i]
                    if a_val[i]>1e-5
                        if capa+2<=d2
                            capa+=2
                            res+=2*ph[i]
                            delta_2_val[i]=2
                            delta_2_val_ph[ph_i]+=2
                        else
                            res+=ph[i]*(d2-capa)
                            delta_2_val[i]=d2-capa
                            delta_2_val_ph[ph_i]+=d2-capa
                            capa=d2
                            break
                        end
                    end
                end
            end
        
            if res >= S + 1e-5 # si notre poids est plus grand que S -> ajout de la contrainte
                println(res, " ", S)
                cstr = @build_constraint(sum(a[i]*(p[i]+delta_2_val[i]*ph[i]) for i in 1:n) <= S)
                MOI.submit(m, MOI.LazyConstraint(cb_data), cstr)
                if symmetry=="no_symmetry"
                    wp_val=Dict((p_i,k) => callback_value(cb_data, wp[p_i, k]) for p_i in p_set for k in 1:Kp_set[p_i])
                    wph_val=Dict((ph_i,k) => callback_value(cb_data, wph[ph_i, k]) for ph_i in ph_set for k in 1:Kph_set[ph_i])
                    Kp_set_star=Dict(p_i => k for p_i in p_set for k in 1:Kp_set[p_i] if wp_val[p_i,k]>1e-5)  
                    Kph_set_star=Dict(ph_i => k for ph_i in ph_set for k in 1:Kph_set[ph_i] if wph_val[ph_i,k]>1e-5)             
                    cstrbis = @build_constraint(sum(sum(wp[p_i,k] for k in Kp_set_star[p_i]:Kp_set[p_i]) for p_i in keys(Kp_set_star))
                                                + sum(sum(wph[ph_i,k] for k in Kph_set_star[ph_i]:Kph_set[ph_i]) for ph_i in keys(Kph_set_star))
                                                <= sum(wp_val[p_i, k] for p_i in p_set for k in 1:Kp_set[p_i] if wp_val[p_i,k]>1e-5)
                                                + sum(wph_val[ph_i, k] for ph_i in ph_set for k in 1:Kph_set[ph_i] if wph_val[ph_i,k]>1e-5)-1)
                    MOI.submit(m, MOI.LazyConstraint(cb_data), cstrbis)
                    cstr1b=@build_constraint(sum(wph[ph_i,k]*ph_i*min(k*2, delta_2_val_ph[ph_i]) for ph_i in keys(Kph_set_star) for k in 1:Kph_set[ph_i])
                                                        + sum(p_i*k*wp[p_i,k] for p_i in p_set for k in 1:Kp_set[p_i]) <= S)

                    MOI.submit(m, MOI.LazyConstraint(cb_data), cstr1b)
                end
            else
                flag_23=true
            end
            # else
            # duree du chemin
            """Même idée sur les chemins : on augmente le poids d'un chemin en trouvant les arcs qui augmente le plus le coût"""
            res=sum([x_val[key]*d[key] for key in Keys])
            capa=0
            delta_1_val = Dict(key => 0.0 for key in Keys) # on initialise tous les delta^1 des arêtes à 0
            delta_1_val_d_D = Dict(d_D => 0.0 for d_D in d_D_set)
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
                            delta_1_val_d_D[d_D]+=d_D[2]
                        else
                            res+=d[i,j]*(d1-capa)
                            delta_1_val[i,j]=d1-capa
                            delta_1_val_d_D[d_D]+=d1-capa
                            capa=d1
                            break
                        end
                    end
                end
            end

            if res >=  z_etoile + 1e-5 # si le coût du chemin est plus grand que la valeur objectif actuelle -> ajout de la contrainte
                iter24+=1
                println(res, " ", z_etoile, " ", flag_23)
                cstr2 = @build_constraint(sum(x[key]*d[key]*(1+delta_1_val[key]) for key in Keys) <= z)
                MOI.submit(m, MOI.LazyConstraint(cb_data), cstr2)
                if symmetry=="no_symmetry"
                    y_val=Dict((d_D,k) => callback_value(cb_data, y[d_D, k]) for d_D in d_D_set for k in 1:K_set[d_D])
                    K_set_star=Dict(d_D => k for d_D in d_D_set for k in 1:K_set[d_D] if y_val[d_D,k]>1e-5)
                    if flag_23
                        if res >= best_integer+0.1 && best_integer>=0
                            cstr3 = @build_constraint(sum(sum(y[d_D,k] for k in K_set_star[d_D]:K_set[d_D]) for d_D in keys(K_set_star))<= sum(y_val[d_D, k] for d_D in d_D_set for k in 1:K_set[d_D] if y_val[d_D,k]>1e-5)-1)
                        else
                            cstr3 = @build_constraint(sum(sum(y[d_D,k] for k in K_set_star[d_D]:K_set[d_D]) for d_D in keys(K_set_star))<= sum(y_val[d_D, k] for d_D in d_D_set for k in 1:K_set[d_D] if y_val[d_D,k]>1e-5)-sum(x_val[i,j]-x[i,j] for i in 1:n for j in deltap[i] if x_val[i,j]>1e-5)/sum([x_val[i,j] for i in 1:n for j in deltap[i]]))
                        end
                        MOI.submit(m, MOI.LazyConstraint(cb_data), cstr3)
                    end
                    cstr2b=@build_constraint(sum(y[d_D,k]* d_D[1]*min(k*d_D[2], delta_1_val_d_D[d_D]) for d_D in keys(K_set_star) for k in 1:K_set[d_D])
                                                + sum(d_D[1]*k*y[d_D,k] for d_D in d_D_set for k in 1:K_set[d_D]) <= z)
                    MOI.submit(m, MOI.LazyConstraint(cb_data), cstr2b)
                else
                    flag_24=true
                end
                if flag_23 && flag_24
                    if best_integer>= 0
                        best_integer=min(best_integer, res)
                    else
                        best_integer=res
                    end
                end
            end
        end
    end

    if with_initial_values=="with initial values"
        JuMP.set_start_value.(x, JuMP.Containers.SparseAxisArray(initial_values[1]))
        JuMP.set_start_value.(a, initial_values[2])
    end

    if symmetry=="no_symmetry"
        if with_initial_values=="with initial values"
            logfile_name = "txtFiles/plans_coupantsALG/no_symmetry/with_initial_values/$name_instance.txt"
        else
            logfile_name = "txtFiles/plans_coupantsALG/no_symmetry/with_initial_values/$name_instance.txt"
        end
    else
        if with_initial_values=="with initial values"
            logfile_name = "txtFiles/plans_coupantsALG/symmetry/with_initial_values/$name_instance.txt"
        else
            logfile_name = "txtFiles/plans_coupantsALG/symmetry/with_initial_values/$name_instance.txt"
        end
    end

    # Obtenir le chemin absolu du fichier de journal dans le répertoire actuel
    logfile_path = abspath(logfile_name)
    logfile = open(logfile_path, "w")
    #redirect_stdout(logfile)


    # On précise que le modèle doit utiliser notre fonction de callback
    MOI.set(m, CPLEX.CallbackFunction(), mon_super_callback)
    optimize!(m)

    println(primal_status(m))
    println(termination_status(m))
    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT

    if feasibleSolutionFound
        # Récupération des valeurs d’une variable 
        # println("Value x : ", JuMP.value.(x))
        println("Valeur de l’objectif : ", JuMP.objective_value(m))
        println("Nombre d'itérations : ", iter)
        close(logfile)
        
        return JuMP.value.(x), JuMP.value.(a)
    end
    close(logfile)
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
    #redirect_stdout(logfile)

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
