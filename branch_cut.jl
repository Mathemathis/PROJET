using JuMP
using CPLEX


function branchAndCut(name_instance, symmetry, with_initial_values, timelimit)
    """Résout la méthode des branch and cut (sans callback)"""
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    if with_initial_values=="with initial values"
        Keys=collect(keys(D))
        d_D_set, K_set, A_d_D_set=calcul_d_D_k_set(p, d, D, S, 100000)
        p_set, ph_set, Kp_set, Kph_set, A_p_set, A_ph_set=calcul_p_ph_k_set(p, ph, S)
        initial_values=get_init_sol(name_instance)

        borne_sup=sum([initial_values[1][key]*d[key] for key in Keys])
        capa=0
        last_d_D=Nothing
        for d_D in reverse(d_D_set)
            if capa>d1-1e-5
                last_d_D=d_D
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
        println("Borne sup ", borne_sup, " ", last_d_D)
        n, s, t, S, d1, d2, p, ph, d, D, initial_values, correspondance = simplify_instance(name_instance, borne_sup, initial_values)
    end
    # initialisation des voisinages entrants et sortants pour chaque sommet
    deltap, deltam=initDelta(d, n)
    println(n)

    Keys=collect(keys(D))
    d_D_set, K_set, A_d_D_set=calcul_d_D_k_set(p, d, D, S, borne_sup)
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
        @variable(m, y[d_D in d_D_set, k in 0:K_set[d_D]], Bin)
        @constraint(m,  [d_D in d_D_set], sum(k*y[d_D,k] for k in 0:K_set[d_D]) ==sum(x[i,j] for (i,j) in A_d_D_set[d_D]))
        @constraint(m,  [d_D in d_D_set], sum(y[d_D,k] for k in 0:K_set[d_D]) ==1)

        @constraint(m, sum(d_D[1]*k*y[d_D,k] for d_D in d_D_set for k in 1:K_set[d_D]) <= z)

        # @constraint(m, [d_D in d_D_set, key in A_d_D_set[d_D]; d_D[1]>=last_d_D[1]-200], y[d_D,0] <=1-x[key])
        
        @variable(m, wp[p_i in p_set, k in 0:Kp_set[p_i]], Bin)
        @constraint(m,  [p_i in p_set], sum(k*wp[p_i,k] for k in 0:Kp_set[p_i]) ==sum(a[i] for i in A_p_set[p_i]))
        @constraint(m,  [p_i in p_set], sum(wp[p_i,k] for k in 0:Kp_set[p_i]) ==1)

        @constraint(m, sum(p_i*k*wp[p_i,k] for p_i in p_set for k in 1:Kp_set[p_i]) <= S)

        @variable(m, wph[ph_i in ph_set, k in 0:Kph_set[ph_i]], Bin)
        @constraint(m,  [ph_i in ph_set], sum(k*wph[ph_i,k] for k in 0:Kph_set[ph_i]) ==sum(a[i] for i in A_ph_set[ph_i]))
        @constraint(m,  [ph_i in ph_set], sum(wph[ph_i,k] for k in 0:Kph_set[ph_i]) ==1)

        # @constraint(m, [ph_i in ph_set, i in A_ph_set[ph_i]], wph[ph_i,0] <=1-a[i])
        # @constraint(m, [p_i in p_set, i in A_p_set[p_i]], wp[p_i,0] <=1-a[i])
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
    iter23=0
    iter24=0
    best_integer=-1

    flag_23=true
    flag_24=true

    while flag_23 || flag_24
        iter+=1
        flag_23=false
        flag_24=false

        optimize!(m)

        # solution courante
        x_val = JuMP.value.(x)
        a_val = JuMP.value.(a)
        z_etoile = JuMP.objective_value(m)
        wp_val=JuMP.value.(wp)
        wph_val= JuMP.value.(wph) 
        Kp_set_star=Dict(p_i => k for p_i in p_set for k in 1:Kp_set[p_i] if wp_val[p_i,k]>1e-5)  
        Kph_set_star=Dict(ph_i => k for ph_i in ph_set for k in 1:Kph_set[ph_i] if wph_val[ph_i,k]>1e-5) 
        y_val= JuMP.value.(y)
        K_set_star=Dict(d_D => k for d_D in d_D_set for k in 1:K_set[d_D] if y_val[d_D,k]>1e-5)         
                

        # sous problème poids du chemin
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


        if res >= S + 1e-5
            iter23 += 1
            println("Ajout d'une contrainte de type 23")
            if symmetry=="no_symmetry"
                
                @constraint(m, sum(sum(wp[p_i,k] for k in Kp_set_star[p_i]:Kp_set[p_i]) for p_i in keys(Kp_set_star))
                                            + sum(sum(wph[ph_i,k] for k in Kph_set_star[ph_i]:Kph_set[ph_i]) for ph_i in keys(Kph_set_star))
                                            <= sum(wp_val[p_i, k] for p_i in p_set for k in 1:Kp_set[p_i] if wp_val[p_i,k]>1e-5)
                                            + sum(wph_val[ph_i, k] for ph_i in ph_set for k in 1:Kph_set[ph_i] if wph_val[ph_i,k]>1e-5)-1)

                @constraint(m, sum(wph[ph_i,k]*ph_i*min(k*2, delta_2_val_ph[ph_i]) for ph_i in keys(Kph_set_star) for k in 1:Kph_set[ph_i])
                                                    + sum(p_i*k*wp[p_i,k] for p_i in p_set for k in 1:Kp_set[p_i]) <= S)


                iter23+=1
            end
            @constraint(m,  sum(a[i]*(p[i]+delta_2_val[i]*ph[i]) for i in 1:n) <= S)
        else
            flag_23=true
        end
        
        # sous probleme duree du chemin
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

        if res >=  z_etoile + 1e-5
            iter24+=1
            @constraint(m,  sum(x[key]*d[key]*(1+delta_1_val[key]) for key in Keys) <= z)
            if symmetry=="no_symmetry"
                if flag_23
                    if res >= best_integer+0.1 && best_integer>=0
                        @constraint(m, sum(sum(y[d_D,k] for k in K_set_star[d_D]:K_set[d_D]) for d_D in keys(K_set_star))<= sum(y_val[d_D, k] for d_D in d_D_set for k in 1:K_set[d_D] if y_val[d_D,k]>1e-5)-1)
                    else
                        @constraint(m, sum(sum(y[d_D,k] for k in K_set_star[d_D]:K_set[d_D]) for d_D in keys(K_set_star))<= sum(y_val[d_D, k] for d_D in d_D_set for k in 1:K_set[d_D] if y_val[d_D,k]>1e-5)-sum(x_val[i,j]-x[i,j] for i in 1:n for j in deltap[i] if x_val[i,j]>1e-5)/sum([x_val[i,j] for i in 1:n for j in deltap[i]]))
                    end
                end
                @constraint(m, sum(y[d_D,k]* d_D[1]*min(k*d_D[2], delta_1_val_d_D[d_D]) for d_D in keys(K_set_star) for k in 1:K_set[d_D])
                                                + sum(d_D[1]*k*y[d_D,k] for d_D in d_D_set for k in 1:K_set[d_D]) <= z)
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

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    if feasibleSolutionFound
        a_val=JuMP.value.(a)
        x_val=JuMP.value.(x)

        Nodes=[i for i in 1:n if a_val[i]>1e-5]
        Nodes_tries=[s]
        Current_node=s
        while Current_node!=t
            for i in Nodes
                if i in deltap[Current_node] && x_val[Current_node, i]>0.1
                    push!(Nodes_tries, i)
                    Current_node=i
                    break
                end
            end
        end
        
        return name_instance, computation_time, [correspondance[i] for i in Nodes], JuMP.objective_value(m), objective_bound(m), iter23, iter24
    end

end
