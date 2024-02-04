using JuMP
using CPLEX

function string_to_array(string_array)
    string_array = replace(string_array, "[" => "")
    string_array = replace(string_array, "]" => "")
    sous_chaines = split(string_array, ",")

    return parse.(Int, sous_chaines)
end

function string_d_to_array(string_array)
    string_array = replace(string_array, ";" => "")
    string_array = replace(string_array, "]" => "")
    sous_chaines = split(string_array, " ")

    return parse.(Float64, sous_chaines)
end

function read_file(file)
    if isfile(file)
        myFile = open(file)
        data = readlines(myFile)
        n = parse(Int64, data[1][5:end])
        s = parse(Int64, data[2][5:end])
        t = parse(Int64, data[3][5:end])
        S = parse(Int64, data[4][5:end])
        d1 = parse(Int64, data[5][5:end])
        d2 = parse(Int64, data[6][5:end])

        p= string_to_array(data[7][5:end])
        ph = string_to_array(data[8][5:end])

        d=Dict()
        D=Dict()
        for line in data[10:end]
            array_line=string_d_to_array(line)
            i = Int(array_line[1])
            j = Int(array_line[2])            
            d[i,j]=Int(array_line[3])
            D[i,j]=array_line[4]
        end
        # Fermer le fichier
        close(myFile)
        return n, s, t, S, d1, d2, p, ph, d, D
    end
end

function djikstra(n, s, d)
    deltap=Dict()
    deltam=Dict()
    for i in 1:n
        deltap[i]=[]
        deltam[i]=[]
    end
    for (i,j) in keys(d)
        # println(i,j)
        push!(deltap[i],j)
        push!(deltam[j],i)
    end
    infini=500000 
    distance=[infini for i in 1:n]
    distance[s]=0
    Q=[i for i in 1:n]
    while Q!=[]
        mini=infini
        sommet=-1
        index_min=0
        for (index_j,j) in enumerate(Q)
            if distance[j]<mini
                mini = distance[j]
                sommet=j
                index_min=index_j
            end
        end
        if index_min==0
            Q=[]
        else
            deleteat!(Q, index_min)
            for k in deltap[sommet]
                if distance[k]>distance[sommet]+d[sommet,k]
                    distance[k]=distance[sommet]+d[sommet,k]
                end
            end
        end
    end
    return distance
end

function simplify_instance(name_instance, borne_sup, initial_values)
    n, s, t, S, d1, d2, p, ph, d, D = read_file("./data/$name_instance")
    initial_x=Dict(key  => initial_values[1][key] for key in keys(d))
    initial_a=[initial_values[2][i] for i in 1:n]
    progress=true
    while progress
        println(n, " ", s, " ", size(collect(keys(d))))
        distance_s=djikstra(n, s, d)
        dreverse=Dict((key[2], key[1])=>d[key] for key in keys(d))
        distance_t=djikstra(n, t, dreverse)
        distance_tot=[distance_s[i]+distance_t[i] for i in 1:n]
        println(distance_t[s])
        new_set=[i for i in 1:n if distance_tot[i]<borne_sup]
        link_new_set=[findfirst(item -> item==i, new_set) for i in 1:n]
        initial_x=Dict((link_new_set[key[1]], link_new_set[key[2]])  => initial_x[key] for key in keys(d) if distance_s[key[1]]+d[key]*(1+D[key])+distance_t[key[2]]<borne_sup && distance_tot[key[1]]<borne_sup && distance_tot[key[2]]<borne_sup)
        Dcopy=D
        D=Dict((link_new_set[key[1]], link_new_set[key[2]]) => D[key] for key in keys(d) if distance_s[key[1]]+d[key]*(1+D[key])+distance_t[key[2]]<borne_sup && distance_tot[key[1]]<borne_sup && distance_tot[key[2]]<borne_sup)
        d=Dict((link_new_set[key[1]], link_new_set[key[2]])  => d[key] for key in keys(d) if distance_s[key[1]]+d[key]*(1+Dcopy[key])+distance_t[key[2]]<borne_sup && distance_tot[key[1]]<borne_sup && distance_tot[key[2]]<borne_sup)
        p=[p[i] for i in 1:n if distance_tot[i]<borne_sup]
        ph=[ph[i] for i in 1:n if distance_tot[i]<borne_sup]
        initial_a=[initial_a[i] for i in 1:n if distance_tot[i]<borne_sup]
        s=link_new_set[s]
        t=link_new_set[t]
        if n==size(new_set)[1]
            progress=false
        else
            n=size(new_set)[1]
        end
    end
    return  n, s, t, S, d1, d2, p, ph, d, D, [initial_x, initial_a]
end

function calcul_d_D_k_set(p, d, D,S, borne_sup)
    d_D=[(d[key], D[key]) for key in collect(keys(D))]
    count_element=Dict(d_D .=> 0)
    min_p=minimum(p)
    for element in d_D
        if count_element[element]<min(S/min_p, borne_sup/element[1])
            count_element[element] += 1
        end
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

function calcul_d_k_set(p, d, S, borne_sup)
    d_s=[d[key] for key in collect(keys(D))]
    count_element=Dict(d_s .=> 0)
    min_p=minimum(p)
    for element in d_s
        if count_element[element]<min(S/min_p, borne_sup/element[1])
            count_element[element] += 1
        end
    end
    d_set=sort([key for key in collect(keys(count_element))], lt = (x, y) -> (x< y))
    A_d_set=Dict()
    for d_i in d_set
        A_d_set[d_i]=[]
    end
    for key in collect(keys(d))
        push!(A_d_set[d[key]], key)
    end
    A_d_set
    return d_set , count_element, A_d_set
end

function calcul_p_ph_k_set(p, ph, S)
    p_ph=[(p[i], ph[i]) for i in 1:size(p)[1]]
    count_element_p=Dict(p .=> 0)
    count_element_ph=Dict(ph .=> 0)

    for element in p_ph
        if count_element_p[element[1]]<S/element[1]
            count_element_p[element[1]] += 1
            count_element_ph[element[2]] += 1
        end
    end
    p_set=sort([key for key in collect(keys(count_element_p))], lt = (x, y) -> (x < y))
    ph_set=sort([key for key in collect(keys(count_element_ph))], lt = (x, y) -> (x < y))
    A_p_set=Dict()
    for p_i in p_set
        A_p_set[p_i]=[]
    end
    A_ph_set=Dict()
    for ph_i in ph_set
        A_ph_set[ph_i]=[]
    end
    for i in 1:size(p)[1]
        push!(A_p_set[p[i]], i)
        push!(A_ph_set[ph[i]], i)
    end
    return p_set, ph_set, count_element_p, count_element_ph, A_p_set, A_ph_set
end

function calcul_p_k_set(p, S)
    count_element_p=Dict(p .=> 0)

    for element in p
        if count_element_p[element]<S/element
            count_element_p[element] += 1
        end
    end
    p_set=sort([key for key in collect(keys(count_element_p))], lt = (x, y) -> (x < y))
    A_p_set=Dict()
    for p_i in p_set
        A_p_set[p_i]=[]
    end
    for i in 1:size(p)[1]
        push!(A_p_set[p[i]], i)
    end
    return p_set, count_element_p, A_p_set
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