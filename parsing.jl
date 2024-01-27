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