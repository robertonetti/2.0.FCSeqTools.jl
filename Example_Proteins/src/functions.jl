function read_fasta(filename, threshold)
    out = do_number_matrix_prot(do_letter_matrix(filename), threshold)
    f = open(filename)
    lines = readlines(f)
    lines = lines[lines.!=""]
    out_re = []
    out_E = []
    for line in lines
        if line[1] == '>'
            num, re_E = split(line, " re_")
            re, E = split(re_E, " E_")
            push!(out_re, parse(Float64, re))
            push!(out_E, parse(Float64, E))
        end
    end
    return out, out_re, out_E
end

# FULL MODEL ENERGIES
function full_model_energy(q, fields, couplings, MSA, L_MSA)
    full_energies = []
    for i in 1:L_MSA
        seq = MSA[i,:]' 
        freq = freq_single_point(seq, q, 0.0) 
        fij = fij_two_point(seq, q, 0.0)
        push!(full_energies, - sum(fij .* couplings) - sum(freq .* fields))
    end
    return full_energies
end

# PROFILE MODEL ENERGIES
function profile_energy(q, fields, MSA, L_MSA)
    energies = []
    for i in 1:L_MSA
        seq = MSA[i,:]'
        freq = freq_single_point(seq, q, 0.0); 
        push!(energies, - sum(freq .* fields))
    end
    return energies
end


# FULL MODEL SINGLE SITE ENTROPIES
function full_model_site_entropy(q, ref_seq, fields, couplings, pseudo_count, threshold)
    L = length(ref_seq)
    S_sites = zeros(L)
    for i ∈ 1:L
        energies_i = zeros(q)
        for a_i ∈ 1:q
            seq = copy(ref_seq)
            seq[i] = a_i
            freq = freq_reweighted(seq', q, pseudo_count, threshold)  
            fij = fij_reweighted(seq', q, pseudo_count, threshold)
            energies_i[a_i] = - sum(fij .* couplings) - sum(freq .* fields)
        end
        p_cond_i = exp.(energies_i) ./ sum(exp.(energies_i))
        S_sites[i] = - sum( p_cond_i .* log.(p_cond_i))
    end
    return S_sites
end

function mutation_MSA(q, seq)
    MSA_mut1 = []
    push!(MSA_mut1, seq)
    c = 0
    for i in seq
        c += 1
        if i != q
            for a in 1:q - 1
                if a != i 
                    new_seq = copy(seq)
                    new_seq[c] = a
                    push!(MSA_mut1, new_seq)
                end
            end
        end
    end

    L_MSA = length(MSA_mut1)
    L_prot = length(MSA_mut1[1])

    MSA_mut = Array{Int64}(undef, L_MSA, L_prot)
    for i in 1:L_MSA
        for j in 1:L_prot
        MSA_mut[i, j] = MSA_mut1[i][j]
        end

    end

    return MSA_mut
end

# DATO UN DMS MSA CALCOLA LE ENERGIE CON I PARAMETRI FORNITI
function compute_dms_ave_energies(q, MSA, wt_seq, fields, couplings)
    L = length(MSA[1,:])
    energies = zeros(L)
    count = zeros(L)

    freq = freq_single_point(wt_seq', q, 0.0) 
    fij = fij_two_point(wt_seq', q, 0.0)
    wt_ene = - sum(fij .* couplings) - sum(freq .* fields)

    for m in 1:length(MSA[:, 1])
        seq = MSA[m,:]
        i = 1
        flag = false
        while flag == false
            if seq[i] != wt_seq[i] 
                flag = true
            else
                i += 1
            end
        end
        freq = freq_single_point(seq', q, 0.0) 
        fij = fij_two_point(seq', q, 0.0)
        energies[i] += - sum(fij .* couplings) - sum(freq .* fields)
        count[i] += 1
    end
    ave_energies = []

    for i in 1:L
        if energies[i] != 0.0
            push!(ave_energies, energies[i] / count[i])
        end
    end
    return ave_energies
end

# DATO L'MSA DMS DI MARTIN, CALCOLA LE ENERGIE MEDIE DEL LORO MODELLO
function martin_compute_dms_ave_energies(q, MSA, wt_seq, MSA_E)
    L = length(MSA[1,:])
    energies = zeros(L)
    count = zeros(L)
    for m in 1:length(MSA[:, 1])
        seq = MSA[m,:]
        i = 1
        flag = false
        while flag == false
            if seq[i] != wt_seq[i] 
                flag = true
            else
                i += 1
            end
        end
        freq = freq_single_point(seq', q, 0.0) 
        fij = fij_two_point(seq', q, 0.0)
        energies[i] += MSA_E[m]
        count[i] += 1
    end
    ave_energies = []

    for i in 1:L
        if energies[i] != 0.0
            push!(ave_energies, energies[i] / count[i])
        end
    end
    return ave_energies
end


# SAVE MODEL
function print_model_to_file_prot(natural_sequences,Jij,h,filename)
    open(filename, "w") do f
        for i in 1:length(natural_sequences[1,:])
            for j in i+1:length(natural_sequences[1,:])
                for k in 1:21    
                    if k==1
                        k2=1        
                    elseif k==2
                        k2=2
                    elseif k==3
                        k2=3
                    elseif k==4
                        k2=4
                    elseif k==5
                        k2=5
                    elseif k==6
                        k2=6
                    elseif k==7
                        k2=7
                    elseif k==8
                        k2=8
                    elseif k==9
                        k2=9
                    elseif k==10
                        k2=10
                    elseif k==11
                        k2=11
                    elseif k==12
                        k2=12
                    elseif k==13
                        k2=13
                    elseif k==14
                        k2=14
                    elseif k==15
                        k2=15
                    elseif k==16
                        k2=16
                    elseif k==17
                        k2=17
                    elseif k==18
                        k2=18
                    elseif k==19
                        k2=19
                    elseif k==20
                        k2=20
                    elseif k==21
                        k2=0
                    end
                        
                    for l in 1:21
                        if l==1
                            l2=1
                        elseif l==2
                            l2=2
                        elseif l==3
                            l2=3
                        elseif l==4
                            l2=4
                        elseif l==5
                            l2=5
                        elseif l==6
                            l2=6
                        elseif l==7
                            l2=7
                        elseif l==8
                            l2=8
                        elseif l==9
                            l2=9
                        elseif l==10
                            l2=10
                        elseif l==11
                            l2=11
                        elseif l==12
                            l2=12
                        elseif l==13
                            l2=13
                        elseif l==14
                            l2=14
                        elseif l==15
                            l2=15
                        elseif l==16
                            l2=16
                        elseif l==17
                            l2=17
                        elseif l==18
                            l2=18
                        elseif l==19
                            l2=19
                        elseif l==20
                            l2=20
                        elseif l==21
                            l2=0
                        end
                        opo=Jij[i,j,21*(k-1)+l]
                        write(f,"\nJ $(i-1) $(j-1) $(k2) $(l2) $opo" )
                    end
                end
            end
        end
        for i in 1:length(natural_sequences[1,:])
            for j in 1:21          
                if j==1
                    j2=1
                elseif j==2
                    j2=2
                elseif j==3
                    j2=3
                elseif j==4
                    j2=4
                elseif j==5
                    j2=5
                elseif j==6
                    j2=6
                elseif j==7
                    j2=7
                elseif j==8
                    j2=8
                elseif j==9
                    j2=9
                elseif j==10
                    j2=10
                elseif j==11
                    j2=11
                elseif j==12
                    j2=12
                elseif j==13
                    j2=13
                elseif j==14
                    j2=14
                elseif j==15
                    j2=15
                elseif j==16
                    j2=16
                elseif j==17
                    j2=17
                elseif j==18
                    j2=18
                elseif j==19
                    j2=19
                elseif j==20
                    j2=20
                elseif j==21
                    j2=0
                end
                opo=h[21*(i-1)+j]
                write(f,"\nh $(i-1) $(j2) $opo" )
            end
        end
    end
end


function count_common_Jij(Jij_1, Jij_2)
    Jij_common1, Jij_common2 = zeros(size(Jij_1)), zeros(size(Jij_1))
    L, qq = size(Jij_1)[1], size(Jij_1)[3] 
    n_ij, n_ij_ab = 0, 0
    n_01, n_10 = 0, 0
    n_00 = 0
    flag = false
    for j in 1:L
        for i in 1:j
            for idx in 1:qq
                if Jij_1[i, j, idx] != 0.0 && Jij_2[i, j, idx] != 0.0
                    n_ij_ab += 1
                    Jij_common1[i, j, idx], Jij_common2[i, j, idx] = Jij_1[i, j, idx], Jij_2[i, j, idx]
                    flag = true
                elseif Jij_1[i, j, idx] == 0.0 && Jij_2[i, j, idx] != 0.0
                    n_01 += 1
                elseif Jij_1[i, j, idx] != 0.0 && Jij_2[i, j, idx] == 0.0
                    n_10 += 1
                elseif Jij_1[i, j, idx] == 0.0 && Jij_2[i, j, idx] == 0.0
                    n_00 += 1
                end
            end
            if flag == true
                n_ij += 1
                flag = false
            end
        end
    end 
    return n_ij, n_ij_ab, n_01, n_10, n_00, Jij_common1, Jij_common2
end

function energy_space_connectivity(q, gen_seq, h, J, contact_list, site_degree)
    L = length(gen_seq[1,:])
    M = size(gen_seq,1)
    ΔE = zeros(M * L * (q-1))
    for n_seq ∈ 1:M
        seq = gen_seq[n_seq,:]
        for i ∈ 1:L 
            a = seq[i]
            remaining_amino = range(1,q)[range(1,q).!= a]
            idx = 0
            for a1 in remaining_amino
                idx += 1
                δe = h[q * (i - 1) + a] - h[q * (i - 1) + a1]
                for j in contact_list[1: site_degree[i], i]
                    b = seq[j]
                    δe += + J[i,j, q * (a - 1) + b] - J[i,j, q * (a1 - 1) + b] 
                end
                index = (n_seq - 1)*L*(q - 1) + (i - 1)*(q - 1) + idx
                ΔE[index] = δe 
            end
        end
    end
    return ΔE
end

function seq_energy(q, fields, couplings, seq)
    freq = freq_single_point(seq', q, 0.0) 
    fij = fij_two_point(seq', q, 0.0)
return - sum(fij .* couplings) - sum(freq .* fields)
end


function read_example_output(file_path)
    f = open(file_path)
    lines = readlines(f)[2 : end]
    n = length(lines) - 2
    edges, elements, scores = Vector{Array{Int64,1}}(undef, n), Vector{Array{Int64,1}}(undef, n), zeros(n)
    idx = 0
    
    for line in lines[1:n]
        idx += 1
        #println(idx," ", line)
        a, b = split(line, "]   (")
        ed = a[2:end]
        el, sc = split(b, " )")
        sc = split(sc, "{")[2]
        edges[idx] = parse.(Int64, split(ed,"  "))
        elements[idx] = parse.(Int64, split(el," "))
        scores[idx] = parse.(Float64, split(sc, "}")[1])
    end
    return edges, elements, scores
end

function compare_edges(edges1, elements1, edges2, elements2)
    n = min(length(edges1), length(edges2))
    L, q = 96, 21
    x_elements, y_elements = falses(L, L, q*q), falses(L, L, q*q)
    x_edges, y_edges = falses(L, L), falses(L, L) 
    n_edges, n_elements = zeros(n), zeros(n)
    n_edges_x, n_edges_y = zeros(n), zeros(n)
    n_elements_x, n_elements_y = zeros(n), zeros(n)

    for i ∈ 1:n
        edge_x, edge_y = edges1[i], edges2[i]
        x_edges[edge_x[1], edge_x[2]] = true
        y_edges[edge_y[1], edge_y[2]] = true
        x_elements[edge_x[1], edge_x[2], elements1[i]] .= true 
        y_elements[edge_y[1], edge_y[2], elements2[i]] .= true
        n_edges[i] = sum(x_edges .* y_edges) 
        n_elements[i] = sum(x_elements .* y_elements) 
        n_edges_x[i], n_edges_y[i] = sum(x_edges), sum(y_edges)
        n_elements_x[i], n_elements_y[i] = sum(x_elements), sum(y_elements)
    end
    return n_edges, n_elements, n_edges_x, n_edges_y, n_elements_x, n_elements_y
end

function read_model(filepath)
    file = open(filepath, "r")
    contenuto = read(file, String)[1:end-1]
    close(file)
    contenuto = split(contenuto, "\n")
    #println(contenuto[end])
    L, q = 96, 21
    J, h = zeros(L, L, q*q), zeros(L * q);

    for element ∈ contenuto
        #println(element)
        if string(element[1]) == "J"
            usl, i, j, a, b, value = split(element, " ")
            i, j, a, b = parse(Int64,i) + 1, parse(Int64,j) + 1, parse(Int64,a) + 1, parse(Int64,b) + 1
            #println(i, j, a, b)
            J[i, j, (a-1)*q + b] = parse(Float64,value)
        elseif string(element[1]) == "h"
            usl, i, a, value = split(element, " ")
            i, a = parse(Int64,i) + 1, parse(Int64,a) + 1
            h[ q*(i - 1) + a] = parse(Float64,value)
        end
    end
    return J, h
end