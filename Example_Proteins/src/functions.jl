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
    full_energies = zeros(L_MSA)
    for i in 1:L_MSA
        seq = MSA[i,:]' 
        freq = freq_single_point(seq, q, 0.0) 
        fij = fij_two_point(seq, q, 0.0)
        full_energies[i] = - sum(fij .* couplings) - sum(freq .* fields)
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
function full_model_site_entropy(q, ref_seq, h, J, pseudo_count, contact_list, site_degree)
    L = length(ref_seq)
    S_sites = zeros(L)

    # compute energy of the reference sequence
    freq =  freq_single_point(ref_seq', q, pseudo_count)  
    fij = fij_two_point(ref_seq', q, pseudo_count)
    ref_energy = - sum(fij .* J) - sum(freq .* h)

    for i ∈ 1:L
        energies_i = zeros(q-1)
        a = ref_seq[i]
        remaining_amino = range(1,q)[range(1,q).!= a]
        idx = 0
        for a_new ∈ remaining_amino
            idx += 1
            δe = h[q * (i - 1) + a] - h[q * (i - 1) + a_new]
            for j in contact_list[1: site_degree[i], i]
                b = ref_seq[j]
                δe += + J[i,j, q * (a - 1) + b] - J[i,j, q * (a_new - 1) + b] 
            end
            energies_i[idx] = ref_energy + δe
        end
        p_cond_i = exp.(energies_i) ./ sum(exp.(energies_i))
        S_sites[i] = - sum( p_cond_i .* log.(p_cond_i))
    end
    return S_sites
end

function make_DMS_and_compute_energy(q, ref_seq, h, J, contact_list, site_degree)
    DMS = []
    DMS_energies = []
    L = length(ref_seq)
    # compute energy of the reference sequence
    freq =  freq_single_point(ref_seq', q, 0)  
    fij = fij_two_point(ref_seq', q, 0)
    ref_energy = - sum(fij .* J) - sum(freq .* h)
    # push the WT
    push!(DMS, seq)
    push!(DMS_energies, ref_energy)

    for i in 1:L
        a = ref_seq[i]
        if a != q # to be sure it is not a gap
            remaining_amino = range(1,q-1)[range(1,q-1).!= a]
            for a_new in remaining_amino
                new_seq = copy(ref_seq)
                new_seq[i] = a_new
                push!(DMS, new_seq)
                # compute energy difference
                δe = h[q * (i - 1) + a] - h[q * (i - 1) + a_new]
                for j in contact_list[1: site_degree[i], i]
                    b = ref_seq[j]
                    δe += + J[i,j, q * (a - 1) + b] - J[i,j, q * (a_new - 1) + b] 
                end
                push!(DMS_energies, ref_energy + δe)
            end
        end
    end

    L_MSA = length(DMS)
    L_prot = length(DMS[1])

    DMS_new = Array{Int64}(undef, L_MSA, L_prot)
    for i in 1:L_MSA
        for j in 1:L_prot
            DMS_new[i, j] = DMS[i][j]
        end
    end
    return MSA_mut, DMS_energies
end

# DATO UN DMS MSA CALCOLA LE ENERGIE CON I PARAMETRI FORNITI
function compute_dms_ave_energies(q, MSA, wt_seq, h, J, contact_list, site_degree)
    L, M = length(MSA[1,:]), length(MSA[:, 1])
     
    energies = []
    ave_energies1 = zeros(L)
    count = zeros(L)

    # compute energy of the reference sequence
    freq = freq_single_point(wt_seq', q, 0.0) 
    fij = fij_two_point(wt_seq', q, 0.0)
    wt_energy = - sum(fij .* J) - sum(freq .* h)

    for m ∈ 1:M
        seq = MSA[m,:]
        # find the mutated position 
        i, flag = 1, false
        while flag == false
            seq[i] != wt_seq[i] ? flag = true : i += 1
        end
        # compute energy difference
        a, a_new = wt_seq[i], seq[i]
        δe = h[q * (i - 1) + a] - h[q * (i - 1) + a_new]
        for j in contact_list[1: site_degree[i], i]
            b = wt_seq[j]
            δe += + J[i,j, q * (a - 1) + b] - J[i,j, q * (a_new - 1) + b] 
        end
        ave_energies1[i] += wt_energy + δe
        push!(energies, wt_energy + δe)
        count[i] += 1
    end
    ave_energies = []
    for i in 1:L
        if energies[i] != 0.0
            push!(ave_energies, ave_energies1[i] / count[i])
        end
    end
    return energies, ave_energies
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
            seq[i] != wt_seq[i] ? flag = true : i += 1
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
    
    for j in 1:L
        for i in 1:j - 1
            for idx in 1:qq
                if Jij_1[i, j, idx] != 0.0 && Jij_2[i, j, idx] != 0.0
                    n_ij_ab += 1
                    Jij_common1[i, j, idx], Jij_common2[i, j, idx] = Jij_1[i, j, idx], Jij_2[i, j, idx]
                elseif Jij_1[i, j, idx] == 0.0 && Jij_2[i, j, idx] != 0.0
                    n_01 += 1
                elseif Jij_1[i, j, idx] != 0.0 && Jij_2[i, j, idx] == 0.0
                    n_10 += 1
                elseif Jij_1[i, j, idx] == 0.0 && Jij_2[i, j, idx] == 0.0
                    n_00 += 1
                end
            end
        end
    end 
    n_ij = count(!iszero, sum(abs.(Jij_1[:,:,i]) for i in 1:size(Jij_1, 3)) .* sum(abs.(Jij_2[:,:,i]) for i in 1:size(Jij_2, 3)))
    return n_ij, n_ij_ab, n_01, n_10, n_00, Jij_common1, Jij_common2
end

function count_common_Jij_three_models(Jij_1, Jij_2, Jij_3)
    Jij_common1, Jij_common2, Jij_common3 = zeros(size(Jij_1)), zeros(size(Jij_1)), zeros(size(Jij_1))
    Jij_12, Jij_13, Jij_23 =  zeros(size(Jij_1)), zeros(size(Jij_1)), zeros(size(Jij_1))
    Jij_123 = zeros(size(Jij_1))
    L, qq = size(Jij_1)[1], size(Jij_1)[3] 

    n_000, n_100, n_010, n_001, n_110, n_101, n_011, n_111 = 0, 0, 0, 0, 0, 0, 0, 0, 0
    n_12, n_13, n_23 = 0, 0, 0
    n_ij = 0
    
    



    for j in 1:L
        for i in 1:j - 1
            for idx in 1:qq
                if Jij_1[i, j, idx] != 0.0 && Jij_2[i, j, idx] != 0.0 && Jij_3[i, j, idx] != 0.0
                    n_111 += 1
                    Jij_common1[i, j, idx], Jij_common2[i, j, idx], Jij_common3[i, j, idx] = Jij_1[i, j, idx], Jij_2[i, j, idx], Jij_3[i, j, idx]
                    Jij_123[i, j, idx] = 1
                    n_12 += 1
                    n_13 += 1
                    n_23 += 1
                elseif Jij_1[i, j, idx] != 0.0 && Jij_2[i, j, idx] == 0.0 && Jij_3[i, j, idx] == 0.0
                    n_100 += 1
                elseif Jij_1[i, j, idx] == 0.0 && Jij_2[i, j, idx] != 0.0 && Jij_3[i, j, idx] == 0.0
                    n_010 += 1
                elseif Jij_1[i, j, idx] == 0.0 && Jij_2[i, j, idx] == 0.0 && Jij_3[i, j, idx] != 0.0
                    n_001 += 1
                elseif Jij_1[i, j, idx] != 0.0 && Jij_2[i, j, idx] != 0.0 && Jij_3[i, j, idx] == 0.0
                    n_110 += 1
                    Jij_12[i, j, idx] = 1
                    n_12 += 1
                elseif Jij_1[i, j, idx] != 0.0 && Jij_2[i, j, idx] == 0.0 && Jij_3[i, j, idx] != 0.0
                    Jij_13[i, j, idx] = 1
                    n_101 += 1
                    n_13 += 1
                elseif Jij_1[i, j, idx] == 0.0 && Jij_2[i, j, idx] != 0.0 && Jij_3[i, j, idx] != 0.0
                    Jij_23[i, j, idx] = 1
                    n_011 += 1
                    n_23 += 1
                elseif Jij_1[i, j, idx] == 0.0 && Jij_2[i, j, idx] == 0.0 && Jij_3[i, j, idx] == 0.0
                    n_000 += 1
                end
            end
        end
    end 
    n_ij = count(!iszero, sum(abs.(Jij_1[:,:,i]) for i in 1:size(Jij_1, 3)) .* sum(abs.(Jij_2[:,:,i]) for i in 1:size(Jij_2, 3)) .* sum(abs.(Jij_3[:,:,i]) for i in 1:size(Jij_3, 3)))
    n_ij_12 = count(!iszero, sum(abs.(Jij_1[:,:,i]) for i in 1:size(Jij_1, 3)) .* sum(abs.(Jij_2[:,:,i]) for i in 1:size(Jij_2, 3)))
    n_ij_13 = count(!iszero, sum(abs.(Jij_1[:,:,i]) for i in 1:size(Jij_1, 3)) .* sum(abs.(Jij_3[:,:,i]) for i in 1:size(Jij_3, 3)))
    n_ij_23 = count(!iszero, sum(abs.(Jij_2[:,:,i]) for i in 1:size(Jij_2, 3)) .* sum(abs.(Jij_3[:,:,i]) for i in 1:size(Jij_3, 3)))
    return n_000, n_100, n_010, n_001, n_110, n_101, n_011, n_111, n_12, n_13, n_23, n_ij, n_ij_12, n_ij_13, n_ij_23,  Jij_common1, Jij_common2, Jij_common3, Jij_12, Jij_13, Jij_23, Jij_123
end



function energy_space_connectivity(q, gen_seq, h, J, contact_list, site_degree)
    L = length(gen_seq[1,:])
    M = size(gen_seq,1)
    ΔE = zeros(M * L * (q-1))
    ΔE_min = zeros(M * L)
    ΔE_mean = zeros(M * L)
    p_tra_mean = zeros(M * L)
    idx_min = 0
    for n_seq ∈ 1:M
        seq = gen_seq[n_seq,:]
        for i ∈ 1:L 
            a = seq[i]
            remaining_amino = range(1,q)[range(1,q).!= a]
            idx = 0
            idx_min += 1
            δe_min = Inf
            for a1 in remaining_amino
                idx += 1
                δe = h[q * (i - 1) + a] - h[q * (i - 1) + a1]
                for j in contact_list[1: site_degree[i], i]
                    b = seq[j]
                    δe += + J[i,j, q * (a - 1) + b] - J[i,j, q * (a1 - 1) + b] 
                end
                index = (n_seq - 1)*L*(q - 1) + (i - 1)*(q - 1) + idx
                ΔE[index] = δe 
                ΔE_mean[idx_min] += δe
                p_tra_mean[idx_min] += exp(-δe)
                if δe_min >= δe 
                    δe_min = δe
                end
            end
            ΔE_mean[idx_min] / length(remaining_amino)
            p_tra_mean[idx_min]/ length(remaining_amino)
            ΔE_min[idx_min] = δe_min
        end
    end
    return ΔE, ΔE_min, ΔE_mean, p_tra_mean
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

function compare_edges(edges1, elements1, edges2, elements2, L)
    n = min(length(edges1), length(edges2))
    q = 96
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
            i, j, a, b = parse(Int64,i) + 1, parse(Int64,j) + 1, parse(Int64,a), parse(Int64,b)
            if a == 0
                a = 21
            elseif b == 0
                b = 21
            end
            #println(i, j, a, b)
            J[i, j, (a-1)*q + b] = parse(Float64,value)
        elseif string(element[1]) == "h"
            usl, i, a, value = split(element, " ")
            i, a = parse(Int64,i) + 1, parse(Int64,a)
            if a == 0
                a = 21
            end
            h[q*(i - 1) + a] = parse(Float64,value)
        end
    end
    return J, h
end

function distance_from_ref_seq(ref_seq, MSA)
    M, L = size(MSA,1), size(MSA,2)
    println(M)
    distance = zeros(M)
    seq_identity = zeros(M)
    for i ∈ 1:M
        common = sum(ref_seq .== nat_MSA[i,:])
        distance[i] = L - common
        seq_identity[i] = common / L 
    end
    return distance, seq_identity
end    

function sequence_identity(MSA_nat, MSA)
    M_nat, M = size(MSA_nat,1), size(MSA, 1)
    L = size(MSA_nat, 2)
    if L != size(MSA, 2)
        error("The two MSAs have different protein length")
    end
    distance = zeros(M)
    seq_identity = zeros(M)
    for m ∈ 1:M
        min = Inf
        for m_nat ∈ 1:M_nat
            val = L - sum(MSA[m, :] .== MSA_nat[m_nat,:])
            #println(val)
            if min >= val
                #println("ok", m, " ", m_nat)
                min = val
            end
        end
        distance[m] = min
        seq_identity[m] = 1 - distance[m]/L
    end
    return distance, seq_identity
end


function sequence_identity_NAT(MSA_nat)
    M_nat, M = size(MSA_nat,1), size(MSA_nat, 1)
    L = size(MSA_nat, 2)
    distance = zeros(M)
    seq_identity = zeros(M)
    for m ∈ 1:M
        min = Inf
        for m_nat ∈ 1:M_nat
            val = L - sum(MSA_nat[m, :] .== MSA_nat[m_nat,:])
            if min >= val && val != 0 
                min = val
            end
        end
        distance[m] = min
        seq_identity[m] = 1 - distance[m]/L
    end
    return distance, seq_identity
end

function counting_different_pairings(MSA, nat_MSA, q)
    pseudo_count = 0.001
    threshold = Float32(pseudo_count /(q*q))
    different_pairs = zeros(size(MSA, 1))
    nat_fij = fij_two_point(nat_MSA, q, pseudo_count)   
    #println("natural sequences fij --- size: ",size(nat_fij), " non-zero elements: ", count(!iszero, nat_fij), " threshold elements: ", count(==(threshold), nat_fij))

    for i ∈ 1:size(MSA,1)
        seq = MSA[i,:]

        fij = fij_two_point(seq', q, 0)
        #println("generated sequences fij --- size: ",size(fij), " non-zero elements: ", count(!iszero, fij), " threshold elements: ", count(==(threshold), fij))
        
        fij =  nat_fij .* fij
        #println("product sequences fij --- size: ",size(fij), " non-zero elements: ", count(!iszero, fij), " threshold elements: ", count(==(threshold), fij))

        different_pairs[i] = count(==(threshold), fij)
    end
    return different_pairs
end