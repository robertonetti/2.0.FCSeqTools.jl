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
        #println("i : ", i)
        energies_i = zeros(q)
        for a_i ∈ 1:q
            seq = copy(ref_seq)
            seq[i] = a_i
            freq = freq_reweighted(seq', q, pseudo_count, threshold)  
            fij = fij_reweighted(seq', q, pseudo_count, threshold)
            energies_i[a_i] = - sum(fij .* couplings) - sum(freq .* fields)
        end

        #println("energies_i = ", energies_i)

        p_cond_i = exp.(energies_i) ./ sum(exp.(energies_i))

        #println("p_cond_i = ", p_cond_i)

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