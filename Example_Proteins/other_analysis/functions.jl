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