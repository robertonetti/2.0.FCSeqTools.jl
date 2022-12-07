# Functions

## Read & Write 


- `do_letter_matrix(filename)`   
Reads a Fasta file and converts it in a letter matrix format

- `do_number_matrix_rna(letter_matrix,threshold)`  
Takes an RNA letter matrix and converts it in a number matrix. Sequences with a % of alignment gaps greater than "threshold" are excluded

- `do_number_matrix_prot(letter_matrix,threshold)`    
Takes a protein letter matrix and converts it in a number matrix. Sequences with a % of alignment gaps greater than "threshold" are excluded

- `print_fasta_to_file_rna(number_matrix,filename,name)`   
Takes an RNA number matrix and prints it in a file named "filename" in Fasta format. The sequences are labelled with ">i_name"

- `print_fasta_rna(number_matrix,name)`  
Takes an RNA number matrix and prints it in Fasta format. The sequences are labelled with ">i_name"

- `one_hot_encode(number_matrix,q)`   
One hot encodes a number matrix with q states ( $q=5$ for RNA and $q=21$ for proteins)


## Statistics


- `freq_single_point(number_matrix,q,pseudo_count)`  
Computes the site frequencies for a number matrix with $q$ states per site. The frequencies of all the configurations that never appear in the data are set to $\frac{pseudo\textunderscore count}{q}$. The output is a vector and the frequency of state $s$ on site $i$ is the $[q\cdot(i-1)+s]$ element of the vector

- `fij_two_point(number_matrix,q,pseudo_count) `  
Computes the two-point site frequencies for a number matrix with $q$ states per site. The frequencies of all the configurations that never appear in the data are set to $\frac{pseudo\textunderscore count}{q^2}$. The output is a matrix and the two-point frequency of state $s$ on site $i$ and state $r$ on site $j$ (with $i$ < $j$) is the $[i,j,q\cdot (s-1) + r]$ element of the matrix

- `correlation_two_point(number_matrix,q,pseudo_count)`  
Computes the two-point correlations for a number matrix with $q$ states per site. the "pseudo_count" is applied to the frequencies needed to the calculation of the correlations. The output is a vector and the two-point correlation of state $s$ on site $i$ and state $r$ on site $j$ (with $i$ < $j$) is the $\bigg[q^2\bigg( L\cdot (i-1) + \frac{i\cdot(i-1)}{2}\bigg)+q^2\cdot(j-1)+q\cdot (s-1) + r\bigg]$ element of the vector

- `weight_of_sequences(number_matrix,threshold)`  
Computes the relative weight of the sequencies. The weight of the sequencies that have a lot of similars in the data-set is reduced. Two sequencies are considered similar if they share a percentage of sites greater than "threshold"

- `freq_reweighted(number_matrix,q,pseudo_count,threshold) `  
Computes the reweighted site frequencies for a number matrix with $q$ states per site. The frequencies of all the configurations that never appear in the data are set to $\frac{pseudo\textunderscore count}{q}$. The output is a vector and the frequency of state $s$ on site $i$ is the $[q\cdot(i-1)+s]$ element of the vector. See `weight_of_sequences` for the "threshold" parameter


- `fij_reweighted(number_matrix,q,pseudo_count,threshold) `  
Computes the reweighted two-point site frequencies for a number matrix with $q$ states per site. The frequencies of all the configurations that never appear in the data are set to $\frac{pseudo\textunderscore count}{q^2}$. The output is a matrix and the two-point frequency of state $s$ on site $i$ and state $r$ on site $j$ (with $i$ < $j$) is the $[i,j,q\cdot (s-1) + r]$ element of the matrix. See `weight_of_sequences` for the "threshold" parameter

- `correlation_reweighted(number_matrix,q,pseudo_count,threshold)  `  
Computes the reweighted two-point correlations for a number matrix with $q$ states per site. the "pseudo_count" is applied to the frequencies needed to the calculation of the correlations. The output is a vector and the two-point correlation of state $s$ on site $i$ and state $r$ on site $j$ (with $i$ < $j$) is the $\bigg[q^2\bigg( L\cdot (i-1) + \frac{i\cdot(i-1)}{2}\bigg)+q^2\cdot(j-1)+q\cdot (s-1) + r\bigg]$ element of the vector. See `weight_of_sequences` for the "threshold" parameter


- `eff_size_family(number_matrix,threshold)`  
Computes the effective size of the dataset taking in account sequence similarity. See `weight_of_sequences` for the "threshold" parameter

- `site_entropy_vector(matrix,q,pseudo_count,threshold)`  
Computes the site entropies for a number matrix with $q$ states per site. the "pseudo_count" is applied to the frequencies needed to the calculation of the entropies. See `weight_of_sequences` for the "threshold" parameter

- `profile_model_entropy(matrix,q,pseudo_count,threshold)`  
Computes the profile model entropy for a number matrix with $q$ states per site. the "pseudo_count" is applied to the frequencies needed to the calculation of the entropy. See `weight_of_sequences` for the "threshold" parameter

- `max_kl_divergence(fij,pij)`  
Takes as input two two-point frequencies matrices and outputs  the $argmax_{ij} D_{kl}( f_{ij} || p_{ij} )$ and the $max_{ij} D_{kl}( f_{ij} || p_{ij} )$


## Generative models & Sampling 


- `profile_model_generation(threshold,q,pseudo_count,number,number_matrix)`  
Predict conserved RNA## Executable Programsfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff - - 
- `gibbs_sampling(q,h_local,Jij,sequences,site_degree,contact_list,sweeps)`  
Predict conserved RNA## Executable Programsfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff -
- `E_A_A(q,pseudo_count,number,number_matrix,filename)`  
Predict conserved RNA## Executable Programsfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff -
- `rna_cm_model_generation(threshold,pseudo_count,number,number_matrix,ss_contact_matrix)	`  
Predict conserved RNA## Executable Programsfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff 
- `reweighted_sample(number_matrix,number,threshold)`  
Predict conserved RNA## Executable Programsfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff 


## Utils


- `ss_matrix_to_dot_bracket(ss_contact_matrix)`  
Takes a secondary structure contact matrix ( 1 if contact 0 if not ) and outputs it in the dot bracket format

- `dot_bracket_to_ss_matrix(dot_bracket_ss)`  
Takes a secondary structure in the dot bracket format and outputs the contact list and secondary structure contact matrix ( 1 if contact 0 if not )


