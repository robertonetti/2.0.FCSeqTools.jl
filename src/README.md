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
Computes the site entropies for a number matrix with $q$ states per site. The "pseudo_count" is applied to the frequencies needed to the calculation of the entropies. See `weight_of_sequences` for the "threshold" parameter

- `profile_model_entropy(matrix,q,pseudo_count,threshold)`  
Computes the profile model entropy for a number matrix with $q$ states per site. The "pseudo_count" is applied to the frequencies needed to the calculation of the entropy. See `weight_of_sequences` for the "threshold" parameter

- `max_kl_divergence(fij,pij)`  
Takes as input two two-point frequencies matrices and outputs  the $\underset{ij}{argmax} \. D_{kl}( f_{ij} || p_{ij} )$ and the $\underset{ij}{max}\. D_{kl}( f_{ij} || p_{ij} )$


## Generative models & Sampling 


- `profile_model_generation(threshold,q,pseudo_count,number,number_matrix)`  
Learns a profile model from a q-states number matrix and samples "number" sequences from it.  The "pseudo_count" and "threshold" are applied to the frequencies needed for the computation of the model parameters. See `weight_of_sequences` for the "threshold" parameter

- `gibbs_sampling(q,h_local,Jij,sequences,site_degree,contact_list,sweeps)`  
Performs N="sweeps" Gibbs MonteCarlo sweeps on a a q-states number matrix specified by "sequences". The "Jij", "h_local", "contact_list" and "site_degree" parameters are used to specify the model and are compatible with the outputs of the `E_A_A` function  

- `E_A_A(q,pseudo_count,number,number_matrix,filename)`  
Applies the Edge Activation Algorithm to generate "number" sequences learning the statistics from the input q-states number matrix. The "pseudo_count" parameter and a similarity threshold of 0.8 are used to learn the statistics. The algorithm is applied for 10000 iteration or until it reaches a model with a two-point Pearson score greater than 0.95. It outputs in order:  "score_vector" that keeps track of the Pearson score each 15 iteration, "likelihood_gain_vector" that keeps track of the likelihood gain each 15 iterations, "sequences" that are the sequences generated by the last model of the iterations, "Jij_couplings" that is the matrix of  the interaction coupling of the model, the coupling of state $s$ on site $i$ and state $r$ on site $j$ (with $i$ < $j$) is the $[i,j,q\cdot (s-1) + r]$ element of the matrix, "h_local" is the vector of local fields and the field of state $s$ on site $i$ is the $[q\cdot(i-1)+s]$ element of the vector, "contact_list" is a list and its i-th elements contains all the sites that interact with the i-th one, "site_degree" is a vector that contains the interaction degree of each site "edge_list" is a matrix with two columns and the two elements of the i-th row are the interaction pair added in the i-th iteration. The progress of the algorithm is both printed and stored in a file named "filemane"

- `rna_cm_model_generation(threshold,pseudo_count,number,number_matrix,ss_contact_matrix)	`  
Learns an RNA covariance model from a number matrix and its associated secondary structure contact matrix ( 1 if contact 0 if not ) and samples "number" sequences from it.  The "pseudo_count" and "threshold" are applied to the frequencies needed for the computation of the model parameters. See `weight_of_sequences` for the "threshold" parameter

- `reweighted_sample(number_matrix,number,threshold)`  
it outputs a sample of "number" sequences from the number matrix that respects the reweighted statistics. See `weight_of_sequences` for the "threshold" parameter




## Utils


- `ss_matrix_to_dot_bracket(ss_contact_matrix)`  
Takes a secondary structure contact matrix ( 1 if contact 0 if not ) and outputs it in the dot bracket format

- `dot_bracket_to_ss_matrix(dot_bracket_ss)`  
Takes a secondary structure in the dot bracket format and outputs the structure contact list and the structure contact matrix ( 1 if contact 0 if not )


