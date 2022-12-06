module FCSeqTools



using Distances
using StatsBase
using Random
using Statistics



function do_letter_matrix(filename)
f=open(filename)
    lines=readlines(f)
    letter_matrix=Array{Char,2}(undef,Int64(1),Int64(0))
    for i in 1:Int64(length(lines))
      if lines[i][1]=='>'
        j=1
        temp=join(lines[i+j])
        while j!=-1
            if i+j+1 <= length(lines)
        	if lines[i+j+1][1]!='>'
                temp=temp*join(lines[i+j+1])
                j=j+1
            else j=-1
            end
            else j=-1
            end
        end
	if i==1
		letter_matrix=hcat(letter_matrix,reshape(collect(temp),1,length(collect(temp))))
	else
        letter_matrix=vcat(letter_matrix,reshape(collect(temp),1,length(collect(temp))))
    end
    end
    end
  return letter_matrix
end



function do_number_matrix_rna(letter_matrix,threshold)
    n_columns=length(letter_matrix[1,:])
    n_rows=length(letter_matrix[:,1])
    number_matrix=zeros(Int8,n_rows,n_columns)
    for i in 1:n_rows
        for j in 1:n_columns
            if letter_matrix[i,j]=='A'
                number_matrix[i,j]=1
            elseif letter_matrix[i,j]=='C'
                number_matrix[i,j]=2
            elseif letter_matrix[i,j]=='G'
                number_matrix[i,j]=3
            elseif letter_matrix[i,j]=='U' || letter_matrix[i,j]=='T' 
                number_matrix[i,j]=4
            elseif letter_matrix[i,j]=='-'
                number_matrix[i,j]=5
            end
        end
    end
    i=1
    while i<=length(number_matrix[:,1])
	   if 0 in number_matrix[i,:]
		number_matrix=number_matrix[setdiff(1:end, i), :]
	   else
    		i=i+1
    	   end
    end  
    i=1 
    while i<=length(number_matrix[:,1])
	   if length(number_matrix[i,:][number_matrix[i,:].==5])/n_columns>=threshold
		number_matrix=number_matrix[setdiff(1:end, i), :]
	   else
	    	i=i+1
	   end
    end
    return number_matrix
end



function do_number_matrix_prot(letter_matrix,threshold)
    columns=length(letter_matrix[1,:])
    rows=length(letter_matrix[:,1])
    number_matrix=zeros(Int8,rows,columns)
    for i in 1:rows
        for j in 1:columns
            if letter_matrix[i,j]=='A'
                number_matrix[i,j]=1
            elseif letter_matrix[i,j]=='C'
                number_matrix[i,j]=2
            elseif letter_matrix[i,j]=='D'
                number_matrix[i,j]=3
            elseif letter_matrix[i,j]=='E'
                number_matrix[i,j]=4
            elseif letter_matrix[i,j]=='F'
                number_matrix[i,j]=5
            elseif letter_matrix[i,j]=='G'
                number_matrix[i,j]=6
            elseif letter_matrix[i,j]=='H'
                number_matrix[i,j]=7
            elseif letter_matrix[i,j]=='I'
                number_matrix[i,j]=8
            elseif letter_matrix[i,j]=='K'
                number_matrix[i,j]=9
            elseif letter_matrix[i,j]=='L'
                number_matrix[i,j]=10
            elseif letter_matrix[i,j]=='M'
                number_matrix[i,j]=11
            elseif letter_matrix[i,j]=='N'
                number_matrix[i,j]=12
            elseif letter_matrix[i,j]=='P'
                number_matrix[i,j]=13
            elseif letter_matrix[i,j]=='Q'
                number_matrix[i,j]=14
            elseif letter_matrix[i,j]=='R'
                number_matrix[i,j]=15
            elseif letter_matrix[i,j]=='S'
                number_matrix[i,j]=16
            elseif letter_matrix[i,j]=='T'
                number_matrix[i,j]=17
            elseif letter_matrix[i,j]=='V'
                number_matrix[i,j]=18
            elseif letter_matrix[i,j]=='W'
                number_matrix[i,j]=19
            elseif letter_matrix[i,j]=='Y'
                number_matrix[i,j]=20
            elseif letter_matrix[i,j]=='-'
                number_matrix[i,j]=21
            end
        end
    end
    i=1
    while i<=length(number_matrix[:,1])
	   if 0 in number_matrix[i,:]
		number_matrix=number_matrix[setdiff(1:end, i), :]
	   else
    		i=i+1
    	   end
    end  
    i=1 
    while i<=length(number_matrix[:,1])
	   if length(number_matrix[i,:][number_matrix[i,:].==21])/columns>=threshold
		number_matrix=number_matrix[setdiff(1:end, i), :]
	   else
	    	i=i+1
	   end
    end
    return number_matrix
end



function one_hot_encode(number_matrix,q)
    one_hot_encoded_data=zeros(Int8,length(number_matrix[:,1]),(q-1)*length(number_matrix[1,:]))
    for i in 1:length(number_matrix[:,1])
        for j in 1:length(number_matrix[1,:])
            if number_matrix[i,j]!=q
            one_hot_encoded_data[i,(q-1)*(j-1)+(number_matrix[i,j])]=1
            end
        end
    end          
    return one_hot_encoded_data
end



function freq_single_point(number_matrix,q,pseudo_count)
    frequencies=zeros(Float32,q*length(number_matrix[1,:]))
    for j in 1:length(number_matrix[:,1])
    	for i in 1:length(number_matrix[1,:])
    		frequencies[(i-1)*q+number_matrix[j,i]]+=1
   	end
    end
    frequencies=frequencies/length(number_matrix[:,1]) 
    return (1-Float32(pseudo_count))*frequencies.+Float32(pseudo_count/q)
end



function fij_two_point(number_matrix,q,pseudo_count)   
    fij=zeros(Float32,length(number_matrix[1,:]),length(number_matrix[1,:]),q^2)
    for i1 in 1:length(number_matrix[:,1])
        for i2 in 1: length(number_matrix[1,:])
            for i3 in i2+1:length(number_matrix[1,:])
                fij[i2,i3,(number_matrix[i1,i2]-1)*q+number_matrix[i1,i3]]+=1
            end
        end
    end
    fij=fij/length(number_matrix[:,1])
    return  (1-Float32(pseudo_count))*fij.+Float32((pseudo_count)/(q^2))
end



function correlation_two_point(number_matrix,q,pseudo_count)  
    L=length(number_matrix[1,:]) 
    fij=fij_two_point(number_matrix,q,pseudo_count)
    frequencies=freq_single_point(number_matrix,q,pseudo_count)
    correlation=zeros(Float32,q*q*Int64(L*(L-1)/2))
    counter=1
    for i in 1:length(number_matrix[1,:])
        for j in i+1:length(number_matrix[1,:])
            for k1 in 1:q
                for k2 in 1:q
                correlation[counter]=fij[i,j,(k1-1)*q+k2]-(frequencies[(i-1)*q+k1]*frequencies[(j-1)*q+k2])
                counter+=1
                end
            end
        end
    end
  return  correlation
end



function weight_of_sequences(number_matrix,threshold)
    N=length(number_matrix[:,1])
    L=length(number_matrix[1,:])
    weights=ones(Float64,N)
    for i in 1:N
        for j in i+1:N
            if (L-hamming(number_matrix[i,:],number_matrix[j,:]))/L>=threshold
                weights[i]+=1
                weights[j]+=1
            end
        end
    end
return weights.^(-1)
end



function fij_reweighted(number_matrix,q,pseudo_count,threshold) 
    weight=weight_of_sequences(number_matrix,threshold)
    fij=zeros(Float32,length(number_matrix[1,:]),length(number_matrix[1,:]),q^2)
    for i1 in 1:length(number_matrix[:,1])
        for i2 in 1: length(number_matrix[1,:])
            for i3 in i2+1:length(number_matrix[1,:])
                fij[i2,i3,(number_matrix[i1,i2]-1)*q+number_matrix[i1,i3]]+=weight[i1]
            end
        end
    end
    fij=fij/sum(weight)
    return  (1-Float32(pseudo_count))*fij.+Float32((pseudo_count)/(q^2))
end



function freq_reweighted(number_matrix,q,pseudo_count,threshold)
    weight=weight_of_sequences(number_matrix,threshold)
    frequencies=zeros(Float32,q*length(number_matrix[1,:]))
    for j in 1:length(number_matrix[:,1])
    	for i in 1:length(number_matrix[1,:])
    		frequencies[(i-1)*q+number_matrix[j,i]]+=weight[j]
   	end
    end
    frequencies=frequencies/sum(weight) 
    return (1-Float32(pseudo_count))*frequencies.+Float32(pseudo_count/q)
end



function correlation_reweighted(number_matrix,q,pseudo_count,threshold)  
    L=length(number_matrix[1,:]) 
    fij=fij_reweighted(number_matrix,q,pseudo_count,threshold) 
    frequencies=freq_reweighted(number_matrix,q,pseudo_count,threshold)
    correlation=zeros(Float32,q*q*Int64(L*(L-1)/2))
    counter=1
    for i in 1:length(number_matrix[1,:])
        for j in i+1:length(number_matrix[1,:])
            for k1 in 1:q
                for k2 in 1:q
                correlation[counter]=fij[i,j,(k1-1)*q+k2]-(frequencies[(i-1)*q+k1]*frequencies[(j-1)*q+k2])
                counter+=1
                end
            end
        end
    end
  return  correlation
end



function eff_size_family(number_matrix,threshold)
	return sum(weight_of_sequences(number_matrix,threshold))
end



function max_kl_divergence(fij,pij)
    L=length(fij[1,:,1])
    kl_matrix=zeros(Float32,L,L)
    for i in 1:L
        for j in i+1:L
            kl_matrix[i,j]=Distances.kl_divergence(fij[i,j,:],pij[i,j,:])
        end
    end
    return  argmax(kl_matrix),maximum(kl_matrix)
end



function print_fasta_to_file_rna(number_matrix,filename,name)
	open(filename, "w") do f
	for i in 1:length(number_matrix[:,1])
		if i==1
		    write(f,">1_",name," \n")
		else
		    write(f,"\n>$(i)_",name," \n")
		end
		for j in 1:length(number_matrix[1,:])
		    if number_matrix[i,j]==1 
		        write(f,"A")
		    elseif number_matrix[i,j]==2
		        write(f,"C")
		    elseif number_matrix[i,j]==3
		        write(f,"G")
		    elseif number_matrix[i,j]==4
		        write(f,"U")
		    elseif number_matrix[i,j]==5
		        write(f,"-")
		    end
		end
	end	    
	end
end



function print_fasta_rna(number_matrix,name)
    for i in 1:length(number_matrix[:,1])
        if i==1
            println(">1_",name)
        else
            println("\n>$(i)_",name)
        end
        for j in 1:length(number_matrix[1,:])
            if number_matrix[i,j]==1 
                print("A")
            elseif number_matrix[i,j]==2
                print("C")
            elseif number_matrix[i,j]==3
                print("G")
            elseif number_matrix[i,j]==4
                print("U")
            elseif number_matrix[i,j]==5
                print("-")
            end
        end
    end
end



function site_entropy_vector(matrix,q,pseudo_count,threshold)
	frequency=freq_reweighted(matrix,q,pseudo_count,threshold)
	entropy_vector=zeros(Float32,length(matrix[1,:]))
	for i in 1:length(matrix[1,:])
	    temp=0
	    for j in 1:q
		temp=temp+frequency[q*(i-1)+j]*log(frequency[q*(i-1)+j])
	    end
	    entropy_vector[i]=temp
	end
	entropy_vector[isnan.(entropy_vector)].= 0
	return -entropy_vector
end



function profile_model_entropy(matrix,q,pseudo_count,threshold)
	return sum(site_entropy_vector(matrix,q,pseudo_count,threshold))
end



function profile_model_generation(threshold,q,pseudo_count,number,number_matrix)
    freq=freq_reweighted(number_matrix,q,pseudo_count,threshold)
    sequences=zeros(Int8,number,length(number_matrix[1,:]))
    for j in 1:length(number_matrix[1,:])
        for i in 1:number
            sequences[i,j]=sample(1:q,Weights(freq[q*(j-1)+1:q*(j-1)+q]))
        end
    end
    return sequences
end



function gibbs_sampling(q,h_local,Jij,sequences,site_degree,contact_list,sweeps)
    w=zeros(Float32,q) 
    w2=zeros(Float32,q) 
    rows=length(sequences[:,1])
    columns=length(sequences[1,:])
        for i3 in 1:sweeps*columns
            if i3%columns==0
                i2=columns
            else
                i2=i3%columns
            end
            w2=h_local[q*(i2-1)+1:q*(i2-1)+q]  
            for i1 in 1:rows 
            copy!(w,w2)
                      
                for i5 in contact_list[1:site_degree[i2],i2]
                for i4 in 1:q
                    if i2<i5
                            w[i4]+=Jij[i2,i5,q*(i4-1)+sequences[i1,i5]]
                    elseif i2>i5
                            w[i4]+=Jij[i5,i2,q*(sequences[i1,i5]-1)+i4]
                    end
                end
             end
            sequences[i1,i2]=sample(1:q,Weights(exp.(w)))
        end
    end
    return sequences
end



function E_A_A(q,pseudo_count,number,number_matrix,filename)
    edge_list=zeros(Int64,0,2)
    n_edges=0
    n_fully_connected_edges=Int64(length(number_matrix[1,:])*(length(number_matrix[1,:])-1)*0.5)
    contact_list=zeros(Int64,length(number_matrix[1,:]),length(number_matrix[1,:]))
    site_degree=zeros(Int64,length(number_matrix[1,:]))
    likelihood_gain_vector=Float32[]
    Jij_couplings=zeros(Float32,length(number_matrix[1,:]),length(number_matrix[1,:]),q*q) 
    h_local=log.(freq_reweighted(number_matrix,q,pseudo_count,0.8))
    sequences=zeros(Int8,number,length(number_matrix[1,:]))
    fij_target=fij_reweighted(number_matrix,q,pseudo_count,0.8)
    cij_target=correlation_reweighted(number_matrix,q,0,0.8) 
    score_vector=Float32[]
    contact_matrix=zeros(Int8,length(number_matrix[1,:]),length(number_matrix[1,:]))
    log_z=Float32(0)
    println("Fully connected model has ",n_fully_connected_edges," edges and a score around ~ 0.95")
    open(filename, "w") do f  
    write(f,"Fully connected model has ","$(n_fully_connected_edges)"," edges and a score around ~ 0.95")  
    for i in 1:10
    	 flush(stdout)   
    	 flush(f)                           
         sequences=gibbs_sampling(q,h_local,Jij_couplings,sequences,site_degree,contact_list,5)
         pij_training=fij_two_point(sequences[1:number-2000,:],q,pseudo_count)
         pij_lgz=fij_two_point(sequences[number-1999:end,:],q,0)  
         if (i-1)%15==0 && i!=1
            cij_model=correlation_two_point(sequences,q,0)  
            score=cor(cij_target,cij_model)    
            score_vector=push!(score_vector,score)            
            print("   Score = ",round(score;digits=3)) 
            write(f,"   Score = ",round(score;digits=3))
            energy1=-sum(fij_two_point(sequences,q,0).*Jij_couplings)-sum(freq_single_point(sequences,q,0).*h_local)               
            print("    <E> = ",round(energy1;digits=2), "    log(Z) = ",round(log_z;digits=2)) 
            print("   S = ",round(log_z+energy1;digits=2))
	    write(f,"    <E> = ","$(round(energy1;digits=2))", "    log(Z) = ","$(round(log_z;digits=2))")  
            write(f,"   S = ","$(round(log_z+energy1;digits=2))")                 
            if score >=Float32(0.95)
            	println("\n \nThe selceted model has ",n_edges," edges and a score = $(round(score;digits=2))")
            	write(f,"\n \nThe selceted model has ","$(n_edges)"," edges and a score = $(round(score;digits=2)) \n")
            	return score_vector ,likelihood_gain_vector, sequences,Jij_couplings,h_local,contact_list,site_degree,edge_list
            end        
         end     
         added_edge,likelihood_gain=max_kl_divergence(fij_target,pij_training)       
         likelihood_gain_vector=push!(likelihood_gain_vector,likelihood_gain)
         print("\n[", added_edge[1] , "  ", added_edge[2], "]  iter: $i" ) 
         write(f,"\n[", "$(added_edge[1])" , "  ", "$(added_edge[2])", "]  iter: $i" ) 
         if contact_matrix[added_edge[1],added_edge[2]]==0
            n_edges+=1
            site_degree[added_edge[1]]+=1
	    site_degree[added_edge[2]]+=1   
	    contact_list[site_degree[added_edge[1]],added_edge[1]]=added_edge[2]
	    contact_list[site_degree[added_edge[2]],added_edge[2]]=added_edge[1]
	    contact_matrix[added_edge[1],added_edge[2]]=1
	    edge_list=vcat(edge_list,[added_edge[1],added_edge[2]]')
        end
        print("   edges: ",n_edges, "   ",round(((n_edges)/n_fully_connected_edges)*100,digits=2) ,"%"    )
        write(f,"   edges: ","$(n_edges)", "   ","$(round(((n_edges)/n_fully_connected_edges)*100,digits=2))" ,"%"    ) 
        log_z+=log(sum((fij_target[added_edge[1],added_edge[2],:]./(pij_training[added_edge[1],added_edge[2],:])).*      (pij_lgz[added_edge[1],added_edge[2],:])))    
        Jij_update=log.(fij_target[added_edge[1],added_edge[2],:]./(pij_training[added_edge[1],added_edge[2],:]))
        Jij_couplings[added_edge[1],added_edge[2],:]+=Jij_update     
    end  
    end
    return score_vector,likelihood_gain_vector, sequences,Jij_couplings,h_local,contact_list,site_degree,edge_list
end



function reweighted_sample(number_matrix,number,threshold)
    weights=weight_of_sequences(number_matrix,threshold)
    reweighted_matrix=zeros(Int8,number,length(number_matrix[1,:]))
    idxs=sample(1:length(number_matrix[:,1]) , Weights(weights),number)
    for i in 1:number
        reweighted_matrix[i,:]=number_matrix[idxs[i],:]
    end
return reweighted_matrix
end



function rna_cm_model_generation(threshold,pseudo_count,number,number_matrix,ss_contact_matrix)	
    sec_proxy_list=findall(!iszero, ss_contact_matrix)
    proxy_idx_1 = getindex.(sec_proxy_list, 1)
    proxy_idx_2 = getindex.(sec_proxy_list, 2)
    sec_contacts=hcat(proxy_idx_1,proxy_idx_2)
    fij=fij_reweighted(number_matrix,5,pseudo_count,threshold)
    sequences=profile_model_generation(threshold,5,number,number_matrix)
    for i in 1:number
        for j in 1:length(sec_contacts[:,1])
            w=fij[sec_contacts[j,1],sec_contacts[j,2],5*((sequences[i,sec_contacts[j,1]])-1)+1:5*((sequences[i,sec_contacts[j,1]])-1)+5]
            sequences[i,sec_contacts[j,2]]=sample(1:5,Weights(w))
        end
    end
    return sequences
end



function ss_matrix_to_dot_bracket(ss_contact_matrix)
	sec_proxy_list=findall(!iszero, ss_contact_matrix)
        proxy_idx_1 = getindex.(sec_proxy_list, 1)
        proxy_idx_2 = getindex.(sec_proxy_list, 2)
        ss_contact_list=hcat(proxy_idx_1,proxy_idx_2)      
 dot_bracket=[]
 for i in 1:length(ss_contact_matrix[1,:])
 	dot_bracket=push!(dot_bracket,".")
 end
 for i in 1:length(ss_contact_list[:,1])
 	dot_bracket[ss_contact_list[i,1]]="("
 	dot_bracket[ss_contact_list[i,2]]=")"
 end
 return join(dot_bracket)
 end

 	
 	
function dot_bracket_to_ss_matrix(dot_bracket_ss)
	 dot_bracket_ss=collect((join(dot_bracket_ss)))
	 for i in 1:length(dot_bracket_ss)
		    if dot_bracket_ss[i]=='{' || dot_bracket_ss[i]=='<' || dot_bracket_ss[i]=='[' || dot_bracket_ss[i]=='(' 
			dot_bracket_ss[i]='('
		    elseif dot_bracket_ss[i]=='}' || dot_bracket_ss[i]=='>' || dot_bracket_ss[i]==']' || dot_bracket_ss[i]==')'
			dot_bracket_ss[i]=')'
		    else     
			dot_bracket_ss[i]='.'
		    end
		end
		len=length(dot_bracket_ss)
		contact_matrix=zeros(Int8,len,len)  
		number_of_contacts=length(dot_bracket_ss[dot_bracket_ss.=='('])
		contact_list=zeros(number_of_contacts,2)                                                 
		ss_proxy=zeros(Int8,len)  
		for i in 1:len
		    if dot_bracket_ss[i]=='('
			ss_proxy[i]=1
		    elseif  dot_bracket_ss[i]==')'
			ss_proxy[i]=2
		    else
			ss_proxy[i]==0
		    end
		end
		contact_matrix=zeros(Int8,len,len)  
		contact_list=zeros(Int64,number_of_contacts,2)  
		idx=1
		for i in 1:len
		    k=0
		    if ss_proxy[i]==2
			ss_proxy[i]=0
			for j in 1:i
			    if k==0
			    if j<i
				if ss_proxy[i-j]==1
				    contact_matrix[i-j,i]=1
				    contact_list[idx,:]=[i-j,i]
				    idx+=1
				    ss_proxy[i-j]=0
				    k=1 
				end
			    end
			    end
			    
			end
		    end
		end
    return contact_list, contact_matrix 
end	
 


end # module
