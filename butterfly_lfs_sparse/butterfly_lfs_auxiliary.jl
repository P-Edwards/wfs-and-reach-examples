using HomologyInferenceWithWeakFeatureSize, NearestNeighbors, LinearAlgebra, Distances, Ripserer, PersistenceDiagrams, ProgressBars

function compute_normal(point,tree,data)
	nearest,dist = knn(tree,point,2)
	nearest_point = data[:,nearest[1]]
	estimated_tangent = nearest_point - point
	estimated_normal = (-estimated_tangent[2],estimated_tangent[1])
end

function angle_between_vectors(vec1,vec2)
	vec1_norm = norm(vec1)
	vec2_norm = norm(vec2)
	dot_product = dot(vec1,vec2)
	the_cos = dot_product/(vec1_norm*vec2_norm)
	if the_cos > 1
		the_cos = 1
	end
	if the_cos < -1
		the_cos = -1
	end
	acos(the_cos)
end

function reconstructed_norm(first_point,axis_point)
	second_point = 2*axis_point - first_point
	return norm(first_point - second_point)
end

function construct_lean_axis(points;beta=pi/5,search_tree=false)
	if search_tree == false
		search_tree = KDTree(points)
	end
	lean_medial_axis = []
	number_of_points = size(points,2)
	normals = [compute_normal(points[:,j],search_tree,points) for j in 1:number_of_points]
	for i in 1:(number_of_points-1)
		if mod(i,500) == 0
			println("Point "*string(i)*" out of "*string(number_of_points-1))
		end
		first_point = points[:,i]
		this_point_axis = []		
		for j in (i+1):number_of_points
			second_point = points[:,j]
			direction = second_point - first_point
			first_angle = angle_between_vectors(direction,normals[i])	
			second_angle = angle_between_vectors(direction,normals[j])
			if max(first_angle,second_angle) <= pi/2 - beta
				v = (first_point+second_point)/2
				push!(this_point_axis,v)
			end			
		end
		if length(this_point_axis) == 0
			continue
		end
		this_point_axis = reduce(vcat,transpose.(this_point_axis))'
		indices,dists = nn(search_tree,this_point_axis)
		c_beta = (1/3)*tan(beta/2)
		this_point_axis_final = [this_point_axis[:,j] for j in 1:size(this_point_axis,2) if dists[j] > c_beta*reconstructed_norm(first_point,this_point_axis[:,j])]		
		if length(this_point_axis_final) > 0
			distances = [reconstructed_norm(first_point,axis_point) for axis_point in this_point_axis_final]				
			min_value,min_index = findmin(distances)
			lean_medial_axis = push!(lean_medial_axis,this_point_axis_final[min_index])
		end
	end
	return lean_medial_axis
end

function distance_to_lean_medial_axis(points;beta=pi/5,mu=nothing,lean_axis=nothing)
	if mu == nothing
		mu = (1/26)*(cos(2*beta)/(1+cos(2*beta)))
	end	
	if lean_axis == nothing
		search_tree = KDTree(points;leafsize=50)
		lean_medial_axis = construct_lean_axis(points;beta=beta,search_tree=search_tree)
		lean_medial_axis = reduce(vcat,transpose.(lean_medial_axis))'
		lean_search_tree = KDTree(lean_medial_axis)
	else
		lean_search_tree = lean_axis
	end
	output = function(point)
		_,dists = nn(lean_search_tree,point)
		return dists*mu
	end
	return output, lean_search_tree
end

function modified_matrix(input_sample,lfs)
    plain_distance_matrix = pairwise(euclidean,[tup[k] for tup in input_sample, k in 1:2]',dims=2)
    big_lfs_matrix = [1 for _ in 1:length(lfs)] .* lfs'
    big_lfs_matrix = big_lfs_matrix + big_lfs_matrix'
    return plain_distance_matrix./big_lfs_matrix
end

function modified_lean_matrix(points, lean_axis_tree)
	indexes, distances = nn(lean_axis_tree,reduce(vcat,transpose.(points))')
	return modified_matrix(points,distances)
end


function construct_grid_points(initial,final,number)
    step_size = (final-initial)/(number-1)
    return [initial + i*step_size for i in 0:(number-1)]
end

function sparsification_tester(start_sample,starting_subsample_density,ending_subsample_density,sparsification_function_family;number_of_discretization_points=10,homology_up_to_degree=1,sample_processing_function=identity,norm_func=norm,has_index=false)
    range_of_proportions = construct_grid_points(starting_subsample_density,ending_subsample_density,number_of_discretization_points)    
    sample_lengths = []
    persistence_diagrams_over_range = Array{Any}(undef,length(range_of_proportions))
    i = 1    
    for proportion in range_of_proportions
        subsample = subsample_with_function(start_sample,sparsification_function_family(proportion),norm_func;has_index=has_index)
        push!(sample_lengths,length(subsample))
        pd = ripserer(sample_processing_function(subsample),dim_max=homology_up_to_degree)    
        persistence_diagrams_over_range[i] = pd 
        i += 1        
    end    
    return Dict("ticks"=>range_of_proportions,"lengths"=> sample_lengths,"persistence_diagrams"=>persistence_diagrams_over_range)
end

function pd_score(pd,degree=1) 
    current_diagram = pd[degree+1]
    most_persistent = current_diagram[end][2] - current_diagram[end][1]
    if length(current_diagram) > 1
        second_most_persistent = current_diagram[end-1][2]-current_diagram[end-1][1]
    else 
        second_most_persistent = 0 
    end
    return most_persistent - second_most_persistent
end

function scores_from_results(results)
	homology_inference_scores = [pd_score(diag) for diag in results["persistence_diagrams"]]
	max_score = maximum(homology_inference_scores)
	homology_inference_scores = [score/max_score for score in homology_inference_scores]
end

function wasserstein_score(pd,base,degree=1)
	current_diagram = pd[degree+1]
	base_diagram = base[degree+1]
	return Wasserstein()(current_diagram,base_diagram)
end

function wasserstein_scores_from_results(results,degree=1)
	base = results["persistence_diagrams"][1]
	return [wasserstein_score(diag,base,degree) for diag in ProgressBar(results["persistence_diagrams"])]
end

function dist_to_empty(pd,degree=1)
	current_diagram = pd[degree+1]
	base_persistence_scores = [((point[2] - point[1])/2)^2 for point in current_diagram]
	return sqrt(sum(base_persistence_scores))
end

function dist_to_empty_from_results(results,degree=1)
	return [dist_to_empty(diag,degree) for diag in ProgressBar(results["persistence_diagrams"])]
end


function number_of_points_from_results(results)
	return results["lengths"]
end

function every_other(ticks,modulus=2)
	if mod(length(ticks),modulus) == 0
	    adjustment = -1
	else
	    adjustment = 0
	end
	return [round(ticks[Int(modulus*i)+1],digits=2) for i in 0:(length(ticks)/modulus + adjustment)]	
end

