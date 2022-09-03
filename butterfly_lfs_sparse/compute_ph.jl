include("butterfly_lfs_auxiliary.jl")
using JLD
original_sample_file_path = "paper_data/butterfly_points_lfs_sparse.csv"
lfs_file_path = "paper_data/local_feature_sizes.csv"
lean_axis_path = "paper_data/lean_medial_axis.csv"


wfs=0.251
constant_wfs_family = function(proportion)
    return HomologyInferenceWithWeakFeatureSize.constant_func_curry(wfs*proportion)
end
standard_wfs_results = sparsification_tester(wfs_homology_inference_results["sample"],0.15,1.0,constant_wfs_family;number_of_discretization_points=20)

@save "standard_wfs_results.jld" standard_wfs_results

# This sample was at uniform 2e-4 density 
using CSV, DataFrames
original_sample = CSV.read(original_sample_file_path,DataFrame,header=false)
original_sample = values.(eachrow(original_sample))
original_sample = [[original_sample[i][1],original_sample[i][2],i] for i in 1:length(original_sample)]
local_feature_sizes = CSV.read(lfs_file_path,DataFrame,header=false)
local_feature_sizes = values.(eachcol(local_feature_sizes))
local_feature_sizes = [entry[1] for entry in local_feature_sizes]
mu = 0.0045
function lfs_return(point) 
    return local_feature_sizes[Int(point[3])]
end
function mu_return(point)
    return mu*lfs_return(point)
end
function modified_norm(point) 
    return norm((point[1],point[2]))
end
sparsified_sample_with_indices = subsample_with_function(original_sample,mu_return,modified_norm;has_index=true)
lfs_family = function(proportion)
    return function(point)
        proportion*lfs_return(point)
    end
end
lfs_transformation = function(points) 
    these_points_lfs = [lfs_return(point) for point in points]
    return modified_matrix(points,these_points_lfs)
end
lfs_baseline_results = sparsification_tester(sparsified_sample_with_indices,mu,0.25,lfs_family;number_of_discretization_points=3,norm_func=modified_norm,has_index=true)

@save "lfs_baseline_results.jld" lfs_baseline_results


lean_medial_axis = CSV.read(lean_axis_path,DataFrame,header=false)
lean_medial_axis = values.(eachrow(lean_medial_axis))
lean_axis_tree = KDTree([tup[k] for k in 1:2,tup in lean_medial_axis])
lean_axis_family = function(proportion)
    dist_func, _ = distance_to_lean_medial_axis([];mu=proportion,lean_axis=lean_axis_tree)
    return dist_func
end
lean_transform = function(points) return modified_lean_matrix(points,lean_axis_tree) end
lean_wfs_results = sparsification_tester([[sample[1], sample[2]] for sample in sparsified_sample_with_indices],mu,1.0,lean_axis_family;number_of_discretization_points=5,sample_processing_function=lean_transform)

@save "lean_wfs_results.jld" lean_wfs_results