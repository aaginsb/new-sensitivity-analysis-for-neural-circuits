Code used in ``Mechanisms for dysregulation of excitatory-inhibitory balance underlying allodynia in dorsal horn neural subcircuits'' by Ginsberg et al. 2024.
	

The three main files in this repository are:

1. generate_sample_space.m
	 * This script implements an algorithm for sampling uniformly-in-space the solution to a hierarchical system of nonlinear inequalities
		 * This algorithm was originally described in  Ginsberg et al. 2024. 
		 * This specific implementation samples uniformly-in-space the system of inequalities defining the allowable parameter space for the ``simple'' circuit, described in Ginsberg et al. 2024.
		 * To run this particular implementation, the user will need to compute maximum and minimum values for each coupling strength (g_abi, g_ie, and g_abe) underlying the circuit, or to use the provided values (see ``g_simple_overall_maxs.mat'' and ``g_simple_overall_mins.mat'').
		
		
   * Files required to run this script: 
		
	   * g_simple_overall_maxs.mat
		* g_simple_overall_mins.mat
		* activation_method.m
		* activation_method_inhibitory.m
		* extrema_of_bad_function.m
		* master_parameter_file.m
		
	* Files produced by this script:
		
	  * sample_points_simple.mat (points from ``simple_points.mat'' that have been de-normalized)
	  * simple_points.mat (a uniform-in-space sample of the normalized allowable parameter space)
		
	
2. data_analysis.m
	
		* This script does limited data analysis of the sampled points from the allowable parameter space for the simple circuit from Ginsberg et al. 2024
		
			* This script produces a 3-D scatter plot of the sampled points.
			* This script also produces plots of voltage vs time and firing-rate vs time for both the inhibitory population and the excitatory population constituting the simple circuit.
		
		* Files required to run this script:
		
			* g_simple_overall_maxs.mat
			* g_simple_overall_mins.mat
			* sample_points_simple.mat
			* simple_points.mat
			
			* activation_method.m
			* make_smooth_average_fr.m
			* master_parameter_file.m
			* poisson_spike_gen.m
			* voltage_model_simple.m
		
		* Files produced by this script:
		
			* N/A
		
	
23. find_shortest_paths_to_allodynia_surface_in_normalized_space.m
	
		* This script finds the shortest path to the allodynia surface from each sampled point in the allowable parameter space for the simple subcircuit
		
			* This implements the algorithm originally described in Ginsberg et al. 2024.
			* This script also makes a parallel-plot representation of the shortest paths.
		
		* Files required to run this script:
		
			* g_simple_overall_maxs.mat
			* g_simple_overall_mins.mat
			* sample_points_simple.mat
			* simple_points.mat
			* activation_method_inhibitory.m
			* master_parameter_file.m
		
		* Files produced by this script:
		
			* displacement_to_allodynia_surface_data_normalized.mat (the displacements from the points in ``sample_points'' corresponding to the shortest paths to the allodynia surface in normalized parameter space)
			* distance_to_allodynia_surface_data_normalized.mat (the corresponding distances to the allodynia surface)
			* nearest_point_on_allodynia_surface_data_normalized.mat (the corresponding nearest points on the allodynia surface)
		
	
