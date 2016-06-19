#--------------- basic bilateral normal filter
#p1_1
#./build-nmake/bin/filter -input_file E:/workspace/geometry/mesh_filter/result/final/p1/5input.stl -output_file E:/workspace/geometry/mesh_filter/result/final/p1/5output1.vtk -filter_type bnf_fe -st 4 -ut 15 -bs 0.1  -rn 1 -need_flip_edge 1 -need_smooth 1 -smooth_threshold 0.936
#./build-nmake/bin/filter -input_file E:/workspace/geometry/mesh_filter/result/final/p1/5input.stl -output_file E:/workspace/geometry/mesh_filter/result/final/p1/5output2.vtk -filter_type bnf_fre -st 4 -ut 15 -bs 0.1 -need_flip_edge 1 -need_smooth 1 -normal_threshold 1.0 -smooth_threshold 0.936 -max_coeff 1.1 -min_coeff 1.0

#p1_2
#./build-nmake/bin/filter -input_file E:/workspace/geometry/mesh_filter/result/final/p1/2input.stl -output_file E:/workspace/geometry/mesh_filter/result/final/p1/2output1.vtk -filter_type bnf_fe -st 4 -ut 15 -bs 0.1  -rn 1 -need_flip_edge 1 -need_smooth 1 -smooth_threshold 0.936
#./build-nmake/bin/filter -input_file E:/workspace/geometry/mesh_filter/result/final/p1/2input.stl -output_file E:/workspace/geometry/mesh_filter/result/final/p1/2output2.vtk -filter_type bnf_fre -st 4 -ut 15 -bs 0.1 -need_flip_edge 1 -need_smooth 1 -normal_threshold 1.0 -smooth_threshold 0.936 -max_coeff 1.1 -min_coeff 1.0

#p_2
#./build-nmake/bin/filter -input_file E:/workspace/geometry/mesh_filter/result/final/p2/12input.stl -output_file E:/workspace/geometry/mesh_filter/result/final/p2/12output1.vtk -filter_type bnf_fe -st 6 -ut 15 -bs 0.1  -rn 3 -need_flip_edge 1 -need_smooth 1 -smooth_threshold 0.936

#rbmls
#./build-nmake/bin/filter -input_file E:/workspace/geometry/mesh_filter/result/final/rbmls/2input.stl -output_file E:/workspace/geometry/mesh_filter/result/final/rbmls/2output_bnf.stl -filter_type bnf_fe -st 4 -ut 15 -bs 0.1  -rn 2 -need_flip_edge 1 -need_smooth 1 -smooth_threshold 0.936
#./build-nmake/bin/filter -input_file E:/workspace/geometry/mesh_filter/result/final/rbmls/5input.stl -output_file E:/workspace/geometry/mesh_filter/result/final/rbmls/5output_bnf.stl -filter_type bnf_fe -st 4 -ut 15 -bs 0.1  -rn 2 -need_flip_edge 1 -need_smooth 1 -smooth_threshold 0.936


#performance
#./build-nmake/bin/filter -input_file E:/workspace/geometry/mesh_filter/result/final/performance/12input.stl -output_file E:/workspace/geometry/mesh_filter/result/final/performance/12output.vtk -filter_type bnf_fe -st 4 -ut 15 -bs 0.1 -need_flip_edge 1 -need_smooth 1 -smooth_threshold 0.936 -rn 2
#./build-nmake/bin/filter -input_file E:/workspace/geometry/mesh_filter/result/final/performance/2input.stl -output_file E:/workspace/geometry/mesh_filter/result/final/performance/2output.vtk -filter_type bnf_fre -st 4 -ut 15 -bs 0.1 -need_flip_edge 1 -need_smooth 1 -normal_threshold 0.9 -smooth_threshold 0.936 -r_max 3 -r_min 2
#./build-nmake/bin/filter -input_file E:/workspace/geometry/mesh_filter/result/final/performance/2input.stl -output_file E:/workspace/geometry/mesh_filter/result/final/performance/2output.vtk -filter_type bnf_fre -st 4 -ut 15 -bs 0.1 -need_flip_edge 1 -need_smooth 1 -normal_threshold 0.9 -smooth_threshold 0.936 -r_max 3 -r_min 2

./build-nmake/bin/filter -input_file E:/workspace/geometry/mesh_filter/result/input/prj_1mesh.stl -output_file E:/workspace/geometry/mesh_filter/result/input/prj_1mesh.vtk -filter_type bnf_fe -st 4 -ut 15 -bs 0.1 -need_flip_edge 1 -need_smooth 1 -smooth_threshold 0.936 -rn 2