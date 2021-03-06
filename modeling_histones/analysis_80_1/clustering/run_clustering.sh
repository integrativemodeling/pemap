module load cuda/9.0.176

export top_dir=/home/ignacia/Research/pE-MAP/mod_hist_80/
export analys_dir=$top_dir/analys
export mod_dir=$analys_dir/GSMs_cl1/
export name=HIS_hist_80

cp $top_dir/density.txt $mod_dir
cp $mod_dir/sample_A.txt $mod_dir/Scores_A.txt
cp $mod_dir/sample_B.txt $mod_dir/Scores_B.txt

ls -lta $analys_dir/GSMs_cl1/sample_A | awk '{print $9}' | grep 'rmf3' > $analys_dir/selected_models_A_cluster1_random.dat
ls -lta $analys_dir/GSMs_cl1/sample_B | awk '{print $9}' | grep 'rmf3' > $analys_dir/selected_models_B_cluster1_random.dat

nohup ~/SOFTW/imp_november/setup_environment.sh python ~/SOFTW/imp-sampcon/pyext/src/Master_Sampling_Exhaustiveness_Analysis.py --sysname $name --path $mod_dir --mode cuda --cores 5 --align --density density.txt --gridsize 3.0 --gnuplot --scoreA Scores_A.txt --scoreB Scores_B.txt  --rmfs_A $analys_dir/selected_models_A_cluster1_random.dat --rmfs_B $analys_dir/selected_models_B_cluster1_random.dat --subsample 20000 > clustering.log &
