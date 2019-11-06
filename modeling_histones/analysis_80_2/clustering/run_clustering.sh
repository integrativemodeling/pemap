module load cuda/9.0.176

export cl=0
export top_dir=/home/ignacia/Research/pE-MAP/mod_hist_80_2/
export analys_dir=$top_dir/analys
export mod_dir=$analys_dir/GSMs_cl$cl/
export name=HIS_hist

cp /home/ignacia/Research/pE-MAP/mod_hist_80/density.txt $mod_dir
cp $mod_dir/sample_A.txt $mod_dir/Scores_A.txt
cp $mod_dir/sample_B.txt $mod_dir/Scores_B.txt

ls -lta $analys_dir/GSMs_cl${cl}/sample_A | awk '{print $9}' | grep 'rmf3' > $analys_dir/selected_models_A_cluster${cl}_random.dat
ls -lta $analys_dir/GSMs_cl${cl}/sample_B | awk '{print $9}' | grep 'rmf3' > $analys_dir/selected_models_B_cluster${cl}_random.dat

nohup ~/SOFTW/imp_november/setup_environment.sh python ~/SOFTW/imp-sampcon/pyext/src/Master_Sampling_Exhaustiveness_Analysis.py --sysname $name --path $mod_dir --mode cuda --align --density density.txt --gridsize 3.0 --gnuplot --scoreA Scores_A.txt --scoreB Scores_B.txt  --rmfs_A $analys_dir/selected_models_A_cluster${cl}_random.dat --rmfs_B $analys_dir/selected_models_B_cluster${cl}_random.dat --subsample 20000 > clustering.log &
