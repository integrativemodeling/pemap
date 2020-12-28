module load python/pyrmsd
module load imp-fast
module load python2/scipy

export cl=5
export top_dir=/home/ignacia/Research/pE-MAP/mod_bac_rnap_new_v2/
export analys_dir=$top_dir/analys
export mod_dir=$analys_dir/GSMs_cl$cl/
export name=BacRNAP_pEMAP

cp $top_dir/density.txt $mod_dir
cp $mod_dir/sample_A.txt $mod_dir/Scores_A.txt
cp $mod_dir/sample_B.txt $mod_dir/Scores_B.txt

ls -lta $analys_dir/GSMs_cl${cl}/sample_A | awk '{print $9}' | grep 'rmf3' > $analys_dir/selected_models_A_cluster${cl}.dat
ls -lta $analys_dir/GSMs_cl${cl}/sample_B | awk '{print $9}' | grep 'rmf3' > $analys_dir/selected_models_B_cluster${cl}.dat

nohup python2.7 ~/SOFTW/imp-sampcon-new/pyext/src/Master_Sampling_Exhaustiveness_Analysis.py --sysname $name --path $mod_dir --align --mode cpu_omp --density density.txt --gridsize 1.0 --gnuplot --scoreA Scores_A.txt --scoreB Scores_B.txt  --rmfs_A $analys_dir/selected_models_A_cluster${cl}.dat --rmfs_B $analys_dir/selected_models_B_cluster${cl}.dat > clustering.log &
