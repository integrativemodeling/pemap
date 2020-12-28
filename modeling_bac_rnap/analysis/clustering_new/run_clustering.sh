module load python/pyrmsd
module load imp-fast

export cl=-1
export top_dir=/home/ignacia/Research/pE-MAP/mod_BacRNAP_PEMAP/
export analys_dir=$top_dir/analys
export mod_dir=$analys_dir/GSMs_cl$cl/
export name=BacRNAP_pEMAP

cp $top_dir/density.txt $mod_dir
cp $mod_dir/sample_A.txt $mod_dir/Scores_A.txt
cp $mod_dir/sample_B.txt $mod_dir/Scores_B.txt

ls -lta $analys_dir/GSMs_cl${cl}/sample_A | awk '{print $9}' | grep 'rmf3' > $analys_dir/selected_models_A_cluster${cl}.dat
ls -lta $analys_dir/GSMs_cl${cl}/sample_B | awk '{print $9}' | grep 'rmf3' > $analys_dir/selected_models_B_cluster${cl}.dat

nohup python ~/SOFTW/imp-sampcon-new/pyext/src/Master_Sampling_Exhaustiveness_Analysis.py --sysname $name --path $mod_dir --mode cuda --align --density density.txt --gridsize 7.5 --gnuplot --scoreA Scores_A.txt --scoreB Scores_B.txt  --rmfs_A $analys_dir/selected_models_A_cluster${cl}.dat --rmfs_B $analys_dir/selected_models_B_cluster${cl}.dat --subsample 20000 > clustering.log &
