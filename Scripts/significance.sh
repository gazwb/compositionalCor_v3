INPUT_PATH="/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/PSO_LES/raw_dat_PSO_LES"
OUTPUT_PATH="/home/gaz/MAARS_p2/Scripts/compositionalCor_v2/PSO_LES/SparCC"
SPARCC_PATH="/home/gaz/Documents/SparCC"
ITER=5

cd $SPARCC_PATH
count=0
for f in $OUTPUT_PATH/Resamplings/boot_*; do
  echo SparCC.py $f -i 5 --cor_file=$OUTPUT_PATH/Bootstraps/sim_cor_$f &
done
