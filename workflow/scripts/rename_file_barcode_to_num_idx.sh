# search the directory of interest and rename any filename with barcode to numerical index

dir_to_search=$1
barcode=$2
num_idx=$3

echo "processing ---- $num_idx"

target_files=$(find $dir_to_search -name "*${barcode}*" -type f)
for i in $target_files
do
    #out_file_name=$(basename $i)
    mv $i ${i/$barcode/$num_idx}
done