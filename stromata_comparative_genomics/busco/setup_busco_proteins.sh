# setup all the busco runs with the busco docker 

#!usr/bin/bash

for i in `cat genome_tag_list.txt`
  do
    mkdir results/$i

    # setup configuration file
    cp myconfig.ini results/$i
    sed -i "s/\/path\/to\/input_file.fna/proteins\/TEMP.fna/g" results/$i/myconfig.ini  
    sed -i "s/\/path\/to\/output_folder/results\/TEMP/g" results/$i/myconfig.ini
    sed -i "s/BUSCO_run/busco_TEMP/g" results/$i/myconfig.ini
    sed -i "s/mode = genome/mode = proteins/g" results/$i/myconfig.ini
    sed -i "s/cpu = 16/cpu = 12/g" results/$i/myconfig.ini
    sed -i "s/lineage_dataset = bacteria/lineage_dataset = eukaryota_odb10/g" results/$i/myconfig.ini
    sed -i "s/TEMP/$i/" results/$i/myconfig.ini

    # setup docker run script
    
    echo "docker run -u $(id -u) -v $(pwd):/results/$i -w /results/$i/ ezlabgva/busco:v4.1.1_cv1 busco -i proteins/$i.fna --config=results/$i/myconfig.ini --out $i --mode proteins --lineage eukaryota_odb10 -f" >> run_busco_docker.sh 
    #echo "docker run -u $(id -u) -v $(pwd):/results/$i -w /results/$i/ ezlabgva/busco:v4.1.1_cv1 busco --config=results/$i/myconfig.ini --lineage eukaryota_odb10 -f" >> run_busco_docker.sh

done

# summarise results (note requires a local busco download).
# wget "https://gitlab.com/ezlab/busco/-/commit/b8cfb1e1892a4bb08a3c0ee34ddcad1f75d87307"
# and unpack
echo "mkdir BUSCO_summaries" >> run_busco_docker.sh
echo "for i in `cat genome_tag_list.txt`; do cp $i/short_summary.specific*.txt  BUSCO_summaries; done" >> run_busco_docker.sh
echo "python busco-4.1.2/scripts/generate_plot.py -wd BUSCO_summaries/" >> run_busco_docker.sh


