To submit test.R1.fastq file  

for SAMPLE in `ls rawdata/test.R1.fastq | cut -f1 -d. | cut -f2 -d/`;  

do  

    sbatch submitcluster.sh $SAMPLE;   

    echo $SAMPLE;  
 
done
