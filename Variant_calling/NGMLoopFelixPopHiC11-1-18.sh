while read i
do

echo "Processing sample $i ..." 
ngm -r /ohta/joanna.rifkin/Genomes/Hi-C/hastate_28Sep2018_nbKXS.fasta -b -q /ohta/felix.beaudry/rawSequence/genomic/GBSpops/40014_$i\_il.fastq.gz -t 16 -p -o /ohta/joanna.rifkin/Alignments/NGM/Hastatulus/HiC/PopRadSeq/RawBam/$i.bam

done < /ohta/joanna.rifkin/Alignments/NGM/Hastatulus/ChicagoAssembly/PopRadSeq/FelixRADSeqPopList.txt
