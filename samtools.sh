cd /bi/group/compneur/work/kiselev/manuel/alignments/

for i in `ls *.bam`
do
	# echo $i
	echo "samtools view -h -o $i.sam $i" | qsub -V -cwd
done
