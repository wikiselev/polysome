
cd /bi/group/compneur/work/kiselev/manuel/alignments

for i in `ls *.sam`
do
	# echo $i
	echo "htseq-count -s no $i /bi/group/compneur/work/kiselev/manuel/scripts/Mus_musculus.GRCm38.78.gtf > /bi/group/compneur/work/kiselev/manuel/htseq-count/$i.count" | qsub -V -cwd
done
