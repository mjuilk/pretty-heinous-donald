mkdir /nfs/data/doad5844/rer/bam/rr4/test1
for bam in $(ls *.bam | sort -V | head -n 1000)
do
	cp ${bam} /nfs/data/doad5844/rer/bam/rr4/test1/${bam}
	echo "${bam} HAS BEEN COPIED TO THE DIRECTORY"
done

mkdir /nfs/data/doad5844/rer/bam/rr4/test2
for bam in $(ls *.bam | sort -V | tail -n +1001)
do
	cp ${bam} /nfs/data/doad5844/rer/bam/rr4/test2/${bam}
	echo "${bam} HAS BEEN COPIED TO THE DIRECTORY"
done 

