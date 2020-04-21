~/ClaraGenomicsAnalysis/build/cudamapper/cudamapper -w 5 -F 0.001 -t $1 $2 $2 > cuda.1.paf  # mandatory output name
~/raven/build/bin/raven -t $1 -p 0 $2
~/ClaraGenomicsAnalysis/build/cudamapper/cudamapper -w 5 -F 0.001 -t $1 uncontained.fasta uncontained.fasta > cuda.2.paf  # mandatory output name
~/ClaraGenomicsAnalysis/build/cudamapper/cudamapper -w 5 -F 0.00001 -t $1 contained.fasta uncontained.fasta > cuda.3.paf  # mandatory output name
~/raven/build/bin/raven -t $1 -p 0 $2 --resume > out
