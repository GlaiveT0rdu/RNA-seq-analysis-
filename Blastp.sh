
#Tell BLAST that the Arabidopsis sequences are a protein database.
makeblastdb -in /data/home/mathieu.bernardini/Rna-seq/Database/Araport11_pep_20220914.fasta -dbtype prot

#Call BLAST to do the search for arabidopsis 
blastp -query /data/home/mathieu.bernardini/Rna-seq/database/RNA_seq_G3.fasta -db /data/home/mathieu.bernardini/Rna-seq/database/Araport11_pep_20220914.fasta -out /data/home/mathieu.bernardini/Rna-seq/blast/arabido/blast_turnip_arabidoG3.txt

grep -e 'Query=' -A 6 /data/home/mathieu.bernardini/Rna-seq/blast/arabido/blast_turnip_arabidoG3.txt > /data/home/mathieu.bernardini/Rna-seq/blast/arabido/blast_reduced_turnip_arabidoG3.txt
grep -e 'Query=' /data/home/mathieu.bernardini/Rna-seq/blast/arabido/blast_reduced_turnip_arabidoG3.txt > /data/home/mathieu.bernardini/Rna-seq/blast/arabido/query_G3column.txt
grep -e ' AT\|***** No hits found *****'  /data/home/mathieu.bernardini/Rna-seq/blast/arabido/blast_reduced_turnip_arabidoG3.txt > /data/home/mathieu.bernardini/Rna-seq/blast/arabido/gene_column_G3.txt
