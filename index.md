
# eQTL Analysis

[Download Script](https://github.com/WentaoCai/eQTL-analysis/) 

### Mapping Reads

You can choose either HISAT2 or [STAR](https://github.com/WentaoCai/RNA-seq/wiki)

    hisat2 -p 4 --dta -x Genome -1 "$i"_1.clean.fq.gz -2 "$i"_2.clean.fq.gz -S "$i".sam
    samtools sort -@ 10 -o "$i".bam "$i".sam

### Get TPM matrix file

    mkdir Expression
    stringtie -p 10 -e -B -G Bos_taurus.ARS-UCD1.2.96.chr.gtf -o "$i".gtf -A ./Expression/"$i".tsv "$i".bam
    cd Expression
    for i in `cat ../Sample_name_List.txt`
    do
    awk '{print $1"\t"$NF}' "$i".tsv|sed '1d'|sort -k1,1>"$i".tpm
    echo -e "ID\t"$i""|cat - "$i".tpm >"$i".tpm.txt
    rm "$i".tpm
    done
    paste *.tpm.txt|awk '{printf("%s\t",$1);for(i=2;i<=NF;i+=2){printf("%s\t",$i)};print ""}'>All.gene.expression.tpm.txt
    cd ..


### Peer comfound factor

 Genes were selected based on expression thresholds of ≥0.1 TPM in ≥20% of samples
 
 Expression values for each gene should be inverse normal transformed across samples
 
 



### Cis-eQTL mapping




### trans-eQTL mapping








### Support or Contact

Any question please [contact me](https://github.com/WentaoCai)
