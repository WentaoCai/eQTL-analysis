
# eQTL Analysis


### Mapping Reads

You can choose either HISAT2 or [STAR](https://github.com/WentaoCai/RNA-seq/wiki)

    hisat2 -p 4 --dta -x Genome -1 "$i"_1.clean.fq.gz -2 "$i"_2.clean.fq.gz -S "$i".sam
    samtools sort -@ 10 -o "$i".bam "$i".sam

### Get FPKM matrix file

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



### Cis-eQTL mapping
#####[RNA-seq](https://github.com/WentaoCai/RNA-seq/wiki)



### trans-eQTL mapping




### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
