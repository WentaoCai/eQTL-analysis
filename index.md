
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

 Genes were selected based on expression thresholds of ≥0.1 TPM in ≥20% of samples and expression values for each gene should be inverse normal transformed across samples.
 
    library("peer")
    library("preprocessCore")
    TPM_tmp0<-read.table("120.gene.expression.tpm.txt",header=T)
    row.names(TPM_tmp0)<-TPM_tmp0$ID
    TPM_tmp0<-TPM_tmp0[,-1]
    #expr=subset(TPM_tmp0,rownames(TPM_tmp0)%in%ccm_pca$sample.id)
    #expr_matrix00=as.data.frame(t(TPM_tmp0))
    #expr_matrix00<-expr_matrix00[,order(factor(colnames(expr_matrix00),levels=ccm_pca$sample.id))]
    count_0.1<-rowSums(TPM_tmp0>0.1)
    nsamples<-227
    expr_matrix00<-TPM_tmp0[count_0.1>=(0.2*nsamples),]
    expr_matrix00_qn<-normalize.quantiles(as.matrix(expr_matrix00)) #normalize.quantiles needs column a chip (sample);row is the gene. normalize between samples.


    INT<-function(x){
      qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
    }
    expr_matrix00_qn_ind<-t(apply(expr_matrix00_qn,MARGIN=1,FUN=INT))#apply to each row, each row represents one gene, observed values for all the samples. scale across samples.
    colnames(expr_matrix00_qn_ind)<-colnames(expr_matrix00)
    rownames(expr_matrix00_qn_ind)<-rownames(expr_matrix00)
    expr_peer<-expr_matrix00_qn_ind
    expr_matrix00_qn_ind<-as.data.frame(expr_matrix00_qn_ind)
    expr_matrix00_qn_ind$id<-row.names(expr_matrix00)
    model=PEER()
    PEER_setPhenoMean(model,as.matrix(t(expr_peer)))
    dim(PEER_getPhenoMean(model))
    PEER_setNk(model,30)
    PEER_setNmax_iterations(model,1000)
    PEER_getNk(model)
    PEER_update(model)
    factors=PEER_getX(model)
    rownames(factors)<-colnames(expr_matrix00)
    colnames(factors)<-c(1:30)
    Alpha = PEER_getAlpha(model)
    pdf(paste0(30,".r_demo_covs.pdf"),width=8,height=8)
    plot(1.0 / Alpha,xlab="Factors", ylab="Factor relevance", main="")
    dev.off()
    factors
    write.table(factors ,"120.peer.covariance.txt",sep="\t",row.names=T,quote =FALSE)
    write.table(expr_matrix00_qn_ind ,"120.expression.qn_ind.txt",sep="\t",row.names=T,quote =FALSE)

 
### Cis-eQTL mapping

#### Permutation    
    #!/usr/bin/bash
    for j in $(seq 1 50);
    do
    fastQTL.static --vcf 227HD.imputated_QC.vcf.gz --bed 227.expression.qn_ind.29.sort1.txt.gz --cov 227covariate.txt.gz --permute 1000 10000 --normal --out ./Permutation/227imputated.permutation_withage.chunk${j}.txt.gz --chunk $j 50&
    done
    wait


#### Nominal    
    #!/usr/bin/bash
    for j in $(seq 1 50);
    do
    fastQTL.static --vcf 227HD.imputated_QC.vcf.gz --bed 227.expression.qn_ind.29.sort1.txt.gz --cov 227covariate.txt.gz  --normal --out ./Nominal/227imputated.nominals.withage.chunk${j}.txt.gz --chunk $j 50&
    done
    wait



### trans-eQTL mapping








### Support or Contact

Any question please [contact me](https://github.com/WentaoCai)
