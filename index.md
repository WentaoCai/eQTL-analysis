
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
    fastQTL.static --vcf 120.imputated_QC.vcf.gz --bed 120.expression.qn_ind.29.sort1.txt.gz --cov 120covariate.txt.gz --permute 1000 10000 --normal --out ./Permutation/120.permutation.chunk${j}.txt.gz --chunk $j 50&
    done
    wait
    
    cat 120.permutation.chunk*.txt.gz>120.permutation.txt.gz

#### Nominal    
    #!/usr/bin/bash
    for j in $(seq 1 50);
    do
    fastQTL.static --vcf 120.imputated_QC.vcf.gz --bed 120.expression.qn_ind.29.sort1.txt.gz --cov 120covariate.txt.gz  --normal --out ./Nominal/120.nominals.chunk${j}.txt.gz --chunk $j 50&
    done
    wait
    
    cat 120.nominals.chunk*.txt.gz>120.nominals.txt.gz

#### Combine Permutation and Nominal to get the final cis-eQTL results

    data = read.table("120.permutation.txt.gz", header=FALSE, stringsAsFactors=FALSE)
    colnames(data) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval","effects", "ppval", "bpval")
    data$bh = p.adjust(data$bpval, method="fdr")
    write.table(data[which(data$bh <= 0.05), ], "120.permutations.fdr.txt", quote=F, row.names=F, col.names=T)
    data = data[which(!is.na(data[, 10])),]
    set0 = data[which(data$bh <= 0.05),] 
    set1 = data[which(data$bh > 0.05),]
    pthreshold = (sort(set1$bpval)[1] - sort(-1.0 * set0$bpval)[1]) / 2
    data$nthresholds = qbeta(pthreshold, data$shape1, data$shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
    nom=read.table("120.nominals.txt.gz",header=F,stringsAsFactors=F)
    colnames(nom)<-c("pid","sid","dst","pval")
    nom$threshold<-data$nthresholds[match(nom$pid,data$pid)]
    write.table(nom[nom$pval<nom$threshold,],"120.nominals.sigs.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

### trans-eQTL mapping
    library("MatrixEQTL")
    useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS
    output_file_name = tempfile();
    pvOutputThreshold = 1e-2;
    errorCovariance = numeric();

    output_file_name_cis = tempfile();
    output_file_name_tra = tempfile();

    pvOutputThreshold_cis = 0.05;
    pvOutputThreshold_tra = 0.0001;
    cisDist = 1e6;

    snpspos<-read.table(snp_pos, header = F, stringsAsFactors = FALSE);
    genepos = read.table(gene_pos, header = F, stringsAsFactors = FALSE);

    snps = SlicedData$new();
    snps
    snps$fileDelimiter = "\t";      # the TAB character
    snps$fileOmitCharacters = "NA"; # denote missing values;
    snps$fileSkipRows = 1;          # one row of column labels
    snps$fileSkipColumns = 1;       # one column of row labels
    snps$fileSliceSize = 200000;      # read file in pieces of 2,000 rows
    snps$LoadFile( Genotype);

    gene = SlicedData$new();
    gene$fileDelimiter = "\t";      # the TAB character
    gene$fileOmitCharacters = "NA"; # denote missing values;
    gene$fileSkipRows = 1;          # one row of column labels
    gene$fileSkipColumns = 1;       # one column of row labels
    gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    gene$LoadFile(Phnotype);

    ## Load covariates

    cvrt = SlicedData$new();
    cvrt$fileDelimiter = "\t";      # the TAB character
    cvrt$fileOmitCharacters = "NA"; # denote missing values;
    cvrt$fileSkipRows = 1;          # one row of column labels
    cvrt$fileSkipColumns = 1;       # one column of row labels
    if(length(Covariate)>0) {
      cvrt$LoadFile(Covariate);
    }

    me = Matrix_eQTL_main(
      snps = snps,
      gene = gene,
      cvrt = cvrt,
      output_file_name     = output_file_name_tra,
      pvOutputThreshold     = pvOutputThreshold_tra,
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = TRUE,
      output_file_name.cis = output_file_name_cis,
      pvOutputThreshold.cis = pvOutputThreshold_cis,
      snpspos = snpspos,
      genepos = genepos,
      cisDist = cisDist,
      pvalue.hist = "qqplot",
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE);

    unlink(output_file_name_tra);
    sum(me$trans$eqtls$FDR< 0.05, na.rm=TRUE)
    TransSig <- subset(me$trans$eqtls, FDR < 0.05)
    TransOrdered <- TransSig[order(TransSig$FDR),]
    write.csv(TransOrdered,paste(Phnotype,".trans.csv",sep = ""))


### Support or Contact

Any question please [contact me](https://github.com/WentaoCai)
