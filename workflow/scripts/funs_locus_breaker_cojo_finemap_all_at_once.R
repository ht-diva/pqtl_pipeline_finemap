# Load packages
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(R.utils))
suppressMessages(library(corrplot))
suppressMessages(library(coloc))
suppressMessages(library(bigsnpr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(stringi))
suppressMessages(library(stringr))
suppressMessages(library(patchwork))
suppressMessages(library(reshape2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(igraph))
suppressMessages(library(purrr))
suppressMessages(library(tidyr))
suppressMessages(library(plyr))
suppressMessages(library(Gviz))
suppressMessages(library(EnsDb.Hsapiens.v75))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))





### locus.breaker
locus.breaker.p <- function(
    res,
    p.sig = 5e-08,
    p.limit = 1e-05,
    hole.size = 250000,
    p.label = "P",
    chr.label = "CHR",
    pos.label = "BP"){

  res <- as.data.frame(res)
  res = res[order(as.numeric(res[, chr.label]), as.numeric(res[,pos.label])), ]
  res = res[which(res[, p.label] > p.limit), ]
  trait.res = c()

  for(j in unique(res[,chr.label])) {
    res.chr = res[which(res[, chr.label] == j), ]
    if (nrow(res.chr) > 1) {
      holes = res.chr[, pos.label][-1] - res.chr[, pos.label][-length(res.chr[,pos.label])]
      gaps = which(holes > hole.size)
      if (length(gaps) > 0) {
        for (k in 1:(length(gaps) + 1)) {
          if (k == 1) {
            res.loc = res.chr[1:(gaps[k]), ]
          }
          else if (k == (length(gaps) + 1)) {
            res.loc = res.chr[(gaps[k - 1] + 1):nrow(res.chr),
            ]
          } else {
            res.loc = res.chr[(gaps[k - 1] + 1):(gaps[k]),
            ]
          }
          if (max(res.loc[, p.label]) > p.sig) {
            start.pos = max(res.loc[, pos.label], na.rm = T)
            end.pos = min(res.loc[, pos.label], na.rm = T)
            chr = j
            best.snp = res.loc[which.max(res.loc[, p.label]),
            ]
            line.res = c(chr, start.pos, end.pos, unlist(best.snp))
            trait.res = rbind(trait.res, line.res)
          }
        }
      } else {
        res.loc = res.chr
        if (min(res.loc[, p.label]) > p.sig) {
          start.pos = max(res.loc[, pos.l.abel], na.rm = T)
          end.pos = min(res.loc[, pos.label], na.rm = T)
          chr = j
          best.snp = res.loc[which.max(res.loc[, p.label]),
          ]
          line.res = c(chr, start.pos, end.pos, unlist(best.snp))
          trait.res = rbind(trait.res, line.res)
        }
      }
    }
    else if (nrow(res.chr) == 1) {
      res.loc = res.chr
      if (max(res.loc[, p.label]) > p.sig) {
        start.pos = max(res.loc[, pos.label], na.rm = T)
        end.pos = min(res.loc[, pos.label], na.rm = T)
        chr = j
        best.snp = res.loc[which.max(res.loc[, p.label]),
        ]
        line.res = c(chr, start.pos, end.pos, unlist(best.snp))
        trait.res = rbind(trait.res, line.res)
      }
    }
  }
  if(!is.null(trait.res)){
    trait.res = as.data.frame(trait.res, stringsAsFactors = FALSE)
    trait.res = trait.res[, -(which(names(trait.res) == chr.label))]
    names(trait.res)[1:3] = c("chr", "start", "end")
    rownames(trait.res) <- NULL
  }
  return(trait.res)
}





### locus.breaker: works with -log10p
locus.breaker <- function(
    res,
    p.sig     = -log10(5e-08),
    p.limit   = -log10(1e-06),
    hole.size = 250000,
    p.label   = "LOG10P",
    chr.label = "CHROM",
    pos.label = "GENPOS"){


  res <- as.data.frame(res)
  res = res[order(as.numeric(res[, chr.label]), as.numeric(res[,pos.label])), ]
  res = res[which(res[, p.label] > p.limit), ]
  trait.res = c()


  for(j in unique(res[,chr.label])) {
    res.chr = res[which(res[, chr.label] == j), ]
    if (nrow(res.chr) > 1) {
      holes = res.chr[, pos.label][-1] - res.chr[, pos.label][-length(res.chr[,pos.label])]
      gaps = which(holes > hole.size)
      if (length(gaps) > 0) {
        for (k in 1:(length(gaps) + 1)) {
          if (k == 1) {
            res.loc = res.chr[1:(gaps[k]), ]
          }
          else if (k == (length(gaps) + 1)) {
            res.loc = res.chr[(gaps[k - 1] + 1):nrow(res.chr),
            ]
          } else {
            res.loc = res.chr[(gaps[k - 1] + 1):(gaps[k]),
            ]
          }
          if (max(res.loc[, p.label]) > p.sig) {
            start.pos = min(res.loc[, pos.label], na.rm = T)
            end.pos = max(res.loc[, pos.label], na.rm = T)
            chr = j
            best.snp = res.loc[which.max(res.loc[, p.label]),
            ]
            line.res = c(chr, start.pos, end.pos, unlist(best.snp))
            trait.res = rbind(trait.res, line.res)
          }
        }
      } else {
        res.loc = res.chr
        if (max(res.loc[, p.label]) > p.sig) {
          start.pos = min(res.loc[, pos.label], na.rm = T)
          end.pos = max(res.loc[, pos.label], na.rm = T)
          chr = j
          best.snp = res.loc[which.max(res.loc[, p.label]),
          ]
          line.res = c(chr, start.pos, end.pos, unlist(best.snp))
          trait.res = rbind(trait.res, line.res)
        }
      }
    }
    else if (nrow(res.chr) == 1) {
      res.loc = res.chr
      if (max(res.loc[, p.label]) > p.sig) {
        start.pos = min(res.loc[, pos.label], na.rm = T)
        end.pos = max(res.loc[, pos.label], na.rm = T)
        chr = j
        best.snp = res.loc[which.max(res.loc[, p.label]),
        ]
        line.res = c(chr, start.pos, end.pos, unlist(best.snp))
        trait.res = rbind(trait.res, line.res)
      }
    }
  }
  if(!is.null(trait.res)){
    trait.res = as.data.frame(trait.res, stringsAsFactors = FALSE)
    trait.res = trait.res[, -(which(names(trait.res) == chr.label))]
    names(trait.res)[1:3] = c("chr", "start", "end")
    rownames(trait.res) <- NULL
  }
  return(trait.res)
}

# Extract with tabix and load-in
  system(paste0(opt$tabix_bin, " -h -D ", mappa," -R ", study_id, "_chr", chr_tabix, "_snplist.tsv > ", study_id, "_chr", chr_tabix, "_map.tsv"))


### cojo.ht ###
### Performs --cojo-slct first to identify all independent SNPs and --cojo-cond then to condition upon identified SNPs
cojo.ht=function(D=dataset_gwas
                , locus_chr=opt$chr
                , locus_start=opt$start
                , locus_end=opt$end
                , p.thresh=1e-4
                , plink.bin= "/ssu/gassu/software/plink/2.00_20211217/plink2"
                , plink.mem=opt$plink2_mem
                , plink.threads=opt$plink2_threads
                , gcta.bin="/ssu/gassu/software/GCTA/1.94.0beta/gcta64"
                , bfile="/processing_data/shared_datasets/ukbiobank/genotypes/LD_reference/p01_output/ukbb_all_30000_random_unrelated_white_british"
                , maf.thresh=1e-4
                , chr.label=opt$chr_label
                , pos.label=opt$pos_label
                , ea.label=opt$ea_label
                , oa.label=opt$oa_label
                , eaf.label=opt$eaf_label
                , se.label=opt$se_labelsb
                , beta.label=opt$beta_label
                , n.label=opt$n_label
                , p.label=opt$p_label
                ){
    
    chr.label <- sym(opt$chr_label)
    pos.label <- sym(opt$pos_label)
    ea.label  <- sym(opt$ea_label)
    oa.label  <- sym(opt$oa_label)
    eaf.label <- sym(opt$eaf_label)
    se.label  <- sym(opt$se_label)
    beta.label <- sym(opt$beta_label)
    n.label   <- sym(opt$n_label)
    p.label   <- sym(opt$p_label)
    
    random.number=stri_rand_strings(n=1, length=20, pattern = "[A-Za-z0-9]")
    
### Produce two snp.lists: 1) all SNPs to compute allele frequency, 2) only snps included in the locus
    write(D$SNP, ncol=1,file=paste0(random.number,".snp.list"))
    write(D %>% filter(!!chr.label==locus_chr, !!pos.label >= locus_start, !!pos.label <= locus_end) %>% pull(SNP), ncol=1,file=paste0(random.number,"_locus_only.snp.list"))

# Compute allele frequency with Plink
    #system(paste0(plink.bin," --bfile ",bfile, locus_chr, " --extract ",random.number,".snp.list --maf ", maf.thresh, " --make-bed --geno-counts --threads ", plink.threads, " --memory ", plink.mem, " 'require'  --out ", random.number))
    system(paste0(plink.bin," --bfile ",bfile, locus_chr, " --extract ",random.number,".snp.list --make-bed --geno-counts --threads ", plink.threads, " --memory ", plink.mem, " 'require'  --out ", random.number))

    freqs <- fread(paste0(random.number,".gcount"))
    freqs$FreqREF=(freqs$HOM_REF_CT*2+freqs$HET_REF_ALT_CTS)/(2*(rowSums(freqs[,c("HOM_REF_CT", "HET_REF_ALT_CTS", "TWO_ALT_GENO_CTS")])))  #### Why doing all this when plink can directly calculate it with --frq?

# Assign allele frequency from the LD reference
    D <- D %>%
      left_join(freqs %>% dplyr::select(ID,FreqREF,REF), by=c("SNP"="ID")) %>%
      mutate(FREQ=ifelse(REF==!!ea.label, FreqREF, (1-FreqREF))) %>% # NO need to compare the alleles -- if use FREQ or FreqREF rather than !!eaf.label, we get empty bC, bC_SE and pC in COJO *_step1.cma.cojo
      #dplyr::select("SNP","ALLELE0","ALLELE1","FREQ","BETA","SE","LOG10P","N", any_of(c("snp_map","type","sdY","s")))
      dplyr::select("SNP",!!ea.label,!!oa.label,FREQ,!!beta.label,!!se.label,!!p.label,!!n.label, any_of(c("snp_map","type","sdY","s")))
  fwrite(D,file=paste0(random.number,"_sum.txt"), row.names=F,quote=F,sep="\t", na=NA)
  cat("\n\nMerge with LD reference...done.\n\n")

# step1 determine independent snps
  system(paste0(gcta.bin," --bfile ", random.number, " --cojo-p ", p.thresh, " --extract ", random.number, "_locus_only.snp.list --cojo-file ", random.number, "_sum.txt --cojo-slct --out ", random.number, "_step1"))
  # system(paste0(gcta.bin," --bfile ", bfile, locus_chr, " --cojo-p ", p.thresh, " --maf ", maf.thresh, " --extract ", random.number, "_locus_only.snp.list --cojo-file ", random.number, "_sum.txt --cojo-slct --out ", random.number, "_step1"))

  if(file.exists(paste0(random.number,"_step1.jma.cojo"))){
    dataset.list=list()
    ind.snp=fread(paste0(random.number,"_step1.jma.cojo")) %>%
      left_join(D %>% dplyr::select(SNP,any_of(c("snp_map","type","sdY", "s"))), by="SNP")

    dataset.list$ind.snps <- data.frame(matrix(ncol = ncol(ind.snp), nrow = 0))
    colnames(dataset.list$ind.snps) <- colnames(ind.snp)
    dataset.list$results=list()

    if(nrow(ind.snp)>1){
      for(i in 1:nrow(ind.snp)){

        write(ind.snp$SNP[-i],ncol=1,file=paste0(random.number,"_independent.snp"))
        print(ind.snp$SNP[-i])

        #system(paste0(gcta.bin," --bfile ",random.number, " --maf ", maf.thresh, " --extract ",random.number,"_locus_only.snp.list --cojo-file ",random.number,"_sum.txt --cojo-cond ",random.number,"_independent.snp --out ",random.number,"_step2"))
        system(paste0(gcta.bin," --bfile ",random.number, " --extract ",random.number,"_locus_only.snp.list --cojo-file ",random.number,"_sum.txt --cojo-cond ",random.number,"_independent.snp --out ",random.number,"_step2"))

        #### STOP ANALYSIS FOR THAT TOP SNP IN CASE OF COLLINEARITY
        if(!file.exists(paste0(random.number,"_step2.cma.cojo"))){
          cat(paste0("\n****WARNING: COJO has encountered a collinearty problem. Affected SNP will be removed from following analysis****\n\n"))
        } else {
          # Re-add type and sdY/s info, and map SNPs!
          step2.res <- fread(paste0(random.number, "_step2.cma.cojo"), data.table=FALSE) %>%
            left_join(D %>% dplyr::select(SNP, any_of(c("snp_map","type","sdY", "s"))), by="SNP") %>%
            dplyr::mutate(cojo_snp=ind.snp$SNP[i])
          # Add SNPs to the ind.snps dataframe
          dataset.list$ind.snps <- rbind(dataset.list$ind.snps, ind.snp[i,])
          # Add conditioned gwas to the results list
          dataset.list$results[[i]]=step2.res
          names(dataset.list$results)[i]=ind.snp$snp_map[i]
          system(paste0("rm ",random.number,"_step2.cma.cojo"))
        }
      }
    } else {

      ### NB: COJO here is performed ONLY for formatting sakes - No need to condition if only one signal is found!!

      write(ind.snp$SNP,ncol=1,file=paste0(random.number,"_independent.snp"))
      #system(paste0(gcta.bin," --bfile ",random.number," --cojo-p ",p.thresh, " --maf ", maf.thresh, " --extract ",random.number,"_locus_only.snp.list --cojo-file ",random.number,"_sum.txt --cojo-cond ",random.number,"_independent.snp --out ",random.number,"_step2"))
      system(paste0(gcta.bin," --bfile ",random.number," --cojo-p ",p.thresh, " --extract ",random.number,"_locus_only.snp.list --cojo-file ",random.number,"_sum.txt --cojo-cond ",random.number,"_independent.snp --out ",random.number,"_step2"))

      step2.res <- fread(paste0(random.number, "_step2.cma.cojo"), data.table=FALSE) %>%
        left_join(D %>% dplyr::select(SNP,!!ea.label, any_of(c("snp_map","type", "sdY", "s"))), by=c("SNP", "refA"=opt$ea_label))

      #### Add back top SNP, removed from the data frame with the conditioning step
      step2.res <- rbind.fill(
        step2.res,
        ind.snp %>% dplyr::select(-bJ,-bJ_se,-pJ,-LD_r)
      )
      step2.res$cojo_snp <- ind.snp$SNP
      step2.res$bC <- step2.res$b
      step2.res$bC_se <- step2.res$se
      step2.res$pC <- step2.res$p

      dataset.list$ind.snps <- rbind(dataset.list$ind.snps, ind.snp)
      dataset.list$results[[1]]=step2.res
      names(dataset.list$results)[1]=ind.snp$snp_map[1]
    }
    # Remove results df possibly empty (in case of collinearity issue)
    dataset.list$results <- dataset.list$results %>% discard(is.null)
  }

  # make directory to save the files
  #locus <- paste0(locus_chr, "_", locus_start, "_", locus_end)
  #system(paste0("mkdir ", locus))
  #system(paste0("mv *",random.number,"* ", locus))

  system(paste0("rm *",random.number,"*"))
  if(exists("dataset.list")){return(dataset.list)}
}




####

finemap.cojo <- function(D, cs_threshold=0.99){
  cojo_snp <- unique(D$cojo_snp)
# Format input
    D <- D %>%
      dplyr::mutate(varbeta=bC_se^2) %>%
      dplyr::select("snp_map","Chr","bp","bC","varbeta","n","pC","freq","type",any_of(c("sdY","s"))) %>%
      rename("snp"="snp_map","chr"="Chr","position"="bp","beta"="bC","N"="n","pvalues"="pC","MAF"="freq")

  D <- as.list(na.omit(D)) ### move to list and keep unique value of "type" otherwise ANNOYING ERROR!
  D$type <- unique(D$type)
  #if(D$type=="cc"){D$s <- unique(D$s)}else{D$sdY <- unique(D$sdY)}
  D$sdY <- unique(D$sdY)

# Finemap
  fine.res <- coloc::finemap.abf(D) %>%
    mutate(cojo_snp=cojo_snp, bC=c(D$beta, NA)) %>%
    arrange(desc(SNP.PP)) %>%
    mutate(cred.set = cumsum(SNP.PP)) %>%
# Add cojo_hit info, to merge with loci table later
    dplyr::select(snp, position, bC, lABF., cred.set, cojo_snp) %>%
    rename("lABF"="lABF.") %>%
    dplyr::filter(snp!="null")

# Identify SNPs part of the credible set (as specified by cs_threshold)
  w <- which(fine.res$cred.set > cs_threshold)[1]
  cs <- fine.res %>%
#    slice(1:w) %>%
    mutate(is_cs=c(rep(TRUE, w), rep(FALSE, (nrow(fine.res)-w)))) %>%
    select(-cred.set)
  return(cs)
}


my_theme <- function(...){
  theme(
  legend.position = c(.15, .95),
  axis.title.x = element_blank(),
  axis.title = element_text(size = 14, face = 2),
  axis.text =  element_text(size = 12, face = 2)
  )
}

### plot.cojo.ht ###
plot.cojo.ht=function(cojo.ht.obj){

  if(nrow(cojo.ht.obj$ind.snps)>1){

    whole.dataset=c()
    for(i in 1:nrow(cojo.ht.obj$ind.snps)){

      tmp=cojo.ht.obj$results[[i]]
      tmp$signal=cojo.ht.obj$ind.snps$snp_map[i]
      whole.dataset=rbind(whole.dataset,tmp)
    }

    p1 <- ggplot(cojo.ht.obj$results[[i]], aes(x=bp,y=-log10(p))) +
      geom_point(alpha=0.6,size=3)+
      theme_classic() + my_theme() +
      geom_point(data=cojo.ht.obj$ind.snps,aes(x=bp,y=-log10(p),fill=snp_map),size=6,shape=23) +
      guides(fill=guide_legend(title="SNP"))

    p2 <- ggplot(whole.dataset,aes(x=bp,y=-log10(pC),color=signal)) +
      facet_grid(signal~.) +
      geom_point(alpha=0.8,size=3) +
      theme_classic() + my_theme() +
      ggtitle("Conditioned results")

    p3 <- p1/p2 + plot_layout(heights = c(1, nrow(cojo.ht.obj$ind.snps)+0.2))

  } else {

    p3 <- ggplot(cojo.ht.obj$results[[1]], aes(x=bp,y=-log10(p))) +
      geom_point(alpha=0.6,size=3)+
      theme_classic() + my_theme() +
      geom_point(data=cojo.ht.obj$ind.snps,aes(x=bp,y=-log10(p),fill=snp_map),size=6,shape=23)
  }
  (p3)
}



### hcolo.cojo.ht ###
hcolo.cojo.ht=function(df1 = conditional.dataset1,
                       df2 = conditional.dataset2,
                       p1=1e-4,
                       p2=1e-4,
                       p12=1e-5
                       ){

  df1 <- df1 %>% rename("lABF.df1"="lABF")
  df2 <- df2 %>% rename("lABF.df2"="lABF")

  p1 <- coloc:::adjust_prior(p1, nrow(df1), "1")
  p2 <- coloc:::adjust_prior(p2, nrow(df2), "2")

  merged.df <- merge(df1, df2, by = "snp")
  p12 <- coloc:::adjust_prior(p12, nrow(merged.df), "12")

  if(!nrow(merged.df))
    stop("dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified")

  merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + lABF.df2)
## add SNP.PP.H4 - post prob that each SNP is THE causal variant for a shared signal
  my.denom.log.abf <- coloc:::logsum(merged.df$internal.sum.lABF)
  merged.df$SNP.PP.H4 <- exp(merged.df$internal.sum.lABF - my.denom.log.abf)

  pp.abf <- coloc:::combine.abf(merged.df$lABF.df1, merged.df$lABF.df2, p1, p2, p12)
  common.snps <- nrow(merged.df)
  results <- c(nsnps=common.snps, pp.abf)

  colo.res <- list(summary=results, results=merged.df, priors=c(p1=p1,p2=p2,p12=p12))
  class(colo.res) <- c("coloc_abf", class(colo.res))

## Save coloc summary
  colo.sum <- data.frame(t(colo.res$summary))

## Save coloc result by SNP
  colo.full_res <- colo.res$results %>% dplyr::select(snp,lABF.df1,lABF.df2,SNP.PP.H4)

## Organise all in a list ( composed of summary + results)
  coloc.final <- list(summary=colo.sum, results=colo.full_res)
  return(coloc.final)
}
