# Load packages
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
#suppressMessages(library(R.utils))
# suppressMessages(library(corrplot))
suppressMessages(library(coloc))
#suppressMessages(library(bigsnpr))
#suppressMessages(library(ggplot2))
#suppressMessages(library(cowplot))
suppressMessages(library(stringi))
#suppressMessages(library(stringr))
suppressMessages(library(patchwork))
# suppressMessages(library(reshape2))
# suppressMessages(library(RColorBrewer))
# suppressMessages(library(igraph))
#suppressMessages(library(purrr))
#suppressMessages(library(tidyr))
#suppressMessages(library(plyr))
# suppressMessages(library(Gviz))
# suppressMessages(library(EnsDb.Hsapiens.v75))
# suppressMessages(library(Matrix))
#suppressMessages(library(dplyr))
suppressMessages(library(Rmpfr))


## Source function for locus.breaker: works with -log10p
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
            res.loc = res.chr[(gaps[k - 1] + 1):nrow(res.chr), ]
          } else {
            res.loc = res.chr[(gaps[k - 1] + 1):(gaps[k]), ]
          }
          if (max(res.loc[, p.label]) > p.sig) {
            start.pos = min(res.loc[, pos.label], na.rm = T)
            end.pos   = max(res.loc[, pos.label], na.rm = T)
            chr = j
            best.snp  = res.loc[which.max(res.loc[, p.label]), ]
            line.res  = c(chr, start.pos, end.pos, unlist(best.snp))
            trait.res = rbind(trait.res, line.res)
          }
        }
      } else {
        res.loc = res.chr
        if (max(res.loc[, p.label]) > p.sig) {
          start.pos = min(res.loc[, pos.label], na.rm = T)
          end.pos   = max(res.loc[, pos.label], na.rm = T)
          chr = j
          best.snp  = res.loc[which.max(res.loc[, p.label]), ]
          line.res  = c(chr, start.pos, end.pos, unlist(best.snp))
          trait.res = rbind(trait.res, line.res)
        }
      }
    }
    else if (nrow(res.chr) == 1) {
      res.loc = res.chr
      if (max(res.loc[, p.label]) > p.sig) {
        start.pos = min(res.loc[, pos.label], na.rm = T)
        end.pos   = max(res.loc[, pos.label], na.rm = T)
        chr = j
        best.snp  = res.loc[which.max(res.loc[, p.label]), ]
        line.res  = c(chr, start.pos, end.pos, unlist(best.snp))
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




# Function to handle missing values once computing p-value from cojo results
safe_pnorm <- function(b, se, p=FALSE) {
  
  # Ensure the vectors are of the same length
  if(length(b) != length(se)) {
    stop("Beta and SE must be of the same length")
  }
  
  k  <- length(b)
  b  <- as.numeric(b)
  se <- as.numeric(se)
  
  # Initialize result vector with NA values
  result <- rep(NA, k)
  
    # Identify non-missing and non-zero indices
  i <- which(!is.na(b) & !is.na(se) & se != 0)
  
  # compute z-score and take absolute, raise digits with mpfr, apply pnorm for non-missing values
  z_score <- b[i] / se[i]
  z_mpfr <- Rmpfr::mpfr(- abs(z_score), 120)
  p_mpfr <- 2 * pnorm(z_mpfr)
  mlog10p <- - log10(p_mpfr)
  
  # print p-value in character format and mlog10p in numeric
  if(p==TRUE){
    # reformat to mpfr character, then to numeric (don't set digits for MLOG10P)
    mlog10p_mpfr <- Rmpfr::formatMpfr(p_mpfr, scientific = TRUE, digits = 6)
    result[i] <- mlog10p_mpfr
  } else {
    mlog10p_mpfr <- Rmpfr::formatMpfr(mlog10p, scientific = TRUE)
    result[i] <- as.numeric(mlog10p_mpfr)
  }

  return(result)
}


### cojo.ht ###
### Performs --cojo-slct first to identify all independent SNPs and --cojo-cond then to condition upon identified SNPs
cojo.ht=function(D=dataset_gwas
                , locus_chr=opt$chr
                , locus_start=opt$start
                , locus_end=opt$end
                , p.cojo=1e-4
                , p.jumper=1e-4
                , plink.bin= "/ssu/gassu/software/plink/2.00_20211217/plink2"
                , plink.mem=opt$plink2_mem
                , plink.threads=opt$plink2_threads
                , gcta.bin="/ssu/gassu/software/GCTA/1.94.0beta/gcta64"
                , bfile="/processing_data/shared_datasets/ukbiobank/genotypes/LD_reference/p01_output/ukbb_all_30000_random_unrelated_white_british"
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

# Compute allele frequency with Plink -- we removed "--maf ", maf.thresh, from original script
    #system(paste0(plink.bin," --bfile ",bfile, locus_chr, " --extract ",random.number,".snp.list --make-bed --geno-counts --threads ", plink.threads, " --memory ", plink.mem, " 'require'  --out ", random.number))

    #freqs <- fread(paste0(random.number,".gcount"))
    #freqs$FreqREF=(freqs$HOM_REF_CT*2+freqs$HET_REF_ALT_CTS)/(2*(rowSums(freqs[,c("HOM_REF_CT", "HET_REF_ALT_CTS", "TWO_ALT_GENO_CTS")])))  #### Why doing all this when plink can directly calculate it with --frq?

# Assign allele frequency from the LD reference
    D <- D %>%
      filter(!!chr.label==locus_chr, !!pos.label >= locus_start, !!pos.label <= locus_end) %>% #bounding the GWAS to variants only falling at the locus solves the problem of phenotipc variance = 0 by GCTA-cojo
      dplyr::mutate_at(vars(SNP), as.character) %>% # to ensure class of joint column is the same
      #left_join(freqs %>% dplyr::select(ID,FreqREF,REF), by=c("SNP"="ID")) %>%
      #mutate(FREQ=ifelse(REF==!!oa.label, FreqREF, (1-FreqREF))) %>% # Need to compare the alleles to avoid getting only one independent SNP per locus, as a results of zero phenotypic variance.
      dplyr::select("SNP",!!ea.label,!!oa.label,!!eaf.label,!!beta.label,!!se.label,!!p.label,!!n.label, any_of(c("type", "N_GWAS")))
  fwrite(D,file=paste0(random.number,"_sum.txt"), row.names=F,quote=F,sep="\t", na=NA)
  cat("\n\nMerge with LD reference...done.\n\n")

# step1 determine independent snps -- we removed "--maf ", maf.thresh, from original script
  system(paste0(gcta.bin," --bfile ", bfile, locus_chr, " --cojo-p ", p.cojo, " --extract ", random.number, "_locus_only.snp.list --cojo-file ", random.number, "_sum.txt --cojo-slct --out ", random.number, "_step1"))

  if(file.exists(paste0(random.number,"_step1.jma.cojo"))){
    dataset.list=list()
    ind.snp=data.table::fread(paste0(random.number,"_step1.jma.cojo")) %>%
      dplyr::mutate(
        SNP = as.character(SNP),   # ensure class of joint column is the same
        p_org = p,                 # keep original p  from GCTA
        pJ_org = pJ,               # keep original pJ from GCTA
        p  = safe_pnorm(b, se, p = TRUE),     # compute unconditional p-value
        pJ = safe_pnorm(bJ, bJ_se, p = TRUE), # compute joint p-value
        mlog10p  = safe_pnorm(b, se),         # compute unconditional MLOG10P
        mlog10pJ = safe_pnorm(bJ, bJ_se)      # compute joint MLOG10P
      ) %>%
      dplyr::relocate(Chr:freq, freq_geno, b:p, p_org, mlog10p, n:pJ, pJ_org, mlog10pJ) %>%  # tidying columns order
      filter(mlog10p > - log10(p.jumper)) %>%   # avoid including any non-significant independent variant to conditional model  
      left_join(D %>% dplyr::select(SNP,any_of(c("type", "N_GWAS", opt$p_label))), by="SNP")

    dataset.list$ind.snps <- data.frame(matrix(ncol = ncol(ind.snp), nrow = 0))
    colnames(dataset.list$ind.snps) <- colnames(ind.snp)
    dataset.list$results=list()

    if(nrow(ind.snp)>1){
      for(i in 1:nrow(ind.snp)){

        write(ind.snp$SNP[-i],ncol=1,file=paste0(random.number,"_independent.snp"))
        print(ind.snp$SNP[-i])

        #Removing "--maf ", maf.thresh, from original script
        system(paste0(gcta.bin," --bfile ",bfile, locus_chr, " --extract ",random.number,"_locus_only.snp.list --cojo-file ",random.number,"_sum.txt --cojo-cond ",random.number,"_independent.snp --out ",random.number,"_step2"))

        #### STOP ANALYSIS FOR THAT TOP SNP IN CASE OF COLLINEARITY
        if(!file.exists(paste0(random.number,"_step2.cma.cojo"))){
          cat(paste0("\n****WARNING: COJO has encountered a collinearty problem. Affected SNP will be removed from following analysis****\n\n"))
        } else {
          # Re-add type and sdY/s info, and map SNPs!
          step2.res <- fread(paste0(random.number, "_step2.cma.cojo"), data.table=FALSE) %>%
            dplyr::mutate(
              SNP = as.character(SNP),   # ensure class of joint column is the same
              p_org = p,                 # keep original p  from GCTA
              pC_org = pC,               # keep original pC from GCTA
              p  = safe_pnorm(b, se, p = TRUE),     # compute unconditional p-value
              pC = safe_pnorm(bC, bC_se, p = TRUE), # compute joint p-value
              mlog10p  = safe_pnorm(b, se),         # compute unconditional MLOG10P
              mlog10pC = safe_pnorm(bC, bC_se)      # compute joint MLOG10P
            ) %>%
            dplyr::relocate(Chr:freq, freq_geno, b:p, p_org, mlog10p, n:pC, pC_org, mlog10pC) %>%  # tidying columns order
            left_join(D %>% dplyr::select(SNP, any_of(c("type", "N_GWAS", opt$p_label))), by="SNP") %>%
            dplyr::mutate(cojo_snp=ind.snp$SNP[i])
          # Add SNPs to the ind.snps dataframe
          dataset.list$ind.snps <- rbind(dataset.list$ind.snps, ind.snp[i,])
          # Add conditioned gwas to the results list
          dataset.list$results[[i]]=step2.res
          names(dataset.list$results)[i]=ind.snp$SNP[i]
          system(paste0("rm ",random.number,"_step2.cma.cojo"))
        }
      }
    } else {

      ### NB: COJO here is performed ONLY for formatting sakes - No need to condition if only one signal is found!!

      write(ind.snp$SNP,ncol=1,file=paste0(random.number,"_independent.snp"))
      #Removing "--maf ", maf.thresh, from original script
      system(paste0(gcta.bin," --bfile ",bfile, locus_chr," --cojo-p ",p.cojo, " --extract ",random.number,"_locus_only.snp.list --cojo-file ",random.number,"_sum.txt --cojo-cond ",random.number,"_independent.snp --out ",random.number,"_step2"))

      step2.res <- fread(paste0(random.number, "_step2.cma.cojo"), data.table=FALSE) %>%
        dplyr::mutate(
          SNP  = as.character(SNP),  # ensure class of joint column is the same
          refA = as.character(refA), # ensure class of joint column is the same
          p_org = p,                 # keep original p  from GCTA
          p  = safe_pnorm(b, se, p = TRUE),     # compute unconditional p-value
          mlog10p  = safe_pnorm(b, se),         # compute unconditional MLOG10P
          ) %>%
        dplyr::relocate(Chr:freq, freq_geno, b:p, p_org, mlog10p) %>%  # tidying columns order
        left_join(D %>% dplyr::select(SNP,!!ea.label, any_of(c("type", "N_GWAS", opt$p_label))), by=c("SNP", "refA"=opt$ea_label))

      #### Add back top SNP, removed from the data frame with the conditioning step
      step2.res <- plyr::rbind.fill(
        step2.res,
        ind.snp %>% dplyr::select(-bJ,-bJ_se,-pJ,-LD_r,-pJ_org,-mlog10pJ)
      )

      step2.res$cojo_snp <- ind.snp$SNP
      step2.res$bC <- step2.res$b
      step2.res$bC_se <- step2.res$se
      step2.res$pC <- step2.res$p
      step2.res$pC_org <- step2.res$p_org
      step2.res$mlog10pC <- step2.res$mlog10p

      dataset.list$ind.snps <- rbind(dataset.list$ind.snps, ind.snp)
      dataset.list$results[[1]]=step2.res
      names(dataset.list$results)[1]=ind.snp$SNP[1]
    }
    # Remove results df possibly empty (in case of collinearity issue)
    dataset.list$results <- dataset.list$results %>% discard(is.null)
  }

  # make directory to save the files
  # locus <- paste0(locus_chr, "_", locus_start, "_", locus_end)
  # system(paste0("mkdir ", locus))
  # system(paste0("mv *",random.number,"* ", locus))

  system(paste0("rm *",random.number,"*"))
  if(exists("dataset.list")){return(dataset.list)}
}



#### finemap.abf from coloc but modified to not include a prior/null in finemapping
finemap.abf_NO_PRIOR <- function(dataset) {
  
  coloc::check_dataset(dataset,"")  # Check all required input is provided  
  df <- coloc::process.dataset(d=dataset, suffix="") # Compute lABF for each SNP 
  
  # Scale  
  my.denom.log.abf <- coloc:::logsum(df$lABF)
  df$SNP.PP <- exp(df$lABF- my.denom.log.abf)
  
  return(df)
}


#### finemap.cojo

finemap.cojo <- function(D, cs_threshold=0.99){
  cojo_snp <- unique(D$cojo_snp)
  # Format input
  D <- D %>%
    dplyr::mutate(
      N = N_GWAS,
      varbeta = bC_se^2,
      varbeta_uncond = se^2,
      MAF = if_else(freq > 0.5, 1 - freq, freq),
      sdY = coloc:::sdY.est(varbeta_uncond, MAF, N)
    ) %>% # to ensure having required column for finemap.abf_NO_PRIOR()
    dplyr::select(Chr, bp, SNP, freq, freq_geno, MAF, b, se, varbeta_uncond, p, mlog10p, bC, bC_se, varbeta, pC, mlog10pC, type, n, N, sdY) %>%
    dplyr::rename(snp=SNP, chr=Chr, position=bp, beta=bC, pvalues=pC) # rename as coloc::process.dataset() requires

  D_list <- as.list(na.omit(D)) ### move to list and keep unique value of "type" otherwise ANNOYING ERROR!
  D_list$type <- unique(D_list$type)
  #if(D$type=="cc"){D$s <- unique(D$s)}else{D$sdY <- unique(D$sdY)}
  D$sdY <- unique(D$sdY)

# Finemap
  fine.res <- finemap.abf_NO_PRIOR(D_list) %>%
    dplyr::left_join(D %>% dplyr::select(chr, position, snp, freq, freq_geno, MAF, b, se, varbeta_uncond, p, mlog10p, beta, bC_se, varbeta, pvalues, mlog10pC, type, n, N, sdY), join_by("snp")) %>%
    dplyr::rename(bC=beta, pC=pvalues, "lABF"="lABF.") %>%
    arrange(desc(SNP.PP)) %>%
    dplyr::mutate(cred.set = cumsum(SNP.PP), cojo_snp=cojo_snp) %>%
    dplyr::select(snp, chr, freq, freq_geno, MAF, b, se, varbeta_uncond, p, mlog10p, bC, bC_se, varbeta, pC, mlog10pC, type, n, N, sdY, lABF, SNP.PP, cred.set, cojo_snp) # Add cojo_hit info, to merge with loci table later

  # Identify SNPs part of the credible set (as specified by cs_threshold)
  w <- which(fine.res$cred.set > cs_threshold)[1]
  cs <- fine.res %>%
    dplyr::mutate(is_cs = c(rep(TRUE, w), rep(FALSE, (nrow(fine.res)-w)))) %>%
    dplyr::select(- cred.set)
  
  return(cs)
}


my_theme <- function(...){
  theme(
  legend.position = c(.15, .65),
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
      tmp$signal=cojo.ht.obj$ind.snps$SNP[i]
      whole.dataset=rbind(whole.dataset,tmp)
      
      # set break points for y-axis
      max_p <- max(whole.dataset %>% select(mlog10p))
      max_y <- max_p + (max_p / 10)
      y_dot <- pretty(range(whole.dataset %>% select(mlog10p)), n = 10)
      
    }
    
    p1 <- ggplot(cojo.ht.obj$results[[i]], aes(x=bp,y=mlog10p)) +
      geom_point(alpha=0.6,size=3)+
      scale_y_continuous(breaks = y_dot, limits = c(0, max_p)) +
      theme_classic() + my_theme() +
      geom_point(data=cojo.ht.obj$ind.snps,aes(x=bp,y=mlog10p,fill=SNP),size=6,shape=23) +
      guides(fill=guide_legend(title="SNP"))
    
    p2 <- ggplot(whole.dataset,aes(x=bp,y=mlog10pC,color=signal)) +
      facet_grid(signal~.) +
      geom_point(alpha=0.8,size=3) +
      theme_classic() + my_theme() +
      ggtitle("Conditioned results")
    
    p3 <- p1/p2 + patchwork::plot_layout(heights = c(1, nrow(cojo.ht.obj$ind.snps)+0.2))
    
  } else {
    
    # set break points for y-axis
    max_p <- max(cojo.ht.obj$results[[1]] %>% select(mlog10p))
    max_y <- max_p + (max_p / 10)
    y_dot <- pretty(range(cojo.ht.obj$results[[1]] %>% select(mlog10p)), n = 10)
    
    p3 <- ggplot(cojo.ht.obj$results[[1]], aes(x=bp,y=mlog10p)) +
      geom_point(alpha=0.6,size=3)+
      scale_y_continuous(breaks = y_dot, limits = c(0, max_p)) +
      theme_classic() + my_theme() +
      geom_point(data=cojo.ht.obj$ind.snps,aes(x=bp,y=mlog10p,fill=SNP),size=6,shape=23)
  }
  return(p3)
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

# function to rename columns, change calss to list, and take unique values as suggested here:
# https://cran.r-project.org/web/packages/coloc/vignettes/a02_data.html
prepare4coloc <- function(data){
  temp  <- data %>% dplyr::rename(beta=bC, pvalues=pC)
  odata <- as.list(na.omit(temp))
  odata$type <- unique(odata$type)
  odata$sdY <- unique(odata$sdY)
  return(odata)
}
