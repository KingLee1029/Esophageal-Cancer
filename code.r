library(TwoSampleMR)
library(MRPRESSO)
library(GagnonMR)
library(ieugwasr)
library(ggplot2)
library(stringr)
library(readxl)
Gwasid <- "finngen_R10_C3_OESOPHAGUS_EXALLC"
DirsOutFile <- "./OutComeTable"
DireQTL <- "./TableMR/"
DirBedFile <- "./EUR"
DirPlink <- "./plink_linux_x86_64_20230116"
CutP1 <- 1e-5
CutP2 <- 5e-5
SensiLoo <- 1;HeteroIVW <- 1
LimitF <- 1;CutF <- 10
clumpLD = TRUE;clump_kb = 10000;clump_r2 = 0.001;clump_p1 = 1
source(paste0(DirPlink,"/ld_clump_localX.R"))
TotalOut <- data.table::fread(paste0(DirsOutFile,"/",Gwasid,".txt.gz"), header = TRUE, sep = "\t")$SNP


CandGen <- as.data.frame(read_xlsx("NIHMS80906-supplement-Table_S1.xlsx",sheet = 1,col_names = T))
gene_human <- read.table("./mart_export.txt",sep="\t",header=T,check.names=F);gene_human <- gene_human[gene_human$`Gene name` != "",]
for (i in 1:nrow(CandGen)) {
 if(is.na(CandGen$hgnc_names[i])){
   if(nrow(gene_human[gene_human$`Gene stable ID` == CandGen$ensembl_gene_id[i],]) != 0){
     CandGen$hgnc_names[i] <- gene_human[gene_human$`Gene stable ID` == CandGen$ensembl_gene_id[i],"Gene name"]
   }
 }
}
gene_human <- CandGen[,c("ensembl_gene_id", "hgnc_names")];colnames(gene_human) <- c("Gene stable ID", "Gene name")
CandGen <- unique(CandGen$hgnc_names);length(CandGen)
outTab0 <- outTab1 <- outTab2 <- KeyGen <- NULL
TotalFiles <- list.files(DireQTL);length(TotalFiles)
TotalFiles <- TotalFiles[str_split_fixed(TotalFiles,"\\-|\\.",4)[,3] %in% gene_human[gene_human$`Gene name` %in% CandGen,]$`Gene stable ID`];length(TotalFiles)
Panel <- paste0("Results_",gsub("-","_",Gwasid));dir.create(Panel)

FileList <- lapply(TotalFiles,function(f){
  tmp <- read_exposure_data(filename = paste0(DireQTL,"/",f), sep = ",",
                            clump = F,
                            phenotype_col = "phenotype",
                            chr_col = "chr",
                            pos_col = "pos",
                            snp_col = "ID",
                            beta_col = "beta",
                            se_col = "se",
                            eaf_col = "eaf",
                            effect_allele_col = "effect allele",
                            other_allele_col = "other allele",
                            pval_col = "pval",
                            samplesize_col = "samplesize")
  tmp$Symbol <- gene_human[gene_human$`Gene stable ID` == str_split_fixed(f,"\\-|\\.",4)[,3],]$`Gene name`
  tmp$id.exposure <- str_split_fixed(f,"\\-|\\.",4)[,3]
  tmp$Source <- f
  return(tmp)
})
combind_expd <- do.call(rbind,FileList)
expdX <- combind_expd[which(combind_expd$pval.exposure < CutP1),]
expdX <- expdX[expdX$Source %in% names(table(expdX[expdX$SNP %in% TotalOut,]$Source)[table(expdX[expdX$SNP %in% TotalOut,]$Source) > 1]),]
Sys.sleep(1)

if(nrow(expdX) > 0){
  
  if(clumpLD){
    SetNamesNum <- c(grep("pval.exposure",colnames(expdX)), grep("SNP",colnames(expdX)))
    colnames(expdX)[SetNamesNum] <- c("pval", "rsid")
    expdList <- lapply(unique(expdX$Source),function(x){
      expTmp <- expdX[expdX$Source == x,]
      ld_clump_localX(dat = expTmp, clump_kb = clump_kb, clump_r2 = clump_r2,
                      clump_p = clump_p1, bfile = DirBedFile,
                      plink_bin = paste0(DirPlink,"/plink"))
    })
    Sys.sleep(1)
    expd0 <- do.call(rbind,expdList)
    colnames(expd0)[SetNamesNum] <- c("pval.exposure", "SNP")
  }else{
    expd0 <- expdX
  }
  expd0 <- expd0[expd0$Source %in% names(table(expd0[expd0$SNP %in% TotalOut,]$Source)[table(expd0[expd0$SNP %in% TotalOut,]$Source) > 1]),]
  Sys.sleep(1)
  
  if(nrow(expd0) > 0){
    
    expd0$F.value <- (expd0$beta.exposure)^2 / (expd0$se.exposure)^2
    if(LimitF == 1){
      expd0 <- expd0[which(expd0$F.value > CutF),]
      expd0 <- expd0[expd0$Source %in% names(table(expd0[expd0$SNP %in% TotalOut,]$Source)[table(expd0[expd0$SNP %in% TotalOut,]$Source) > 1]),]
      Sys.sleep(1)
    }
    
    if(nrow(expd0) > 0){
      
      outTab0 <- rbind(outTab0, expd0)
      outd0 <- read_outcome_data(snps = unique(expd0$SNP),
                                 filename = paste0(DirsOutFile,"/",Gwasid,".txt.gz"),
                                 sep = "\t",
                                 snp_col = "SNP",
                                 chr_col = "chr.outcome",
                                 pos_col = "pos.outcome",
                                 effect_allele_col = "effect_allele.outcome",
                                 other_allele_col = "other_allele.outcome",
                                 beta_col = "beta.outcome",
                                 se_col = "se.outcome",
                                 pval_col = "pval.outcome",
                                 eaf_col = "eaf.outcome",
                                 samplesize_col = "samplesize.outcome",
                                 id_col = "id.outcome",
                                 phenotype_col = "outcome")
      outd0 <- unique(outd0[which(outd0$pval.outcome > CutP2),])
      if(nrow(outd0) > 0){
        
        expd0 <- expd0[expd0$SNP %in% outd0$SNP,]
        expd0 <- expd0[expd0$Source %in% names(table(expd0[expd0$SNP %in% TotalOut,]$Source)[table(expd0[expd0$SNP %in% TotalOut,]$Source) > 1]),]
        for (i in intersect(gene_human[gene_human$`Gene stable ID` %in% expd0$id.exposure,]$`Gene name`, CandGen)) {

          expd1 <- expd0[expd0$id.exposure %in% gene_human[gene_human$`Gene name` == i,]$`Gene stable ID`,]
          outd1 <- outd0[outd0$SNP %in% expd1$SNP,]
          if(nrow(outd1) > 1){
            
            dat1 <- TwoSampleMR::harmonise_data(expd1, outd1)
            outTab1 <- merge(outTab1, dat1, all = T)
            if(nrow(outTab1) == 0){
              outTab1 <- merge(outTab1, dat1, all = T)
            }
            if(nrow(dat1) > 1){
              
              res <- TwoSampleMR::mr(dat1, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
              if(nrow(res) > 0){
                if(min(res$pval) < 0.05){
                  dir.create(paste0(Panel,"/",i))
                  
                  p <- mr_scatter_plot(res, dat1)
                  pdf(paste0(Panel,"/",i,"/PointPlot.pdf"),width = 8,height = 8)
                  print(p)
                  invisible(dev.off())
                  
                  
                  res_OR <- generate_odds_ratios(res)
                  write.table(res_OR,file = paste0(Panel,"/",i,"/odds_ratios.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
                  
                
                  res_ht <- mr_heterogeneity(dat1)
                  write.table(res_ht,file = paste0(Panel,"/",i,"/heterogeneity.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
                  
                
                  dat_PRS <- dat1[rowSums(is.na(dat1[,c("beta.outcome", "beta.exposure", "se.outcome", "se.exposure")])) == 0,]
                  if(nrow(dat_PRS) > 3){
                    res_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome",
                                            SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_PRS,
                                            NbDistribution = 1000,  SignifThreshold = 0.05, seed = 123456)
                    b <- gsub("<","",res_presso$`MR-PRESSO results`$`Global Test`$Pvalue)
                    save(res_presso,file = paste0(Panel,"/",i,"/presso_",b,"_.rda"))
                  }
                  
                 
                  res_pl <- mr_pleiotropy_test(dat1)
                  write.table(res_pl,file = paste0(Panel,"/",i,"/pleiotropy.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
                  
                 
                  res_loo <- mr_leaveoneout(dat1)
                  if(max(table(res_loo$SNP)) > 1){
                    res_loo <- res_loo[order(res_loo$SNP, res_loo$p, decreasing = T),]
                    res_loo <- res_loo[!duplicated(res_loo$SNP, fromLast = T),]
                  }
                  p <- mr_leaveoneout_plot(res_loo)
                  pdf(paste0(Panel,"/",i,"/LeaveoneoutPlot.pdf"),width = 8,height = ((nrow(res_loo)/10) + 3))
                  print(p)
                  invisible(dev.off())
                  write.table(res_loo,file = paste0(Panel,"/",i,"/leaveoneout.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
                  
                 
                  res_single <- mr_singlesnp(dat1)
                  if(max(table(res_single$SNP)) > 1){
                    res_single <- res_single[order(res_single$SNP, res_single$p, decreasing = T),]
                    res_single <- res_single[!duplicated(res_single$SNP, fromLast = T),]
                  }
                  p <- mr_forest_plot(res_single)
                  pdf(paste0(Panel,"/",i,"/ForestPlot.pdf"),width = 6.5,height = ((nrow(res_single)/10) + 3))
                  print(p)
                  invisible(dev.off())
                  p <- mr_funnel_plot(res_single)
                  pdf(paste0(Panel,"/",i,"/FunnelPlot.pdf"),width = 8,height = 8)
                  print(p)
                  invisible(dev.off())
                  p <- mr_density_plot(res_single, res)
                  pdf(paste0(Panel,"/",i,"/DensityPlot.pdf"),width = 8,height = 5.5)
                  print(p)
                  invisible(dev.off())
                  write.table(res_single,file = paste0(Panel,"/",i,"/singlesnp.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
                }
                
                res$Symbol <- i
                outTab2 <- rbind(outTab2, res)
                if(dir.exists(paste0(Panel,"/",i))){
                  
                  if(res[res$method == "Inverse variance weighted",]$pval < 0.05){
                    
                    KeyGen <- c(KeyGen, i)
                   
                    if(SensiLoo == 1){
                      if(unique(res$nsnp) > 1){
                        if(length(setdiff(res_loo$SNP, "All")) == 0){
                          repeat{
                            a <- file.rename(paste0(Panel,"/",i),paste0(Panel,"/lo_",i))
                            if(a){
                              break
                            }
                          }
                          KeyGen <- setdiff(KeyGen, i)
                        }
                      }
                    }
                    
                    if(HeteroIVW == 1){
                      if(unique(res$nsnp) > 1){
                        if(res_ht[res_ht$method == "Inverse variance weighted",]$Q_pval < 0.05){
                          repeat{
                            a <- file.rename(paste0(Panel,"/",i),paste0(Panel,"/ht_",i))
                            if(a){
                              break
                            }
                          }
                          KeyGen <- setdiff(KeyGen, i)
                        }
                      }
                    }
                  }else{
                    repeat{
                      a <- file.rename(paste0(Panel,"/",i),paste0(Panel,"/pv_",i))
                      if(a){
                        break
                      }
                    }
                  }
                }
              }
            }
            Sys.sleep(1)
            rm(list = intersect(ls(),c("expd1", "outd1", "res", "dat1", "res_presso", "res_OR", "res_ht", "res_pl", "res_loo", "res_single", "p")))
          }
        }
      }
    }
  }
}
Sys.sleep(1)
write.table(outTab0,file = paste0("./",Panel,"/0.table.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
write.table(outTab1,file = paste0("./",Panel,"/1.table.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
write.table(outTab2,file = paste0("./",Panel,"/2.table.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

KeyGen <- unique(KeyGen);length(KeyGen)
tmp <- outTab2[outTab2$Symbol %in% KeyGen & outTab2$method == "Inverse variance weighted",]
if(nrow(tmp[which(tmp$pval < 0.05),]) > 1){
  cat("IVW pval < 0.05")
  KeyGen <- tmp[which(tmp$pval < 0.05),]$Symbol
}
length(KeyGen)
save(KeyGen,file = "./KeyGen.rda")





library(TwoSampleMR)
library(MRPRESSO)
library(GagnonMR)
library(ieugwasr)
library(ggplot2)
library(stringr)
library(readxl)
Gwasid <- "GCST90018841"
DirsOutFile <- "/OutComeTable"
DireQTL <- "/TableMR/"
DirBedFile <- "/EUR"
DirPlink <- "/plink_win64"
CutP1 <- 1e-5
CutP2 <- 5e-5
SensiLoo <- 1;HeteroIVW <- 1
LimitF <- 1;CutF <- 10
clumpLD = TRUE;clump_kb = 10000;clump_r2 = 0.001;clump_p1 = 1
source(paste0(DirPlink,"/ld_clump_localX.R"))
TotalOut <- data.table::fread(paste0(DirsOutFile,"/",Gwasid,".txt.gz"), header = TRUE, sep = "\t")$SNP


CandGen <- read.table("./table.OR.txt",sep="\t",header=T,check.names=F);nrow(CandGen)
gene_human <- CandGen[,c("id.exposure", "Symbol")];colnames(gene_human) <- c("Gene stable ID", "Gene name");nrow(gene_human)
CandGen <- unique(CandGen$Symbol);length(CandGen)
outTab0 <- outTab1 <- outTab2 <- KeyGen <- NULL
TotalFiles <- list.files(DireQTL);length(TotalFiles)
TotalFiles <- TotalFiles[str_split_fixed(TotalFiles,"\\-|\\.",4)[,3] %in% gene_human[gene_human$`Gene name` %in% CandGen,]$`Gene stable ID`];length(TotalFiles)
Panel <- paste0("Results_",gsub("-","_",Gwasid));dir.create(Panel)

FileList <- lapply(TotalFiles,function(f){
  if(nrow(gene_human[gene_human$`Gene stable ID` == str_split_fixed(f,"\\-|\\.",4)[,3],]) > 0){
    tmp <- read_exposure_data(filename = paste0(DireQTL,"/",f), sep = ",",
                              clump = F,
                              phenotype_col = "phenotype",
                              chr_col = "chr",
                              pos_col = "pos",
                              snp_col = "ID",
                              beta_col = "beta",
                              se_col = "se",
                              eaf_col = "eaf",
                              effect_allele_col = "effect allele",
                              other_allele_col = "other allele",
                              pval_col = "pval",
                              samplesize_col = "samplesize")
    if(nrow(gene_human[gene_human$`Gene stable ID` == str_split_fixed(f,"\\-|\\.",4)[,3],]) > 0){
      tmp <- tmp[which(tmp$pval.exposure < CutP1),]
      if(length(intersect(tmp$SNP, TotalOut)) > 0){
        tmp$Symbol <- gene_human[gene_human$`Gene stable ID` == str_split_fixed(f,"\\-|\\.",4)[,3],]$`Gene name`
        tmp$id.exposure <- str_split_fixed(f,"\\-|\\.",4)[,3]
        tmp$Source <- f
        return(tmp)
      }
    }
  }
})
combind_expd <- do.call(rbind,FileList)
expdX <- combind_expd[which(combind_expd$pval.exposure < CutP1),]
expdX <- expdX[expdX$Source %in% names(table(expdX[expdX$SNP %in% TotalOut,]$Source)[table(expdX[expdX$SNP %in% TotalOut,]$Source) > 1]),]
Sys.sleep(1)

if(nrow(expdX) > 0){
  
  if(clumpLD){
    SetNamesNum <- c(grep("pval.exposure",colnames(expdX)), grep("SNP",colnames(expdX)))
    colnames(expdX)[SetNamesNum] <- c("pval", "rsid")
    expdList <- lapply(unique(expdX$Source),function(x){
      expTmp <- expdX[expdX$Source == x,]
      ld_clump_localX(dat = expTmp, clump_kb = clump_kb, clump_r2 = clump_r2,
                      clump_p = clump_p1, bfile = DirBedFile,
                      plink_bin = paste0(DirPlink,"/plink.exe"))
    })
    Sys.sleep(1)
    expd0 <- do.call(rbind,expdList)
    colnames(expd0)[SetNamesNum] <- c("pval.exposure", "SNP")
  }else{
    expd0 <- expdX
  }
  expd0 <- expd0[expd0$Source %in% names(table(expd0[expd0$SNP %in% TotalOut,]$Source)[table(expd0[expd0$SNP %in% TotalOut,]$Source) > 1]),]
  Sys.sleep(1)
  
  if(nrow(expd0) > 0){
    
    expd0$F.value <- (expd0$beta.exposure)^2 / (expd0$se.exposure)^2
    if(LimitF == 1){
      expd0 <- expd0[which(expd0$F.value > CutF),]
      expd0 <- expd0[expd0$Source %in% names(table(expd0[expd0$SNP %in% TotalOut,]$Source)[table(expd0[expd0$SNP %in% TotalOut,]$Source) > 1]),]
      Sys.sleep(1)
    }
    
    if(nrow(expd0) > 0){
      
      outTab0 <- rbind(outTab0, expd0)
      outd0 <- read_outcome_data(snps = unique(expd0$SNP),
                                 filename = paste0(DirsOutFile,"/",Gwasid,".txt.gz"),
                                 sep = "\t",
                                 snp_col = "SNP",
                                 chr_col = "chr.outcome",
                                 pos_col = "pos.outcome",
                                 effect_allele_col = "effect_allele.outcome",
                                 other_allele_col = "other_allele.outcome",
                                 beta_col = "beta.outcome",
                                 se_col = "se.outcome",
                                 pval_col = "pval.outcome",
                                 eaf_col = "eaf.outcome",
                                 samplesize_col = "samplesize.outcome",
                                 id_col = "id.outcome",
                                 phenotype_col = "outcome")
      outd0 <- unique(outd0[which(outd0$pval.outcome > CutP2),])
      if(nrow(outd0) > 0){
        
        expd0 <- expd0[expd0$SNP %in% outd0$SNP,]
        expd0 <- expd0[expd0$Source %in% names(table(expd0[expd0$SNP %in% TotalOut,]$Source)[table(expd0[expd0$SNP %in% TotalOut,]$Source) > 1]),]
        for (i in sort(unique(gene_human[gene_human$`Gene stable ID` %in% expd0$id.exposure & gene_human$`Gene name` != "",]$`Gene name`))) {
          
          expd1 <- expd0[expd0$id.exposure %in% gene_human[gene_human$`Gene name` == i,]$`Gene stable ID`,]
          outd1 <- outd0[outd0$SNP %in% expd1$SNP,]
          if(nrow(outd1) > 1){
            
            dat1 <- TwoSampleMR::harmonise_data(expd1, outd1)
            outTab1 <- merge(outTab1, dat1, all = T)
            if(nrow(outTab1) == 0){
              outTab1 <- merge(outTab1, dat1, all = T)
            }
            if(nrow(dat1) > 1){
              
              res <- TwoSampleMR::mr(dat1, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
              if(nrow(res) > 0){
                if(min(res$pval) < 0.05){
                  dir.create(paste0(Panel,"/",i))
                  
                  p <- mr_scatter_plot(res, dat1)
                  pdf(paste0(Panel,"/",i,"/PointPlot.pdf"),width = 8,height = 8)
                  print(p)
                  invisible(dev.off())
                  
                  
                  res_OR <- generate_odds_ratios(res)
                  write.table(res_OR,file = paste0(Panel,"/",i,"/odds_ratios.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
                  
                 
                  res_ht <- mr_heterogeneity(dat1)
                  write.table(res_ht,file = paste0(Panel,"/",i,"/heterogeneity.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
                  
               
                  dat_PRS <- dat1[rowSums(is.na(dat1[,c("beta.outcome", "beta.exposure", "se.outcome", "se.exposure")])) == 0,]
                  if(nrow(dat_PRS) > 3){
                    res_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome",
                                            SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_PRS,
                                            NbDistribution = 1000,  SignifThreshold = 0.05, seed = 123456)
                    b <- gsub("<","",res_presso$`MR-PRESSO results`$`Global Test`$Pvalue)
                    save(res_presso,file = paste0(Panel,"/",i,"/presso_",b,"_.rda"))
                  }
                  
                 
                  res_pl <- mr_pleiotropy_test(dat1)
                  write.table(res_pl,file = paste0(Panel,"/",i,"/pleiotropy.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
                  
                 
                  res_loo <- mr_leaveoneout(dat1)
                  if(max(table(res_loo$SNP)) > 1){
                    res_loo <- res_loo[order(res_loo$SNP, res_loo$p, decreasing = T),]
                    res_loo <- res_loo[!duplicated(res_loo$SNP, fromLast = T),]
                  }
                  p <- mr_leaveoneout_plot(res_loo)
                  pdf(paste0(Panel,"/",i,"/LeaveoneoutPlot.pdf"),width = 8,height = ((nrow(res_loo)/10) + 3))
                  print(p)
                  invisible(dev.off())
                  write.table(res_loo,file = paste0(Panel,"/",i,"/leaveoneout.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
                  
                  
                  res_single <- mr_singlesnp(dat1)
                  if(max(table(res_single$SNP)) > 1){
                    res_single <- res_single[order(res_single$SNP, res_single$p, decreasing = T),]
                    res_single <- res_single[!duplicated(res_single$SNP, fromLast = T),]
                  }
                  p <- mr_forest_plot(res_single)
                  pdf(paste0(Panel,"/",i,"/ForestPlot.pdf"),width = 6.5,height = ((nrow(res_single)/10) + 3))
                  print(p)
                  invisible(dev.off())
                  p <- mr_funnel_plot(res_single)
                  pdf(paste0(Panel,"/",i,"/FunnelPlot.pdf"),width = 8,height = 8)
                  print(p)
                  invisible(dev.off())
                  p <- mr_density_plot(res_single, res)
                  pdf(paste0(Panel,"/",i,"/DensityPlot.pdf"),width = 8,height = 5.5)
                  print(p)
                  invisible(dev.off())
                  write.table(res_single,file = paste0(Panel,"/",i,"/singlesnp.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
                }
                
                res$Symbol <- i
                outTab2 <- rbind(outTab2, res)
                if(dir.exists(paste0(Panel,"/",i))){
                  
                  if(res[res$method == "Inverse variance weighted",]$pval < 0.05){
                    
                    KeyGen <- c(KeyGen, i)
                    
                    if(SensiLoo == 1){
                      if(unique(res$nsnp) > 1){
                        if(length(setdiff(res_loo$SNP, "All")) == 0){
                          repeat{
                            a <- file.rename(paste0(Panel,"/",i),paste0(Panel,"/lo_",i))
                            if(a){
                              break
                            }
                          }
                          KeyGen <- setdiff(KeyGen, i)
                        }
                      }
                    }
                    
                    if(HeteroIVW == 1){
                      if(unique(res$nsnp) > 1){
                        if(res_ht[res_ht$method == "Inverse variance weighted",]$Q_pval < 0.05){
                          repeat{
                            a <- file.rename(paste0(Panel,"/",i),paste0(Panel,"/ht_",i))
                            if(a){
                              break
                            }
                          }
                          KeyGen <- setdiff(KeyGen, i)
                        }
                      }
                    }
                  }else{
                    repeat{
                      a <- file.rename(paste0(Panel,"/",i),paste0(Panel,"/pv_",i))
                      if(a){
                        break
                      }
                    }
                  }
                }
              }
            }
            Sys.sleep(1)
            rm(list = intersect(ls(),c("expd1", "outd1", "res", "dat1", "res_presso", "res_OR", "res_ht", "res_pl", "res_loo", "res_single", "p")))
          }
        }
      }
    }
  }
}
Sys.sleep(1)
write.table(outTab0,file = paste0("./",Panel,"/0.table.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
write.table(outTab1,file = paste0("./",Panel,"/1.table.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
write.table(outTab2,file = paste0("./",Panel,"/2.table.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

KeyGen <- unique(KeyGen);length(KeyGen)
tmp <- outTab2[outTab2$Symbol %in% KeyGen & outTab2$method == "Inverse variance weighted",]
if(nrow(tmp[which(tmp$pval < 0.05),]) > 1){
  cat("IVW pval < 0.05")
  KeyGen <- tmp[which(tmp$pval < 0.05),]$Symbol
}
length(KeyGen)
save(KeyGen,file = "./KeyGen.rda")




library(locuscomparer)
library(TwoSampleMR)
library(GagnonMR)
library(ieugwasr)
library(ggplot2)
library(coloc)
load("./KeyGen.rda")
Gwasid <- "finngen_R10_C3_OESOPHAGUS_EXALLC"
DirsOutFile <- "./SourceData/"
DirSMR <- "./smr-1.3.1-win-x86_64/"
DireQTL <- "./eQTLGen/"
GwasCaseNum <- 619;GwasControlNum <- 314193
ExpCutP <- 1e-5


(CandGen <- unique(KeyGen));length(CandGen)
gene_human <- read.table("./table.OR.txt",sep="\t",header=T,check.names=F);nrow(gene_human)
gene_human <- gene_human[,c("id.exposure", "Symbol")];colnames(gene_human) <- c("Gene stable ID", "Gene name");nrow(gene_human)

TotalOutCome <- read.table(paste0(DirsOutFile,Gwasid,".gz"), sep="\t", header=T, check.names=F, comment.char = "")
colnames(TotalOutCome)[1] <- c("Chr")
TotalOutCome <- TotalOutCome[TotalOutCome$rsids != "",c("rsids", "Chr", "pos", "alt", "ref", "pval", "beta", "sebeta", "af_alt")]
colnames(TotalOutCome) <- c("rsid", "Chr", "BP", "effectAllele", "otherAllele", "p", "beta", "se", "eaf")
head(TotalOutCome)
outTab <- list()
for (i in CandGen) {

  if(nrow(gene_human[gene_human$`Gene name` == i,]) != 0){
    dir.create(i)
    if(length(unique(gene_human[gene_human$`Gene name` == i,]$`Gene stable ID`)) > 1){
      tmp <- read.table(paste0("./",i,".txt"),sep="\t",header=T,check.names=F,quote = "")
      gene_human[gene_human$`Gene name` == i,]$`Gene stable ID` <- unique(tmp$id.exposure)
    }
    shell(paste0(DirSMR,"smr-1.3.1-win.exe --beqtl-summary ",
                 DireQTL,"eQTLGen --query ", ExpCutP," --gene ",
                 unique(gene_human[gene_human$`Gene name` == i,]$`Gene stable ID`)," --out ./",i,"/myquery"))
    Sys.sleep(1)
    if(file.exists(paste0(i,"/myquery.txt"))){
      
      eqtl <- read.table(paste0(i,"/myquery.txt"),sep="\t",header=T,check.names=F)
      eqtl$maf <- ifelse(eqtl$Freq > 0.5, (1-eqtl$Freq), eqtl$Freq)
      head(eqtl)
      write.table(eqtl,row.names = F,file=paste0(i,"/eqtl.txt"),sep="\t",quote=F)
      write.table(eqtl[,c("SNP", "p")],row.names = F,file=paste0(i,"/eqtl_in.txt"),sep="\t",quote=F)
      gwas <- TotalOutCome[TotalOutCome$rsid %in% eqtl$SNP,]
      gwas$maf <- ifelse(gwas$eaf > 0.5, (1-gwas$eaf), gwas$eaf)
      gwas$varbeta <- gwas$se^2
      gwas <- gwas[gwas$maf > 0.3,]
      gwas <- gwas[order(gwas$rsid, abs(gwas$beta), decreasing = T),]
      gwas <- gwas[!duplicated(gwas$rsid, fromLast = F),]
     
      write.table(gwas,row.names = F,file=paste0(i,"/gwas.txt"),sep="\t",quote=F)
      write.table(gwas[,c("rsid","p")],row.names = F,file=paste0(i,"/gwas_in.txt"),sep="\t",quote=F)
      
      if(nrow(gwas) > 0){
        
        p <- locuscompare(in_fn1=paste0("./",i,"/gwas_in.txt"), in_fn2=paste0("./",i,"/eqtl_in.txt"), 
                          title1="GWAS", title2="eQTL", marker_col1= "rsid", 
                          pval_col1="p", marker_col2="SNP", pval_col2="p")
        pdf(file=paste0("./",i,"/locuscompare.pdf"),width=7,height=6)
        print(p)
        dev.off()
        
        colnames(eqtl)[1] <- "rsid"
        if(min(eqtl$p) <= 0){
          eqtl[which(eqtl$p <= 0),]$p <- sort(unique(eqtl$p[!(eqtl$p <= 0)]), decreasing = F)[1] / 10
        }
        if(min(gwas$p) <= 0){
          gwas[which(gwas$p > 1),]$p <- sort(unique(gwas$p[!(gwas$p <= 0)]), decreasing = F)[1] / 10
        }
        inputDat <- merge(eqtl, gwas, by = "rsid", all = F, suffixes = c("_eqtl", "_gwas"))
       
        Res <- coloc.abf(dataset1 = list(pvalues = inputDat$p_eqtl, type = "quant", N = 31864, snp = inputDat$rsid, MAF=inputDat$maf_eqtl),
                         dataset2 = list(pvalues = inputDat$p_gwas, type = "cc", s=GwasCaseNum/(GwasCaseNum+GwasControlNum), N = (GwasCaseNum+GwasControlNum), snp = inputDat$rsid, MAF=inputDat$maf_gwas))
        
        tmp <- Res$results
        write.table(tmp[order(tmp$SNP.PP.H4, decreasing = T),],row.names = F,col.names = T,file=paste0(i,"/Colocalization.",GwasCaseNum,".",GwasControlNum,".txt"),sep="\t",quote=F)
        write.table(as.data.frame(Res$summary),row.names = T,col.names = F,file=paste0(i,"/Colocalization.summary.txt"),sep="\t",quote=F)
        write.table(as.data.frame(Res$priors),row.names = T,col.names = F,file=paste0(i,"/Colocalization.priors.txt"),sep="\t",quote=F)
        outTab[[i]] <- Res
      }
    }
  }
}


CutNum <- 0.75
KeyGen <- NULL
ResTab <- NULL
for (i in names(outTab)) {
  
  a <- sort(outTab[[i]]$results$SNP.PP.H4, decreasing = T)
  if(a[1] > CutNum){
    
    file.copy(paste0("./",i,"/locuscompare.pdf"), paste0("./",i,".pdf"), overwrite = T)
    tmp <- outTab[[i]]$results
    tmp <- tmp[which(tmp$SNP.PP.H4 > CutNum),]
    tmp$Symbol <- i
    ResTab <- rbind(ResTab, tmp)
    KeyGen <- c(KeyGen, i)
  }
}
length(KeyGen)
save(KeyGen,file = "./KeyGen.rda")
write.table(ResTab,row.names = F,col.names = T,file=paste0("./Colocalization.result.xls"),sep="\t",quote=F)



library(cowplot)
library(clusterProfiler)
library(enrichplot)
library(plyr)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) r





selectedGeneID <- c("IGLV2-11")                             


mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

gsym.fc <- read.table("input.txt", header = T)
dim(gsym.fc)

gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")


gsym.fc.id <- merge(gsym.fc, gsym.id, by="SYMBOL", all=F)

gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T),]

id.fc <- gsym.fc.id.sorted$logFC
names(id.fc) <- gsym.fc.id.sorted$ENTREZID

kk <- gseKEGG(id.fc, organism = "hsa")

kk.gsym <- setReadable(kk, 'org.Hs.eg.db', 
                       'ENTREZID')

sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T),]



write.csv(sortkk,"gsea_output.csv", quote = F, row.names = F)               



geneSetID <- c("hsa04514", "hsa04062", "hsa04970")


for (i in geneSetID) {
  gseaplot(kk, i)
  myGeneList <- enrichplot:::gsInfo(kk, i)
  row.names(myGeneList) <- gsym.fc$gsym
  myGeneList$id <- gsym.fc$ENTREZID 
  write.csv(myGeneList, paste0("gsea_genelist_", i, "_group1.csv"))
}

x <- kk
geneList <- position <- NULL 


gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
gsdata$gsym <- rep(gsym.fc.id.sorted$SYMBOL,3)                                 


p.res <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
  geom_line(aes_(y = ~runningScore, color= ~Description), size=1) +
  scale_color_manual(values = mycol) +
  

  geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) +
  ylab("Enrichment\n Score") +
  
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  
  theme(legend.position = "top", legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent")) +
  
  theme(axis.text.y=element_text(size = 12, face = "bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))


rel_heights <- c(1.5, .5, 1.5) 

i <- 0
for (term in unique(gsdata$Description)) {
  idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
  gsdata[idx, "ymin"] <- i
  gsdata[idx, "ymax"] <- i + 1
  i <- i + 1
}


p2 <- ggplot(gsdata, aes_(x = ~x)) +
  geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
  xlab(NULL) + ylab(NULL) + 
  scale_color_manual(values = mycol) + 
  
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  
  theme(legend.position = "none",
        plot.margin = margin(t=-.1, b=0,unit="cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank()) +
  scale_y_continuous(expand=c(0,0))


df2 <- p.res$data
df2$y <- p.res$data$geneList[df2$x]
df2$gsym <- p.res$data$gsym[df2$x]


selectgenes <- data.frame(gsym = selectedGeneID)
selectgenes <- merge(selectgenes, df2, by = "gsym")
selectgenes <- selectgenes[selectgenes$position == 1,]
head(selectgenes)


p.pos <- ggplot(selectgenes, aes(x, y, fill = Description, color = Description, label = gsym)) + 
  geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0), 
               color = "grey") +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = mycol, guide=FALSE) + 
  scale_color_manual(values = mycol, guide=FALSE) + 
  
  geom_hline(yintercept = 0, lty = 2, lwd = 0.2) +
  ylab("Ranked list\n metric") +
  xlab("Rank in ordered dataset") +
  
  theme_bw() +
  theme(axis.text.y=element_text(size = 12, face = "bold"),
        panel.grid = element_blank()) +
  
geom_text_repel(data = selectgenes, 
                show.legend = FALSE, 
                direction = "x", 
                ylim = c(2, NA),
                angle = 90, 
                size = 2.5, box.padding = unit(0.35, "lines"), 
                point.padding = unit(0.3, "lines")) +
  theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))



plotlist <- list(p.res, p2, p.pos)
n <- length(plotlist)
plotlist[[n]] <- plotlist[[n]] +
  theme(axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.text.x = element_text(size = 12, face = "bold"))

plot_grid(plotlist = plotlist, ncol = 1, align="v", rel_heights = rel_heights)

ggsave("GSEA_multi_pathways.pdf", width=6, height=5)


library(clusterProfiler)
library(GOplot)
library(tidyverse)
library(data.table)
library(ggraph)
library(tidygraph)
library(dplyr)


sortkk <- kk.gsym[kk.gsym@result$Description %like% "Cell adhesion molecules" |
                    kk.gsym@result$Description %like% "Chemokine signaling pathway" |
                    kk.gsym@result$Description %like% "Salivary secretion",]  

go <- data.frame(Category = "KEGG",
                 ID = sortkk$ID,
                 Term = sortkk$Description, 
                 Genes = gsub("/", ", ", sortkk$core_enrichment), 
                 adj_pval = sortkk$p.adjust)


genelist <- data.frame(ID = gsym.fc.id$SYMBOL, logFC = gsym.fc.id$logFC)
circ <- circle_dat(go, genelist)
head(circ)



write.csv(circ[,c(3,5,6)],"very_easy_input.csv", quote = F, row.names = F)


df <- read.csv("very_easy_input.csv")
head(df)
source(file = "gather_graph_node.R")
source(file = "gather_graph_edge.R")
nodes <- gather_graph_node(df, index = c("term", "genes"), value = "logFC", root="all")
edges <- gather_graph_edge(df, index = c("term", "genes"), root = "all")
nodes <- nodes %>% mutate_at(c("node.level","node.branch"),as.character)
head(nodes, 10)




graph <- tbl_graph(nodes, edges)


gc <- ggraph(graph, layout = 'dendrogram', circular = TRUE) + 
  geom_edge_diagonal(aes(color = node1.node.branch,
                         filter=node1.node.level!="all"), 
                     alpha = 1/3,edge_width=1) + 
  geom_node_point(aes(size = node.size, 
                      color = node.branch,
                      filter=node.level!="all"), alpha = 1/3) + 
  scale_size(range = c(0.5,8)) + 
  theme(legend.position = "none") + 
  


  scale_edge_color_brewer(palette = "Set1") + 
  scale_color_brewer(palette = "Set1") +
  
  geom_node_text(
    aes(
      x = 1.048 * x, 
      y = 1.048 * y, ?
      label = node.short_name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = leaf,
      color = node.branch
    ),
    size = 3, hjust = 'outward') +
  

  geom_node_text(
    aes(label=node.short_name,
        filter = !leaf & (node.level != "all"),
        color = node.branch),
    fontface="bold",
    size=3,
    family="sans"
  ) + 
  theme(panel.background = element_rect(fill = NA)) +
  coord_cartesian(xlim=c(-1.3,1.3),ylim = c(-1.3,1.3)) 

gc


ggsave("ccgraph_color.pdf", width = 14, height = 14)



gc1 <- ggraph(graph, layout = 'dendrogram', circular = TRUE) + 

  geom_edge_diagonal(aes(color = node1.node.branch,
                         filter=node1.node.level!="all"), 
                     alpha = 0.5, 
                     edge_width=2.5) + 
  scale_edge_color_manual(values = c("#61C3ED","red","purple","darkgreen")) + 
  

  geom_node_point(aes(size = node.size,
                      filter=node.level!="all"), 
                  color = "#61C3ED") + 
  scale_size(range = c(0.5,3)) + 
  theme(legend.position = "none") + 
  

  geom_node_text(
    aes(
      x = 1.05 * x, 
      y = 1.05 * y, 
      label = node.short_name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = leaf
    ),
    color="black", 
    size = 3, hjust = 'outward') +
  

  geom_node_text(
    aes(label=node.short_name,
        filter = !leaf & (node.level != "all")
    ),
    color="black", 
    fontface="bold",
    size=6,
    family="sans"
  ) + 
  theme(panel.background = element_rect(fill = NA)) + 
  coord_cartesian(xlim=c(-1.3,1.3),ylim = c(-1.3,1.3))

gc1


ggsave("2.Ccgraph.pdf",width = 14,height = 14)





source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "symbol.txt", perm=100, QN=TRUE)



library(RcisTarget)
library(visNetwork)
library(reshape2)
library(stringr)
library(ggplot2)
library(ggpubr)
library(AUCell)
library(limma)
library(e1071)
library(Hmisc)
library(dplyr)
library(aplot)
library(doRNG)
library(doMC)
library(DT)

GseGroup <- read.table("group.txt",header=T,sep="\t",check.names=F)

compare_method = "anova"

ciberRes=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)

ciber <- ciberRes[,setdiff(colnames(ciberRes), c("P-value", "Correlation", "RMSE"))]
ciber$Acc <- rownames(ciber)
if(sort(unique(GseGroup$Tissue))[1] == "Normal"){
  datGroup <- GseGroup[order(GseGroup$Tissue),]
}else{
  datGroup <- GseGroup[order(GseGroup$Tissue, decreasing = T),]
}
datGroup$Acc <- factor(datGroup$Acc, levels = datGroup$Acc)
data <- dplyr::inner_join(datGroup,ciber,by="Acc")
data_p <- melt(data, id.vars = colnames(datGroup))
head(data_p)
data_p$Acc <- factor(data_p$Acc, levels = datGroup$Acc)

datGroup$p <- "Group"
p2 <- ggplot(datGroup,aes(Acc, p, fill=Tissue)) +
  geom_tile() +
  scale_fill_manual(values=c("#1CFA04", "#C705FF")) +
  scale_y_discrete(position = "right") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(text = element_text(size = 15)) +
  theme(axis.text.x = element_blank())+
  labs(fill = "Group")


sample_names <- datGroup$Acc 

p1 <- ggplot(data_p, aes(x = Acc, y = value, fill=variable)) +
  geom_bar(stat="identity", position = "fill", width = 0.5) +
  geom_col(position = 'fill', width = 0.6) +
  guides(fill=guide_legend(title = NULL)) +
  ylab("Relative Percent") + xlab("") +
  theme_bw() +
  theme(axis.ticks.length=unit(0.5,'cm')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
  theme(axis.text=element_text(size = 15)) +
  scale_x_discrete(labels = sample_names) + 
  theme(axis.text.y=element_text(colour = "black", vjust=0,size = 5)) +
  theme(axis.title =element_text(size = 20)) +
  theme(text = element_text(size = 15)) +
  scale_y_continuous(expand=c(0,0.05))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 10))

p <- p1 %>% insert_top(p2, height = 0.05)
ggsave(p, filename = "Immune infiltration.pdf", width = (nrow(datGroup)/55)+15, height = 6,limitsize = FALSE)



rt=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)

library(corrplot)
pdf("corHeatmap.pdf",height=15,width=15)              
corrplot(corr=cor(rt),
         method = "color",
         order = "hclust",
         tl.col="black",
         addCoef.col = "black",
         number.cex = 1,
         col=colorRampPalette(c("blue", "white", "red"))(50),
)
dev.off()







library(ggpubr)

pFilter=0.99

rt=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)    

data=rt



Type=read.table("group.txt",sep="\t",check.names=F,row.names=1,header=F)
Type=Type[row.names(data),]
colnames(Type)=c("cluster","Subtype")

outTab=data.frame()
data=cbind(data,Type)


for(i in colnames(data[,1:(ncol(data)-2)]))
{
  rt1=data[,c(i,"Subtype")]
  colnames(rt1)=c("expression","Subtype")
  ksTest<-kruskal.test(expression ~ Subtype, data = rt1)
  pValue=ksTest$p.value
  if(pValue<pFilter){
    outTab=rbind(outTab,cbind(rt1,gene=i))
    print(pValue)
  }
}
write.table(outTab,file="data.txt",sep="\t",row.names=F,quote=F)

data=read.table("data.txt",sep="\t",header=T,check.names=F)       
data$Subtype=factor(data$Subtype, levels=c("Normal","Tumor"))
p=ggboxplot(data, x="gene", y="expression",color = "grey",fill = "Subtype",
            ylab="Expression",
            xlab="",
            palette =c("skyblue","pink") )
p=p+rotate_x_text(45)
p
pdf(file="boxplot.pdf",width=12,height=4)                          
p+stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif",method="kruskal.test")
dev.off()

library(dplyr)
library(ggplot2)
data %>% 
  filter(Subtype %in% c("Normal","Tumor")) %>% 
  ggplot(aes(x= gene, y= expression, fill = Subtype, color = Subtype))+
  geom_boxplot(alpha=0.3)+
  scale_fill_manual(name= "Subtype", values = c("deepskyblue", "hotpink"))+
  scale_color_manual(name = "Subtype", values = c("dodgerblue", "plum3"))+
  theme_bw()+labs(x="", y="Expression")+
  theme(axis.text.x = element_text( vjust = 1,size = 12, hjust = 1,colour = "black"),legend.position="top")+
  rotate_x_text(45)+stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif",method="wilcox")






library(ggplot2)
library(ggpubr)
library(SimDesign)
library(cowplot)
library(dplyr)
library(GSVA)
library(limma)
library(stringr)

jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")


expr <- read.table("symbol.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
expr<-log2(expr+1)
gene <- read.table("genelist.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = F)
ciber <- read.table("CIBERSORT-Results.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

for (i in gene$V1) {
  message(paste0("analysis of ",i," starts..."))
  subexpr <- as.numeric(expr[i,])
  names(subexpr) <- colnames(expr)
  
  lsam <- names(subexpr[subexpr < median(subexpr)])
  hsam <- names(subexpr[subexpr >= median(subexpr)])
  
 
  dat <- as.numeric(expr[i,]); names(dat) <- colnames(expr)
  comsam <- intersect(names(dat), rownames(ciber))
  tmp1 <- dat[comsam]
  tmp2 <- ciber[comsam,]
  
  var <- colnames(ciber)
  data <- data.frame(var)
  for (j in 1:length(var)){
    test <- cor.test(as.numeric(tmp2[,j]),tmp1,method = "spearman") 
    data[j,2] <- test$estimate                                            
    data[j,3] <- test$p.value
  }
  names(data) <- c("symbol","correlation","pvalue")
  data <- as.data.frame(na.omit(data))
  data %>% 
    filter(pvalue <0.05) %>%  
    ggplot(aes(correlation,forcats::fct_reorder(symbol,correlation))) +
    geom_segment(aes(xend=0,yend=symbol)) +
    geom_point(aes(col=pvalue,size=abs(correlation))) +
    scale_colour_gradientn(colours=c("#7fc97f","#984ea3")) +
    scale_size_continuous(range =c(2,8))  +
    theme_minimal() +
    ylab(NULL)
  ggsave(paste0("correlation between cibersort and expression of ", i,".pdf"),width = 8,height = 6)
}


a5 <- read.table("nomogram.txt",header = T,sep = "\t", quote = "",fill = T)
head(a5)
dim(a5)
a5[,"LRP1"]=log2(a5[,"LRP1"]+1)
a5[,"RORC"]=log2(a5[,"RORC"]+1)
a5[,"IGLV2.11"]=log2(a5[,"IGLV2.11"]+1)
a5 <- a5[,-1]
dim(a5)
head(a5)
a5$died <- a5$fustat==1
a5$futime <- as.numeric(a5$futime)
a5$fustat <- as.numeric(a5$fustat)
a5$Gender <- as.numeric(a5$Gender)
a5$Stage <- as.numeric(a5$Stage)
a5$T <- as.numeric(a5$T)
a5$M <- as.numeric(a5$M)
a5$N <- as.numeric(a5$N)
a5$LRP1 <- as.numeric(a5$LRP1)
a5$RORC <- as.numeric(a5$RORC)
a5$IGLV2.11 <- as.numeric(a5$IGLV2.11)
data6 <- a5
head(data6)
library(rms)




dd<-datadist(a5)
options(datadist="dd")
options(na.action="na.delete")
summary(data6$futime)

coxpbc<-cph(formula = Surv(futime,died) ~ Gender +  Stage + T + M + N+ LRP1+ RORC+ IGLV2.11,data=a5,x=T,y=T,surv =T,na.action=na.delete)

print(coxpbc)

surv<-Survival(coxpbc) 
surv3<-function(x) surv(365,x)                    
surv4<-function(x) surv(1095,x)                   

x<-nomogram(coxpbc,fun = list(surv3,surv4),lp=T,
            funlabel = c('1-year survival Probability','3-year survival Probability'),
            maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))

pdf("nomogram_classical.pdf",width = 12,height = 10)
plot(x, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()


f5<-cph(formula = Surv(futime,died) ~ Gender +  Stage + T + M + N+ LRP1+ RORC+ IGLV2.11,data=a5,x=T,y=T,surv = T,na.action=na.delete,time.inc = 365) 


cal5<-calibrate(f5, cmethod="KM", method="boot",u=365,m=70,B=1000)

pdf("calibration_1y.pdf",width = 8,height = 8)
plot(cal5,
     lwd = 2,
     lty = 1,
     errbar.col = c("#2166AC"),
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal5[,c('mean.predicted',"KM")], 
      type = 'b',
      lwd = 2, 
      pch = 16, 
      col = c("#2166AC")) 
mtext("")
box(lwd = 1) 
abline(0,1,lty = 3, 
       lwd = 2, 
       col = c("#224444")
) 
dev.off()



f8<-cph(formula = Surv(futime,died) ~ Gender +  Stage + T + M + N+ LRP1+ RORC+ IGLV2.11,data=a5,x=T,y=T,surv = T,na.action=na.delete,time.inc = 1095) 
cal8<-calibrate(f8, cmethod="KM", method="boot",u=1095,m=70,B=1000)

pdf("calibration_3y.pdf",width = 8,height = 8)
plot(cal8,
     lwd = 2,
     lty = 1,
     errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#B2182B"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal8[,c('mean.predicted',"KM")],
      type= 'b',
      lwd = 2,
      col = c("#B2182B"),
      pch = 16)
mtext("")
box(lwd = 1)
abline(0,1,lty= 3,
       lwd = 2,
       col =c("#224444"))
dev.off()


pdf("calibration_compare.pdf",width = 8,height = 8)
plot(cal5,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", 
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal8,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal8[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", 
       legend = c("1-year","3-year"), 
       col =c("#2166AC","#B2182B"), 
       lwd = 2,
       cex = 1.2,
       bty = "n")
dev.off()




library(timeROC)

ROC.DSST<-timeROC(T=data6$futime,
                  delta=data6$fustat,
                  marker=data6$LRP1,
                  other_markers=as.matrix(data6[,c( "Gender", "Stage", "T", "M", "N", "RORC", "IGLV2.11")]),
                  weighting="cox",cause=1,
                  times=c(365,1095,1825),ROC=TRUE)

cols <- c("#FD0207", "#0122D4", "#26D401")
pdf("nomogram_ROC.pdf",width = 6,height = 6)
plot(ROC.DSST,time=365, col=cols[1], title = "")
plot(ROC.DSST,time=1095,add=TRUE,col=cols[2])
plot(ROC.DSST,time=1825,add=TRUE,col=cols[3])
legend("bottomright", legend=c(paste0("1-year AUC:",round(ROC.DSST$AUC["t=365"], 4)),
                               paste0("3-year AUC:",round(ROC.DSST$AUC["t=1095"], 4)),
                               paste0("5-year AUC:",round(ROC.DSST$AUC["t=1825"], 4))),
       col=cols, lwd=2)
dev.off()



library(rmda)
library(ggDCA)
library(ggplot2)
library(rms)
library(caret)
data(dcaData)
head(dcaData)

set.seed(123)

rt=read.table("nomogram.txt",header=T,sep="\t",check.names=F,row.names=1)
lrm1 <- lrm(fustat ~ LRP1, rt)
lrm2 <- lrm(fustat ~ Gender, rt)
lrm3 <- lrm(fustat ~ Stage, rt)
lrm4 <- lrm(fustat ~ T, rt)
lrm5 <- lrm(fustat ~ M, rt)
lrm6 <- lrm(fustat ~ N, rt)
lrm7 <- lrm(fustat ~ RORC, rt)
lrm8 <- lrm(fustat ~ IGLV2.11, rt)


dca_lrm <- dca(lrm1, lrm2, lrm3, lrm4,lrm5,lrm6, lrm7,lrm8,model.names = c("LRP1", "RORC","IGLV2.11", "Gender", "Stage", "T", "M", "N"
))

ggplot(dca_lrm)



library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)
library(data.table)

logFCfilter=1              
adjPvalFilter=0.05         


pbmc=readRDS("CellAnnoted.cellType_1.rds")
Idents(pbmc)=pbmc$cellType_1


showGenes=c("LRP1","RORC")


pdf(file="03.markerViolin.pdf",width=10,height=10)
VlnPlot(object = pbmc, features = showGenes)
dev.off()


pdf(file="03.markerScatter.pdf",width=8,height=6)
FeaturePlot(object = pbmc, features = showGenes, cols = c("grey", "blue"),raster=FALSE)
dev.off()


p <- DotPlot(pbmc,features = c( 'LRP1', 'RORC'),cols = "RdYlBu"
)+ scale_size_continuous(range = c(0, 10))+ theme(
  panel.border = element_rect(colour = 'black'),
  axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ),
  legend.position = "top",
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.8, "cm"),
);p




library(SeuratData)
library(patchwork)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(Seurat)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)



getScatterplot <- function(object, gene1, gene2, cor.method = "pearson",
                           jitter.num = 0.15, pos = TRUE){
  if (!gene1 %in% rownames(object)) {
    print("gene1 was not found")
    if (!gene2 %in% rownames(object)) {
      print("gene2 was not found")
    }
  }else{
    exp.mat <- GetAssayData(object = object, assay = "RNA") %>% .[c(gene1,gene2),] %>% 
      as.matrix() %>% t() %>% as.data.frame()
    if (pos) {
      if(nrow(exp.mat[which(exp.mat[,1] > 0 & exp.mat[,2] > 0),]) > (nrow(exp.mat)*0.01)){
        exp.mat <- exp.mat[which(exp.mat[,1] > 0 & exp.mat[,2] > 0),]
      }else{
        exp.mat <- exp.mat[which(exp.mat[,1] > 0 | exp.mat[,2] > 0),]
      }
    }
    colnames(exp.mat) <- c("Var1", "Var2")
    plots <- ggplot(data=exp.mat, mapping = aes_string(x = "Var1", y = "Var2")) + 
      geom_smooth(method = 'lm', se = T, color='red', size=1) +
      stat_cor(method = cor.method)+ labs(x=gene1, y=gene2) +
      geom_jitter(width = jitter.num, height = jitter.num, color = "black", size = 1, alpha=1)+
      theme_bw()+
      theme(panel.grid=element_blank(),
            legend.text=element_text(colour= 'black',size=10),
            axis.text= element_text(colour= 'black',size=10),
            axis.line= element_line(colour= 'black'),
            panel.border = element_rect(size = 1, linetype = "solid", colour = "black"),
            panel.background=element_rect(fill="white"))
    return(plots)
  }
}




showGenes <- c("LRP1","RORC")
geneCard <- intersect(c("BRCA2","BRCA1","TP53","ATM","CDH1"),rownames(pbmc))


for (j in geneCard) {
  for (i in showGenes) {
    
    p1 <- FeaturePlot(pbmc, features = c(i, j),
                      blend = TRUE, cols = c("gray80","red", "green"), 
                      pt.size = 0.5, raster = F) +  
      theme(aspect.ratio = 1)
    p2 <- getScatterplot(pbmc, gene1 = j, gene2 = i, 
                         jitter.num = 0.15, pos = TRUE) +
      theme(aspect.ratio = 1)
    p <- CombinePlots(plots = list(p1, p2), ncol = 2, rel_widths = c(4, 1))
    pdf(file=paste0(j," ~ ",i,".pdf"),width=15,height=4)
    print(p)
    dev.off()
  }
}


library(clusterProfiler)
library(magrittr)
library(ggplot2)
library(cowplot)
library(Seurat)
library(AUCell)
library(limma)
library(dplyr)
library(plyr)
library(GSVA)
ShowGen <- c("LRP1","RORC")


pbmc <- readRDS("CellAnnoted.cellType_1.rds")
pbmc[["celltype"]] <- as.data.frame(Idents(pbmc))


geneset <- read.gmt("h.all.v7.5.1.symbols.gmt")
geneset <- split(geneset$gene, geneset$term)
genesetInfo <- read.delim("GenesetInfo.txt", sep = ",")
genesetInfo <- subset(genesetInfo, Classification != "")
levels = c("Immune", "Metabolism", "Signaling", "Proliferation")
genesetInfo$Classification <- factor(genesetInfo$Classification, levels)
genesetInfo <- arrange(genesetInfo, genesetInfo$Classification, genesetInfo$geneset)
genesetInfo$geneset <- factor(genesetInfo$geneset, levels = genesetInfo$geneset)
mat <- GetAssayData(object = pbmc, assay = "RNA", slot = "data")
cells_rankings <- AUCell_buildRankings(mat, nCores=1, plotStats = F)
score <- AUCell_calcAUC(geneset, cells_rankings, nCores = 1, 
                        aucMaxRank = nrow(cells_rankings)*0.05)
score <- getAUC(score)

compare="Hexp-Lexp"
adjust.method = "bonferroni"
DP.list <- lapply(ShowGen, function(celltype){
  
  
  
  DatGroup <- FetchData(pbmc, vars = c("Group", celltype), slot = "data")
  group = setNames(object = ifelse(DatGroup[,1] > median(DatGroup[,1]), "Hexp", "Lexp"),
                   nm = rownames(DatGroup))
  if(length(unique(group)) > 1){
    design = model.matrix(~ 0 + factor(group))
    colnames(design) = levels(factor(group))
    rownames(design) = names(group)
    
    contrast.matrix = makeContrasts(compare, levels = design) # should be Test-Control
    fit = lmFit(score[, names(group)], design)
    fit2 = contrasts.fit(fit, contrast.matrix)
    fit2 = eBayes(fit2) 
    DPs = topTable(fit2, coef=1, n=Inf, adjust.method = adjust.method)
    DPs$Celltypes = celltype
    DPs$Pathway = rownames(DPs)
    return(DPs)
  }})
DPs <- do.call(rbind, DP.list)
write.table(DPs, file = "output_hallmark.txt", sep = "\t", row.names = F, col.names = T, quote = F)
plot.data <- DPs
plot.data <- subset(plot.data, Pathway %in% genesetInfo$geneset)
plot.data$Celltypes <- factor(plot.data$Celltypes)
plot.data$Pathway <- factor(plot.data$Pathway, levels = genesetInfo$geneset)
plot.data$FDR<- cut(plot.data$adj.P.Val, breaks = c(0, 1e-125, 1e-75, 1e-25, 1),
                    include.lowest = T)
plot.data$FDR <- factor(as.character(plot.data$FDR),
                        levels = rev(levels(plot.data$FDR)))
levels(plot.data$Pathway) <- tolower(gsub("HALLMARK_", "", levels(plot.data$Pathway)))

color = c("#4682B4", "#FFFFFF", "#CD2626") 
class.color = c("Immune" = "#D58986", "Metabolism" = "#80554C",
                "Signaling" = "#71AC7A", "Proliferation" = "#E8D4B4") 
p1 <- ggplot(plot.data, aes(x = Pathway, y = Celltypes, color = logFC, size = FDR)) +
  geom_point() +
  scale_color_gradient2(low = color[1], mid = color[2], high = color[3]) +
  geom_hline(yintercept = seq(min(as.numeric(plot.data$Celltypes))-0.5,
                              max(as.numeric(plot.data$Celltypes))+0.5),
             color = "grey80") +
  geom_vline(xintercept = seq(min(as.numeric(plot.data$Pathway))-0.5,
                              max(as.numeric(plot.data$Pathway))+0.5),
             color = "grey80") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_blank(), legend.position = "top")

p2 <- ggplot(genesetInfo, aes(x = geneset, y = 1, fill = Classification)) +
  geom_tile() +
  theme_classic() +
  scale_fill_manual(values = class.color) +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank(), legend.position = "bottom")
p <- plot_grid(p1, p2, ncol = 1, align = 'v', rel_heights = c(10, 2)) # 拼接图像
ggsave(p, filename = paste0("./hallmark.pdf"), width = 12, height = 9)

