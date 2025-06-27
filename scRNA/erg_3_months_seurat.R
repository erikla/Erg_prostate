### RNA processing with Seurat

library(Seurat)
library(SingleCellExperiment)

library(ggplot2)
library(gridExtra)
library(viridis)
library(data.table)
library(magrittr)

# https://github.com/immunogenomics/presto
library(presto)

base.dir <- '/data/projects/sawyers/erg/analysis/scRNA/3months/seurat/'
robj.dir <- file.path(base.dir,'robj')
plot.dir <- file.path(base.dir,'plot')
data.dir <- file.path(base.dir,'data')
spreadsheets.dir <- file.path(base.dir,'spreadsheets')
scrna.object <- file.path(robj.dir,'erg_seurat.RDS')
epi.scrna.object <- file.path(robj.dir,'erg_epi_seurat.RDS')
epi_with_SV.scrna.object <- file.path(robj.dir,'erg_epi_with_SV_seurat.RDS')

source("~/lib/scripts/R/Colors.R")

## All murine specific
lum_mature_markers <- c("Cd24a","Foxa1","Gata3","Krt8","Krt18","Krt19","Nkx3-1")
basal_markers <- c("Itga6","Krt5","Krt14","Krt15","Krt17","Trp63")
club_markers <- c("Lcn2","Wfdc2","Tacstd2","Krt4","Mmp7")
hillock_markers <- c("Krt13","S100a16","S100a14","Cstb","Lypd3")

epithelial_markers <- list('basal'= c('Krt5','Sox4'),
                           'intermediate luminal' = c('Ppp1r1b','Clu','Tacstd2','Krt4'),
                           'differentiated luminal' = c('Nkx3-1','Sbp'))
markers_list <- list(
  'Prostate luminal mature'=lum_mature_markers,
  'Prostate basal'=basal_markers,
  'Prostate L1 secretory luminal'=c('Dpp4','Ly6a','Prom1'),  # DNE Cd133=Prom1 Cd26=Dpp4, Sca1=Ly6a, Trop2=Tacstd2
  'Prostate L2 stem'=c('Ly6a','Psca','Tacstd2'),
  'Prostate L3'=c('Foxi1'),
  'Prostate proliferation'=c('Mki67',"Mcm7","Pcna",'Rest',"Top2a"),
  'AR dependent'=c('Ar','Fkbp5','Hpgd','Fgfr2','Lgr4','Egf','Egfr'),
  # https://www.nature.com/articles/s41419-020-2671-1 (KLF5 downregulation regulates IGF1 promoting tumor invasion)
  # https://www.cell.com/cell-reports/pdf/S2211-1247(18)31841-2.pdf (Klf4 expression stem cell like, loss induces aggressive cancer features)
  'PCa'=c("Igf1","Klf5","Klf4")
)

violin_markers <- list(
  'luminal_markers'=c("Cd24a","Foxa1","Krt8","Krt18","Krt19"),
  'basal_markers'=c("Krt5","Krt14","Krt15","Krt17","Trp63"),
  'club_markers'=c("Lcn2","Wfdc2","Tacstd2","Krt4","Mmp7"),
  'hillock_markers'=c("Krt13","S100a16","S100a14","Cstb","Lypd3")
)

chromatin_modifiers_list <- list(
  'H3K4 methyl' = c('Kmt2a','Kmt2b','Kmt2c','Kmt2d', 'Setd1a','Setd1b'),
  'H3K methyl' = c('Dot1l','Nsd1','Nsd2','Nsd3','Setd2'),
  #'Histone acetyl'
  'p300/CBP family' = c('Ep300', 'Crebbp'),
  'GCN5 family' = c('Nr6a1', 'Kat2b'),
  'MYST family' = c('Kat5', 'Kat6a', 'Kat6b', 'Kat7', 'Kat8')
)

#### Gene signature sets ####
Nelson_response_ANDR_UP_mouse <- c("Slc26a2","Orm1","Orm3","Orm2","Gsr","Aldh1a3","Gucy1a1","Hpgd","Actn1","Azgp1","Cpd","B4galt1","Klk1b8","Klk1b1","Klk1b9","Klk1b11","Klk1b26","Klk1b27","Klk1b21","Klk1b22","Klk1b16","Klk1b24","Klk1b3","Klk1b4","Klk1b5","Klk1","Hmgcs1","Id2","Itgav","Krt8","Krt19","Lifr","Sat1","Snap25","Sord","Uap1","Srp19","Ugdh","Vapa","Cdc14b","Plpp1","Inpp4b","B2m","Fkbp5","Klf4","Acsl3","Idi1","Rab4a","Sms","Homer2","Klk4","Scd2","Scd4","Tpd52","Akap12","H1f0","Maf","Selenop","Lman1","Sgk1","Tmprss2","Rbm10","Abcc4","Tsc22d1","Ndrg1","Nkx3-1","Phyh","Unc13b","Appbp2","Pdlim5","Camkk2","Iqgap2","Adamts1","Adrm1","Gm9774","Abhd2","Ptpn21","Pakap","Dctn3","Ell2","Dnajb9","Spdef","Stk39","Tnfaip8","Herc3","Dhcr24","Sec24d","Pgm3","Pias1","Hsd17b14","Hes6","Slc38a2","Pmepa1","Dbi","Elovl5","Zbtb10","Steap4","1700025G04Rik")
Nelson_response_no_NEG_ANDR_UP_mouse <- c("1700025G04Rik","Abcc4","Abhd2","Acsl3","Actn1","Adamts1","Adrm1","Akap12","Aldh1a3","Appbp2","Azgp1","B2m","B4galt1","Camkk2","Cdc14b","Cpd","Dbi","Dctn3","Dhcr24","Dnajb9","Ell2","Elovl5","Fkbp5","Gm9774","Gsr","Gucy1a1","H1f0","Herc3","Hes6","Hmgcs1","Homer2","Hpgd","Hsd17b14","Id2","Idi1","Inpp4b","Iqgap2","Itgav","Klf4","Klk1","Klk1b1","Klk1b11","Klk1b16","Klk1b21","Klk1b22","Klk1b24","Klk1b26","Klk1b27","Klk1b3","Klk1b4","Klk1b5","Klk1b8","Klk1b9","Klk4","Krt19","Krt8","Lifr","Lman1","Maf","Ndrg1","Nkx3-1","Orm1","Orm2","Orm3","Pakap","Pdlim5","Pgm3","Phyh","Pias1","Plpp1","Pmepa1","Ptpn21","Rab4a","Rbm10","Sat1","Scd2","Scd4","Sec24d","Selenop","Sgk1","Slc26a2","Slc38a2","Sms","Snap25","Sord","Spdef","Srp19","Steap4","Stk39","Tmprss2","Tnfaip8","Tpd52","Tsc22d1","Uap1","Ugdh","Unc13b","Vapa","Zbtb10")
Hieronymus_Androgen_mouse <- c("AA986860","Abcc4","Acsl3","Adam7","Ar","Cd200","Cenpn","Eaf2","Ell2","Fkbp5","Glra2","Gnmt","Herc3","Klk1","Klk1b1",
                               "Klk1b11","Klk1b16","Klk1b21","Klk1b22","Klk1b24","Klk1b26","Klk1b27","Klk1b3","Klk1b4","Klk1b5","Klk1b8","Klk1b9","Maf","Man1a","Mapre2",
                               "Med28","Mphosph9","Nkx3-1","Nnmt","Pip4k2b","Pmepa1","Ptger4","Srp72","Tbc1d9b","Tmprss2","Tnk1","Zbtb10")

ar_targets_list = list(
  'Metabolism' = c('Dhcr24','Idi1','Scd2','Scd4','Sms','Phyh','Hpgd','Aldh1a3','Hmgcs1','Sord','Ugdh','B4galt1','Dbi'), #Ddc
  'Proliferation/differentiation/apoptosis' =  c('Ell2','Maf','Ndrg1'), #,'Cdk8','Hes6','Id2','Ctbp1'
  'Transport/trafficking' = c('Slc26a2','Abcc4','Srp19','Ank','Lman1','Vapa','Sec24d','Fkbp5','Snap25','Azgp1','Appbp2'), #'Scnn1a'
  'Transcription regulation' = c('Klf4','Arid5b'), #Myc MRF2=Arid5b
  'Protease/protease Inhibitor' = c('Klk4','Cpd','Tmprss2'
                                    #'Klk1b8','Klk1b9',
                                    #"Klk1b11","Klk1b27","Klk1b21","Klk1b24","Klk1b3",
                                    #'Klk1b3','Klk1b16','Klk1b26','Klk1b4','Klk1b5', #Serpin1
  ),
  'Signal transduction' = c('Lifr','Inpp4b','Camkk2','Ptpn21','Tpd52','Gm49337','Akap12','Iqgap2','Sdc4','Mertk'), #'Bche', 'Peg3', 'Pik3r3',
  'Structure/motility/adhesion' = c('Adrm1','Krt19','Actn1','Krt8','Gm9774','Itgav','Dctn3'), #'Mmp16','Fn1'
  'Stress response' = c('Orm1','Orm3','Orm2','Dnajb9','Gsr'),
  'Other functions' = c('Rbm10','Uap1','Herc3','B2m') #'St7'
)

INPUT_LABELS = list(
  '1855-EPC-3m' = 'Erg+Pten- 3 months 1',
  '1856-EPC-3m' = 'Erg+Pten- 3 months 2',
  '1923-PC-3m' = 'Pten- 3 months 1',
  '1928-PC-3m' = 'Pten- 3 months 2',
  '1853-EPh-3m' = 'Wildtype 3 months'
)

CONDITION_LABELS <- list(
  '1855-EPC-3m' = 'Erg+Pten-',
  '1856-EPC-3m' = 'Erg+Pten-',
  '1923-PC-3m' = 'Pten-',
  '1928-PC-3m' = 'Pten-',
  '1853-EPh-3m' = 'Wildtype'
)

sample_palette <- c('#00B050','#4251A2','#FFCC00','#E60012','#34BDEF')
cluster_palette <- c(as.vector(pals::alphabet(26)), as.vector(pals::alphabet2(26)))
condition_palette <- c('#272E6A','#D51F26','#4ED956') #,'gray90')
names(condition_palette) <- c('Erg+Pten-','Pten-','Wildtype') #,NA)
condition_palette_continuous <- list('Erg+Pten-'=pals::brewer.blues(20), 
                                     'Pten-'=pals::brewer.reds(20),
                                     'Wildtype'=pals::brewer.greens(20))

bionames_palette_list = sample_palette
names(bionames_palette_list) =  unname(INPUT_LABELS)

cluster_palette_list = cluster_palette
names(cluster_palette_list) =  as.character(seq_along(cluster_palette))

#epiClass.palette <- c("#0086ED","#4FC601","#FF4A46","#BC23FF")
epiClass.palette <- c("#222222","#F3C300","#875692","#008856")
names(epiClass.palette) <- c("Epi_Basal_1","Epi_Luminal_1","Epi_Luminal_2Psca","Epi_Luminal_3Foxi1")

## Epithelial cells as defined with Wouter scores ##
epi_clusters <- c(2,4,9,13,15,17,18,20,21,25,26)


#clusters 5,6,7 => L1
#cluster 2 => Basal
#cluster 1,10 => IM
#cluster 8 => Mixed
#cluster 3,4,9,11 => Tumor-L2

celltype_palette_list <- list("Basal" = "#1C86EE",
                              "IM" = "#D51F26",
                              "L1" = "#FFC125",
                              "L2"= "#4C005C",
                              "Tumor-L2" = "#FF90C9",
                              "Mixed" = "#767171")

#### Step 1. Create the data set ####
mergeMultiomeSeurat <- function(){
  min_cells <- 10 # "min_cells" parameter is the minimum number of cells where the gene must be detected in order to be included in the data scaling
  min_features <- 300
  sequencing.dir <- '/storage/projects/sawyers/erg/scRNA/processed/'
  sc.data <- list.files(path = sequencing.dir, full.names = T)
  sc.data <- grep("WF-2020",sc.data, value = T, invert = T)
  sc.data
  (rna.counts <- list.files(path = sc.data, pattern = '^filtered_feature_bc_matrix.h5$', recursive = T, full.names = T))
  (sc.sample <- basename(sc.data))
  sample.list <- lapply(sc.sample, function(x){
    list('rna'=grep(x, rna.counts, value = T))
  }) %>% set_names(sc.sample)
  
  print("Found the following samples")
  print(unlist(sample.list))
  
  #### 1. Create seurat object ####
  seurat.obj <- lapply(names(sample.list), function(i){
    x <- sample.list[[i]]$rna
    mat <- Read10X_h5(x)
    CreateSeuratObject(counts =  mat, min.cells = min_cells, project=i)
  })
  (names(seurat.obj) <- sapply(seurat.obj, function(x){ unique(levels(x$orig.ident)) }))
  length(seurat.obj)
  #### 2. Merge unnormalized seurat objects ####
  # https://satijalab.org/seurat/articles/merge_vignette.html
  erg.scRNA <- merge(x = seurat.obj[[1]],
                     y = c(seurat.obj[[2]],
                           seurat.obj[[3]],
                           seurat.obj[[4]],
                           seurat.obj[[5]]),
                     # already done above 
                     add.cell.ids = names(seurat.obj),
                     project = 'Erg')
  erg.scRNA
  
  #### 3. Attributes ####
  ## Better for Integration in ArchR
  head(Cells(erg.scRNA))
  table(gsub("_.*","",Cells(erg.scRNA)))
  erg.scRNA <- RenameCells(object = erg.scRNA, new.names = gsub("_","#",Cells(erg.scRNA)) )
  head(Cells(erg.scRNA))
  table(gsub("#.*","",Cells(erg.scRNA)))
  
  ## Bionames
  erg.scRNA@meta.data$bioNames <- sapply(erg.scRNA@meta.data$orig.ident, function(i) INPUT_LABELS[[i]], USE.NAMES = F)
  table(erg.scRNA@meta.data$bioNames)
  ## Condition
  erg.scRNA@meta.data$condition <- sapply(erg.scRNA@meta.data$orig.ident, function(i) CONDITION_LABELS[[i]])
  table(erg.scRNA@meta.data$condition)
  
  ### dsRED and Foxa1 are driven by same promoter in all except empty vector and BL6 ###
  ### Sum dsRED and FoxA1 counts for all but EV ###
  
  table(gsub("#.*","",names(erg.scRNA@assays$RNA@counts['Erg',])))
  table(gsub("#.*","",names(erg.scRNA@assays$RNA@counts['ErgIRESGFP',])))
  
  ## Just record percentage of endogenous Foxa1
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "^Erg$", col.name = "percent_endogenous_erg")
  ## Just record percentage of plasmid
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "ErgIRESGFP", col.name = "percent_ErgIRESGFP")
  
  sapply(unique(erg.scRNA$bioNames), function(x){
    sum(erg.scRNA@assays$RNA@counts['Erg', which(erg.scRNA$bioNames==x) ])})
  
  sapply(unique(erg.scRNA$bioNames), function(x){
    sum(erg.scRNA@assays$RNA@counts['ErgIRESGFP', which(erg.scRNA$bioNames==x) ])})
  
  tapply(erg.scRNA$percent_ErgIRESGFP, erg.scRNA$bioNames, range)
  
  # Replace Foxa1 counts with endogenous and plasmid counts within WT and mutants (exclude EV)
  (idx <-   which(!erg.scRNA@meta.data$bioNames %in% c("Wildtype 3 months")))
  table(gsub("#.*","",names(erg.scRNA@assays$RNA@counts['Erg',idx])))
  sum(erg.scRNA@assays$RNA@counts['Erg',idx])
  sum(erg.scRNA@assays$RNA@counts['ErgIRESGFP',idx])
  
  identical(erg.scRNA@assays$RNA@counts, erg.scRNA@assays$RNA@data)
  
  # R
  (empty.idx <-   which(erg.scRNA@meta.data$bioNames %in% c("Wildtype 3 months")))
  erg.scRNA@assays$RNA@counts['ErgIRESGFP',empty.idx] <- 0
  erg.scRNA@assays$RNA@data['ErgIRESGFP',empty.idx] <- 0
  
  # Replace Foxa1 counts with endogenous and plasmid counts
  # Since we specified EV and BL6 have 0 plasmid counts their totals will remain the same
  # Row replacement is very slow with sparse matrix
  # so we replace the old data using rbind
  replaceErg <- function(obj, slot='counts'){
    
    m <- GetAssayData(object = obj, slot = slot)
    message("Original data [",nrow(m),'x',ncol(m),']')
    rorder <- rownames(m)
    y <- m
    
    x <- t( as.matrix(data.frame("Erg"=m['Erg',] + m['ErgIRESGFP',])) )
    ridx <- which(rownames(y)=='Erg')
    y <- y[-ridx,]
    y <- rbind(y, x)
    y <- y[rorder,]
    if(!identical( rorder, rownames(y) ))
      stop("Problem with rownames")
    
    message("Updated data [",nrow(y),'x',ncol(y),']')
    
    obj <- SetAssayData(object = obj, slot = slot, new.data = y)
  }
  erg.scRNA <- replaceErg(erg.scRNA, slot = 'counts')
  erg.scRNA <- replaceErg(erg.scRNA, slot = 'data')
  identical(erg.scRNA@assays$RNA@counts, erg.scRNA@assays$RNA@data)
  
  tapply(erg.scRNA@assays$RNA@counts['Erg',], erg.scRNA$bioNames, sum)
  
  tapply(erg.scRNA@assays$RNA@counts['ErgIRESGFP',], erg.scRNA$bioNames, sum)
  
  #
  ## Add reporter counts to metadata
  erg.scRNA$ErgIRESGFP_counts <- 0
  erg.scRNA$ErgIRESGFP_counts <- erg.scRNA@assays$RNA@counts['ErgIRESGFP',]
  ## Remove the gene reporter
  message("Removing Erg reporter")
  counts <- GetAssayData(object = erg.scRNA, assay = 'RNA')
  (rem <- which(rownames(counts) == 'ErgIRESGFP'))
  if(length(rem) > 0){
    counts <- counts[-c(rem),]
    erg.scRNA <- subset(erg.scRNA, features = rownames(counts))
  }
  identical(erg.scRNA@assays$RNA@counts, erg.scRNA@assays$RNA@data)
  
  head(colnames(erg.scRNA))
  tail(colnames(erg.scRNA))
  table(gsub('#.*',"",colnames(erg.scRNA)))
  table(erg.scRNA$orig.ident)
  
  #### 3. Add percentages and filter
  # https://satijalab.org/seurat/articles/sctransform_vignette.html
  # store mitochondrial percentage in object meta data
  ribosomal_genes = fread("/home/erik/annotation/genome/tables/mouse_ribosomal_genes.txt", header = F)
  ribosomal_genes <- paste(paste0("^",unname(unlist(ribosomal_genes))), collapse = "|")
  grep("^Mt-|^mt-",rownames(erg.scRNA), value = T)
  grep(ribosomal_genes,rownames(erg.scRNA), ignore.case = T, value = T)
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "^Mt-|^mt-", col.name = "percent_mito")
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = ribosomal_genes, col.name = "percent_ribo")
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "^Malat1", col.name = "percent_malat1")
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "^Foxa1", col.name = "percent_foxa1")
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "^Krt4$", col.name = "percent_krt4")
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "^Krt5$", col.name = "percent_krt5")
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "^Krt8$", col.name = "percent_krt8")
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "^Krt18$", col.name = "percent_krt18")
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "^Krt19$", col.name = "percent_krt19")
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "^Trp63$", col.name = "percent_trp63")
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "^Epcam", col.name = "percent_epcam")
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "^Pten", col.name = "percent_pten")
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "^Erg", col.name = "percent_erg")
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "^Trp53$", col.name = "percent_trp53")
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "^Nkx3-1", col.name = "percent_nkx3.1")
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "^Pbsn$", col.name = "percent_pbsn")
  erg.scRNA <- PercentageFeatureSet(erg.scRNA, pattern = "Tmprss2", col.name = "percent_tmprss2")
  
  range(erg.scRNA$percent_mito, na.rm = T)
  range(erg.scRNA$percent_erg, na.rm = T)
  
  for(i in colnames(erg.scRNA@meta.data)){
    idx <- which(is.nan(erg.scRNA@meta.data[[i]]))
    if(length(idx)>0){
      message(i," has missing %")
      erg.scRNA@meta.data[idx,i] <- 0
    }
  }
  
  #### 4. Basic QC filtering ####
  # Visualize QC metrics as a violin plot
  dtplot <- data.table(erg.scRNA@meta.data)[,list(bioNames, 
                                                  nFeature_RNA, nCount_RNA, 
                                                  percent_malat1,
                                                  percent_mito, percent_ribo, 
                                                  percent_erg, percent_endogenous_erg, 
                                                  percent_krt5, percent_krt8,
                                                  percent_pten, percent_trp53)]
  
  dtplot <- melt(dtplot, id.vars = "bioNames")
  df <- data.frame(variable=c('percent_mito','nFeature_RNA','nCount_RNA'), y=c(15,300,100000))
  
  png(filename = file.path(plot.dir, 'pre_processing_violin.png'), width = 20, height = 12, bg = 'white', res = 150, units = 'in')
  ggplot(dtplot, aes(x=bioNames, y=value, fill=bioNames, label=variable)) +
    geom_violin(trim = F, scale = 'width') +
    scale_fill_manual("", values = bionames_palette_list) +
    xlab("") + ylab("") +
    geom_boxplot(width=.4, fill='white', alpha=0.6, outlier.size = 0.0) +
    geom_hline(data = df, aes(yintercept=y), color='red', linetype='dashed') +
    facet_wrap(~ variable, nrow = 3, scales = 'free_y') + theme_bw(14) +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          panel.grid = element_blank(), panel.background = element_blank(), 
          plot.background = element_blank(), legend.background = element_blank(),
          strip.background = element_blank(), strip.text = element_text(face = 'bold'), )
  dev.off()
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  plot1 <- FeatureScatter(erg.scRNA, feature1 = "nCount_RNA", feature2 = "percent_mito", group.by = 'bioNames') +
    facet_grid(.~ colors) + scale_color_manual("", values = bionames_palette_list) +
    geom_hline(yintercept = df[df$variable=='percent_mito','y'], color='red', linetype='dashed') +
    theme(plot.background = element_blank(), panel.background = element_blank())
  plot2 <- FeatureScatter(erg.scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'bioNames') +
    facet_grid(.~ colors) + scale_color_manual("", values = bionames_palette_list) +
    geom_hline(yintercept = 300, color='red', linetype='dashed') + 
    geom_vline(xintercept = df[df$variable=='nCount_RNA','y'], color='red', linetype='dashed') +
    theme(plot.background = element_blank(), panel.background = element_blank())
  plot3 <- FeatureScatter(erg.scRNA, feature1 = "nFeature_RNA", feature2 = "percent_ribo", group.by = 'bioNames') +
    facet_grid(.~ colors) + scale_color_manual("", values = bionames_palette_list) +
    geom_hline(yintercept = df[df$variable=='percent_ribo','y'], color='red', linetype='dashed') + 
    geom_vline(xintercept = df[df$variable=='nFeature_RNA','y'], color='red', linetype='dashed') +
    theme(plot.background = element_blank(), panel.background = element_blank())
  png(filename = file.path(plot.dir, 'pre_processing_scatter.png'), width = 20, height = 14, bg = 'transparent', res = 200, units = 'in')
  do.call(grid.arrange, c(list(plot1,plot2,plot3), nrow=3))
  dev.off()
  
  
  lowRNA <- WhichCells(erg.scRNA, expression = nFeature_RNA < df[df$variable=='nFeature_RNA','y'])
  
  message("Filtering cells")
  erg.scRNA <- subset(erg.scRNA, cells=setdiff(Cells(erg.scRNA), lowRNA))
  erg.scRNA <- subset(erg.scRNA, subset = nCount_RNA > 200 & 
                        nCount_RNA < df[df$variable=='nCount_RNA','y'] & 
                        percent_mito < df[df$variable=='percent_mito','y'])
  identical(erg.scRNA@assays$RNA@counts, erg.scRNA@assays$RNA@data)
  
  erg.scRNA
  
  counts <- GetAssayData(erg.scRNA, slot="counts", assay="RNA")   
  genes.count.expression <- rowCounts(counts>0)
  names(genes.count.expression) <- rownames(counts)
  which(genes.count.expression < 5)
  genes.percent.expression <- rowMeans(counts>0 )*100   
  
  # Visualize QC metrics as a violin plot
  dtplot <- data.table(erg.scRNA@meta.data)[,list(bioNames, nFeature_RNA, nCount_RNA, 
                                                  percent_malat1,
                                                  percent_mito, percent_ribo, 
                                                  percent_erg, percent_endogenous_erg, 
                                                  percent_krt5, percent_krt8,
                                                  percent_pten, percent_trp53)]
  dtplot <- melt(dtplot, id.vars = "bioNames")
  
  png(filename = file.path(plot.dir, 'post_processing_filter_violin.png'), width = 20, height = 12, bg = 'white', res = 150, units = 'in')
  ggplot(dtplot, aes(x=bioNames, y=value, fill=bioNames, label=variable)) +
    geom_violin(trim = F, scale = 'width') +
    scale_fill_manual("", values = bionames_palette_list) +
    xlab("") + ylab("") +
    geom_boxplot(width=.4, fill='white', alpha=0.6, outlier.size = 0.0) +
    geom_hline(data = df, aes(yintercept=y), color='red', linetype='dashed') +
    facet_wrap(~ variable, nrow = 3, scales = 'free') + theme_bw(14) +
    theme(axis.text.x = element_text(angle=-45, hjust=0),
          panel.grid = element_blank(), panel.background = element_blank(), 
          plot.background = element_blank(), legend.background = element_blank(),
          strip.background = element_blank(), strip.text = element_text(face = 'bold'), )
  dev.off()
  
  plot1 <- FeatureScatter(erg.scRNA, feature1 = "nCount_RNA", feature2 = "percent_mito", group.by = 'bioNames') +
    facet_grid(.~ colors) + scale_color_manual("", values = bionames_palette_list) +
    geom_hline(yintercept = df[df$variable=='percent_mito','y'], color='red', linetype='dashed') +
    theme(plot.background = element_blank(), panel.background = element_blank())
  plot2 <- FeatureScatter(erg.scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'bioNames') +
    facet_grid(.~ colors) + scale_color_manual("", values = bionames_palette_list) +
    geom_hline(yintercept = 300, color='red', linetype='dashed') + 
    geom_vline(xintercept = df[df$variable=='nCount_RNA','y'], color='red', linetype='dashed') +
    theme(plot.background = element_blank(), panel.background = element_blank())
  plot3 <- FeatureScatter(erg.scRNA, feature1 = "nFeature_RNA", feature2 = "percent_ribo", group.by = 'bioNames') +
    facet_grid(.~ colors) + scale_color_manual("", values = bionames_palette_list) +
    geom_hline(yintercept = df[df$variable=='percent_ribo','y'], color='red', linetype='dashed') + 
    geom_vline(xintercept = df[df$variable=='nFeature_RNA','y'], color='red', linetype='dashed') +
    theme(plot.background = element_blank(), panel.background = element_blank())
  png(filename = file.path(plot.dir, 'post_processing_scatter.png'), width = 20, height = 14, bg = 'transparent', res = 200, units = 'in')
  do.call(grid.arrange, c(list(plot1,plot2,plot3), nrow=3))
  dev.off()
  
  dtplot <- data.table(erg.scRNA@meta.data, keep.rownames = 'cells')[,list(cells, bioNames, 
                                                                           nFeature_RNA, nCount_RNA, 
                                                                           percent_epcam,
                                                                           percent_mito, percent_ribo, 
                                                                           percent_erg, percent_endogenous_erg, 
                                                                           percent_ErgIRESGFP,
                                                                           percent_krt5, percent_krt8,
                                                                           percent_pten, percent_trp53)]
  
  features <- c("percent_mito","percent_ribo","percent_epcam","percent_erg","percent_ErgIRESGFP","percent_pten")
  bionames <- unique(dtplot$bioNames)
  pl <- lapply(bionames, function(j){
    toplot <- dtplot[bioNames==j]
    gp <- lapply(features, function(i){
      ggplot(toplot[order(get(i))], aes(x=nCount_RNA, y=nFeature_RNA, color=get(i))) + geom_point() +
        xlab("No. UMI per cell") + ylab("No. Genes per cell") + scale_color_viridis(i, discrete = F) +
        theme_bw(14) + theme(panel.grid = element_line(size = .5), 
                             text = element_text(color='black'), 
                             legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5)) +
        guides(colour = guide_colorbar(title.position = "left"))
    }) %>% set_names(features)
    gp
  }) %>% set_names(bionames)
  
  for(i in names(pl)){
    png(filename = file.path(plot.dir,paste0('scRNA_scatter_',i,'.png')), width = 30, height = 4, units = 'in', res = 150)
    do.call(grid.arrange, c(pl[[i]], nrow=1, top=i))
    dev.off()
  }
  
  ## Cell cycle
  library(Seurat)
  cc.genes
  #https://doi.org/10.1016/j.cell.2015.05.002 - Table S2
  cell_cycle_markers_df = fread("/home/erik/annotation/markers/Human_Cell_Cycle_Macosko_et_al.txt", sep='\t')
  s.genes <- Hmisc::capitalize( tolower(cell_cycle_markers_df$S[cell_cycle_markers_df$S!=""]) )
  g2m.genes <- Hmisc::capitalize( tolower(cell_cycle_markers_df$G2.M[cell_cycle_markers_df$G2.M!=""]) )
  m.genes <- Hmisc::capitalize( tolower(cell_cycle_markers_df$M[cell_cycle_markers_df$M!=""]) )
  
  #### 5. Normalize merged objets using Sctransform ####
  # run sctransform, default is 3000 variable features
  # An object of class Seurat 
  # 17038 features across 35456 samples within 1 assay 
  # Active assay: RNA (17038 features, 0 variable features)
  
  erg.scRNA <- SCTransform(erg.scRNA,
                           vst.flavor="v2",
                           method = "glmGamPoi",
                           vars.to.regress = c("percent_mito","nFeature_RNA","nCount_RNA"),
                           return.only.var.genes = FALSE,
                           verbose = T)
  
  ## Regress out signal from SCT assay
  erg.scRNA <- CellCycleScoring(object = erg.scRNA, s.features = s.genes, g2m.features = g2m.genes, assay = 'SCT')
  ## now that we have SCT normalized cell-cycle scores we run normalization again and regress these out
  erg.scRNA <- SCTransform(erg.scRNA,
                           vst.flavor="v2",
                           assay = 'RNA',
                           new.assay.name = 'SCT' ,
                           method = "glmGamPoi",
                           return.only.var.genes = FALSE,
                           vars.to.regress = c('percent_mito', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))
  erg.scRNA
  # Active assay: SCT (17037 features, 3000 variable features)
  # 1 other assay present: RNA
  
  # find genes expressed in <1 percent of cells
  counts <- GetAssayData(erg.scRNA, slot="counts", assay="RNA")
  genes.percent.expression <- rowMeans(counts>0 )*100
  names(genes.percent.expression[genes.percent.expression<1]) %>% sort
  
  # this gene is removed by SCT
  rem.gene <- setdiff(rownames(erg.scRNA@assays$RNA@counts), rownames(erg.scRNA@assays$SCT@counts))
  if(length(rem.gene)){
    message("Removing following gene(s):",paste(rem.gene, collapse = ","))
    genes.filter <- setdiff(rownames(counts), rem.gene)
    erg.scRNA[["RNA"]] <- CreateAssayObject(counts = erg.scRNA[["RNA"]]@counts[genes.filter , ] )
  }
  
  ## Normalize RNA assay
  erg.scRNA <- NormalizeData(object = erg.scRNA, assay='RNA', normalization.method = "LogNormalize")
  erg.scRNA <- CellCycleScoring(object = erg.scRNA,
                                s.features = s.genes, g2m.features = g2m.genes,
                                assay = 'RNA')
  erg.scRNA <- FindVariableFeatures(erg.scRNA, selection.method = "vst", assay = 'RNA')
  ## Regress out signal from RNA assay
  erg.scRNA <- ScaleData(erg.scRNA,
                         assay = 'RNA',
                         vars.to.regress = c('percent_mito', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'))
  
  #### 6. Attributes ####
  saveRDS(erg.scRNA, file = scrna.object)
  
  #### 7. PCA and UMAP ####
  # erg.scRNA[["SCT"]]@scale.data contains the residuals (normalized values), and is used directly as input to PCA
  # erg.scRNA[["SCT"]]@counts are corrected UMI counts
  # erg.scRNA[["SCT"]]@data are log-normalized corrected UMI counts
  
  # library(future)
  # # max ram allowed here is 1000 * 1024^2 = 1GB
  options('future.globals.maxSize' = 1000 * 1024^2)
  
  erg.scRNA.rna <- RunPCA(erg.scRNA, verbose = FALSE, assay = 'RNA')
  
  # determine dimensionality
  pdf(file = file.path(plot.dir, 'pca_dimension_elbow_regress_RNA.pdf'), width = 6, height = 6)
  print(ElbowPlot(erg.scRNA.rna, ndims = 50))
  dev.off()
  
  png(filename = file.path(plot.dir, 'pca_dimension_heatmap_regress_RNA.png'), width = 8, height = 16, units = 'in', bg = 'transparent', res = 150)
  DimHeatmap(erg.scRNA.rna, dims=1:15, cells = 2000, assays = 'RNA', balanced = T)
  dev.off()
  
  # evaluate umap
  DefaultAssay(erg.scRNA.rna) <- 'RNA'
  erg.scRNA.rna <- RunUMAP(erg.scRNA.rna, dims= 1:30, n.neighbors = 30L,
                           min.dist = 0.3, spread = 1, verbose = FALSE, assay = 'RNA', metric = 'cosine')
  
  g1 <- DimPlot(erg.scRNA.rna, label = F, group.by = 'nCount_RNA') + NoLegend()
  
  g1 <- DimPlot(erg.scRNA.rna, label = F, group.by = 'bioNames', cols = sample_palette) + NoLegend()
  
  table(erg.scRNA.rna$bioNames)
  
  ### SCT
  erg.scRNA.sct <- RunPCA(erg.scRNA, verbose = FALSE, assay = 'SCT') #features=NULL, meaning PCA is run for variable features for the Assay.
  
  # determine dimensionality
  pdf(file = file.path(plot.dir, 'pca_dimension_elbow_regress_SCT.pdf'), width = 6, height = 6)
  print(ElbowPlot(erg.scRNA.sct, ndims = 50))
  dev.off()
  
  png(filename = file.path(plot.dir, 'pca_dimension_heatmap_regress_SCT.png'), width = 8, height = 16, units = 'in',
      bg = 'transparent', res = 150)
  DimHeatmap(erg.scRNA.sct, dims=1:15, cells = 2000, balanced = T)
  dev.off()
  
  # evaluate umap
  DefaultAssay(erg.scRNA.sct) <- 'SCT'
  erg.scRNA.sct <- RunUMAP(erg.scRNA.sct, dims= 1:30, n.neighbors = 30L, seed.use = 42,
                           min.dist = 0.5, spread = 1.2, verbose = FALSE, assay = 'SCT', metric = 'cosine')
  
  g1 <- DimPlot(erg.scRNA.sct, label = F, group.by = 'nCount_RNA') + NoLegend()
  
  g1 <- DimPlot(erg.scRNA.sct, label = F, group.by = 'bioNames', cols = sample_palette) + NoLegend()
  
  g1 <- DimPlot(erg.scRNA.sct, label = F, group.by = 'condition', cols = condition_palette) + NoLegend()
  
  
  #### 8. Choose Assay ####
  DefaultAssay(erg.scRNA)
  erg.scRNA <- RunPCA(erg.scRNA, verbose = FALSE, assay = 'SCT') #features=NULL, meaning PCA is run for variable features for the Assay.

  pdf(file = file.path(plot.dir, 'pca_dimension_elbow_regress.pdf'), width = 6, height = 6)
  print(ElbowPlot(erg.scRNA, ndims = 50))
  dev.off()
  
  png(filename = file.path(plot.dir, 'pca_dimension_heatmap_regress.png'), width = 8, height = 16, units = 'in',
      bg = 'transparent', res = 150)
  DimHeatmap(erg.scRNA, dims=1:15, cells = 2000, balanced = T)
  dev.off()
  
  # evaluate umap
  erg.scRNA <- RunUMAP(erg.scRNA, dims= 1:30, n.neighbors = 30L, seed.use = 42,
                       min.dist = 0.5, spread = 1.2, verbose = FALSE, assay = 'SCT', metric = 'cosine')
  
  g1 <- DimPlot(erg.scRNA, label = F, group.by = 'nCount_RNA') + NoLegend()
  g1$data$nCount_RNA <- as.numeric(g1$data$nCount_RNA)

  g1 <- DimPlot(erg.scRNA, label = F, group.by = 'bioNames', cols = sample_palette) + NoLegend()
  
  g1 <- DimPlot(erg.scRNA, label = F, group.by = 'condition', cols = condition_palette) + NoLegend()
  
  erg.scRNA <- FindNeighbors(erg.scRNA,
                             dims = 1:30,
                             annoy.metric = 'cosine',
                             assay = 'SCT',
                             verbose = FALSE)
  
  
  saveRDS(erg.scRNA, file = scrna.object)
  
  #### 9. Evaluate clustering resolution ####
  library(cluster)
  library(future)
  plan()
  res <- seq(from = 0.1, to = 2.0, by = 0.1)
  erg.scRNA <- FindClusters(erg.scRNA, verbose = FALSE, algorithm = 4, resolution = res, method='igraph')
  
  resm <- grep("SCT_snn_res", colnames(erg.scRNA@meta.data), value = T)
  (resm <- resm[order(gsub("SCT_snn_res.","",resm))])
  
  sils <- mclapply(resm, function(r){
    sil <- silhouetteScore(obj = erg.scRNA, reduction = 'umap', clusters =  r, dims = 1:2)
  }, mc.cores = 2) %>% set_names(resm)
  
  
  um <- lapply(resm, function(r){
    score <- summary(sils[[r]])$avg.width
    g1 <- DimPlot(erg.scRNA, label = TRUE, group.by = r, cols = cluster_palette) + 
      ggtitle(label = paste("Cluster res.:", gsub("SCT_snn_res.","",r),
                            "\nmean silhouette ", formatC(score, digits = 3)))
  }) %>% set_names(resm)
  length(um)
  
  png(filename = file.path(plot.dir, 'silhouette_scores_regress.png'), width = 20, height = 14, units = 'in', res=150, bg = 'transparent')
  do.call(grid.arrange, c(um[1:10], nrow=3, top='resolutions'))
  dev.off()
  
  png(filename = file.path(plot.dir, 'silhouette_scores_regress_2.png'), width = 20, height = 14, units = 'in', res=150, bg = 'transparent')
  do.call(grid.arrange, c(um[11:20], nrow=3, top='resolutions'))
  dev.off()
  
  pdf(file = file.path(plot.dir, 'silhouette_scores_regress.pdf'), width = 10, height = 8, onefile=T)
  for(r in resm){
    g <- um[[r]]
    p1 <- ~plot(sils[[r]], main=paste("resolution:",r), cex.names = par("cex.axis")) # RStudio sometimes does not display silhouette plots correctly
    grid.arrange(cowplot::as_grob(p1), g, nrow=1, widths=c(.35,.65))
  }
  dev.off()
  
  
  #### 10. Pick clustering resolution ####
  cresolution <- 0.9
  erg.scRNA <- FindClusters(erg.scRNA, verbose = FALSE, algorithm = 4, resolution = cresolution, method = 'igraph')
  
  cresolution <- paste0('SCT_snn_res.', cresolution)
  
  erg.scRNA@meta.data$seurat_clusters <- erg.scRNA@meta.data[,cresolution]
  identical(erg.scRNA@meta.data$seurat_clusters, erg.scRNA@meta.data[[cresolution]])
  
  Idents(object = erg.scRNA) <- cresolution
  
  # Plot UMAPS
  g1 <- DimPlot(erg.scRNA, label = TRUE, cols = cluster_palette) + NoLegend()
  
  #### 11. save object ####
  saveRDS(erg.scRNA, file = scrna.object)
  message("Saved seurat objects to: ", scrna.object)
  
}


# Description of Assay and Slots
# https://github.com/satijalab/seurat/wiki/Assay
# RNA: for gene level : D.E., violin plots
# counts = raw counts
# data = log1p(counts) after normalization divide counts by 10K (for instance)
# scale.data = normalized and scaled for dim reduction
#
# SCT: cell level, UMAP clustering
# counts = (corrected UMI counts if all cells were sequenced with same depth) counts
# data = log1p(counts)
# scale.data = pearson residuals

## Which assay for what?
# https://github.com/satijalab/seurat/issues/4082
# RNA@data : D.E., Feature/Violin plots
# SCT@scale.data : cell level, UMAP clustering
#
# https://github.com/satijalab/seurat/issues/2899
# You can use the corrected Pearson residuals (stored in the scale.data slot of the integrated assay) for joint analysis, like PCA and clustering
# You should never use integrated/corrected values when performing DE.
# We typically use the @data slot (log-normalized corrected UMI counts) instead of scale.data (pearson residuals).
# Not suggested but if you do run DE on the scale.data slot (for any assay), the logFC values are not meaningful.
# Typically use:
# scale.data slot for heat maps
# Data slot for Feature/Violin plots.

#### Step 2. Analyses ####
main <- function(){
  #buildMultiomeSeurat()
  
  library(future)
  # max ram allowed here is 1000 * 1024^2 = 1GB
  options('future.globals.maxSize' = 1000 * 2048^2)
  # check the current active plan
  plan()
  # change the current plan to access parallelization
  plan("multiprocess", workers = 4)
  #plan(strategy = 'sequential')
  
  #### Load full scRNA object ####
  erg.scRNA <- readRDS(file=scrna.object)
  
  #### QC plots ####
  FeatureScatter(erg.scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

  gv <- ViolinPlotQC(erg.scRNA)
  
  png(filename = file.path(plot.dir, 'scRNA_counts_violin.png'), units = 'in', width = 18, height = 5, bg = 'transparent', res = 250)
  do.call(grid.arrange, c(gv, nrow=1))
  dev.off()
  
  
  #### Standard UMAPS ####
  # nCount_RNA is the total number of molecules detected within a cell. 
  # nFeature_RNA is the number of genes detected in each cell.
  
  g <- DimPlot(erg.scRNA, label = F, group.by = 'nCount_RNA') + NoLegend()
  
  g0 <- DimPlot(erg.scRNA, label = F, group.by = 'condition', cols = condition_palette) %>% 
    
  g1 <- DimPlot(erg.scRNA, label = F, group.by = 'bioNames', cols = bionames_palette_list)
  
  g2 <- DimPlot(erg.scRNA, label = F, group.by = 'seurat_clusters', cols = cluster_palette)
  
  #### UMAP highlight cells ####
  g0 <- FeaturePlot(erg.scRNA, features = 'percent_epcam', reduction = 'umap', cols = brewer.blues(10))
  

  #### Cell cycle scores ####
  x <- c('G2M.Score', 'S.Score')
  (g2 <- FeaturePlot(erg.scRNA, features = x, pt.size = 0.2, ncol = 1, slot = 'scale.data', blend = T))
  
  ## Cell cycle check
  cell_cycle_markers_df = fread("/home/erik/annotation/markers/Human_Cell_Cycle_Macosko_et_al.txt", sep='\t')
  s.genes <- Hmisc::capitalize( tolower(cell_cycle_markers_df$S[cell_cycle_markers_df$S!=""]) )
  g2m.genes <- Hmisc::capitalize( tolower(cell_cycle_markers_df$G2.M[cell_cycle_markers_df$G2.M!=""]) )
  m.genes <- Hmisc::capitalize( tolower(cell_cycle_markers_df$M[cell_cycle_markers_df$M!=""]) )
  
  ## Plot PCA (same for both assays) ##
  cc <- RunPCA(erg.scRNA, features = c(s.genes, g2m.genes))
  g1 <- DimPlot(cc, reduction = 'pca', group.by = 'Phase', split.by = 'Phase')
  g2 <- DimPlot(cc, reduction = 'pca', group.by = 'Phase')
  png(file = file.path(plot.dir,'sct_cell_cycle.png'), width = 10, height = 4, bg = 'transparent', units = 'in', res = 150)
  grid.arrange(g1, g2, nrow=1, widths=c(.65,.35))
  dev.off()
  
  erg.scRNA <- AddModuleScore(object = erg.scRNA,
                              features = list(m.genes),
                              assay = 'RNA',
                              ctrl = 20,
                              name = 'M.score', # appends a 1
                              search = F) # searches HUGO contaiing human gene names only
  colnames(erg.scRNA@meta.data)[which(colnames(erg.scRNA@meta.data)=='M.score1')] <- 'M.score'
  
  g0 <- DimPlot(erg.scRNA, label = F, group.by = 'Phase', cols = as.vector(pals::glasbey(3)))
  
  
  #### Cell cycle violin ####
  g0 <- VlnPlot(erg.scRNA, features = 'G2M.Score', assay = 'RNA', slot='data',
                pt.size = 0, flip = T, group.by = 'seurat_clusters', stack = F, cols = cluster_palette_list) + xlab("")
  df <- g0$data
  g0 <- ggplot(df, aes(x=ident, y=G2M.Score, fill=ident)) +
    geom_violin(trim = T, scale = 'width', draw_quantiles = c(.25,.5,.75)) +
    scale_fill_manual("", values = cluster_palette_list) +
    xlab("") + ylab("Expression level") + ggtitle("G2M score") +
    theme_classic(14) + theme(plot.title = element_text(face = 'bold', hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  
  g1 <- VlnPlot(erg.scRNA, features = 'S.Score', assay = 'RNA', slot='data',
                pt.size = 0, flip = T, group.by = 'seurat_clusters', stack = F, cols = cluster_palette_list) + xlab("")
  df <- g1$data
  g1 <- ggplot(df, aes(x=ident, y=S.Score, fill=ident)) +
    geom_violin(trim = T, scale = 'width', draw_quantiles = c(.25,.5,.75)) +
    scale_fill_manual("", values = cluster_palette_list) +
    xlab("") + ylab("Expression level") + ggtitle("S score") +
    theme_classic(14) + theme(plot.title = element_text(face = 'bold', hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  pdf(file = file.path(plot.dir,'cell_cycle_violin.pdf'), width = 5, height = 8, bg = 'transparent')
  grid_arrange_list_shared_legend(list(g0, g1), nrow = 2, legend.position = 'right')
  dev.off()
  
  
  ###------------------------------------###
  #### Load reclustered epithelial data ####
  ###------------------------------------###
  ss_epi <- BuildEpiSeurat(obj)
  ss_epi <- readRDS(file=epi.scrna.object)
  unique(ss_epi$seurat_clusters)
  table(ss_epi$seurat_clusters)
  
  DimPlot(ss_epi, group.by = 'seurat_clusters', reduction = 'umap', cols = cluster_palette)
  
  DimPlot(ss_epi, label = F, group.by = 'lum_bas_int_classification', reduction = 'umap', cols = luminal_basal_int_palette)

  #### epi Umaps ####
  # Plot UMAPS
 
  DimPlot(ss_epi, group.by = 'bioNames', cols = bionames_palette_list, reduction = 'umap')
  
  DimPlot(ss_epi, group.by = 'seurat_clusters', cols = cluster_palette, reduction = 'umap')
  
  FeaturePlot(ss_epi, features = 'percent_epcam', reduction = 'umap')
  
  FeaturePlot(ss_epi, features = 'Erg', reduction = 'umap', cols = rocket_palette[-1])
  
  #### Fig 1D ####
  FeaturePlot(ss_epi, features = 'Erg', reduction = 'umap', cols = c("lightgrey", "red"))
  
  FeaturePlot(ss_epi, features = 'percent_ErgIRESGFP', reduction = 'umap')
  
  FeaturePlot(ss_epi, features = 'Tmprss2', reduction = 'umap')
  
  DimPlot(ss_epi, label = F, group.by = 'lum_bas_int_classification', cols = luminal_basal_int_palette, reduction = 'umap')
  
  FeaturePlot(ss_epi, features = 'probability_luminal', cols = luminal_basal_palette_continuous, reduction = 'umap')
  
  DimPlot(ss_epi, group.by = 'condition', cols = condition_palette, reduction = 'umap')
  
  #### Epi barplots ####
  PlotBarVariableInClusters(seuratObj = ss_epi, 
                            orientation = 'vertical',
                            legend.title = "Samples",
                            by.group = 'seurat_clusters',
                            variable = 'bioNames', 
                            palette = bionames_palette_list, 
                            legend.column = 1,
                            text.angle=0)
  
  PlotBarVariableInClusters(seuratObj = ss_epi, 
                            orientation = 'vertical',
                            legend.title = "Genotype",
                            by.group = 'seurat_clusters',
                            variable = 'condition', 
                            palette = condition_palette, 
                            legend.column = 1,
                            text.angle=0)
  
  f <- c("S.Score","G2M.Score")
  g1 <- ViolinPlotEnhance(ss_epi,
                          features = f,
                          assay = 'RNA',
                          slot = 'data',
                          pt.size = 0,
                          flip = T,
                          group.by = 'seurat_clusters',
                          stack = T,
                          fill.by = 'ident', 
                          y_axis_label = "Normalized scores",
                          title = gsub("_"," ","cell cycle scores"),
                          #pal = brewer.reds(20)
                          pal = pals::coolwarm(25),
                          x_axis_text_angle = 90, 
                          strip.text.size = 14,
                          #pal = rev(pals::kovesi.linear_blue_5_95_c73(20)),
  ) + xlab("")
  
  g1 + theme(legend.position = 'none')
  
  #### Fig. 7A ####
  g <- ViolinPlotEnhance(ss_epi,
                         features = unlist(chromatin_modifiers_list),
                         assay = 'RNA',
                         slot = 'data',
                         pt.size = 0,
                         flip = T,
                         group.by = 'cell_type',
                         stack = T,
                         fill.by = 'ident',
                         title = "Chromatin modifiers",
                         pal = pals::coolwarm(25),
                         x_axis_text_angle = 90, 
                         strip.text.size = 14
  ) + xlab("")
  g
  
  # mol cell 2007 myles brown
  # AR collaborating factors FoxA1, GATA2 and Oct1 (Wang et al., 2007) and AR coactivator MED1 (TRAP220)
  f <- 'percent_epcam'
  g <- FeaturePlot(erg.scRNA, features = f)
  
  f <- 'percent_ErgIRESGFP'
  g <- FeaturePlot(erg.scRNA, features = f)
  
  f <- 'ErgIRESGFP_counts'
  g <- FeaturePlot(ss_epi, features = f)
  
  f <- 'Krt5'
  g <- FeaturePlot(ss_epi, features = f)
  
  ps <- lapply(c('Foxa1','Krt8','Krt18',
                 'Krt5','Krt14','Trp63'), function(f){
                   #ps <- lapply(c(markers_list$`Prostate luminal mature`, markers_list$`Prostate basal`), function(f){
                   FeaturePlot(erg.scRNA, features = f) 
                 })
  do.call(grid.arrange, c(ps, nrow=1, top='L/P markers'))
  
  ps <- lapply(c('Nkx3-1','Pbsn','Tmprss2'), function(f){
    FeaturePlot(erg.scRNA, features = f)
  })
  do.call(grid.arrange, c(ps, nrow=1))
  
  # cell division factors found by strings.db search of Cdca3
  cdf <- c('Cdca3','Cdc20','Top2a','Cdca8','Cdk1','Kif4','Birc5','Ccnb2','Nusap1','Pbk','Kif11','Mki67')
  ps <- lapply(cdf, function(f){
    g <- FeaturePlot(erg.scRNA, features = f)
  })
  do.call(grid.arrange, c(ps, nrow=2))
  
  #### Score Wouter's markers ####
  refineWouterL1L2MarkerClassification()
  
  #### UMAP epi cell types ####
  g0 <- DimPlot(ss_epi, group.by = c('wouter_class'), cols=wouter_pal)
  
  #### cell type violin ####
  features <- c('Epi_Basal_1', 'Epi_Luminal_1', 'Epi_Luminal_2Psca')
  g <- lapply(features, function(f){
    
    g1 <- VlnPlot(ss_epi, features = f, assay = 'RNA', slot='data',
                  pt.size = 0, flip = T, group.by = 'condition', stack = F, cols = condition_palette) + xlab("")
    df <- g1$data
    g1 <- ggplot(df, aes(x=ident, y=get(f), fill=ident)) +
      geom_violin(trim = T, scale = 'width', draw_quantiles = c(.25,.5,.75)) +
      scale_fill_manual("", values = condition_palette) +
      xlab("") + ylab("Geneset Score") + ggtitle(f) +
      theme_classic(14) + theme(plot.title = element_text(face = 'bold', hjust = 0.5),
                                axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
  })
  pdf(file = file.path(plot.dir,'EPI_violin_condition.pdf'), width = 10, height = 6, bg = 'transparent')
  grid_arrange_list_shared_legend(g, nrow = 3, legend.position = 'right')
  dev.off()
  
  features <- c('Epi_Basal_1', 'Epi_Luminal_1', 'Epi_Luminal_2Psca')
  g <- lapply(features, function(f){
    
    g1 <- VlnPlot(ss_epi, features = f, assay = 'RNA', slot='data',
                  pt.size = 0, flip = T, group.by = 'seurat_clusters', stack = F, cols = cluster_palette) + xlab("")
    df <- g1$data
    g1 <- ggplot(df, aes(x=ident, y=get(f), fill=ident)) +
      geom_violin(trim = T, scale = 'width', draw_quantiles = c(.25,.5,.75)) +
      scale_fill_manual("", values = cluster_palette) +
      xlab("") + ylab("Geneset Score") + ggtitle(f) +
      theme_classic(22) + theme(plot.title = element_text(face = 'bold', hjust = 0.5),
                                axis.text.x = element_text(angle = 0, vjust = 1, hjust=1), 
                                legend.position = "none")
    
  })
  
  pdf(file = file.path(plot.dir,'EPI_violin_cluster.pdf'), width = 9, height = 8, bg = 'transparent')
  do.call(grid.arrange, c(g, nrow=3))
  dev.off()
  
  
  #### Refine cell type labels ####
  #clusters 5,6,7,12 => L1
  #cluster 2 => Bas
  #cluster 1,10 => IM
  #cluster 8 => Mixed
  #cluster 3,4,9,11 => Tumor-L2
  ss_epi@meta.data$cell_type <- ""
  ss_epi@meta.data[ss_epi@meta.data$seurat_clusters %in% c(5,6,7,12),]$cell_type <- 'L1'
  ss_epi@meta.data[ss_epi@meta.data$seurat_clusters %in% c(2),]$cell_type <- 'Basal'
  ss_epi@meta.data[ss_epi@meta.data$seurat_clusters %in% c(1,10),]$cell_type <- 'IM'
  ss_epi@meta.data[ss_epi@meta.data$seurat_clusters %in% c(8),]$cell_type <- 'Mixed'
  ss_epi@meta.data[ss_epi@meta.data$seurat_clusters %in% c(9),]$cell_type <- 'L2'
  ss_epi@meta.data[ss_epi@meta.data$seurat_clusters %in% c(3,4,11),]$cell_type <- 'Tumor-L2'
  table(ss_epi$seurat_clusters, ss_epi$cell_type)
  
  #### Fig 1B ####
  g0 <- DimPlot(ss_epi, group.by = 'cell_type', ncol = 1, cols = celltype_palette_list)

  #### Fig 1C ####
  pl <- lapply( names(condition_palette), function(x){
    ss_cond <- subset(ss_epi, condition==x)
    g1 <- DimPlot(ss_cond, group.by = c('condition','cell_type'), reduction = 'umap', cols = celltype_palette_list) + ggtitle(x)
    g1[[2]]
  })
  do.call(grid.arrange, c(pl, nrow=1))
  
  #### Fig 1E ####
  f <- c("S.Score","G2M.Score")
  g <- ViolinPlotEnhance(ss_epi,
                          features = f,
                          assay = 'RNA',
                          slot = 'data',
                          pt.size = 0,
                          flip = T,
                          group.by = 'seurat_clusters',
                          stack = T,
                          fill.by = 'ident', 
                          y_axis_label = "Normalized scores",
                          title = gsub("_"," ","cell cycle scores"),
                          #pal = brewer.reds(20)
                          pal = pals::coolwarm(25),
                          x_axis_text_angle = 90, 
                          strip.text.size = 14
  ) + xlab("")
  g
  
  #### Fig 1F ####
  g <- lapply(names(violin_markers), function(x){
    ViolinPlotEnhance(ss_epi,
                      features = violin_markers[[x]],
                      assay = 'RNA',
                      slot = 'data',
                      pt.size = 0,
                      flip = T,
                      group.by = 'seurat_clusters',
                      stack = T,
                      fill.by = 'ident',
                      title = gsub("_"," ",x), 
                      x_axis_text_angle = 0,
                      pal = pals::coolwarm(25),
                      ) + xlab("") + theme(legend.position = 'none')
  })
  do.call(grid.arrange, c(g, nrow=2))
  
  #### Save Seurat object ####
  saveRDS(ss_epi, file = epi.scrna.object)
  
  #### Load Seurat object ####
  ss_epi <- readRDS(file = epi.scrna.object)
  
  
  #---------------------------------------------------------------------#
  #### DE Markers (1-vs-all) ####
  #---------------------------------------------------------------------#
  # bionames
  DEsample_v_others <- wilcoxauc(erg.scRNA, group_by = 'bioNames', seurat_assay = 'RNA', assay = 'data')
  DEsample_v_others <- setDT(DEsample_v_others)
  setnames(DEsample_v_others,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DEsample_v_others[,cluster2 := 'others']
  DEsample_v_others <- DEsample_v_others[order(statistic, decreasing = T)]
  DEsample_v_others[, head(.SD, 10), by=cluster1]
  
  # clusters
  grep("SCT_snn_res", colnames(erg.scRNA@meta.data), value = T)
  Idents(erg.scRNA) #SCT_snn_res_1.0
  DEcluster_v_others <- wilcoxauc(erg.scRNA, group_by = 'seurat_clusters', seurat_assay = 'RNA', assay = 'data')
  DEcluster_v_others <- setDT(DEcluster_v_others)
  setnames(DEcluster_v_others,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DEcluster_v_others[,cluster2 := 'others']
  DEcluster_v_others <- DEcluster_v_others[order(statistic, decreasing = T)]
  View(DEcluster_v_others[, head(.SD, 5), by=cluster1])
  
  # conditions
  DEcondition_v_others <- wilcoxauc(erg.scRNA, group_by = 'condition', seurat_assay = 'RNA', assay = 'data')
  DEcondition_v_others <- setDT(DEcondition_v_others)
  setnames(DEcondition_v_others,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DEcondition_v_others[,cluster2 := 'others']
  DEcondition_v_others <- DEcondition_v_others[order(statistic, decreasing = T)]
  DEcondition_v_others[, head(.SD, 5), by=cluster1]
  
  # lum_bas_int_classification
  DElum_bas_int_classification_v_others <- wilcoxauc(erg.scRNA, group_by = 'lum_bas_int_classification', seurat_assay = 'RNA', assay = 'data')
  DElum_bas_int_classification_v_others <- setDT(DElum_bas_int_classification_v_others)
  setnames(DElum_bas_int_classification_v_others,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DElum_bas_int_classification_v_others[,cluster2 := 'others']
  DElum_bas_int_classification_v_others <- DElum_bas_int_classification_v_others[order(statistic, decreasing = T)]
  DElum_bas_int_classification_v_others[, head(.SD, 10), by=cluster1]
  
  
  ## Add gene annotation and sort outputs
  gene_annot <- fread(file.path('/data/projects/sawyers/erg/analysis/scRNA/3months/seurat/data/GRCm39_gene_annotation.txt'), col.names = c("geneID","gene","description"))
  gene_annot <- gene_annot[, lapply(.SD, paste0, collapse=","), .SDcols='description', by=gene]
  
  DEsample_v_others <- merge(DEsample_v_others, gene_annot, by='gene', all.x=T)
  DEsample_v_others <- DEsample_v_others[order(cluster1, statistic, decreasing = T)]
  #DEsample_v_others[,rank:=avg_logFC*10*exp(1-p_val_adj)]
  DEsample_v_others <- DEsample_v_others[order(rank, decreasing = T)]
  DEsample_v_others[,rank:= avg_logFC*log(statistic)]
  plot(DEsample_v_others[cluster1=='Wildtype 3 months' ,rank],
       DEsample_v_others[cluster1=='Wildtype 3 months' ,auc])
  
  
  DEcluster_v_others <- merge(DEcluster_v_others, gene_annot, by='gene', all.x=T)
  DEcluster_v_others <- DEcluster_v_others[order(cluster1, statistic, decreasing = T)]
  DEcluster_v_others[,rank:= avg_logFC*log(statistic)]
  DEcluster_v_others <- DEcluster_v_others[order(rank, decreasing = T)]
  plot(DEcluster_v_others[cluster1=='1' ,rank],
       DEcluster_v_others[cluster1=='1' ,auc])
  
  DEcondition_v_others <- merge(DEcondition_v_others, gene_annot, by='gene', all.x=T)
  DEcondition_v_others <- DEcondition_v_others[order(cluster1, statistic, decreasing = T)]
  #DEcondition_v_others[,rank:=avg_logFC*10*exp(1-p_val_adj)]
  DEcondition_v_others[,rank:= avg_logFC*log(statistic)]
  DEcondition_v_others <- DEcondition_v_others[order(rank, decreasing = T)]
  plot(DEcondition_v_others[cluster1=='Pten-' ,rank],
       DEcondition_v_others[cluster1=='Pten-' ,auc])
  
  
  DElum_bas_int_classification_v_others <- merge(DElum_bas_int_classification_v_others, gene_annot, by='gene', all.x=T)
  DElum_bas_int_classification_v_others <- DElum_bas_int_classification_v_others[order(cluster1, statistic, decreasing = T)]
  DElum_bas_int_classification_v_others[,rank:=avg_logFC*10*exp(1-p_val_adj)]
  DElum_bas_int_classification_v_others <- DElum_bas_int_classification_v_others[order(rank, decreasing = T)]
  plot(DElum_bas_int_classification_v_others[cluster1=='Intermediate' ,rank],
       DElum_bas_int_classification_v_others[cluster1=='Intermediate' ,auc])
  
  
  DEcell_type_v_others <- merge(DEcell_type_v_others, gene_annot, by='gene', all.x=T)
  DEcell_type_v_others <- DEcell_type_v_others[order(cluster1, statistic, decreasing = T)]
  DEcell_type_v_others[,rank:= avg_logFC*log(statistic)]
  DEcell_type_v_others <- DEcell_type_v_others[order(rank, decreasing = T)]
  plot(DEcell_type_v_others[cluster1=='IM' ,rank],
       DEcell_type_v_others[cluster1=='IM' ,auc])
  
  
  
  #### DE spreadsheet ####
  ## One vs all ##
  # write
  writexl::write_xlsx(x = DEsample_v_others, path = file.path(spreadsheets.dir,'scRNA_sample_vs_others_wilcox_de.xlsx'))
  writexl::write_xlsx(x = DEcluster_v_others, path = file.path(spreadsheets.dir,'scRNA_cluster_vs_others_wilcox_de.xlsx'))
  writexl::write_xlsx(x = DEcondition_v_others, path = file.path(spreadsheets.dir,'scRNA_condition_vs_others_wilcox_de.xlsx'))
  writexl::write_xlsx(x = DElum_bas_int_classification_v_others, path = file.path(spreadsheets.dir,'scRNA_lum_bas_int_classification_v_others_wilcox_de.xlsx'))
  
  # load
  gene_annot <- fread(file.path('/home/erik/annotation/transcriptome/tables/GRCm39_gene_annotation.txt'), col.names = c("geneID","gene","description"))
  gene_annot <- gene_annot[, lapply(.SD, paste0, collapse=","), .SDcols='description', by=gene]
  
  DEsample_v_others <- setDT(readxl::read_xlsx(file.path(spreadsheets.dir,'scRNA_sample_vs_others_wilcox_de.xlsx')))
  DEcluster_v_others <- setDT(readxl::read_xlsx(file.path(spreadsheets.dir,'scRNA_cluster_vs_others_wilcox_de.xlsx')))
  DEcondition_v_others <- setDT(readxl::read_xlsx(file.path(spreadsheets.dir,'scRNA_condition_vs_others_wilcox_de.xlsx')))
  DElum_bas_int_classification_v_others <- setDT(readxl::read_xlsx(file.path(spreadsheets.dir,'scRNA_lum_bas_int_classification_v_others_wilcox_de.xlsx')))
  
  rem <- grep("^Gm[0-9]+$",rownames(ss_epi), value=T)
  rem <- c(rem, grep("^mt-",rownames(ss_epi), value=T, ignore.case = T))
  rem <- c(rem, grep("^[0-9]",rownames(ss_epi), value=T, ignore.case = T))
  ribosomal_genes = fread("/home/erik/annotation/genome/tables/mouse_ribosomal_genes.txt", header = F)
  ribosomal_genes <- paste(paste0("^",unname(unlist(ribosomal_genes))), collapse = "|")
  rem <- c(rem, grep(ribosomal_genes, rownames(ss_epi), value=T))
  
  cairo_pdf(filename = file.path(plot.dir, 'Sample_DEupregulated_barplot.pdf'), width = 5, height = 3, bg = 'transparent')
  # dt <- DEsample_v_others[!gene %in% rem & p_val_adj < 0.01 & avg_logFC > 0.1, unique(gene), by=cluster1][,.N,by=cluster1]
  dt <- DEsample_v_others[p_val_adj < 0.01 & avg_logFC > 0.1, unique(gene), by=cluster1][,.N,by=cluster1]
  ggplot(dt,
         aes(x=reorder(cluster1, N), y=N, fill=cluster1)) +
    geom_col() + xlab("") + ylab("Number of genes") +
    ggtitle(label = "Signifcant upregulated genes", subtitle = "adjusted pval < 0.01 & LogFC > 0.1") +
    scale_fill_manual("",values = sample_palette) + theme_classic(14) +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1),
          panel.grid = element_blank(),
          legend.background = element_blank(),
          legend.position = 'none',
          panel.background = element_blank(),
          plot.title = element_text(hjust=1), plot.subtitle = element_text(hjust=1),
          plot.background = element_blank()) +
    lemon::coord_capped_flip(expand = F)
  dev.off()
  
  #---------------------------------------------------------------------#
  #### DE epi sample (1-vs-all) ####
  #---------------------------------------------------------------------#
  gene_annot <- fread(file.path(data.dir,'GRCm39_gene_annotation.txt'), col.names = c("geneID","gene","description"))
  gene_annot <- gene_annot[, lapply(.SD, paste0, collapse=","), .SDcols='description', by=gene]
  rem <- grep("^Gm[0-9]+$",rownames(ss_epi), value=T)
  rem <- c(rem, grep("^mt-",rownames(ss_epi), value=T, ignore.case = T))
  rem <- c(rem, grep("^[0-9]",rownames(ss_epi), value=T, ignore.case = T))
  ribosomal_genes = fread("/home/erik/annotation/genome/tables/mouse_ribosomal_genes.txt", header = F)
  ribosomal_genes <- paste(paste0("^",unname(unlist(ribosomal_genes))), collapse = "|")
  rem <- c(rem, grep(ribosomal_genes, rownames(ss_epi), value=T))
  
  ##### BioNames #####
  DEsample_v_others_epi <- wilcoxauc(ss_epi, group_by = 'bioNames', seurat_assay = 'RNA', assay = 'data')
  DEsample_v_others_epi <- setDT(DEsample_v_others_epi)
  setnames(DEsample_v_others_epi,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DEsample_v_others_epi[,cluster2 := 'others']
  DEsample_v_others_epi <- DEsample_v_others_epi[order(statistic, decreasing = T)]
  DEsample_v_others_epi[, head(.SD, 10), by=cluster1]
  DEsample_v_others_epi[,rank:=avg_logFC*10*exp(1-p_val_adj)]
  DEsample_v_others_epi <- merge(DEsample_v_others_epi, gene_annot, by='gene', all.x=T)
  DEsample_v_others_epi <- DEsample_v_others_epi[order(rank, decreasing = T)]
  
  writexl::write_xlsx(x = DEsample_v_others_epi, 
                      path = file.path(spreadsheets.dir,'scRNA_episample_vs_others_wilcox_de.xlsx'))
  
  ##### Epithelial clusters #####
  # 2,4,7,14,15,17,18,19,20,21
  ##epierg.scRNA = subset(erg.scRNA, subset = seurat_clusters %in% epi_clusters)
  table(ss_epi$seurat_clusters)
  DEepicluster_v_others <- wilcoxauc(ss_epi, group_by = 'seurat_clusters', seurat_assay = 'RNA', assay = 'data')
  DEepicluster_v_others <- setDT(DEepicluster_v_others)
  setnames(DEepicluster_v_others,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DEepicluster_v_others[,cluster2 := 'others']
  DEepicluster_v_others <- DEepicluster_v_others[order(statistic, decreasing = T)]
  table(DEepicluster_v_others$cluster1)
  
  DEepicluster_v_others <- merge(DEepicluster_v_others, gene_annot, by='gene', all.x=T)
  DEepicluster_v_others <- DEepicluster_v_others[order(cluster1, statistic, decreasing = T)]
  DEepicluster_v_others[,rank:= avg_logFC*log(statistic)]
  DEepicluster_v_others <- DEepicluster_v_others[order(rank, decreasing = T)]
  
  writexl::write_xlsx(x = DEepicluster_v_others, path = file.path(spreadsheets.dir,'scRNA_epicluster_vs_others_wilcox_de.xlsx'))
  DEepicluster_v_others <- setDT(readxl::read_xlsx(file.path(spreadsheets.dir,'scRNA_epicluster_vs_others_wilcox_de.xlsx')))
  
  ##### Celltype #####
  DEcell_type_v_others <- wilcoxauc(ss_epi, group_by = 'cell_type', seurat_assay = 'RNA', assay = 'data')
  DEcell_type_v_others <- setDT(DEcell_type_v_others)
  setnames(DEcell_type_v_others,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DEcell_type_v_others[,cluster2 := 'others']
  DEcell_type_v_others <- DEcell_type_v_others[order(statistic, decreasing = T)]
  DEcell_type_v_others[, head(.SD, 5), by=cluster1]
  
  writexl::write_xlsx(x = DEcell_type_v_others, path = file.path(spreadsheets.dir,'scRNA_cell_type_v_others_wilcox_de.xlsx'))
  
  DEcell_type_v_others <- setDT(readxl::read_xlsx(file.path(spreadsheets.dir,'scRNA_cell_type_v_others_wilcox_de.xlsx')))
  table(DEcell_type_v_others$cluster1)
  DEcell_type_v_others[cluster1=='IM']
  topgenes <- DEcell_type_v_others[cluster1=='IM', head(gene, 20)]
  gu <- lapply(topgenes, function(f){
    FeaturePlot(ss_epi, features = f, cols = pals::plasma(20))
  }) %>% set_names(topgenes)
  do.call(grid.arrange, c(gu, nrow=4))
  
  #---------------------------------------------------------------------#
  #### DE epi Markers (pairwise) ####
  #---------------------------------------------------------------------#
  
  comb <- combn(x = as.integer(unique(ss_epi$seurat_clusters)), m = 2)
  comb
  DEcluster_v_cluster <- lapply(seq_len(ncol(comb)), function(i){
    dt <- wilcoxauc(ss_epi, group_by = 'seurat_clusters', groups_use = c(comb[,i]), seurat_assay = 'RNA', assay = 'data')
    dt <- setDT(dt)
    c1 <- comb[1,i]
    c2 <- comb[2,i]
    setnames(dt,
             old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
             new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
    dt[,cluster2 := ifelse(cluster1==c1, c2, c1)]
    dt <- dt[order(statistic, decreasing = T)]
    dt[,rank:=avg_logFC*10*exp(1-p_val_adj)]
    dt <- merge(dt, gene_annot, by='gene', all.x=T)
    dt[order(rank, decreasing = T)]
  })
  names(DEcluster_v_cluster) <- apply(comb, 2, function(x){ paste0(x[1],"_vs_",x[2]) })
  
  sapply(DEcluster_v_cluster, nrow)
  
  DEcluster_v_cluster$`13_vs_26`[, head(.SD, 10), by=cluster1]
  
  writexl::write_xlsx(x = DEcluster_v_cluster, path = file.path(spreadsheets.dir,'scRNA_epicluster_vs_epicluster_wilcox_de.xlsx'))
  
  DEcluster_v_cluster <- setDT(read_xlsx(path = file.path(spreadsheets.dir,'scRNA_epicluster_vs_epicluster_wilcox_de.xlsx')))
  
  #---------------------------------------------------------------------#
  #### DE epi Markers (condition within cluster) ####
  #---------------------------------------------------------------------#
  # clusters
  gene_annot <- fread(file.path(data.dir,'GRCm39_gene_annotation.txt'), col.names = c("geneID","gene","description"))
  gene_annot <- gene_annot[, lapply(.SD, paste0, collapse=","), .SDcols='description', by=gene]
  DEepi_condition_within_cluster <- lapply(epi_clusters, function(i){
    ss <- subset(erg.scRNA, subset = seurat_clusters == i)
    dt <- setDT(wilcoxauc(ss, group_by = 'condition', seurat_assay = 'RNA', assay = 'data'))
    setnames(dt,
             old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
             new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
    dt[,epi_cluster:=i]
    dt[,cluster2 := 'others']
    dt[,rank:=avg_logFC*10*exp(1-p_val_adj)]
    dt <- merge(dt, gene_annot, by='gene', all.x=T)
    dt[order(rank, decreasing = T)]
  })
  names(DEepi_condition_within_cluster) <- epi_clusters
  
  DEepi_condition_within_cluster$`2`[, head(.SD, 10), by=cluster1]
  
  writexl::write_xlsx(x = DEepi_condition_within_cluster, path = file.path(spreadsheets.dir,'scRNA_condition_within_epicluster_wilcox_de.xlsx'))
  
  rem <- grep("^Gm[0-9]+$",rownames(erg.scRNA), value=T)
  rem <- c(rem, grep("^mt-",rownames(erg.scRNA), value=T, ignore.case = T))
  rem <- c(rem, grep("^[0-9]",rownames(erg.scRNA), value=T, ignore.case = T))
  ribosomal_genes = fread("/home/erik/annotation/genome/tables/mouse_ribosomal_genes.txt", header = F)
  ribosomal_genes <- paste(paste0("^",unname(unlist(ribosomal_genes))), collapse = "|")
  rem <- c(rem, grep(ribosomal_genes, rownames(erg.scRNA), value=T))
  
  #---------------------------------------------------------------------#
  #### DE epi Markers (per each condition one vs others clusters) ####
  #---------------------------------------------------------------------#
  condition <- unique(ss_epi$condition)
  gene_annot <- fread(file.path(data.dir,'GRCm39_gene_annotation.txt'), col.names = c("geneID","gene","description")) %>% 
    .[, lapply(.SD, paste0, collapse=","), .SDcols='description', by=gene]
  DEepi_percondition_cluster_vs_others <- lapply(conditions, function(j){
    ss <- subset(epierg.scRNA, subset = condition == j)
    dt <- wilcoxauc(ss, group_by = 'seurat_clusters', seurat_assay = 'RNA', assay = 'data')
    dt <- setDT(dt)
    setnames(dt,
             old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
             new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
    dt[,cluster2 := 'others']
    dt[,condition:=j]
    dt <- dt[order(statistic, decreasing = T)]
    dt[,rank:=avg_logFC*10*exp(1-p_val_adj)]
    dt <- merge(dt, gene_annot, by='gene', all.x=T)
    dt[order(rank, decreasing = T)]
  })
  names(DEepi_percondition_cluster_vs_others) <- conditions
  
  writexl::write_xlsx(x = DEepi_percondition_cluster_vs_others, path = file.path(spreadsheets.dir,'scRNA_percondition_one_epicluster_vs_others_wilcox_de.xlsx'))
  
  
  #---------------------------------------------------------------------#
  #### DE epi Markers (per each celltype one vs others condition) ####
  # ## i.e. for all IM cells compare Erg+ vs Others
  #---------------------------------------------------------------------#
  celltypes <- unique(ss_epi$cell_type)
  gene_annot <- fread(file.path(data.dir,'GRCm39_gene_annotation.txt'), col.names = c("geneID","gene","description")) %>% 
    .[, lapply(.SD, paste0, collapse=","), .SDcols='description', by=gene]
  DEepi_percelltypes_vs_others <- lapply(celltypes, function(j){
    ss <- subset(ss_epi, subset = cell_type == j)
    dt <- wilcoxauc(ss, group_by = 'condition', seurat_assay = 'RNA', assay = 'data')
    dt <- setDT(dt)
    setnames(dt,
             old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
             new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
    dt[,cluster2 := 'others']
    dt[,cell_type:=j]
    dt <- dt[order(statistic, decreasing = T)]
    dt[,rank:=avg_logFC*10*exp(1-p_val_adj)]
    dt <- merge(dt, gene_annot, by='gene', all.x=T)
    dt[order(rank, decreasing = T)]
  })
  names(DEepi_percelltypes_vs_others) <- celltypes
  View(DEepi_percelltypes_vs_others$IM)
  
  writexl::write_xlsx(x = DEepi_percelltypes_vs_others, path = file.path(spreadsheets.dir,'scRNA_percelltype_one_condition_vs_others_wilcox_de.xlsx'))
  DEepi_percelltypes_vs_others <- read_xlsx(path = file.path(spreadsheets.dir,'scRNA_percelltype_one_condition_vs_others_wilcox_de.xlsx'))
  setDT(DEepi_percelltypes_vs_others)
  table(DEepi_percelltypes_vs_others$cell_type)
  
  #---------------------------------------------------------------------#
  #### DE epi Markers (per each celltype pairwise by condition) ####
  ## i.e. for all IM cells compare Erg+ vs Pten-
  #---------------------------------------------------------------------#
  celltypes <- unique(ss_epi$cell_type)
  gene_annot <- fread(file.path(data.dir,'GRCm39_gene_annotation.txt'), col.names = c("geneID","gene","description")) %>% 
    .[, lapply(.SD, paste0, collapse=","), .SDcols='description', by=gene]
  DEepi_percelltypes_pairwise_condition <- lapply(celltypes, function(j){
    ss <- subset(ss_epi, subset = cell_type == j)
    comb <- combn(x = unique(ss$condition), m = 2)
    DTL <- lapply(seq_len(ncol(comb)), function(i){
      dt <- wilcoxauc(ss, group_by = 'condition', groups_use = c(comb[,i]), seurat_assay = 'RNA', assay = 'data')
      dt <- setDT(dt)
      c1 <- comb[1,i]
      c2 <- comb[2,i]
      setnames(dt,
               old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
               new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
      dt[,cluster2 := ifelse(cluster1==c1, c2, c1)]
      dt[,cell_type:=j]
      dt <- dt[order(statistic, decreasing = T)]
      dt[,rank:=avg_logFC*10*exp(1-p_val_adj)]
      dt <- merge(dt, gene_annot, by='gene', all.x=T)
      dt[order(rank, decreasing = T)]
    }) %>% rbindlist
  })
  names(DEepi_percelltypes_pairwise_condition) <- celltypes
  View(DEepi_percelltypes_pairwise_condition$IM)
  
  celltypes
  writexl::write_xlsx(x = DEepi_percelltypes_pairwise_condition, path = file.path(spreadsheets.dir,'scRNA_percelltype_condition_pairwise_wilcox_de.xlsx'))
  
  
  
  
  #---------------------------------------------------------------------#
  #### DE epi Markers (per each condition pairwise cluster) ####
  #---------------------------------------------------------------------#
  # clusters
  conditions <- unique(epierg.scRNA$condition)
  gene_annot <- fread(file.path(data.dir,'GRCm39_gene_annotation.txt'), 
                      col.names = c("geneID","gene","description")) %>% 
    .[, lapply(.SD, paste0, collapse=","), .SDcols='description', by=gene]
  DEepi_condition_between_cluster <- lapply(conditions, function(j){
    ss <- subset(epierg.scRNA, subset = condition == j)
    ##
    comb <- combn(x = as.integer(unique(ss$seurat_clusters)), m = 2)
    DTL <- lapply(seq_len(ncol(comb)), function(i){
      dt <- wilcoxauc(ss, group_by = 'seurat_clusters', groups_use = c(comb[,i]), seurat_assay = 'RNA', assay = 'data')
      dt <- setDT(dt)
      c1 <- comb[1,i]
      c2 <- comb[2,i]
      setnames(dt,
               old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
               new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
      dt[,cluster2 := ifelse(cluster1==c1, c2, c1)]
      dt[,condition:=j]
      dt <- dt[order(statistic, decreasing = T)]
      dt[,rank:=avg_logFC*10*exp(1-p_val_adj)]
      dt <- merge(dt, gene_annot, by='gene', all.x=T)
      dt[order(rank, decreasing = T)]
    }) %>% rbindlist
    ##
  })
  names(DEepi_condition_between_cluster) <- conditions
  
  DEepi_condition_between_cluster$`Erg+Pten-`[, head(.SD, 10), by=cluster1]
  
  lapply( names(DEepi_condition_between_cluster), function(i){
    cluster1 <- unique(DEepi_condition_between_cluster[[i]]$cluster1)
    for(j in cluster1)
      writexl::write_xlsx(x = DEepi_condition_between_cluster[[i]][cluster1==j],
                          path = file.path(spreadsheets.dir,
                                           paste0('scRNA_epi_condition_',i,'_cluster',j,'_between_cluster_wilcox_de.xlsx')))
  })
  
  rem <- grep("^Gm[0-9]+$",rownames(erg.scRNA), value=T)
  rem <- c(rem, grep("^mt-",rownames(erg.scRNA), value=T, ignore.case = T))
  rem <- c(rem, grep("^[0-9]",rownames(erg.scRNA), value=T, ignore.case = T))
  ribosomal_genes = fread("/home/erik/annotation/genome/tables/mouse_ribosomal_genes.txt", header = F)
  ribosomal_genes <- paste(paste0("^",unname(unlist(ribosomal_genes))), collapse = "|")
  rem <- c(rem, grep(ribosomal_genes, rownames(erg.scRNA), value=T))
  
  #---------------------------------------------------------------------#
  #### DE epi Markers (for each cluster one vs other conditions) ####
  #---------------------------------------------------------------------#
  condition <- unique(ss_epi$condition)
  clusters <- unique(ss_epi$seurat_clusters)
  gene_annot <- fread(file.path(data.dir,'GRCm39_gene_annotation.txt'), col.names = c("geneID","gene","description")) %>% 
    .[, lapply(.SD, paste0, collapse=","), .SDcols='description', by=gene]
  DEepi_percluster_condition_vs_others <- lapply(clusters, function(j){
    ss <- subset(ss_epi, subset = seurat_clusters == j)
    if( length(unique(ss$condition)) > 1 ){
      dt <- wilcoxauc(ss, group_by = 'condition', seurat_assay = 'RNA', assay = 'data')
      dt <- setDT(dt)
      setnames(dt,
               old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
               new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
      dt[,cluster2 := 'others']
      dt[,cluster:=j]
      dt <- dt[order(statistic, decreasing = T)]
      dt[,rank:=avg_logFC*10*exp(1-p_val_adj)]
      dt <- merge(dt, gene_annot, by='gene', all.x=T)
      dt[order(rank, decreasing = T)]
    }
    #
  })
  names(DEepi_percluster_condition_vs_others) <- paste0("cluster_",clusters)
  rem <- which(sapply(DEepi_percluster_condition_vs_others, is.null))
  if(length(rem) > 0)
    DEepi_percluster_condition_vs_others <- DEepi_percluster_condition_vs_others[-rem]
  
  writexl::write_xlsx(x = DEepi_percluster_condition_vs_others, 
                      path = file.path(spreadsheets.dir,'scRNA_per_epicluster_condition_one_vs_others_wilcox_de.xlsx'))
  
  #### DE sample vs EV ####
  gene_annot <- fread(file.path(data.dir,'GRCm39_gene_annotation.txt'), col.names = c("geneID","gene","description"))
  gene_annot <- gene_annot[, lapply(.SD, paste0, collapse=","), .SDcols='description', by=gene]
  # bionames
  DEsample_v_EV <- lapply(grep("^Foxa1", unique(erg.scRNA@meta.data$bioNames), value = T), function(i){
    DT <- wilcoxauc(erg.scRNA, group_by = 'bioNames', groups_use = c(i, 'Empty Vector'), seurat_assay = 'RNA', assay = 'data')
    DT <- setDT(DT)
    setnames(DT,
             old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
             new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
    DT[,cluster2 := ifelse(cluster1=='Empty Vector',i,'Empty Vector')]
    DT[order(statistic, decreasing = T)]
  }) %>% rbindlist()
  DEsample_v_EV <- merge(DEsample_v_EV, gene_annot, by='gene', all.x=T)
  DEsample_v_EV[,rank:=avg_logFC*10*exp(1-p_val_adj)]
  plot(DEsample_v_EV[cluster1=='Foxa1 WildType', rank], DEsample_v_EV[cluster1=='Foxa1 WildType', auc])
  writexl::write_xlsx(x = DEsample_v_EV, path = file.path(spreadsheets.dir,'scRNA_bioNames_vs_EV_wilcox_de.xlsx'))
  
  DEsample_v_EV[p_val_adj < 0.01 & avg_logFC > 0, .N, by=list(cluster1, cluster2)]
  cairo_pdf(filename = file.path(plot.dir, 'Sample_vs_EV_DEupregulated_barplot.pdf'), width = 5, height = 5, bg = 'transparent')
  DT <- DEsample_v_EV[p_val_adj < 0.01, .N, by=list(cluster1, cluster2, sign(avg_logFC))][cluster1 != 'Empty Vector'][order(cluster1)]
  ggplot(DT, aes(x=cluster1, y=N, fill=factor(sign, levels = c(1,-1), labels = c('up','down')))) +
    geom_col(width = .60, position = position_dodge(width = 0.65)) +
    xlab("") + ylab("Number of genes") + ggtitle(label = "Signifcant differential genes vs Empty Vector", subtitle = "adjusted pval < 0.01 & LogFC > 0") +
    #scale_fill_manual("",values = sample_palette) +
    scale_fill_manual("",values = c('red','blue')) +
    theme_classic(14) +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1),
          panel.grid = element_blank(),
          legend.background = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5),
          plot.background = element_blank()) +
    lemon::coord_capped_cart(expand = F)
  dev.off()
  #---------------------------------------------------------------------#
  
  
  #### DE Luminal Basal ####
  ## Add gene annotation and sort outputs
  gene_annot <- fread(file.path(data.dir,'GRCm39_gene_annotation.txt'), col.names = c("geneID","gene","description"))
  gene_annot <- gene_annot[, lapply(.SD, paste0, collapse=","), .SDcols='description', by=gene]
  
  # luminal basal
  DEluminal_v_basal <- wilcoxauc(erg.scRNA, group_by = 'luminal_basal_classification', seurat_assay = 'RNA', assay = 'data')
  DEluminal_v_basal <- setDT(DEluminal_v_basal)
  setnames(DEluminal_v_basal,
           old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
           new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
  DEluminal_v_basal[,cluster2 := 'others']
  DEluminal_v_basal <- DEluminal_v_basal[order(statistic, decreasing = T)]
  DEluminal_v_basal[,rank:=avg_logFC*10*exp(1-p_val_adj)]
  DEluminal_v_basal <- merge(DEluminal_v_basal, gene_annot, by='gene', all.x=T)
  DEluminal_v_basal <- DEluminal_v_basal[order(cluster1, statistic, decreasing = T)]
  
  # luminal vs basal within each sample
  DEluminal_v_basal_per_sample <- lapply(unique(erg.scRNA@meta.data$bioNames), function(i){
    ss <- subset(erg.scRNA, subset = bioNames == i)
    dt <- setDT(wilcoxauc(ss, group_by = 'luminal_basal_classification', seurat_assay = 'RNA', assay = 'data'))
    setnames(dt,
             old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
             new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
    dt[,cluster2 := 'others']
    dt[, bioNames := i]
    dt <- dt[order(statistic, decreasing = T)]
  }) %>% rbindlist
  DEluminal_v_basal_per_sample[,rank:=avg_logFC*10*exp(1-p_val_adj)]
  plot(rank ~ auc, DEluminal_v_basal_per_sample[bioNames=='Foxa1 WildType'])
  # Luminal vs basal within each sample
  DEluminal_v_basal_per_sample <- merge(DEluminal_v_basal_per_sample, gene_annot, by='gene', all.x=T)
  DEluminal_v_basal_per_sample <- DEluminal_v_basal_per_sample[order(cluster1, statistic, decreasing = T)]
  
  
  
  # Write spreadsheets
  writexl::write_xlsx(x = DEluminal_v_basal, path = file.path(spreadsheets.dir,'scRNA_luminal_vs_basal_wilcox_de.xlsx'))
  writexl::write_xlsx(x = DEluminal_v_basal_per_sample, path = file.path(spreadsheets.dir,'scRNA_luminal_vs_basal_per_sample_wilcox_de.xlsx'))
  writexl::write_xlsx(x = DEluminal_v_basal_per_cluster, path = file.path(spreadsheets.dir,'scRNA_luminal_vs_basal_per_cluster_wilcox_de.xlsx'))
  
  # Read spreadsheets
  DEluminal_v_basal <- setDT(readxl::read_xlsx(file.path(spreadsheets.dir,'scRNA_luminal_vs_basal_wilcox_de.xlsx')))
  DEluminal_v_basal_per_sample <- setDT(readxl::read_xlsx(file.path(spreadsheets.dir,'scRNA_luminal_vs_basal_per_sample_wilcox_de.xlsx')))
  DEluminal_v_basal_per_cluster <- setDT(readxl::read_xlsx(file.path(spreadsheets.dir,'scRNA_luminal_vs_basal_per_cluster_wilcox_de.xlsx')))
  
  
  # luminal vs basal within each cluster
  DEluminal_v_basal_per_cluster <- lapply(unique(erg.scRNA@meta.data$seurat_clusters), function(i){
    ss <- subset(erg.scRNA, subset = seurat_clusters == i)
    dt <- setDT(wilcoxauc(ss, group_by = 'luminal_basal_classification', seurat_assay = 'RNA', assay = 'data'))
    setnames(dt,
             old=c("feature","group","logFC","pct_in","pct_out","pval","padj"),
             new = c("gene","cluster1","avg_logFC","percent_cluster1","percent_cluster2","p_val","p_val_adj"))
    dt[,cluster2 := 'others']
    dt[, seurat_cluster := i]
    dt <- dt[order(statistic, decreasing = T)]
  }) %>% rbindlist
  DEluminal_v_basal_per_cluster[,rank:=avg_logFC*10*exp(1-p_val_adj)]
  DEluminal_v_basal_per_cluster <- merge(DEluminal_v_basal_per_cluster, gene_annot, by='gene', all.x=T)
  DEluminal_v_basal_per_cluster <- DEluminal_v_basal_per_cluster[order(cluster1, statistic, decreasing = T)]
  
  ### Score AR target markers ####
  ar_targets_combine <- unique(ar_targets_combine <- unlist(ar_targets_list, use.names = F))
  ss_epi <- AddModuleScore(object = ss_epi,
                           features = list(ar_targets_combine),
                           assay = 'RNA',
                           name = 'ar_targets', # appends a 1
                           search = F) # searches HUGO contains human gene names only
  colnames(ss_epi@meta.data)[which(colnames(ss_epi@meta.data)=='ar_targets1')] <- 'ar_targets'
  
  ss_epi <- AddModuleScore(object = ss_epi,
                           features = list(Nelson_response_ANDR_UP_mouse),
                           assay = 'RNA',
                           name = 'Nelson_response_ANDR_UP_mouse', # appends a 1
                           search = F) # searches HUGO contains human gene names only
  colnames(ss_epi@meta.data)[which(colnames(ss_epi@meta.data)=='Nelson_response_ANDR_UP_mouse1')] <- 'Nelson_response_ANDR_UP_mouse'
  
  ss_epi <- AddModuleScore(object = ss_epi,
                           features = list(Nelson_response_no_NEG_ANDR_UP_mouse),
                           assay = 'RNA',
                           name = 'Nelson_response_no_NEG_ANDR_UP_mouse', # appends a 1
                           search = F) # searches HUGO contains human gene names only
  colnames(ss_epi@meta.data)[which(colnames(ss_epi@meta.data)=='Nelson_response_no_NEG_ANDR_UP_mouse1')] <- 'Nelson_response_no_NEG_ANDR_UP_mouse'
  
  ss_epi <- AddModuleScore(object = ss_epi,
                           features = list(Hieronymus_Androgen_mouse),
                           assay = 'RNA',
                           name = 'Hieronymus_Androgen_mouse', # appends a 1
                           search = F) # searches HUGO contains human gene names only
  colnames(ss_epi@meta.data)[which(colnames(ss_epi@meta.data)=='Hieronymus_Androgen_mouse1')] <- 'Hieronymus_Androgen_mouse'
  
  
  #### AR violin by condition ####
  g0 <- VlnPlot(ss_epi, features = 'Ar', assay = 'RNA', slot='data', log = F,
                pt.size = 0, flip = T, group.by = 'condition', stack = F, cols = condition_palette) + xlab("")
  df <- g0$data
  g0 <- ggplot(df, aes(x=ident, y=Ar, fill=ident)) +
    geom_violin(trim = T, scale = 'width', draw_quantiles = c(.25,.5,.75)) +
    scale_fill_manual("", values = condition_palette) +
    xlab("") + ylab("Expression level") + ggtitle("Androgen Receptor (Ar)") +
    scale_y_continuous(limits = c(0,4), expand = c(0,0)) +
    theme_classic(14) + theme(plot.title = element_text(face = 'bold', hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  g1 <- VlnPlot(ss_epi, features = 'ar_targets', assay = 'RNA', slot='data',
                pt.size = 0, flip = T, group.by = 'condition', stack = F, cols = condition_palette) + xlab("")
  df <- g1$data
  g1 <- ggplot(df, aes(x=ident, y=ar_targets, fill=ident)) +
    geom_violin(trim = T, scale = 'width', draw_quantiles = c(.25,.5,.75)) +
    scale_fill_manual("", values = condition_palette) +
    xlab("") + ylab("Expression level") + ggtitle("Ar targets") +
    theme_classic(14) + theme(plot.title = element_text(face = 'bold', hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  pdf(file = file.path(plot.dir,'Ar_violin_epithelial.pdf'), width = 4, height = 6, bg = 'transparent')
  grid_arrange_list_shared_legend(list(g0, g1), nrow = 2, legend.position = 'right')
  dev.off()
  
  #### AR violin by cell type ####
  g0 <- VlnPlot(ss_epi, features = 'Ar', assay = 'RNA', slot='data', log = F,
                pt.size = 0, flip = T, group.by = 'cell_type', stack = F)
  df <- g0$data
  g0 <- ggplot(df, aes(x=ident, y=Ar, fill=ident)) +
    geom_violin(trim = T, scale = 'width', draw_quantiles = c(.25,.5,.75)) +
    scale_fill_manual("", values = celltype_palette_list) +
    xlab("") + ylab("Expression level") + ggtitle("Androgen Receptor (Ar)") +
    scale_y_continuous(limits = c(0,4), expand = c(0,0)) +
    theme_classic(14) + theme(plot.title = element_text(face = 'bold', hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  g1 <- VlnPlot(ss_epi, features = 'ar_targets', assay = 'RNA', slot='data',
                pt.size = 0, flip = T, group.by = 'cell_type', stack = F) + xlab("")
  df <- g1$data
  g1 <- ggplot(df, aes(x=ident, y=ar_targets, fill=ident)) +
    geom_violin(trim = T, scale = 'width', draw_quantiles = c(.25,.5,.75)) +
    scale_fill_manual("", values = celltype_palette_list) +
    xlab("") + ylab("Expression level") + ggtitle("Ar targets") +
    theme_classic(14) + theme(plot.title = element_text(face = 'bold', hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  g2 <- VlnPlot(ss_epi, features = 'Nelson_response_ANDR_UP_mouse', assay = 'RNA', slot='data',
                pt.size = 0, flip = T, group.by = 'cell_type', stack = F) + xlab("")
  df <- g2$data
  g2 <- ggplot(df, aes(x=ident, y=Nelson_response_ANDR_UP_mouse, fill=ident)) +
    geom_violin(trim = T, scale = 'width', draw_quantiles = c(.25,.5,.75)) +
    scale_fill_manual("", values = celltype_palette_list) +
    xlab("") + ylab("Expression level") + ggtitle("Nelson_response_ANDR_UP_mouse") +
    theme_classic(14) + theme(plot.title = element_text(face = 'bold', hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  g3 <- VlnPlot(ss_epi, features = 'Hieronymus_Androgen_mouse', assay = 'RNA', slot='data',
                pt.size = 0, flip = T, group.by = 'cell_type', stack = F) + xlab("")
  df <- g3$data
  g3 <- ggplot(df, aes(x=ident, y=Hieronymus_Androgen_mouse, fill=ident)) +
    geom_violin(trim = T, scale = 'width', draw_quantiles = c(.25,.5,.75)) +
    scale_fill_manual("", values = celltype_palette_list) +
    xlab("") + ylab("Expression level") + ggtitle("Hieronymus_Androgen_mouse") +
    theme_classic(14) + theme(plot.title = element_text(face = 'bold', hjust = 0.5), 
                              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  pdf(file = file.path(plot.dir,'Ar_violin_epithelial_celltype.pdf'), width = 5, height = 13, bg = 'transparent')
  grid_arrange_list_shared_legend(list(g0, g1, g2, g3), nrow = 4, legend.position = 'right')
  dev.off()
  
  
  g2 <- VlnPlot(ss_epi, features = 'Ar', assay = 'RNA', slot='data',
                pt.size = 0, flip = T, group.by = 'condition', stack = F, cols = condition_palette)
  df <- g2$data
  g2 <- ggplot(df, aes(x=ident, y=Ar, fill=ident)) +
    geom_violin(trim = T, scale = 'width', draw_quantiles = c(.25,.5,.75)) +
    scale_fill_manual("", values = condition_palette) +
    xlab("") + ylab("Expression level") + ggtitle("Ar") +
    theme_classic(14) + theme(plot.title = element_text(face = 'bold', hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  g3 <- VlnPlot(ss_epi, features = 'ar_targets', assay = 'RNA', slot='data',
                pt.size = 0, flip = T, group.by = 'condition', stack = F, cols = condition_palette)
  df <- g3$data
  g3 <- ggplot(df, aes(x=ident, y=ar_targets, fill=ident)) +
    geom_violin(trim = T, scale = 'width', draw_quantiles = c(.25,.5,.75)) +
    scale_fill_manual("", values = condition_palette) +
    xlab("") + ylab("Expression level") + ggtitle("Ar targets") +
    theme_classic(14) + theme(plot.title = element_text(face = 'bold', hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  pdf(file = file.path(plot.dir,'Ar_violin_epithelial.pdf'), width = 9, height = 9, bg = 'transparent')
  do.call(grid.arrange, c(list(g0, g1, g2, g3), nrow = 2))
  dev.off()
  
  
  colors <- colorRampPalette(c("royalblue4","royalblue2","deepskyblue",
                               "grey90","orangered3","red3","darkred"), space="rgb")(1024);
  FeaturePlot(ss_epi, features = 'ar_targets', ncol = 1, slot = 'data', cols = colors)
  
  FeaturePlot(ss_epi, features = 'Ar', reduction = 'umap')
  
 
}

#### Helpers ####


#### Returns expression data matrix by ident (i.e. total expression by cluster) ####
TotalExpression <- function (object, assays = NULL, features = NULL, return.seurat = FALSE,
                             add.ident = NULL, slot = "data", use.scale = FALSE, use.counts = FALSE,
                             verbose = TRUE, ...)
{
  library(rlang)
  Seurat:::CheckDots(..., fxns = "CreateSeuratObject")
  if (use.scale) {
    .Deprecated(msg = "'use.scale' is a deprecated argument, please use the 'slot' argument instead")
    slot <- "scale.data"
  }
  if (use.counts) {
    .Deprecated(msg = "'use.counts' is a deprecated argument, please use the 'slot' argument instead")
    if (use.scale) {
      warning("Both 'use.scale' and 'use.counts' were set; using counts",
              call. = FALSE, immediate. = TRUE)
    }
    slot <- "counts"
  }
  fxn.total <- switch(EXPR = slot, data = function(x) {
    rowSums(x = expm1(x = x))
  }, rowSums)
  object.assays <- Seurat:::FilterObjects(object = object, classes.keep = "Assay")
  assays <- assays %||% object.assays
  ident.orig <- Idents(object = object)
  orig.levels <- levels(x = Idents(object = object))
  ident.new <- c()
  if (!all(assays %in% object.assays)) {
    assays <- assays[assays %in% object.assays]
    if (length(assays) == 0) {
      stop("None of the requested assays are present in the object")
    }
    else {
      warning("Requested assays that do not exist in object. Proceeding with existing assays only.")
    }
  }
  if (!is.null(x = add.ident)) {
    new.data <- Seurat:::FetchData(object = object, vars = add.ident)
    new.ident <- paste(Idents(object)[rownames(x = new.data)],
                       new.data[, 1], sep = "_")
    Idents(object, cells = rownames(new.data)) <- new.ident
  }
  data.return <- list()
  for (i in 1:length(x = assays)) {
    data.use <- GetAssayData(object = object, assay = assays[i],
                             slot = slot)
    features.assay <- features
    if (length(x = intersect(x = features, y = rownames(x = data.use))) <
        1) {
      features.assay <- rownames(x = data.use)
    }
    data.all <- list(data.frame(row.names = features.assay))
    for (j in levels(x = Idents(object))) {
      temp.cells <- WhichCells(object = object, idents = j)
      features.assay <- unique(x = intersect(x = features.assay,
                                             y = rownames(x = data.use)))
      if (length(x = temp.cells) == 1) {
        data.temp <- (data.use[features.assay, temp.cells])
        if (slot == "data") {
          data.temp <- expm1(x = data.temp)
        }
      }
      if (length(x = temp.cells) > 1) {
        data.temp <- fxn.total(data.use[features.assay,
                                        temp.cells, drop = FALSE])
      }
      data.all[[j]] <- data.temp
      if (verbose) {
        message(paste("Finished totalling", assays[i],
                      "for cluster", j))
      }
      if (i == 1) {
        ident.new <- c(ident.new, as.character(x = ident.orig[temp.cells[1]]))
      }
    }
    names(x = ident.new) <- levels(x = Idents(object))
    data.return[[i]] <- do.call(cbind, data.all)
    names(x = data.return)[i] <- assays[[i]]
  }
  if (return.seurat) {
    toRet <- CreateSeuratObject(counts = data.return[[1]],
                                project = "Total", assay = names(x = data.return)[1],
                                ...)
    toRet <- SetAssayData(object = toRet, assay = names(x = data.return)[1],
                          slot = "data", new.data = log1p(x = as.matrix(x = data.return[[1]])))
    if (length(x = data.return) > 1) {
      for (i in 2:length(x = data.return)) {
        toRet[[names(x = data.return)[i]]] <- CreateAssayObject(counts = data.return[[i]])
        toRet <- SetAssayData(object = toRet, assay = names(x = data.return)[i],
                              slot = "data", new.data = log1p(x = as.matrix(x = data.return[[i]])))
      }
    }
    if (DefaultAssay(object = object) %in% names(x = data.return)) {
      DefaultAssay(object = toRet) <- DefaultAssay(object = object)
    }
    Idents(toRet, cells = colnames(x = toRet)) <- ident.new[colnames(x = toRet)]
    Idents(object = toRet) <- factor(x = Idents(object = toRet),
                                     levels = as.character(x = orig.levels), ordered = TRUE)
    toRet <- ScaleData(object = toRet, verbose = verbose)
    return(toRet)
  }
  else {
    return(data.return)
  }
}


silhouetteScore <- function(obj, reduction='umap', clusters='seurat_clusters', dims=c(1:30)){
  require(cluster)
  
  dist.matrix <- dist(x = Embeddings(object = obj[[reduction]])[, dims])
  clusters <- obj@meta.data[[clusters]]
  sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
}

unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
  registerDoSEQ()
  gc()
}

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

LuminalBasalClassification <- function(){
  library(MASS)
  library(plotly)
  library(ggrastr)
  
  #### Score CD24 CD49 from bulk ####
  ## Wouter's bulkRNA
  # w_bulk_rna <- setDT(readxl::read_xlsx('/home/erik/projects/sawyers/foxa1/analysis/multiome/spreadsheets/DESEQ_condition_cd49f_vs_cd24.xlsx'))
  # w_bulk_rna[,source:='wouter']
  # ## Abbas's bulkRNA
  # a_bulk_rna <- setDT(readxl::read_xlsx('/home/erik/projects/sawyers/foxa1/analysis/multiome/spreadsheets/Abbas_DESEQcondition_cd49f_vs_cd24a.xlsx'))
  # a_bulk_rna[,source:='abbas']
  # # standalone EV
  # scRNA_DE <- setDT(readxl::read_xlsx('/home/erik/projects/sawyers/foxa1/analysis/multiome/RNA_standalone/spreadsheets/scRNA_cluster_vs_others_wilcox_de.xlsx'))
  # scRNA_DE[,source:='standalone']
  # intersect(w_bulk_rna[padj < 0.01 & log2FoldChange > 1, gene],
  #           bulk_rna[padj < 0.01 & log2FoldChange > 1, gene]) %>% sort
  
  erg.scRNA$luminal.Scores <- NULL
  erg.scRNA$basal.Scores <- NULL
  (features = list(lum_mature_markers, basal_markers))
  
  sl <- lapply(unique(erg.scRNA$bioNames), function(i){
    sset <- subset(x = erg.scRNA, subset = (bioNames == i))
    sset <- AddModuleScore(object = sset,
                           features = features,
                           assay = 'RNA',
                           nbin = 24,
                           ctrl = 100,
                           seed = 2,
                           name = 'luminal_v_basal.Score',
                           force = T,
                           search = F)
    colnames(sset@meta.data)[which(colnames(sset@meta.data)=='luminal_v_basal.Score1')] <- 'luminal.Scores'
    colnames(sset@meta.data)[which(colnames(sset@meta.data)=='luminal_v_basal.Score2')] <- 'basal.Scores'
    sset
  })
  names(sl) <- unique(erg.scRNA$bioNames)
  
  
  ## Gather basal_luminal scores for each sample
  dt <- lapply(names(sl), function(i){
    o <- sl[[i]]
    data.table(cells=rownames(o@meta.data), 
               percent_epcam=o@meta.data$percent_epcam,
               o@meta.data[,grep("luminal.Scores|basal.Scores|bioNames", colnames(o@meta.data))])
  }) %>% rbindlist()
  
  dt$luminal_basal <- 'Luminal'
  dt[basal.Scores > luminal.Scores,]$luminal_basal <- 'Basal'
  
  if(identical(Cells(erg.scRNA), dt$cells)){
    erg.scRNA$luminal.Scores <- dt$luminal.Scores
    erg.scRNA$basal.Scores <- dt$basal.Scores
    erg.scRNA$luminal_basal <- dt$luminal_basal
  }
  colnames(erg.scRNA@meta.data)
  
  FeaturePlot(ss_epi, features = c('luminal.Scores'), cols = rev(pals::brewer.rdbu(11)))
  
  FeaturePlot(ss_epi, features = c('basal.Scores'), cols = rev(pals::brewer.rdbu(11)))
  
  DimPlot(ss_epi, group.by = c('luminal_basal'), cols = c("dodgerblue","goldenrod"))
  
  
  #### Define epithelial cells and clusters ####
  dt <- data.table(erg.scRNA@meta.data, keep.rownames = 'cells')
  dt[, quantile(percent_epcam, c(0.25)), by=seurat_clusters]
  
  ds <- merge(dt[percent_epcam <=0, .(NoEpcam=.N), by=seurat_clusters],
              dt[percent_epcam > 0, .(Epcam=.N), by=seurat_clusters], 
              by='seurat_clusters', all.x=T, all.y=T)
  ds[is.na(NoEpcam), NoEpcam:=0]
  ds[is.na(Epcam), Epcam:=0]
  ds[,FracEpcam := Epcam/(Epcam+NoEpcam)]
  
  png(filename = file.path(plot.dir,"fracEpcamBar.png"), width = 8, height = 4.5, units = 'in', res = 150)
  ggplot(ds, aes(x=reorder(seurat_clusters, -FracEpcam), y=FracEpcam, fill=FracEpcam)) + xlab("") +
    geom_col() + theme_bw(14) + scale_fill_fermenter(n.breaks = 10, palette = "PuOr")
  dev.off()
  
  png(filename = file.path(plot.dir,"fracEpcamBox.png"), width = 7, height = 4, units = 'in', res = 150)
  ggplot(dt[, percent_epcam, by=seurat_clusters], 
         aes(x=reorder(seurat_clusters, -percent_epcam), y=percent_epcam)) + xlab("") +
    geom_boxplot(fill='red')  + theme_bw(14) #+ scale_y_continuous(breaks = seq(0, 2, by=0.1))
  dev.off()
  
  
  ## All cells - density
  gl <- lapply(unique(dt$bioNames), function(i){
    ggplot(dt[bioNames==i], aes(x=luminal.Scores, y=basal.Scores)) +
      geom_point((aes(color=luminal_basal))) +
      ggtitle(i) +
      geom_density_2d(data = dt, bins=10,contour_var = "density", col='black', size=.1) +
      geom_density_2d(bins=10,contour_var = "density", col='blue') +
      theme_bw()
  })
  g <- ggplot(dt, aes(x=luminal.Scores, y=basal.Scores)) +
    geom_density_2d(bins=20,contour_var = "density") +
    ggtitle("Density of all cells") +
    theme_bw()
  gl <- append(gl, list(g))
  
  png(filename = file.path(plot.dir,'luminal_basal_sample_scores_markers.png'), width = 10, height = 8, units = 'in', bg = 'transparent', res = 150)
  grid_arrange_list_shared_legend(gl, nrow = 3, ncol = 3, legend.position = 'top')
  dev.off()
  
  ## EPCAM + cells - density
  dt <- data.table(ss_epi@meta.data, keep.rownames = 'cells')
  
  gl <- lapply(unique(dt$bioNames), function(i){
    ggplot(dt[bioNames==i & percent_epcam > 0], aes(x=luminal.Scores, y=basal.Scores)) +
      geom_point((aes(color=luminal_basal))) +
      ggtitle(i) +
      geom_density_2d(data = dt, bins=10,contour_var = "density", col='black', size=.1) +
      geom_density_2d(bins=10,contour_var = "density", col='blue') +
      theme_bw()
  })
  g <- ggplot(dt[percent_epcam > 0], aes(x=luminal.Scores, y=basal.Scores)) +
    #geom_point((aes(color=luminal_basal))) +
    geom_density_2d(bins=20,contour_var = "density") +
    ggtitle("Density of data set") +
    theme_bw()
  gl <- append(gl, list(g))
  
  png(filename = file.path(plot.dir,'luminal_basal_sample_scores_markers_epcam_cells.png'), width = 10, height = 8, units = 'in', bg = 'transparent', res = 150)
  grid_arrange_list_shared_legend(gl, nrow = 3, ncol = 3, legend.position = 'top')
  dev.off()
  
  ## Inspection of density
  library(MASS)
  library(plotly)
  dt <- data.table(erg.scRNA@meta.data)
  dt[,cells:= rownames(erg.scRNA@meta.data)]
  class <- factor(dt$luminal_basal, labels=unique(dt$luminal_basal))
  
  X <- dt[,list(luminal.Scores, basal.Scores)]
  rownames(X) <- dt$cells
  
  den3d <- kde2d(X$luminal.Scores, X$basal.Scores)
  persp(den3d, box=FALSE)
  
  plot_ly(x=den3d$x, y=den3d$y, z=den3d$z) %>%
    add_surface() %>%
    layout(scene = list(xaxis=list(title = 'luminal'),
                        yaxis=list(title = 'basal'),
                        zaxis=list(title = 'Density')))
  
  #### Kmeans ####
  library(factoextra)
  fv <- fviz_nbclust(x = X, FUNcluster = kmeans, method = "silhouette")
  fv # shows optimal number of clusters is 2.
  
  png(file.path(plot.dir, 'KMeansOptimalClusters.png'), width = 6, height = 8, units = 'in', res = 150, bg = 'white')
  # create an apporpriate viewport.  Modify the dimensions and coordinates as needed
  vp.Bottom <- viewport(height=unit(.45, "npc"), width=unit(1, "npc"), 
                        just=c("left","top"), 
                        y=0.45, x=0)
  # plot your base graphics 
  par(mfrow=c(2,2))
  (cl2 <- kmeans(x = X, centers = 2))
  plot(X, col=cl2$cluster)
  points(cl2$centers, pch=8, cex=2, col='white') #col=1:2,
  
  (cl3 <- kmeans(x = X, centers = 3))
  plot(X, col=cl3$cluster)
  points(cl3$centers, pch=8, cex=2, col='white') #col=1:2,
  # plot the ggplot using the print command
  print(fv, vp=vp.Bottom)
  dev.off()
  
  # totss = total variance in data
  cl2$totss
  cl3$totss
  # size of each cluster
  cl2$size
  cl3$size
  fvcl <- fviz_cluster(object = cl3, data = X)
  fvcl
  
  #### MVGaussian ####
  library(mclust)
  #### All cells ####
  mod <- densityMclust(X)
  summary(mod)
  plot(mod, what = "BIC")
  plot(mod, what = "density")
  
  BIC <- mclustBIC(X)
  modBIC <- Mclust(X, x=BIC)
  png(filename = file.path(plot.dir,'BIC_luminal_basal_mclust_diagnostics.png'), width = 10, height = 10, units = 'in', res = 150)
  par(mfrow=c(2,2))
  plot(modBIC, what = "BIC", main=T)
  plot(modBIC, what = "classification", main=T, col=pals::alphabet(22))
  plot(modBIC, what = "density", main=T)
  plot(modBIC, what = "uncertainty", main=T, col=pals::alphabet(22))
  dev.off()
  
  
  
  # ------------------------------------#
  #### Epi cells ####
  #### dte is data from epithelial cells
  dte <- data.table(ss_epi@meta.data)
  dte[,cells:= rownames(ss_epi@meta.data)]
  Xe <- dte[,list(luminal.Scores, basal.Scores)]
  rownames(Xe) <- dte$cells
  
  den3d <- kde2d(Xe$luminal.Scores, Xe$basal.Scores)
  persp(den3d, box=FALSE)
  
  plot_ly(x=den3d$x, y=den3d$y, z=den3d$z) %>%
    add_surface() %>%
    layout(scene = list(xaxis=list(title = 'luminal'),
                        yaxis=list(title = 'basal'),
                        zaxis=list(title = 'Density')), 
           title="EPCAM+ cells")
  
  #### MVGaussian ####
  library(mclust)
  #### Epithelial cells as defined previously with Wouter scores ####
  erg.scRNA@meta.data$seurat_clusters_class <- 'Unknown'
  
  erg.scRNA@meta.data[rownames(erg.scRNA@meta.data) %in% dte$cells,'seurat_clusters_class'] <- 'epithelial'
  
  ## Inspection of density
  dt <- data.table(erg.scRNA@meta.data[,grep("percent_epcam|luminal.Scores|basal.Scores|luminal_basal|bioNames|seurat_clusters", 
                                             colnames(erg.scRNA@meta.data))], keep.rownames = "cells")
  table(dt$bioNames)
  Class <- factor(dte$luminal_basal, labels=unique(dte$luminal_basal))
  table(Class)
  
  table(dte$bioNames)
  Xe <- dte[,luminal.Scores, basal.Scores]
  rownames(Xe) <- dte$cells
  
  ## modBIC ##
  BIC <- mclustBIC(Xe)
  modBIC <- Mclust(Xe, x=BIC)
  summary(modBIC)
  table(Class, modBIC$classification)
  # evaluate clustering solution (0=random, 1=upper boundary, higher values are better)
  adjustedRandIndex(Class, modBIC$classification)
  
  png(filename = file.path(plot.dir,'BIC_epi_luminal_basal_mclust_diagnostics.png'), width = 10, height = 10, units = 'in', res = 150)
  par(mfrow=c(2,2))
  plot(modBIC, what = "BIC", main=T)
  plot(modBIC, what = "classification", main=T, col=pals::alphabet(22))
  plot(modBIC, what = "density", main=T)
  plot(modBIC, what = "uncertainty", main=T, col=pals::alphabet(22))
  dev.off()
  
  drmod <- MclustDR(modBIC, lambda = 1)
  summary(drmod)
  plot(drmod, what = "contour")
  # On the same subspace we can also plot the uncertainty boundaries corresponding to the MAP classification:
  plot(drmod, what = "boundaries", ngrid = 200)
  #and then add a circle around the misclassified observations
  miscl <- classError(Class, modBIC$classification)$misclassified
  points(drmod$dir[miscl,], pch = 1, cex = 2)
  
  
  ## mod2 ##
  mod2 <- Mclust(Xe, G = 2, prior = priorControl())
  summary(mod2)
  table(Class, mod2$classification)
  mclust2Dplot(Xe, classification = mod2$classification, parameters = mod2$parameters, colors = pals::alphabet(26))
  # evaluate clustering solution (0=random, 1=upper boundary, higher values are better)
  adjustedRandIndex(Class, mod2$classification) # 0.7897
  
  
  png(filename = file.path(plot.dir,'MOD2_epi_luminal_basal_mclust_diagnostics.png'), width = 10, height = 10, units = 'in', res = 150)
  par(mfrow=c(2,2))
  plot(mod2, what = "BIC", main=T)
  plot(mod2, what = "classification", main=T, sub=paste(mod2$modelName,"model"), col=pals::alphabet(22))
  plot(mod2, what = "density", main=T, sub=paste(mod2$modelName,"model"), )
  plot(mod2, what = "uncertainty", main=T, sub=paste(mod2$modelName,"model"), col=pals::alphabet(22))
  dev.off()
  
  drmod <- MclustDR(mod2, lambda = 1)
  summary(drmod)
  png(filename = file.path(plot.dir,'MOD2_epi_luminal_basal_mclust_MclustDR.png'), width = 8, height = 8, units = 'in', res = 150)
  # On the same subspace we can also plot the uncertainty boundaries corresponding to the MAP classification:
  #plot(drmod, what = "boundaries", ngrid = 200)
  plot(drmod, what = "contour")
  # On the same subspace we can also plot the uncertainty boundaries corresponding to the MAP classification:
  plot(drmod, what = "boundaries", ngrid = 200)
  #and then add a circle around the misclassified observations
  miscl <- classError(Class, mod2$classification)$misclassified
  points(drmod$dir[miscl,], pch = 1, cex = 1)
  title(main=paste("Adjusted Rand Index:", round(adjustedRandIndex(Class, mod2$classification), 4)))
  dev.off()
  
  
  ## mod3 ##
  mod3 <- Mclust(Xe, G = 3, prior = priorControl())
  #mod3 <- Mclust(Xe, G = 3, modelNames = "EEV") #VEV, EEV?, EVV? VVV
  summary(mod3)
  table(Class, mod3$classification)
  mclust2Dplot(Xe, classification = mod3$classification, parameters = mod3$parameters, colors = pals::alphabet(26))
  # evaluate clustering solution (0=random, 1=upper boundary, higher values are better)
  adjustedRandIndex(Class, mod3$classification) # 0.70644
  
  png(filename = file.path(plot.dir,'MOD3_epi_luminal_basal_mclust_diagnostics.png'), width = 10, height = 10, units = 'in', res = 150)
  par(mfrow=c(2,2))
  plot(mod3, what = "BIC", main=T)
  plot(mod3, what = "classification", main=T, sub=paste(mod3$modelName,"model"), col=pals::alphabet(22))
  plot(mod3, what = "density", main=T, sub=paste(mod3$modelName,"model"))
  plot(mod3, what = "uncertainty", main=T, sub=paste(mod3$modelName,"model"), col=pals::alphabet(22))
  dev.off()
  
  drmod <- MclustDR(mod3, lambda = 1)
  summary(drmod)
  png(filename = file.path(plot.dir,'MOD3_epi_luminal_basal_mclust_MclustDR.png'), width = 8, height = 8, units = 'in', res = 150)
  # On the same subspace we can also plot the uncertainty boundaries corresponding to the MAP classification:
  #plot(drmod, what = "boundaries", ngrid = 200)
  plot(drmod, what = "contour")
  # On the same subspace we can also plot the uncertainty boundaries corresponding to the MAP classification:
  plot(drmod, what = "boundaries", ngrid = 200)
  #and then add a circle around the misclassified observations
  miscl <- classError(Class, mod3$classification)$misclassified
  points(drmod$dir[miscl,], pch = 1, cex = 1)
  title(main=paste("Adjusted Rand Index:", round(adjustedRandIndex(Class, mod3$classification), 4)))
  dev.off()
  
  ICL3 <- mclustICL(Xe)
  png(filename = file.path(plot.dir,'MOD3_epi_luminal_basal_mclust_ICL.png'), width = 6, height = 6, units = 'in', res = 150)  
  plot(ICL3)
  dev.off()
  
  
  #### Cluster uncertainty ####
  uncerPlot(mod2$z)
  sort(mod2$uncertainty, decreasing = TRUE) %>% head(20)
  
  ## Should not be the same
  identical(Xe$luminal.Scores, unname(erg.scRNA$luminal.Scores))
  identical(Xe$basal.Scores, unname(erg.scRNA$basal.Scores))
  
  # mod2$z contains conditional probabilities for 2 classes 1 and 2.
  # These are the certainty of classification
  # apply(mod2$z, 1, max) are the certainties of classification
  # 1-apply(mod2$z, 1, max) are the uncertainties of classification
  mod2$z
  dta <- data.table(class=mod2$classification, certain=(1-mod2$uncertainty))
  ggplot(dta, aes(x=as.factor(class), y=certain, fill=as.factor(class))) + geom_violin() + theme_bw()
  
  quantile(mod2$uncertainty, probs=seq(0,1, by=.25))
  # 0%          25%          50%          75%         100% 
  # 5.329071e-15 2.096478e-05 8.691648e-04 2.170055e-02 6.591430e-01 
  #           0%          25%          50%          75%         100%
  # 2.365796e-09 5.635173e-03 5.360654e-02 2.176453e-01 4.999247e-01
  # Here, all data is less than .659
  # Here, 3/4 of the data are smaller than .0217
  # Here, 1/2 of the data are smaller than .000869
  
  # erg.scRNA$luminal_basal_classification <- ifelse(mod2$classification==1, 'Luminal',
  #                                                       ifelse(mod2$classification==2,'Basal','Unknown'))
  
  #### Set cell classification ####
  erg.scRNA$luminal_basal_classification <- 'Unknown'
  erg.scRNA$luminal_basal_classification_uncertainty <- 1
  dt <- data.table(erg.scRNA@meta.data, keep.rownames = 'cells')
  setkey(dt, cells)
  
  dtclass2 <- dt[names(mod2$classification),]
  identical(dtclass2$cells, names(mod2$classification))
  table(dtclass2$luminal_basal, mod2$classification)
  basal_idx <- 2
  lum_idx <- 1
  mix_idx <- 3
  
  epi_cells <- names(mod2$classification)
  lum_cells <- names(which(mod2$classification==lum_idx))
  basal_cells <- names(which(mod2$classification==basal_idx))
  erg.scRNA@meta.data[lum_cells,]$luminal_basal_classification <- 'Luminal'
  erg.scRNA@meta.data[basal_cells,]$luminal_basal_classification <- 'Basal'
  
  erg.scRNA@meta.data[names(mod2$uncertainty),]$luminal_basal_classification_uncertainty <- mod2$uncertainty
  erg.scRNA$luminal_basal_classification_certainty <- 1-mod2$uncertainty
  
  erg.scRNA$probability_luminal <- 0
  erg.scRNA$probability_basal <- 0
  erg.scRNA@meta.data[epi_cells,]$probability_luminal <- mod2$z[epi_cells,lum_idx]
  erg.scRNA@meta.data[epi_cells,]$probability_basal <- mod2$z[epi_cells,basal_idx]
  erg.scRNA$probability_other <- 1-(erg.scRNA$probability_luminal + erg.scRNA$probability_basal)
  
  table(erg.scRNA$luminal_basal_classification)
  
  
  # density estimate
  dens <- densityMclust(Xe)
  png(filename = file.path(plot.dir,'MOD2_epi_luminal_basal_mclust_density.png'), width = 6, height = 12, units = 'in', res = 150)
  par(mfrow=c(3,1))
  plot(dens, what = 'density', data = Xe, grid = 200, points.cex = 0.5, drawlabels = FALSE)
  plot(dens, what = 'density', type = 'image', col = 'steelblue', grid = 200)
  plot(dens, what = 'density', type = 'persp', theta = -25, phi = 20)
  dev.off()
  
  
  ###--------------------------------###
  #### DA using WT as train ####
  nrow(dte)
  dttrain <- dte[cells %in% names(mod2$classification),]
  identical(dttrain$cells, names(mod2$classification))
  dttrain <- dttrain[bioNames=='Wildtype 3 months']
  Xtrain <- dttrain[,list(basal.Scores, luminal.Scores)]
  rownames(Xtrain) <- dttrain[, cells]
  
  # determine class for EV
  mod2 <- Mclust(Xtrain, G = 2, prior = priorControl())
  summary(mod2)
  plot(mod2, what = "classification", main=T, sub=paste(mod2$modelName,"model"), col=pals::alphabet(22))
  plot(mod2, what = "density", main=T, sub=paste(mod2$modelName,"model"), )
  
  Class <- factor(dttrain$luminal_basal, labels=unique(dttrain$luminal_basal))
  table(Class, mod2$classification)
  Class.train <- ifelse(mod2$classification==basal_idx,"Basal","Luminal")
  
  dttest <- dte[bioNames!='Wildtype 3 months']
  Xtest <- dttest[,list(basal.Scores, luminal.Scores)]
  rownames(Xtest) <- dttest[, cells]
  Class.test <- dttest$luminal_basal
  
  lbClustDA <- MclustDA(data = Xtrain, G = 1, class = Class.train, modelNames = mod2$modelName)
  summary(lbClustDA, parameters=T)
  summary(lbClustDA, newdata = Xtest, newclass = Class.test)
  
  plot(lbClustDA, what = 'classification', main=T)
  plot(lbClustDA, what = "scatterplot", main=T)
  plot(lbClustDA, what = "classification", newdata = Xtest, main=T)
  plot(lbClustDA, what = "scatterplot", newdata = Xtest, main=T)
  plot(lbClustDA, what = "train&test", newdata = Xtest)
  
  png(filename = file.path(plot.dir,'ClustDA_luminal_basal_mclust_diagnostics.png'), width = 10, height = 12, units = 'in', res = 150)
  par(mfrow=c(3,2))
  plot(mod2, what = "classification", main=T, sub=paste("Training EV",mod2$modelName,"model"), col=pals::alphabet(22))
  plot(mod2, what = "density", main=T, sub=paste("Training EV",mod2$modelName,"model"), )
  plot(lbClustDA, what = "classification", main=T, col=pals::alphabet(22))
  plot(lbClustDA, what = "scatterplot", main=T )
  plot(lbClustDA, what = "error", main = T)
  plot(lbClustDA, what = "train&test", newdata = Xtest, main=T)
  dev.off()
  
  predictTrain <- predict(lbClustDA)
  predictTrain
  
  predictTest <- predict(lbClustDA, Xtest)
  predictTest
  
  dtrain <- data.table(cells=names(mod2$classification),
                       classification=mod2$classification, 
                       basal_certainty=c(mod2$z[,1]),
                       luminal_certainty=c(mod2$z[,2]))
  
  dtest <- data.table(cells=rownames(Xtest),
                      classification=c(as.character(predictTest$classification)),
                      basal_certainty=predictTest$z[,"Basal"],
                      luminal_certainty=predictTest$z[,"Luminal"])
  
  dtMclustDA <- data.table(cells=c(rownames(Xtrain), rownames(Xtest)),
                           classification=c(as.character(predictTrain$classification),
                                            as.character(predictTest$classification)),
                           basal_certainty=c(predictTrain$z[,"Basal"], predictTest$z[,"Basal"]),
                           luminal_certainty=c(predictTrain$z[,"Luminal"], predictTest$z[,"Luminal"]))
  
  head(dtMclustDA)
  plot(dtMclustDA$basal_certainty)
  plot(dtMclustDA$luminal_certainty)
  
  drmod <- MclustDR(lbClustDA, lambda = 1)
  summary(drmod)
  
  dtMclustDA <- dtMclustDA[data.table(cells = Cells(erg.scRNA)), on = "cells"] # cells = Cells(erg.scRNA), on="cells"]
  identical(dtMclustDA$cells, Cells(erg.scRNA))
  dtMclustDA[is.na(classification), classification:='Unknown']
  dtMclustDA[is.na(basal_certainty), basal_certainty:=NA]
  dtMclustDA[is.na(luminal_certainty), luminal_certainty:=NA]
  
  erg.scRNA$train_test_luminal_basal_classification <- dtMclustDA$classification
  erg.scRNA$train_test_probability_luminal <- dtMclustDA$luminal_certainty
  erg.scRNA$train_test_probability_basal <- dtMclustDA$basal_certainty
  
  #### UMAP luminal basal cell probabilities ####
  g0 <- FeaturePlot(erg.scRNA, features = c('probability_luminal'),
                    pt.size = 0.5, ncol = 1, slot = 'data')
  g0$data <- g0$data[order(g0$data$probability_luminal, decreasing = T),]
  nonepi <- rownames(erg.scRNA@meta.data[erg.scRNA@meta.data$seurat_clusters_class=='Unknown',])
  g0$data[nonepi,]$probability_luminal <- NA
  (g0 <- UmapPlotEnhance(obj = erg.scRNA, 
                         p1 = g0,
                         label.size = 5,
                         pt.size = .8, 
                         break.unit = .2,
                         title = 'Epithelial Luminal probability', 
                         legend.title = "Probability cell is luminal",
                         pal=luminal_basal_palette_continuous, 
                         resolution = 300,
                         filename.prefix = 'Erg_luminal_probability_Umap_markers_epi_reclust'))
  
  gc <- DimPlot(erg.scRNA, group.by = c('seurat_clusters'))
  gc$data$seurat_clusters <- as.integer(gc$data$seurat_clusters) #make continuous for continuous layering
  (gc <- UmapPlotEnhance(obj = erg.scRNA, 
                         p1 = gc,
                         addLabels = T,
                         label.size = 4,
                         pt.size = .8, 
                         title = 'Probabillity Luminal', 
                         pal=NA,
                         resolution = 300))
  png(filename = file.path(plot.dir,"Erg_luminal_probability_Umap_clusters_epi_reclust.png"), width = 5, height = 4, res = 300, units = 'in', bg = 'transparent')
  g0 + gc$layers
  dev.off()
  
  
  #### Umap LB probabilities train/test ####
  g0 <- FeaturePlot(erg.scRNA, features = c('train_test_probability_luminal'),
                    pt.size = 0.5, ncol = 1, slot = 'data')
  g0$data <- g0$data[order(g0$data$train_test_probability_luminal, decreasing = T),]
  nonepi <- rownames(erg.scRNA@meta.data[erg.scRNA@meta.data$seurat_clusters_class=='Unknown',])
  g0$data[nonepi,]$train_test_probability_luminal <- NA
  (g0 <- UmapPlotEnhance(obj = erg.scRNA, 
                         p1 = g0,
                         label.size = 5,
                         pt.size = .8, 
                         break.unit = .2,
                         title = 'Epithelial Luminal probability', 
                         legend.title = "Probability cell is luminal",
                         pal=luminal_basal_palette_continuous, 
                         resolution = 300,
                         filename.prefix = 'Erg_luminal_probability_Umap_markers_epi_train_test'))
  
  gc <- DimPlot(erg.scRNA, group.by = c('seurat_clusters'))
  gc$data$seurat_clusters <- as.integer(gc$data$seurat_clusters) #make continuous for continuous layering
  (gc <- UmapPlotEnhance(obj = erg.scRNA, 
                         p1 = gc,
                         addLabels = T,
                         label.size = 4,
                         pt.size = .8, 
                         title = 'Probabillity Luminal', 
                         pal=NA,
                         resolution = 300))
  png(filename = file.path(plot.dir,"Erg_luminal_probability_Umap_clusters_epi_train_test.png"), width = 5, height = 4, res = 300, units = 'in', bg = 'transparent')
  g0 + gc$layers
  dev.off()
  
  #### Do we choose the fit of all epi data or the test training for probabilities and cell type assignment? ####
  
  ## Look at scatter plots for each
  #### Fits ####
  #### 2d density basal, luminal scores ####
  require(ggrastr)
  
  g0 <- ggplot(dte, aes(x=basal.Scores, y=luminal.Scores)) +
    geom_point_rast((aes(color=luminal_basal_classification)), shape=1, size=.8, alpha=.8, stroke=1) +
    geom_density_2d_filled(bins=10, linewidth=0, alpha=.5, color='black') +
    #scale_fill_manual("Density", values = c("white",pals::ocean.amp(10))) +
    #scale_fill_manual("Density", values = c("white",rev(pals::kovesi.linear_grey_10_95_c0(10)))) +
    scale_fill_manual("Density", values = c(NA, pals::ocean.amp(10)), na.value = NA) +
    scale_color_manual("Classification", values = luminal_basal_palette) +
    ggtitle("Epithelial cells") +
    scale_x_continuous(limits = c(-1.5,3)) + 
    scale_y_continuous(limits = c(-1.5,3)) + 
    geom_hline(yintercept = 0, linewidth=.1) + 
    geom_vline(xintercept = 0, linewidth=.1) + 
    theme_bw(12) +
    theme(panel.grid = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          legend.background = element_blank(),
          legend.key.size = unit(24, units='points'),
          legend.text = element_text(size = 12)) +
    guides(colour = guide_legend(override.aes = list(size=4, alpha=1)),
           fill = guide_legend(override.aes = list(size=.5, alpha=1)))
  png(filename = file.path(plot.dir,'Luminal_Basal_MVG_Epi_classification_eval_reclust.png'), width = 7.5, height = 6, units = 'in', res=200, bg = 'transparent')
  print(g0)
  dev.off()
  
  
  gl <- lapply(unique(dte$bioNames), function(i){
    ggplot(dte[bioNames==i], aes(x=basal.Scores, y=luminal.Scores)) +
      geom_point_rast((aes(color=luminal_basal_classification)), shape=1, size=.8, alpha=.8, stroke=1) +
      geom_density_2d_filled(bins=10, linewidth=0, alpha=.5, color='black') +
      #scale_fill_manual("Density", values = c("white",pals::ocean.amp(10))) +
      #scale_fill_manual("Density", values = c("white",rev(pals::kovesi.linear_grey_10_95_c0(10)))) +
      scale_fill_manual("Density", values = c(NA, pals::ocean.amp(10)), na.value = NA) +
      scale_color_manual("Classification", values = luminal_basal_palette) +
      ggtitle(label = i, subtitle = "Epithelial cells") +
      scale_x_continuous(limits = c(-1.5,3)) + 
      scale_y_continuous(limits = c(-1.5,3)) + 
      geom_hline(yintercept = 0, linewidth=.1) + 
      geom_vline(xintercept = 0, linewidth=.1) + 
      theme_bw(12) +
      theme(panel.grid = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank(),
            legend.background = element_blank(),
            legend.key.size = unit(24, units='points'),
            legend.text = element_text(size = 12)) +
      guides(colour = guide_legend(override.aes = list(size=4, alpha=1)),
             fill = guide_legend(override.aes = list(size=.5, alpha=1)))
  })
  png(filename = file.path(plot.dir,'Luminal_Basal_MVG_classification_by_sample_markers_eval_reclust.png'), width = 12, height = 8, units = 'in', res=200, bg = 'transparent')
  grid_arrange_list_shared_legend(gl, ncol=3, legend.position = 'right')
  dev.off()
  
  #### 2d density Luminal basal probabilities ####
  gl <- lapply(unique(dte$bioNames), function(i){
    ggplot(dte[bioNames==i], aes(x=basal.Scores, y=luminal.Scores)) +
      geom_point_rast((aes(color=probability_luminal)), shape=1, size=.8, alpha=.5, stroke=1) +
      geom_density_2d(bins=10, linewidth=.5, alpha=.5, color='black') +
      #scale_fill_manual("Density", values = c("white",pals::ocean.amp(10))) +
      scale_fill_manual("Density", values = c(NA, pals::ocean.amp(10)), na.value = NA) +
      #scale_color_manual("Classification", values = luminal_basal_palette_continuous) +
      #scale_color_binned(type = "viridis", n.breaks=10) +
      scale_color_stepsn("Prob. cell is luminal",colours = luminal_basal_palette_continuous, n.breaks=10) +
      ggtitle(label = i, subtitle = "Epithelial cells") +
      scale_x_continuous(limits = c(-1.5,3)) + 
      scale_y_continuous(limits = c(-1.5,3)) + 
      geom_hline(yintercept = 0, linewidth=.1) + 
      geom_vline(xintercept = 0, linewidth=.1) + 
      theme_bw(12) +
      theme(panel.grid = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank(),
            legend.background = element_blank(),
            legend.key.size = unit(24, units='points'),
            legend.text = element_text(size = 12)) +
      guides(#colour = guide_legend(override.aes = list(size=10, alpha=1)),
        fill = guide_legend(override.aes = list(size=.5, alpha=1)))
  })
  png(filename = file.path(plot.dir,'Luminal_Basal_MVG_classification_by_sample_prob_reclust.png'), width = 12, height = 8, units = 'in', res=200, bg = 'transparent')
  grid_arrange_list_shared_legend(gl, ncol=3, legend.position = 'right')
  dev.off()
  
  unique(dte$condition)
  condition_order <- c("Wildtype", "Pten-", "Erg+Pten-")
  gl <- lapply(condition_order, function(i){
    ggplot(dte[condition==i], aes(x=basal.Scores, y=luminal.Scores)) +
      geom_point_rast((aes(color=probability_luminal)), shape=1, size=.8, alpha=.5, stroke=1) +
      geom_density_2d(bins=10, linewidth=.5, alpha=.5, color='black') +
      #scale_fill_manual("Density", values = c("white",pals::ocean.amp(10))) +
      scale_fill_manual("Density", values = c(NA, pals::ocean.amp(10)), na.value = NA) +
      #scale_color_manual("Classification", values = luminal_basal_palette_continuous) +
      #scale_color_binned(type = "viridis", n.breaks=10) +
      scale_color_stepsn("Prob. cell is luminal",colours = luminal_basal_palette_continuous, n.breaks=10) +
      ggtitle(label = i, subtitle = "Epithelial cells") +
      scale_x_continuous(limits = c(-1.5,3)) + 
      scale_y_continuous(limits = c(-1.5,3)) + 
      geom_hline(yintercept = 0, linewidth=.1) + 
      geom_vline(xintercept = 0, linewidth=.1) + 
      coord_fixed() +
      theme_bw(12) +
      theme(panel.grid = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank(),
            legend.background = element_blank(),
            legend.key.size = unit(24, units='points'),
            legend.text = element_text(size = 12)) +
      guides(#colour = guide_legend(override.aes = list(size=10, alpha=1)),
        fill = guide_legend(override.aes = list(size=.5, alpha=1))) 
    #lemon::coord_capped_cart(expand = F)
  })
  png(filename = file.path(plot.dir,'Luminal_Basal_MVG_classification_by_condition_prob_reclust.png'), width = 12, height = 5, units = 'in', res=200, bg = 'transparent')
  grid_arrange_list_shared_legend(gl, ncol=3, legend.position = 'right')
  dev.off()
  
  
  gl <- lapply(sort(unique(dte$seurat_clusters)), function(i){
    ggplot(dte[seurat_clusters==i], aes(x=basal.Scores, y=luminal.Scores)) +
      geom_point_rast((aes(color=probability_luminal)), shape=1, size=.8, alpha=.5, stroke=1) +
      geom_density_2d(bins=10, linewidth=.5, alpha=.5, color='black') +
      scale_fill_manual("Density", values = c(NA, pals::ocean.amp(10)), na.value = NA) +
      scale_color_stepsn("Prob. cell is luminal",colours = luminal_basal_palette_continuous, n.breaks=10) +
      ggtitle(label = i, subtitle = "Epithelial cells") +
      scale_x_continuous(limits = c(-1.5,3)) + 
      scale_y_continuous(limits = c(-1.5,3)) + 
      geom_hline(yintercept = 0, linewidth=.1) + 
      geom_vline(xintercept = 0, linewidth=.1) + 
      theme_bw(12) + 
      #facet_wrap( ~condition) +
      theme(panel.grid = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank(),
            legend.background = element_blank(),
            legend.key.size = unit(24, units='points'),
            legend.text = element_text(size = 12)) +
      guides(#colour = guide_legend(override.aes = list(size=10, alpha=1)),
        fill = guide_legend(override.aes = list(size=.5, alpha=1)))
  })
  png(filename = file.path(plot.dir,'Luminal_Basal_MVG_classification_by_cluster_prob_reclust.png'), width = 12, height = 12, units = 'in', res=200, bg = 'transparent')
  grid_arrange_list_shared_legend(gl, ncol=3, legend.position = 'right')
  dev.off()


}


refineWouterL1L2MarkerClassificationEpiCells <- function(){
  library(MASS)
  library(plotly)
  
  # get sheets
  wouter_markers <- '/home/erik/projects/sawyers/foxa1/analysis/multiome_integrate/spreadsheets/Table S8 DEG Intact Ranksum global all cells.xlsx'
  (sheets <- readxl::excel_sheets(wouter_markers))
  data <- lapply(sheets, function(x){
    setDT(readxl::read_xlsx(wouter_markers, sheet = x))
  }) %>% set_names(paste0("",sheets))
  
  (sigs <- sheets[-1])
  # sigs <- c("Epi_Basal_1","Epi_Luminal_1","Epi_Luminal_2Psca")
  data <- lapply(sigs, function(x){
    # sort
    ds <- data[[x]][Sig_filter==T][order(AU_ROC, decreasing = T)][order(BHFDR_Q, decreasing = F)]
    ds[1:50]
  }) %>% set_names(sigs)
  
  for(i in grep("Epi", colnames(ss_epi@meta.data), ignore.case = T, value = T)){
    ss_epi[[i]] <- NULL
  }
  
  for(i in sigs){
    ss_epi <- AddModuleScore(object = ss_epi,
                             features = list(data[[i]]$GeneSymbol),
                             assay = 'RNA',
                             ctrl = 100,
                             name = i, # appends a 1
                             search = F) # searches HUGO containing human gene names only
  }
  grep("Epi", colnames(ss_epi@meta.data), ignore.case = T, value = T)
  
  grep(paste(sigs, collapse = '|'), colnames(ss_epi@meta.data), ignore.case = T, value = T)
  
  for(i in sigs)
    colnames(ss_epi@meta.data)[which(colnames(ss_epi@meta.data)==paste0(i,'1'))] <- i
  
  grep(paste(sigs, collapse = '|'), colnames(ss_epi@meta.data), ignore.case = T, value = T)
  
  dt <- data.table(ss_epi@meta.data[,grep(paste(sigs, collapse = '|'), colnames(ss_epi@meta.data))])
  
  ss_epi@meta.data$wouter_class <- colnames(dt)[max.col(dt)]
  
  dt <- data.table(ss_epi@meta.data[,grep(paste("^Epi_", collapse = '|'), colnames(ss_epi@meta.data))])
  
  dt <- dt[,list(Epi_Basal_1, Epi_Luminal_1, Epi_Luminal_2Psca)]
  
  ss_epi@meta.data$epi_class <- colnames(dt)[max.col(dt)]
  
  #### Umap ####
  g <- FeaturePlot(ss_epi, features = c('wouter_class'),
                   pt.size = 0.2, ncol = 1, slot = 'scale.data') %>%
    UmapPlotEnhance(obj = ss_epi, p1 = ., pt.size = .8,
                    scale.color = F, pal = wouter_pal, legend.title.rotation = 0, title = "",
                    legend.title = 'Wouter class', break.unit = 0.1)
  png(filename = file.path(plot.dir,'wouter_class_epi_umap.png'), width = 10, height = 5, units = 'in', res = 150, bg = 'transparent')
  g
  dev.off()
  
  glabels <- data.frame(ss_epi@reductions$umap@cell.embeddings, 
                        seurat_clusters=ss_epi$seurat_clusters)
  glabels$seurat_clusters <- as.character(glabels$seurat_clusters)
  glabels <- UmapPlotEnhance(obj = ss_epi, 
                             p1 = glabels, 
                             label.size = 5, 
                             pt.size = 0.2, addLabels = T, pal=NA)
  png(filename = file.path(plot.dir,'wouter_class_cluster_epi_umap.png'), width = 10, height = 5, units = 'in', res = 150, bg = 'transparent')
  g + glabels$layers
  dev.off()
  
  
  
  #### Optional ####
  ## Look at density
  dt <- as.data.table(ss_epi@meta.data, keep.rownames = 'cells')
  X <- as.matrix(dt[,list(Epi_Basal_1, Epi_Luminal_1, Epi_Luminal_2Psca)])
  class <- dt$epi_class
  
  #### Kmeans ####
  library(factoextra)
  fc <- fviz_nbclust(x = X, FUNcluster = kmeans)
  
  (cl <- kmeans(x = X, centers = 3))
  plot(X, col=cl$cluster)
  points(cl$centers, col=1:2, pch=8, cex=2)
  print(cl)
  # totss = total variance in data
  cl$totss
  # size of each cluster
  cl$size
  fviz_cluster(object = cl, data = X)
  
  
  #### MVGaussian ####
  library(mclust)
  mod <- densityMclust(X)
  summary(mod)
  plot(mod, what = "BIC")
  plot(mod, what = "density")
  
  BIC <- mclustBIC(X)
  plot(BIC)
  mod1 <- Mclust(X, x=BIC)
  plot(mod1, what = "classification")
  
  pdf(file = file.path(plot.dir,'EPI-class_diagnostics_epi.pdf'), width = 6, height = 6, onefile = T)
  plot(mod, what = "BIC")
  plot(mod, what = "density")
  plot(mod1, what = "classification")
  dev.off()
  
  
  #### Cluster uncertainty ####
  library(mclust)
  mod3 <- Mclust(X, G = 3)
  summary(mod3)
  # shows 1 is blue (basal) and 2 is yellow (luminal)
  pdf(file = file.path(plot.dir,'EPI-class_mclust_diagnostics.pdf'), width = 6, height = 6, onefile = T)
  plot(mod3, what = "BIC", main=T)
  plot(mod3, what = "classification", main=T, col=pals::alphabet())
  plot(mod3, what = "density", main=T)
  plot(mod3, what = "uncertainty", main=T, col=pals::alphabet())
  dev.off()
  
  uncerPlot(mod3$z)
  sort(mod3$uncertainty, decreasing = TRUE) %>% head(20)
  table(class, mod3$classification)
  
  ebasal_idx <- 1
  eluminal2_idx <- 2
  eluminal1_idx <- 3
  ss_epi$epi_class_classification <- ifelse(mod3$classification==ebasal_idx,'Epi_Basal_1',
                                            ifelse(mod3$classification==eluminal1_idx,'Epi_Luminal_1','Epi_Luminal_2Psca'))
  ss_epi$epi_class_classification_uncertainty <- mod3$uncertainty
  ss_epi$Epi_Basal_1_certainty <- mod3$z[,ebasal_idx]
  ss_epi$Epi_Luminal_1_certainty <- mod3$z[,eluminal1_idx]
  ss_epi$Epi_Luminal_2_certainty <- mod3$z[,eluminal2_idx]
  ss_epi$epi_class_classification_uncertainty <- mod3$uncertainty
  ss_epi$epi_class_classification_certainty <- 1-mod3$uncertainty
  
  
  ## by sample
  dt <- data.table(ss_epi@meta.data)
  xlim <- range(dt$Epi_Basal_1)
  ylim <- range(dt$Epi_Luminal_1)
  gp <- lapply(sort(unique(dt$bioNames)), function(x){
    idx <- which(dt$bioNames==x)
    set <- dt[idx]
    ggplot(dt, aes(x=Epi_Luminal_1, y=Epi_Basal_1)) +
      geom_point_rast(data=dt, (aes(color=epi_class_classification)), size=.8, alpha=0.2) +
      # show rings
      geom_density_2d(bins=15, linewidth=.2, color='white') +
      # fill density
      geom_density_2d_filled(data = set, contour_var = "ndensity", color=NA, alpha=0.7) +
      # fill colors
      #scale_fill_manual(values = rev(pals::kovesi.linear_kryw_5_100_c64(10))) +
      scale_color_manual("Classification", values = c("dodgerblue","chartreuse3","orangered1"), guide = 'legend') +
      scale_fill_manual("Density", values = c(NA, pals::brewer.bupu(9)), na.value = NA) +
      #scale_fill_manual("Density", values = c(NA, pals::brewer.reds(9)), na.value = NA) +
      ggtitle(x) +
      theme_bw(12) +
      theme(panel.grid = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank(),
            legend.background = element_blank(),
            legend.key.size = unit(16, units='points'),
            legend.text = element_text(size = 12), legend.position = 'right') +
      guides(fill = guide_legend(override.aes = list(size=4, alpha=1), reverse = TRUE),
             col = guide_legend(override.aes = list(size=5, alpha = 1, shape=20)))
  })
  png(filename = file.path(plot.dir,'Epi-class_MVG_classification_by_sample_markers_epi.png'), width = 11, height = 8, units = 'in', res=200, bg = 'transparent')
  grid_arrange_list_shared_legend(gp,nrow = 2, legend.position='right')
  dev.off()
  
  ## by celltype
  dt <- data.table(ss_epi@meta.data)
  xlim <- range(dt$basal.Scores)
  ylim <- range(dt$luminal.Scores)
  gp <- lapply(sort(unique(dt$cell_type)), function(x){
    idx <- which(dt$cell_type==x)
    set <- dt[idx]
    ggplot(dt, aes(x=basal.Scores, y=luminal.Scores)) +
      geom_point_rast(data=dt, (aes(color=epi_class_classification)), size=.8, alpha=0.2) +
      # show rings
      geom_density_2d(bins=15, linewidth=.2, color='white') +
      # fill density
      geom_density_2d_filled(data = set, contour_var = "ndensity", color=NA, alpha=0.7) +
      # fill colors
      #scale_fill_manual(values = rev(pals::kovesi.linear_kryw_5_100_c64(10))) +
      scale_color_manual("Classification", values = c("dodgerblue","chartreuse3","orangered1"), guide = 'legend') +
      scale_fill_manual("Density", values = c(NA, pals::brewer.bupu(9)), na.value = NA) +
      #scale_fill_manual("Density", values = c(NA, pals::brewer.reds(9)), na.value = NA) +
      ggtitle(x) +
      theme_bw(12) +
      theme(panel.grid = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank(),
            legend.background = element_blank(),
            legend.key.size = unit(16, units='points'),
            legend.text = element_text(size = 12), legend.position = 'right') +
      guides(fill = guide_legend(override.aes = list(size=4, alpha=1), reverse = TRUE),
             col = guide_legend(override.aes = list(size=5, alpha = 1, shape=20)))
  })
  png(filename = file.path(plot.dir,'Epi-class_MVG_classification_by_cell_type_markers_epi.png'), width = 12, height = 6, units = 'in', res=200, bg = 'transparent')
  grid_arrange_list_shared_legend(gp,nrow = 2, legend.position='right')
  dev.off()
  
  
  dt <- data.table(erg.scRNA@meta.data[,grep("Epi_Basal|Epi_Luminal|epi_class", colnames(erg.scRNA@meta.data))])
  ggplot(dt, aes(x=Epi_Basal_1, y=Epi_Luminal_1)) +
    geom_point(aes(color=epi_class)) +
    geom_density_2d(bins=40, size=.2, color='black') +
    theme_bw()
  
  class <- dt$epi_class
  X <- dt[,grep("Epi_Basal|Epi_Luminal", colnames(dt)), with=F]
  
  
  epiClass.palette <- c("dodgerblue","chartreuse3","orangered1")
  names(epiClass.palette) <- c("Epi_Basal_1","Epi_Luminal_1","Epi_Luminal_2Psca")
  # Original
  g1 <- FeaturePlot(erg.scRNA, features = c('epi_class'),
                    pt.size = 0.2, ncol = 1, slot = 'data') %>%
    UmapPlotEnhance(obj = erg.scRNA, p1 = .,
                    label.size = 5,
                    pal = epiClass.palette,
                    axis.label = "", pt.size = .4,
                    legend.title = 'Epi class score')
  # classified
  g2 <- FeaturePlot(erg.scRNA, features = c('epi_class_classification'),
                    pt.size = 0.2, ncol = 1, slot = 'scale.data') %>%
    UmapPlotEnhance(obj = erg.scRNA, p1 = .,
                    label.size = 5,
                    pal = epiClass.palette,
                    axis.label = "", pt.size = .4,
                    legend.title = 'Epi class classification')
  
  png(filename = file.path(plot.dir,'Epi_Class_Classification_Umap.png'),
      width = 7, height = 7.5, units = 'in',
      bg = 'transparent', res = 200)
  grid.arrange(g1, g2, nrow=2)
  dev.off()
  
  # classification uncertainty
  g3 <- FeaturePlot(erg.scRNA, features = c('epi_class_classification_uncertainty'),
                    pt.size = 0.2, ncol = 1, slot = 'scale.data') %>%
    UmapPlotEnhance(obj = erg.scRNA, p1 = ., pt.size = .8,
                    scale.color = F, pal = pals::brewer.orrd(10),
                    legend.title = 'Epi class uncertainty', break.unit = 0.1,
                    filename.prefix = 'Epi_class_classification_uncertainty_umap')
  g3
  
  dt <- data.table(erg.scRNA@meta.data[,c('bioNames', 'seurat_clusters', 'epi_class_classification_uncertainty', 'epi_class_classification')])
  position=position_dodge(1)
  pdf(file = file.path(plot.dir, 'epi_class_classification_metrics.pdf'), width = 6, height = 4, onefile = T)
  ggplot(dt, aes(x=epi_class_classification_uncertainty,
                 color=epi_class_classification,
                 fill=epi_class_classification)) +
    geom_histogram(aes(y=..density..), position='identity', alpha=0.1, binwidth=.02) + #,  binwidth=.01) +
    geom_density(alpha=0.2) +
    xlab("classification uncertainty") +
    #scale_fill_manual("", values = luminal_basal_palette) +
    #scale_color_manual("", values = luminal_basal_palette) +
    theme_bw(14) + theme(panel.grid = element_blank()) + lemon::coord_capped_cart(expand = F)
  
  ggplot(dt, aes(x=epi_class_classification_uncertainty, color=bioNames)) +
    #geom_histogram(aes(y=..density..), position='identity', alpha=0.5) + #,  binwidth=.01) +
    geom_density(alpha=0.2) +
    #scale_fill_manual("", values = sample_palette) +
    scale_color_manual("", values = sample_palette) +
    theme_bw(9) + theme(panel.grid = element_blank()) + lemon::coord_capped_cart(expand = F)
  
  ggplot(dt, aes(y=epi_class_classification_uncertainty, x=seurat_clusters, color=seurat_clusters)) +
    geom_boxplot(alpha=0.2) +
    scale_color_manual("", values = cluster_palette, guide=guide_legend(ncol = 2)) +
    theme_bw(9) + theme(panel.grid = element_blank()) + lemon::coord_capped_cart(ylim = c(0,1), expand = T)
  dev.off()
  
  
  ## Density of EPI-class
  UmapPlotDensity(obj = erg.scRNA, feature = 'epi_class_classification',
                  filename.prefix = 'EPI-classDensityPlots', nrow = 1, pt.size = 0.8)
  
  ### 2d Density ###
  g0 <- ggplot(dte, aes(x=probability_basal, y=probability_luminal)) +
    geom_point_rast((aes(color=luminal_basal_classification)), shape=1, size=.8, alpha=.8, stroke=1) +
    geom_density_2d_filled(bins=50, size=0, alpha=.5, color='black') +
    #scale_fill_manual("Density", values = c("white",pals::ocean.amp(10))) +
    #scale_fill_manual("Density", values = c("white",rev(pals::kovesi.linear_grey_10_95_c0(10)))) +
    scale_fill_manual("Density", values = c(NA, pals::ocean.amp(50)), na.value = NA) +
    scale_color_manual("Classification", values = luminal_basal_palette) +
    ggtitle("Epithelial cells") +
    scale_x_continuous(limits = c(0,1)) + 
    scale_y_continuous(limits = c(0,1)) + 
    geom_hline(yintercept = 0, size=.1) + 
    geom_vline(xintercept = 0, size=.1) + 
    theme_bw(12) + facet_wrap(~ bioNames) +
    theme(panel.grid = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          legend.background = element_blank(),
          legend.key.size = unit(24, units='points'),
          legend.text = element_text(size = 12)) +
    guides(colour = guide_legend(override.aes = list(size=4, alpha=1)),
           fill = guide_legend(override.aes = list(size=.5, alpha=1))) +
    lemon::coord_capped_cart(expand = F)
  png(filename = file.path(plot.dir,'Luminal_Basal_MVG_Epi_classification_eval.png'), width = 7.5, height = 6, units = 'in', res=200, bg = 'transparent')
  print(g0)
  dev.off()
  
  ### Bar plot ####
  PlotBarVariableInClusters(seuratObj = erg.scRNA, variable = 'wouter_class', palette = wouter_pal,
                            orientation = "vertical", width_multiplier = 1.8, filename = "Wouter_barplot")
  
  
  # Barplot cluster percentage by each and mark cluster as luminal or basal -like
  dt <- setDT(erg.scRNA@meta.data[,c('bioNames','epi_class_classification')])[
    ,sum:=.N, by=list(bioNames)][
      ,sum_group := .N, by = c("bioNames", "epi_class_classification")][, prop := sum_group/sum] %>% unique
  dtsum <- unique(dt[,list(bioNames, sum)])
  dt <- dt[order(sum, decreasing = F)]
  dt[,bioNames:=factor(x = bioNames, levels = unique(bioNames), ordered = T)]
  g0 <- ggplot(dt, aes(x=bioNames, y=prop)) +
    geom_bar(stat='identity', aes(fill=epi_class_classification)) +
    geom_hline(yintercept = .50, linetype='dashed', size=.5, color='gray30') +
    geom_text(data = dtsum, aes( label = paste0("n=",sum), y=1+.1, x=bioNames ), size=3) +
    scale_y_continuous(breaks = seq(0,1,by=.25), labels = scales::percent(seq(0,1,by=.25))) +
    scale_fill_manual("Cell type",values = epiClass.palette) +
    xlab("Bionames") + ylab("Cell % within Sample") +
    coord_flip(ylim = c(-0.01, 1.0), expand = F, clip='off') +
    theme_bw(12) + theme(plot.background = element_blank(),
                         panel.grid = element_blank(),
                         legend.background = element_blank(),
                         legend.box.margin=margin(10,4,10,20))
  cairo_pdf(filename = file.path(plot.dir,'bioNames_proportion_EPI-class_markers.pdf'), width = 5.8, height = 1.8, bg = 'transparent')
  print(g0)
  dev.off()
  
  # Barplot cluster percentage by each and mark cluster as luminal or basal -like
  dt <- setDT(erg.scRNA@meta.data[,c('seurat_clusters','epi_class_classification')])[
    ,sum:=.N, by=seurat_clusters][
      ,sum_group := .N, by = c("seurat_clusters", "epi_class_classification")][, prop := sum_group/sum] %>% unique
  dtsum <- dt[,list(seurat_clusters, sum, prop=1)]
  g1 <- ggplot(dt, aes(x=reorder(seurat_clusters, sum), y=prop)) +
    geom_bar(stat='identity', aes(fill=epi_class_classification)) +
    geom_hline(yintercept = .50, linetype='dashed', size=.5, color='gray30') +
    geom_text(data = dtsum, aes( label = paste0("n=",sum), y=prop, x=seurat_clusters ), size=3, hjust=-0.1) +
    scale_y_continuous(breaks = seq(0,1,by=.25), labels = scales::percent(seq(0,1,by=.25))) +
    scale_fill_manual("Cell type",values = epiClass.palette) +
    xlab("Cluster") + ylab("Cell % within cluster") +
    coord_flip(ylim = c(-0.01, 1.0), expand = F, clip='off') +
    theme_bw(12) + theme(plot.background = element_blank(),
                         panel.grid = element_blank(),
                         legend.background = element_blank(),
                         legend.box.margin=margin(10,4,10,20))
  cairo_pdf(filename = file.path(plot.dir,'cluster_proportion_EPI-class_markers.pdf'), width = 4, height = 5, bg = 'transparent')
  print(g1)
  dev.off()
  
  ### Save Seurat object ####
  saveRDS(erg.scRNA, file = scrna.object)
  
}

refineWouterL1L2MarkerClassificationAllCells <- function(){
  library(MASS)
  library(plotly)
  
  # get sheets
  wouter_markers <- '/home/erik/projects/sawyers/foxa1/analysis/multiome_integrate/spreadsheets/Table S8 DEG Intact Ranksum global all cells.xlsx'
  (sheets <- readxl::excel_sheets(wouter_markers))
  data <- lapply(sheets, function(x){
    setDT(readxl::read_xlsx(wouter_markers, sheet = x))
  }) %>% set_names(paste0("",sheets))
  
  (sigs <- sheets[-1])
  # sigs <- c("Epi_Basal_1","Epi_Luminal_1","Epi_Luminal_2Psca")
  data <- lapply(sigs, function(x){
    # sort
    ds <- data[[x]][Sig_filter==T][order(AU_ROC, decreasing = T)][order(BHFDR_Q, decreasing = F)]
    ds[1:50]
  }) %>% set_names(sigs)
  
  #for(i in grep("Epi", colnames(erg.scRNA@meta.data), ignore.case = T, value = T)){
  #  erg.scRNA[[i]] <- NULL
  #}
  
  for(i in sigs){
    ss_epi <- AddModuleScore(object = erg.scRNA,
                             features = list(data[[i]]$GeneSymbol),
                             assay = 'RNA',
                             ctrl = 100,
                             name = i, # appends a 1
                             search = F) # searches HUGO containing human gene names only
  }
  grep("Epi", colnames(erg.scRNA@meta.data), ignore.case = T, value = T)
  
  grep(paste(sigs, collapse = '|'), colnames(erg.scRNA@meta.data), ignore.case = T, value = T)
  
  for(i in sigs)
    colnames(erg.scRNA@meta.data)[which(colnames(erg.scRNA@meta.data)==paste0(i,'1'))] <- i
  
  grep(paste(sigs, collapse = '|'), colnames(erg.scRNA@meta.data), ignore.case = T, value = T)
  
  dt <- data.table(erg.scRNA@meta.data[,grep(paste(sigs, collapse = '|'), colnames(erg.scRNA@meta.data))])
  
  erg.scRNA@meta.data$wouter_class <- colnames(dt)[max.col(dt)]
  
  #### Umap ####
  g <- FeaturePlot(erg.scRNA, features = c('wouter_class'),
                   pt.size = 0.2, ncol = 1, slot = 'scale.data') %>%
    UmapPlotEnhance(obj = erg.scRNA, p1 = ., pt.size = .8,
                    scale.color = F, pal = wouter_pal, legend.title.rotation = 0, title = "",
                    legend.title = 'Wouter class', break.unit = 0.1)
  png(filename = file.path(plot.dir,'wouter_class_umap.png'), width = 10, height = 5, units = 'in', res = 150, bg = 'transparent')
  g
  dev.off()
  
  glabels <- data.frame(erg.scRNA@reductions$umap@cell.embeddings, seurat_clusters=erg.scRNA$seurat_clusters)
  glabels$seurat_clusters <- as.character(glabels$seurat_clusters)
  glabels <- UmapPlotEnhance(obj = erg.scRNA, p1 = glabels, 
                             label.size = 5, 
                             pt.size = 0.2, addLabels = T, pal=NA)
  png(filename = file.path(plot.dir,'wouter_class_cluster_umap.png'), width = 10, height = 5, units = 'in', res = 150, bg = 'transparent')
  g + glabels$layers
  dev.off()
  
  #### Optional ####
  ## Look at density
  den3d <- kde2d(X$Epi_Basal_1, X$Epi_Luminal_1, X$Epi_Luminal_2Psca)
  persp(den3d, box=FALSE)
  
  plot_ly(x=den3d$x, y=den3d$y, z=den3d$z) %>%
    add_surface() %>%
    layout(scene = list(xaxis=list(title = 'Epi-Basal1'),
                        yaxis=list(title = 'Epi-Luminal1'),
                        zaxis=list(title = 'Epi-Luminal2')))
  
  #### Kmeans ####
  library(factoextra)
  fc <- fviz_nbclust(x = X, FUNcluster = kmeans)
  
  (cl <- kmeans(x = X, centers = 3))
  plot(X, col=cl$cluster)
  points(cl$centers, col=1:2, pch=8, cex=2)
  print(cl)
  # totss = total variance in data
  cl$totss
  # size of each cluster
  cl$size
  fviz_cluster(object = cl, data = X)
  
  
  #### MVGaussian ####
  library(mclust)
  mod <- densityMclust(X)
  summary(mod)
  plot(mod, what = "BIC")
  plot(mod, what = "density")
  
  BIC <- mclustBIC(X)
  mod1 <- Mclust(X, x=BIC)
  plot(mod1, what = "classification")
  table(class, mod1$classification)
  
  pdf(file = file.path(plot.dir,'EPI-class_diagnostics.pdf'), width = 6, height = 6, onefile = T)
  plot(mod, what = "BIC")
  plot(mod, what = "density")
  plot(mod1, what = "classification")
  dev.off()
  
  
  #### Cluster uncertainty ####
  library(mclust)
  mod3 <- Mclust(X, G = 3, prior = priorControl())
  summary(mod3)
  # shows 1 is blue (basal) and 2 is yellow (luminal)
  
  table(class, mod3$classification)
  lum_idx <- 1
  pdf(file = file.path(plot.dir,'EPI-class_mclust_diagnostics.pdf'), width = 6, height = 6, onefile = T)
  plot(mod3, what = "BIC", main=T)
  plot(mod3, what = "classification", main=T, col=pals::alphabet())
  plot(mod3, what = "density", main=T)
  plot(mod3, what = "uncertainty", main=T, col=pals::alphabet())
  dev.off()
  
  uncerPlot(mod2$z)
  sort(mod2$uncertainty, decreasing = TRUE) %>% head(20)
  
  
  erg.scRNA$epi_class_classification <- ifelse(mod3$classification==1,'Epi_Luminal_1',
                                               ifelse(mod3$classification==3,'Epi_Basal_1','Epi_Luminal_2Psca'))
  erg.scRNA$epi_class_classification_uncertainty <- mod3$uncertainty
  
  ## by sample
  dt <- data.table(erg.scRNA@meta.data)
  xlim <- range(dt$Epi_Basal_1)
  ylim <- range(dt$Epi_Luminal_1)
  gp <- lapply(sort(unique(dt$bioNames)), function(x){
    idx <- which(dt$bioNames==x)
    set <- dt[idx]
    ggplot(dt, aes(x=Epi_Luminal_1, y=Epi_Basal_1)) +
      geom_point_rast(data=dt, (aes(color=epi_class_classification)), size=.8, alpha=0.2) +
      # show rings
      geom_density_2d(bins=15, size=.2, color='white') +
      # fill density
      geom_density_2d_filled(data = set, contour_var = "ndensity", color=NA, alpha=0.7) +
      # fill colors
      #scale_fill_manual(values = rev(pals::kovesi.linear_kryw_5_100_c64(10))) +
      scale_color_manual("Classification", values = c("dodgerblue","chartreuse3","orangered1"), guide = 'legend') +
      scale_fill_manual("Density", values = c(NA, pals::brewer.bupu(9)), na.value = NA) +
      #scale_fill_manual("Density", values = c(NA, pals::brewer.reds(9)), na.value = NA) +
      ggtitle(x) +
      theme_bw(12) +
      theme(panel.grid = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank(),
            legend.background = element_blank(),
            legend.key.size = unit(16, units='points'),
            legend.text = element_text(size = 12), legend.position = 'right') +
      guides(fill = guide_legend(override.aes = list(size=4, alpha=1), reverse = TRUE),
             col = guide_legend(override.aes = list(size=5, alpha = 1, shape=20)))
  })
  png(filename = file.path(plot.dir,'Epi-class_MVG_classification_by_sample_markers.png'), width = 11, height = 8, units = 'in', res=200, bg = 'transparent')
  grid_arrange_list_shared_legend(gp,nrow = 2, legend.position='right')
  dev.off()
  
  
  dt <- data.table(erg.scRNA@meta.data[,grep("Epi_Basal|Epi_Luminal|epi_class", colnames(erg.scRNA@meta.data))])
  ggplot(dt, aes(x=Epi_Basal_1, y=Epi_Luminal_1)) +
    geom_point(aes(color=epi_class)) +
    geom_density_2d(bins=40, size=.2, color='black') +
    theme_bw()
  
  class <- dt$epi_class
  X <- dt[,grep("Epi_Basal|Epi_Luminal", colnames(dt)), with=F]
  
  
  epiClass.palette <- c("dodgerblue","chartreuse3","orangered1")
  names(epiClass.palette) <- c("Epi_Basal_1","Epi_Luminal_1","Epi_Luminal_2Psca")
  # Original
  g1 <- FeaturePlot(erg.scRNA, features = c('epi_class'),
                    pt.size = 0.2, ncol = 1, slot = 'data') %>%
    UmapPlotEnhance(obj = erg.scRNA, p1 = .,
                    label.size = 5,
                    pal = epiClass.palette,
                    axis.label = "", pt.size = .4,
                    legend.title = 'Epi class score')
  # classified
  g2 <- FeaturePlot(erg.scRNA, features = c('epi_class_classification'),
                    pt.size = 0.2, ncol = 1, slot = 'scale.data') %>%
    UmapPlotEnhance(obj = erg.scRNA, p1 = .,
                    label.size = 5,
                    pal = epiClass.palette,
                    axis.label = "", pt.size = .4,
                    legend.title = 'Epi class classification')
  
  png(filename = file.path(plot.dir,'Epi_Class_Classification_Umap.png'),
      width = 7, height = 7.5, units = 'in',
      bg = 'transparent', res = 200)
  grid.arrange(g1, g2, nrow=2)
  dev.off()
  
  # classification uncertainty
  g3 <- FeaturePlot(erg.scRNA, features = c('epi_class_classification_uncertainty'),
                    pt.size = 0.2, ncol = 1, slot = 'scale.data') %>%
    UmapPlotEnhance(obj = erg.scRNA, p1 = ., pt.size = .8,
                    scale.color = F, pal = pals::brewer.orrd(10),
                    legend.title = 'Epi class uncertainty', break.unit = 0.1,
                    filename.prefix = 'Epi_class_classification_uncertainty_umap')
  g3
  
  dt <- data.table(erg.scRNA@meta.data[,c('bioNames', 'seurat_clusters', 'epi_class_classification_uncertainty', 'epi_class_classification')])
  position=position_dodge(1)
  pdf(file = file.path(plot.dir, 'epi_class_classification_metrics.pdf'), width = 6, height = 4, onefile = T)
  ggplot(dt, aes(x=epi_class_classification_uncertainty,
                 color=epi_class_classification,
                 fill=epi_class_classification)) +
    geom_histogram(aes(y=..density..), position='identity', alpha=0.1, binwidth=.02) + #,  binwidth=.01) +
    geom_density(alpha=0.2) +
    xlab("classification uncertainty") +
    #scale_fill_manual("", values = luminal_basal_palette) +
    #scale_color_manual("", values = luminal_basal_palette) +
    theme_bw(14) + theme(panel.grid = element_blank()) + lemon::coord_capped_cart(expand = F)
  
  ggplot(dt, aes(x=epi_class_classification_uncertainty, color=bioNames)) +
    #geom_histogram(aes(y=..density..), position='identity', alpha=0.5) + #,  binwidth=.01) +
    geom_density(alpha=0.2) +
    #scale_fill_manual("", values = sample_palette) +
    scale_color_manual("", values = sample_palette) +
    theme_bw(9) + theme(panel.grid = element_blank()) + lemon::coord_capped_cart(expand = F)
  
  ggplot(dt, aes(y=epi_class_classification_uncertainty, x=seurat_clusters, color=seurat_clusters)) +
    geom_boxplot(alpha=0.2) +
    scale_color_manual("", values = cluster_palette, guide=guide_legend(ncol = 2)) +
    theme_bw(9) + theme(panel.grid = element_blank()) + lemon::coord_capped_cart(ylim = c(0,1), expand = T)
  dev.off()
  
  
  ## Density of EPI-class
  UmapPlotDensity(obj = erg.scRNA, feature = 'epi_class_classification',
                  filename.prefix = 'EPI-classDensityPlots', nrow = 1, pt.size = 0.8)
  
  ### 2d Density ###
  g0 <- ggplot(dte, aes(x=probability_basal, y=probability_luminal)) +
    geom_point_rast((aes(color=luminal_basal_classification)), shape=1, size=.8, alpha=.8, stroke=1) +
    geom_density_2d_filled(bins=50, size=0, alpha=.5, color='black') +
    #scale_fill_manual("Density", values = c("white",pals::ocean.amp(10))) +
    #scale_fill_manual("Density", values = c("white",rev(pals::kovesi.linear_grey_10_95_c0(10)))) +
    scale_fill_manual("Density", values = c(NA, pals::ocean.amp(50)), na.value = NA) +
    scale_color_manual("Classification", values = luminal_basal_palette) +
    ggtitle("Epithelial cells") +
    scale_x_continuous(limits = c(0,1)) + 
    scale_y_continuous(limits = c(0,1)) + 
    geom_hline(yintercept = 0, size=.1) + 
    geom_vline(xintercept = 0, size=.1) + 
    theme_bw(12) + facet_wrap(~ bioNames) +
    theme(panel.grid = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          legend.background = element_blank(),
          legend.key.size = unit(24, units='points'),
          legend.text = element_text(size = 12)) +
    guides(colour = guide_legend(override.aes = list(size=4, alpha=1)),
           fill = guide_legend(override.aes = list(size=.5, alpha=1))) +
    lemon::coord_capped_cart(expand = F)
  png(filename = file.path(plot.dir,'Luminal_Basal_MVG_Epi_classification_eval.png'), width = 7.5, height = 6, units = 'in', res=200, bg = 'transparent')
  print(g0)
  dev.off()
  
  ### Bar plot ####
  PlotBarVariableInClusters(seuratObj = erg.scRNA, variable = 'wouter_class', palette = wouter_pal,
                            orientation = "vertical", width_multiplier = 1.8, filename = "Wouter_barplot")
  
  
  # Barplot cluster percentage by each and mark cluster as luminal or basal -like
  dt <- setDT(erg.scRNA@meta.data[,c('bioNames','epi_class_classification')])[
    ,sum:=.N, by=list(bioNames)][
      ,sum_group := .N, by = c("bioNames", "epi_class_classification")][, prop := sum_group/sum] %>% unique
  dtsum <- unique(dt[,list(bioNames, sum)])
  dt <- dt[order(sum, decreasing = F)]
  dt[,bioNames:=factor(x = bioNames, levels = unique(bioNames), ordered = T)]
  g0 <- ggplot(dt, aes(x=bioNames, y=prop)) +
    geom_bar(stat='identity', aes(fill=epi_class_classification)) +
    geom_hline(yintercept = .50, linetype='dashed', size=.5, color='gray30') +
    geom_text(data = dtsum, aes( label = paste0("n=",sum), y=1+.1, x=bioNames ), size=3) +
    scale_y_continuous(breaks = seq(0,1,by=.25), labels = scales::percent(seq(0,1,by=.25))) +
    scale_fill_manual("Cell type",values = epiClass.palette) +
    xlab("Bionames") + ylab("Cell % within Sample") +
    coord_flip(ylim = c(-0.01, 1.0), expand = F, clip='off') +
    theme_bw(12) + theme(plot.background = element_blank(),
                         panel.grid = element_blank(),
                         legend.background = element_blank(),
                         legend.box.margin=margin(10,4,10,20))
  cairo_pdf(filename = file.path(plot.dir,'bioNames_proportion_EPI-class_markers.pdf'), width = 5.8, height = 1.8, bg = 'transparent')
  print(g0)
  dev.off()
  
  # Barplot cluster percentage by each and mark cluster as luminal or basal -like
  dt <- setDT(erg.scRNA@meta.data[,c('seurat_clusters','epi_class_classification')])[
    ,sum:=.N, by=seurat_clusters][
      ,sum_group := .N, by = c("seurat_clusters", "epi_class_classification")][, prop := sum_group/sum] %>% unique
  dtsum <- dt[,list(seurat_clusters, sum, prop=1)]
  g1 <- ggplot(dt, aes(x=reorder(seurat_clusters, sum), y=prop)) +
    geom_bar(stat='identity', aes(fill=epi_class_classification)) +
    geom_hline(yintercept = .50, linetype='dashed', size=.5, color='gray30') +
    geom_text(data = dtsum, aes( label = paste0("n=",sum), y=prop, x=seurat_clusters ), size=3, hjust=-0.1) +
    scale_y_continuous(breaks = seq(0,1,by=.25), labels = scales::percent(seq(0,1,by=.25))) +
    scale_fill_manual("Cell type",values = epiClass.palette) +
    xlab("Cluster") + ylab("Cell % within cluster") +
    coord_flip(ylim = c(-0.01, 1.0), expand = F, clip='off') +
    theme_bw(12) + theme(plot.background = element_blank(),
                         panel.grid = element_blank(),
                         legend.background = element_blank(),
                         legend.box.margin=margin(10,4,10,20))
  cairo_pdf(filename = file.path(plot.dir,'cluster_proportion_EPI-class_markers.pdf'), width = 4, height = 5, bg = 'transparent')
  print(g1)
  dev.off()
  
  ### Save Seurat object ####
  saveRDS(erg.scRNA, file = scrna.object)
  
}

BuildEpiSeurat <- function(obj){
  # subset for epithelial cluster class
  ss_epi <- subset(obj, seurat_clusters_class == 'epithelial')
  # keep epi wouter
  ss_epi <- subset(ss_epi, wouter_class %in% c("Epi_Basal_1","Epi_Luminal_1","Epi_Luminal_2Psca","Epi_Luminal_3Foxi1"))
  table(ss_epi$seurat_clusters)
  # remove SV
  ##ss_epi <- subset(ss_epi, seurat_clusters != 20) # remove SV
  
  # recluster data
  # evaluate umap
  ss_epi <- RunUMAP(ss_epi, dims= 1:30, n.neighbors = 30, 
                    seed.use = 42,
                    min.dist = 0.3, spread = 1, verbose = FALSE, 
                    assay = 'RNA', metric = 'cosine')
  
  ss_epi@reductions
  
  #### Clusters for EPI ####
  ss_epi <- FindNeighbors(ss_epi,
                          dims = 1:30, 
                          reduction='pca',
                          annoy.metric = 'cosine',
                          assay = 'SCT',
                          verbose = FALSE)
  
  saveRDS(ss_epi, epi.scrna.object)
  
  ## Look at the UMAPS below ##
  
  
  res <- seq(from = 0.1, to = 2.0, by = 0.1)
  ss_epi <- FindClusters(ss_epi, 
                         verbose = FALSE, algorithm = 4, 
                         resolution = res)
  
  res <- grep("SCT_snn_res.", colnames(ss_epi@meta.data), value = T)
  (res <- res[order(gsub("SCT_snn_res.","",res))])
  
  sils <- lapply(res, function(r){
    sil <- silhouetteScore(obj = ss_epi, reduction = 'umap', # new umap generated by me
                           clusters =  r, dims = 1:2)
    #summary(sil)$avg.width
  }) %>% set_names(res)
  
  um <- lapply(res, function(r){
    score <- summary(sils[[r]])$avg.width
    g1 <- DimPlot(ss_epi, label = TRUE, reduction = "umap", group.by = r) + NoLegend()
    g1 <- UmapPlotEnhance(obj = ss_epi, 
                          outline.size = 0.2, 
                          outline.color = 'gray50',
                          p1 = g1, 
                          coord1 = "umap_1",
                          coord2 = "umap_2",
                          label.size = 10,
                          title.size = 12,
                          pal=godsnot_102,
                          title = paste("Cluster res.:", r,
                                        "\nmean silhouette ", formatC(score, digits = 3)))
  }) %>% set_names(res)
  
  
  png(filename = file.path(plot.dir, 'epi_silhouette_scores_regress.png'), width = 24, height = 18, units = 'in', res=300, bg = 'transparent')
  do.call(grid.arrange, c(um, nrow=4, top='Clustering resolutions'))
  dev.off()
  
  pdf(file = file.path(plot.dir, 'epi_silhouette_scores_regress.pdf'), width = 10, height = 8, onefile=T)
  for(r in res){
    g <- um[[r]]
    p1 <- ~plot(sils[[r]], main=paste(r), cex.names = par("cex.axis")) # RStudio sometimes does not display silhouette plots correctly
    grid.arrange(cowplot::as_grob(p1), g, nrow=1, widths=c(.35,.65))
  }
  dev.off()
  
  
  
  # Plot UMAPS
  glabel <- DimPlot(ss_epi, group.by = 'seurat_clusters', reduction = 'umap')
  glabel <- UmapPlotEnhance(obj = ss_epi, p1 = glabel, addLabels = T, pt.size = 1,
                            coord1 = 'umap_1', coord2 = 'umap_2',
                            label.size = 5, pal=NA,
                            axis.label = NA, title="", legend.title = '', 
                            scale.color = F)
  
  g1 <- DimPlot(ss_epi, group.by = 'bioNames', reduction = 'umap')
  g1 <- UmapPlotEnhance(obj = ss_epi, p1 = g1, addLabels = F, pt.size = 1,
                        coord1 = 'umap_1', coord2 = 'umap_2',
                        label.size = 10, pal=bionames_palette_list,
                        axis.label = NA, title="", legend.title = '', 
                        scale.color = F, filename.prefix = 'UMAP_epi_bioNames')
  
  
  #### Pick clustering resolution ####
  cresolution <- 0.3
  (cresolution <- paste0("SCT_snn_res.",cresolution))
  
  ss_epi@meta.data$seurat_clusters <- ss_epi@meta.data[,cresolution]
  identical(ss_epi@meta.data$seurat_clusters, ss_epi@meta.data[[cresolution]])
  table(ss_epi@meta.data$seurat_clusters)
  is.factor(ss_epi$seurat_clusters)
  
  g1 <- DimPlot(ss_epi, label = F, group.by = 'seurat_clusters', reduction = 'umap')
  UmapPlotEnhance(obj = ss_epi, p1 = g1, addLabels = F, 
                  coord1 = 'umap_1', coord2 = 'umap_2',
                  label.size = 10, pal=cluster_palette,
                  axis.label = NA, title="", legend.title = '', 
                  scale.color = F, filename.prefix = 'UMAP_clusters')
  
  
  #### Save EPI object ####
  #ss_epi$seurat_clusters <- factor(x = ss_epi$seurat_clusters, levels = c(1,2,10,3,4,9,11,5,6,7,12,8))
  ss_epi$seurat_clusters <- factor(x = ss_epi$seurat_clusters, levels = c(2, 1, 10, 6, 7, 12, 5, 9, 3, 4, 11, 8))
  
  saveRDS(ss_epi, epi.scrna.object)
  
  #### Load reclustered EPI object ####
  ss_epi <- readRDS(epi.scrna.object)
  
}

GetEpiSeurat <- function(obj){
  ss_epi <- subset(obj, seurat_clusters_class == 'epithelial')
  ss_epi <- subset(ss_epi, seurat_clusters != 20) # remove SV
}

#### Plots ####
UmapPlotEnhance <- function(obj,
                            p1,
                            addLabels=F,
                            axis.label=NA, #"orig.ident",
                            outline.size=.1,
                            outline.color='black', #'grey80',
                            na.outline.size=.05,
                            pal=NA,
                            panel.border.size=0,
                            show.legend=T,
                            label.size=4,
                            title.size=14,
                            scale.color=F,
                            pt.size=0.8,
                            pt.alpha=1,
                            title=NA,
                            legend.title="",
                            legend.title.rotation=90,
                            legend.title.size=10,
                            legend.height.factor=.4,
                            legend.position="right",
                            legend.justification="top",
                            legend.direction="vertical",
                            legend.text.size=12,
                            legend.show.na=T,
                            break.unit=1,
                            fixed_coord=T,
                            panel_width_in=3,
                            vmin=NA,
                            vmax=NA,
                            na.color="gray80",
                            filename.prefix,
                            nrow=NULL,
                            coord1='UMAP_1',
                            coord2='UMAP_2',
                            random_order=F,
                            small_axis_title=T,
                            resolution=150,
                            as.pdf=F){
  # Plot the integration UMAP
  require(shadowtext)
  require(ggrastr)
  require(ggh4x)
  require(gtable)
  
  dt <- NULL
  
  if(is.data.frame(p1))
    dt <- data.table(p1, keep.rownames = 'rn')
  else
    dt <- data.table(p1$data, keep.rownames = 'rn')
  cols <- colnames(dt)
  
  panel_width <- unit(panel_width_in, "in")
  if(is.na(title)){
    title <- data.table::last(cols)
  }
  
  if(!coord1 %in% colnames(dt)){
    coord1 <- colnames(dt)[2]
    coord2 <- colnames(dt)[3]
  }
  
  setnames(dt, old = c(coord1,coord2), new = c("x","y"))
  setnames(dt, data.table::last(colnames(dt)), new = c("color"))
  
  ## This creates 0 valued points not in faceting window
  if(axis.label %in% colnames(obj@meta.data)){
    dt[,axis.label := obj@meta.data[dt$rn,eval(axis.label)]]
    dt <- lapply(unique(dt$axis.label), function(i){
      ds <- copy(dt[axis.label!=i])
      dt_zero <- copy(ds)
      dt_zero$color <- NA
      dt_zero$axis.label <- i
      rbindlist(list(dt_zero, ds))
    }) %>% rbindlist
  }
  
  if(random_order)
    dt <- dt[order(sample(color), na.last = F)]
  else
    dt <- dt[order(color, decreasing = F, na.last = F)]
  
  p1 <- NULL
  if(all(is.na(pal))){
    p1 <- ggplot(dt, aes(x=x, y=y)) #, fill=color
    outline.size <- 0
    
    p1 <- p1 + geom_point_rast(size=pt.size, shape=21, stroke=outline.size, alpha=pt.alpha) +
      theme_bw()
  }
  else{
    if(length(grep("_multiome", dt$color))){
      dt[,color := gsub("_multiome","",color)]
      if(!is.null(names(pal)))
        names(pal) <- gsub("_multiome","",names(pal))
    }
    p1 <- ggplot(dt, aes(x=x, y=y, fill=color))
    
    if(outline.size > 0){
      # this is na outline
      p1 <- p1 + geom_point_rast(size=pt.size, shape=1, stroke=na.outline.size, col=outline.color) # standard stroke #'grey40'
      #  this is NA fill
      p1 <- p1 + geom_point_rast(data=dt[is.na(color)], size=pt.size, shape=21, stroke=0, alpha=pt.alpha)
      
      # this is non NA outline
      p1 <- p1 + geom_point_rast(data=dt[!is.na(color)],
                                 size=pt.size, shape=1, stroke=outline.size, col=outline.color) # standard stroke #'grey40'
    }
    # Non na fill
    p1 <- p1 + geom_point_rast(data=dt[!is.na(color)], size=pt.size, shape=21, stroke=0, alpha=pt.alpha)
  }
  
  p1 <- p1 + xlab( gsub("_","",coord1) ) + ylab( gsub("_","",coord2) ) +
    theme_bw() +
    theme(text = element_text(size = 15),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(size=title.size, hjust = 0.5),
          legend.text = element_text(size = legend.text.size),
          panel.border = element_rect(linewidth = panel.border.size, color='black'),
          strip.background = element_blank(),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          plot.margin = margin(t=0,r=0,b=0,l=0),
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.direction = legend.direction,
          legend.position = legend.position,
          legend.justification = legend.justification,
          legend.title = element_text(angle = legend.title.rotation, size = legend.title.size),
          legend.background = element_blank(), # get rid of legend bg element_rect(fill = 'transparent')
          legend.box.background = element_blank(), # get rid o
          legend.key=element_blank(), 
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(t=0,r=0,b=0,l=0)) #trbl
  if(fixed_coord)
    p1 <- p1 + coord_fixed(clip = 'off')
  
  if(addLabels){
    dt[,label := color]
    median.labels <- dt[,.(x=median(x),y=median(y), color), by=label]
    median.labels <- median.labels[order(label)]
    median.labels <- unique(median.labels)
    p1 <- p1 + geom_shadowtext(data = median.labels,
                               aes(x=x, y=y, label=label, fontface='bold'),
                               color='black',
                               bg.color='white',
                               bg.r=0.1,
                               label.padding=0.01,
                               size=label.size)
  }
  
  #  get dimensions of plot
  gt <- ggplotGrob(p1)
  panel <- gtable_filter(gt, "panel")
  # height is scaled compared to width
  panel_height_in <- panel_width_in*as.numeric(gsub("null","",panel$heights))
  # set the panel height based on width
  panel_height <- unit(panel_height_in, "in") # reduce for title
  
  p1 <- p1 + ggtitle(label = title) + theme(plot.title = element_text(size = title.size))
  
  # Do we want color fill?
  if(all(is.na(pal))){
    ##nothing
  }
  else {
    if(is.numeric(dt$color)){
      message("continuous color")
      if(missing(pal))
        pal <- pals::brewer.purples(10)
      
      my.breaks <- NULL
      
      lim <- range(dt$color, na.rm = T)
      
      message("value limits [",lim[1],":",lim[2],"]")
      
      if(!is.na(vmin))
        #lim[1] <- min(dt[color>=vmin]$color, na.rm = T)
        lim[1] <- min(vmin, lim, na.rm = T)
      
      if(!is.na(vmax))
        lim[2] <- max(vmax, lim, na.rm = T)
      
      if(scale.color){
        lim <- c(-max(abs(lim), na.rm = T), max(abs(lim), na.rm = T))
        numbreaks <- 6
        if(lim[2] > 1)
          break.unit <- ceiling(lim[2]/numbreaks)
        else
          break.unit <- round(lim[2]/numbreaks, digits=2)
        my.breaks <- seq(from=-break.unit*numbreaks,
                         to=break.unit*numbreaks,
                         by = break.unit)
      }
      else if(lim[1] > 0)
        lim <- c(0,max(abs(lim), na.rm = T))
      
      if(is.null(break.unit)){
        break.unit <- 1
        if( identical( range(c(0,1)), lim ) )
          break.unit <- 0.2
      }
      if(break.unit > lim[2]){
        break.unit <- 0.2
      }
      if(is.null(my.breaks)){
        my.breaks <- seq(from=lim[1], to=lim[2], by=break.unit)
      }
      
      lim <- c(floor(1000*lim[1])/1000, ceiling(1000*lim[2])/1000)
      message("scale limits [",lim[1],":",lim[2],"]")
      
      p1 <- p1 + scale_fill_gradientn(legend.title,
                                      guide = "colourbar",
                                      colours = pal,
                                      na.value = na.color,
                                      limits = lim,
                                      breaks = my.breaks)
      
      legend.title.position <- ifelse(legend.title=="","right","left")
      p1 <- p1 + guides(fill = guide_colourbar(title.position = legend.title.position,
                                               ticks = T, 
                                               ticks.linewidth = 1,
                                               ticks.colour = 'black',
                                               frame.colour = 'black',
                                               frame.linewidth = 0.2,
                                               nbin = 1000, # high to set extreme ends
                                               barwidth = .6,
                                               barheight = (panel_height / ifelse(is.null(nrow), 1, nrow))*legend.height.factor))
      
      # accomodate labels
      # p1 <- p1 + theme(plot.margin = margin(t = 0, r = 1, b = 0, l = 0, unit = "cm"),
      #                  legend.box.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      #                  legend.margin=margin(t = 0, r = 1, b = 0, l = 1, unit = "cm"))
      if(!show.legend){
        p1 <- p1 + theme(legend.position = 'none')
        panel_width_in <- panel_width_in+.5
        panel_height_in <- panel_height_in+.5
      }
      else{
        panel_width_in <- panel_width_in+1.25
        panel_height_in <- panel_height_in+.5
      }
    }
    #discrete
    else{
      n <- length(unique(dt$color))
      message("discrete color : ", n)
      
      if(legend.show.na==F & any(is.na(unique(dt$color)))){
        message("Removing NA from legend")
        n <- n-1 
      }
      
      if(missing(pal)){
        message("Palette missing")
        pal <- pals::brewer.accent(n = n)
      }
      
      p1 <- p1 + scale_fill_manual(legend.title, values = pal, na.value = na.color)
      p1 <- p1 + guides(fill = guide_legend(override.aes = list(size=5), ncol=ifelse(n>10,ifelse(n>20,3,2),1)))
      # discrete needs more space than continuous to accomodate labels
      if(!show.legend){
        p1 <- p1 + theme(legend.position = 'none')
        panel_width_in <- panel_width_in+.5
        panel_height_in <- panel_height_in+.5
      }
      else{
        panel_width_in <- panel_width_in+2.5
        panel_height_in <- panel_height_in+.5
      }
    }
  }
  
  if(legend.title!=""){
    p1 <- p1 + theme(legend.title = element_text(angle = legend.title.rotation)) #trbl
  }
  
  p1 <- p1 + force_panelsizes(cols = panel_width, rows = panel_height)
  
  if(small_axis_title){
    xmin <- min(dt[,x])
    xmax <- xmin + abs(xmin)*.2
    ymin <- min(dt[,y])
    ymax <- ymin + abs(ymin)*.2
    message("xmax",xmax)
    message("xmin",xmin)
    
    p1 <- p1 +
      # vertical arrow
      annotate("segment", x=xmin, y=ymin, xend=xmin, yend=ymax, linewidth=.5,
               arrow = arrow(type = "closed", length = unit(0.15, "cm")), color="black") +
      annotate("segment", x=xmin, y=ymin, xend=xmax, yend=ymin, linewidth=.5, 
               arrow = arrow(type = "closed", length = unit(0.15, "cm")), color="black") +
      theme(axis.title.x =  element_text(hjust = 0, vjust = 2.4, size = 12),
            axis.title.y = element_text(hjust = 0, vjust = -1.6, size = 12))
  }
  
  if("axis.label" %in% colnames(dt)){
    if(!is.null(nrow))
      p1 <- p1 + facet_wrap(~ axis.label, nrow = nrow)
    else
      p1 <- p1 + facet_wrap(~ axis.label)
  }
  
  if(!missing(filename.prefix)){
    if(as.pdf){
      f <- file.path(plot.dir,paste0(filename.prefix,".pdf"))
      message("Plotting to:",f)
      pdf(file = f, width = panel_width_in, height = panel_height_in, bg = 'transparent')
      print(p1)
      dev.off()
    } else{
      f <- file.path(plot.dir,paste0(filename.prefix,".png"))
      message("Plotting to:",f)
      png(file = f,
          width = panel_width_in, height = panel_height_in, units = 'in', bg = 'transparent', res = resolution)
      print(p1)
      dev.off()
    }
  }
  p1
}


UmapPlotEnhanceOld <- function(obj,
                               p1,
                               addLabels=F,
                               axis.label=NA, #"orig.ident",
                               outline.size=.2,
                               outline.color='black', #'grey80',
                               pal=NA,
                               panel.border.size=0.2,
                               show.legend=T,
                               label.size=4,
                               title.size=14,
                               scale.color=F,
                               pt.size=0.8,
                               pt.alpha=1,
                               title=NA,
                               legend.title="",
                               legend.title.rotation=90,
                               legend.title.size=10,
                               legend.height.factor=1,
                               legend.position=c(1.3, .85),
                               legend.direction="vertical",
                               legend.show.na=T,
                               break.unit=1,
                               fixed_coord=T,
                               panel_width_in=3,
                               vmin=NA,
                               vmax=NA,
                               na.color="gray80",
                               filename.prefix,
                               nrow=NULL,
                               coord1='UMAP_1',
                               coord2='UMAP_2',
                               small_axis_title=T,
                               resolution=150,
                               as.pdf=F){
  # Plot the integration UMAP
  require(shadowtext)
  require(ggrastr)
  require(ggh4x)
  require(gtable)
  
  dt <- NULL
  
  if(is.data.frame(p1))
    dt <- data.table(p1, keep.rownames = 'rn')
  else
    dt <- data.table(p1$data, keep.rownames = 'rn')
  cols <- colnames(dt)
  
  panel_width <- unit(panel_width_in, "in")
  if(is.na(title)){
    title <- data.table::last(cols)
  }
  
  if(!coord1 %in% colnames(dt)){
    coord1 <- colnames(dt)[2]
    coord2 <- colnames(dt)[3]
  }
  
  setnames(dt, old = c(coord1,coord2), new = c("x","y"))
  setnames(dt, data.table::last(colnames(dt)), new = c("color"))
  
  ## This creates 0 valued points not in faceting window
  if(axis.label %in% colnames(obj@meta.data)){
    dt[,axis.label := obj@meta.data[dt$rn,eval(axis.label)]]
    dt <- lapply(unique(dt$axis.label), function(i){
      ds <- copy(dt[axis.label!=i])
      dt_zero <- copy(ds)
      dt_zero$color <- NA
      dt_zero$axis.label <- i
      rbindlist(list(dt_zero, ds))
    }) %>% rbindlist
  }
  
  dt <- dt[order(color, decreasing = F, na.last = F)]
  
  
  p1 <- NULL
  if(all(is.na(pal))){
    p1 <- ggplot(dt, aes(x=x, y=y)) #, fill=color
    outline.size <- 0
  }
  else{
    if(length(grep("_multiome", dt$color))){
      dt[,color := gsub("_multiome","",color)]
      if(!is.null(names(pal)))
        names(pal) <- gsub("_multiome","",names(pal))
    }
    p1 <- ggplot(dt, aes(x=x, y=y, fill=color))
  }
  
  
  stroke= outline.size
  
  if(outline.size > 0){
    stroke <- outline.size # deintensified stroke
    p1 <- p1 + geom_point_rast(size=pt.size, shape=1, stroke=outline.size, col=outline.color) # standard stroke #'grey40'
  }
  
  p1 <- p1 + geom_point_rast(size=pt.size, shape=21, stroke=0, alpha=pt.alpha) +
    theme_bw()
  
  p1 <- p1 + xlab( gsub("_","",coord1) ) + ylab( gsub("_","",coord2) ) +
    theme(legend.direction = legend.direction,
          legend.position = legend.position,
          text = element_text(size = 15),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(size=title.size, hjust = 0.5),
          legend.text = element_text(size = 12),
          panel.border = element_rect(linewidth = panel.border.size, color='black'),
          strip.background = element_blank(),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          plot.margin = margin(0,0,0,0),
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.title = element_text(angle = legend.title.rotation, size = legend.title.size),
          legend.background = element_blank(), # get rid of legend bg element_rect(fill = 'transparent')
          legend.box.background = element_blank(), # get rid o
          legend.key=element_blank(), 
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(t=0,r=0,b=0,l=0)) #trbl
  if(fixed_coord)
    p1 <- p1 + coord_fixed(clip = 'off')
  
  if(addLabels){
    dt[,label := color]
    median.labels <- dt[,.(x=median(x),y=median(y), color), by=label]
    median.labels <- median.labels[order(label)]
    median.labels <- unique(median.labels)
    p1 <- p1 + geom_shadowtext(data = median.labels,
                               aes(x=x, y=y, label=label, fontface='bold'),
                               color='black',
                               bg.color='white',
                               bg.r=0.1,
                               label.padding=0.01,
                               size=label.size)
  }
  
  #  get dimensions of plot
  gt <- ggplotGrob(p1)
  panel <- gtable_filter(gt, "panel")
  # height is scaled compared to width
  panel_height_in <- panel_width_in*as.numeric(gsub("null","",panel$heights))
  # set the panel height based on width
  panel_height <- unit(panel_height_in, "in") # reduce for title
  
  p1 <- p1 + ggtitle(label = title) + theme(plot.title = element_text(size = title.size))
  
  # Do we want color fill?
  if(all(is.na(pal))){
    ##nothing
  }
  else {
    if(is.numeric(dt$color)){
      message("continuous color")
      if(missing(pal))
        pal <- pals::brewer.purples(10)
      
      lim <- range(dt$color, na.rm = T)
      
      if(!is.na(vmin))
        #lim[1] <- min(dt[color>=vmin]$color, na.rm = T)
        lim[1] <- min(vmin, dt$color, na.rm = T)
      
      if(!is.na(vmax))
        lim[2] <- max(vmax, dt$color, na.rm = T)
      
      if(scale.color)
        lim <- c(-max(abs(dt$color), na.rm = T), max(abs(dt$color), na.rm = T))
      else if(lim[1] > 0)
        lim <- c(0,max(abs(dt$color), na.rm = T))
      
      if(is.null(break.unit)){
        break.unit <- 1
        if( identical( range(c(0,1)), lim ) )
          break.unit <- 0.2
      }
      
      p1 <- p1 + scale_fill_gradientn(legend.title,
                                      guide = "colourbar",
                                      colours = pal,
                                      na.value = na.color,
                                      limits=lim,
                                      breaks=seq(from=floor(lim[1]), to=ceiling(lim[2]), by=break.unit))
      
      legend.title.position <- ifelse(legend.title=="","right","left")
      p1 <- p1 + guides(fill = guide_colourbar(title.position = legend.title.position,
                                               ticks = T,
                                               ticks.linewidth = 1,
                                               ticks.colour = 'black',
                                               frame.colour = 'black',
                                               frame.linewidth = 0.2,
                                               nbin = 1000, # high to set extreme ends
                                               barwidth = .5,
                                               barheight = (panel_height / ifelse(is.null(nrow), 1, nrow))*legend.height.factor))
      
      # accomodate labels
      p1 <- p1 + theme(plot.margin = margin(t = 0, r = 2, b = 0, l = 0, unit = "cm"),
                       legend.box.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
                       legend.margin=margin(t = 0, r = 2, b = 0, l = 1, unit = "cm"))
      if(!show.legend){
        p1 <- p1 + theme(legend.position = 'none')
        panel_width_in <- panel_width_in+.5
        panel_height_in <- panel_height_in+.5
      }
      else{
        panel_width_in <- panel_width_in+1.25
        panel_height_in <- panel_height_in+.5
      }
    }
    #discrete
    else{
      n <- length(unique(dt$color))
      message("discrete color : ", n)
      
      if(legend.show.na==F & any(is.na(unique(dt$color)))){
        message("Removing NA from legend")
        n <- n-1 
      }
      
      if(missing(pal))
        pal <- pals::brewer.accent(n = n)
      
      p1 <- p1 + scale_fill_manual(legend.title, values = pal[1:n], na.value = na.color)
      p1 <- p1 + guides(fill = guide_legend(override.aes = list(size=5), ncol=ifelse(n>10,ifelse(n>20,3,2),1)))
      # discrete needs more space than continuous to accomodate labels
      p1 <- p1 + theme(plot.margin = margin(t = 0, r = 3, b = 0, l = 0, unit = "cm"),
                       legend.box.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
                       legend.margin=margin(t = 0, r = 3, b = 0, l = 1, unit = "cm"))
      if(!show.legend){
        p1 <- p1 + theme(legend.position = 'none')
        panel_width_in <- panel_width_in+.5
        panel_height_in <- panel_height_in+.5
      }
      else{
        panel_width_in <- panel_width_in+1.5
        panel_height_in <- panel_height_in+.5
      }
    }
  }
  
  if(legend.title!=""){
    p1 <- p1 + theme(legend.title = element_text(angle = legend.title.rotation)) #trbl
  }
  
  p1 <- p1 + force_panelsizes(cols = panel_width, rows = panel_height)
  
  if(small_axis_title){
    xmin <- dt[,min(x)]
    ymin <- dt[,min(y)]
    
    p1 <- p1 +
      geom_segment(x=xmin, y=ymin, xend=xmin+2, yend=ymin, linewidth=.5, arrow = arrow(length = unit(.2, "cm"), type = 'closed')) +
      geom_segment(x=xmin, y=ymin, xend=xmin, yend=ymin+2, linewidth=.5, arrow = arrow(length = unit(.2, "cm"), type = 'closed')) +
      theme(axis.title.x =  element_text(hjust = 0, size = 12), 
            axis.title.y = element_text(hjust = 0, size = 12))
    
  }
  
  if("axis.label" %in% colnames(dt)){
    if(!is.null(nrow))
      p1 <- p1 + facet_wrap(~ axis.label, nrow = nrow)
    else
      p1 <- p1 + facet_wrap(~ axis.label)
  }
  
  if(!missing(filename.prefix)){
    if(as.pdf){
      f <- file.path(plot.dir,paste0(filename.prefix,".pdf"))
      message("Plotting to:",f)
      pdf(file = f, width = plot.width, height = plot.height, bg = 'transparent')
      print(p1)
      dev.off()
    } else{
      f <- file.path(plot.dir,paste0(filename.prefix,".png"))
      message("Plotting to:",f)
      png(file = f,
          width = panel_width_in, height = panel_height_in, units = 'in', bg = 'transparent', res = resolution)
      print(p1)
      dev.off()
    }
  }
  p1
}


UmapPlotDensity <- function(obj,
                            feature='Phase',
                            title.size=16,
                            pt.size=0.8,
                            panel_height_in=3,
                            pal=rev(pals::kovesi.linear_kryw_5_100_c64(10)), #pal=pals::viridis(10), #pal=pals::parula(10)
                            filename.prefix,
                            nrow=NULL,
                            as.pdf=F){
  
  
  get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  
  dt <- merge( data.table( obj@reductions$umap@cell.embeddings, keep.rownames = 'cellid' ),
               data.table( obj@meta.data, keep.rownames = 'cellid' ), by='cellid')
  dt[,temp_features := get(feature)]
  
  require(ggh4x)
  p <- lapply(unique(dt[,temp_features]), function(s){
    # select data of interest
    toplot <- dt[temp_features == s]
    if(nrow(toplot) > 1){
      toplot[,density := get_density(x = UMAP_1, y=UMAP_2, n=200)]
      toplot[,density := density/max(density)]
      
      g <- ggplot(toplot, aes(x=UMAP_1, y=UMAP_2, color=density)) +
        geom_point_rast(data = dt, size = pt.size, color = 'gray80') + # include all points and color gray
        geom_point_rast(size = pt.size) +
        scale_color_gradientn('', guide = "colourbar", colours = pal, limits=c(0,1),
                              breaks=seq(from=0, to=1, by=.1)) +
        ggtitle(paste("Density",s)) + xlab("UMAP1") + ylab("UMAP2") +
        theme_bw(12) +
        theme(panel.grid = element_blank(),
              plot.title = element_text(size=title.size, hjust = 0.5),
              plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
              strip.background = element_blank(),
              #panel.border = element_blank(),
              panel.border = element_rect(linewidth = 0.2, color='black'),
              panel.background = element_rect(fill = "transparent"), # bg of the panel
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              legend.background = element_blank(),
              legend.key.size = unit(16, units='points'),
              legend.title = element_text(angle = 90, hjust = 0.5),
              legend.text = element_text(size = 12),
              legend.position = 'right',
              legend.margin=margin(0,0,0,0),
              legend.box.margin=margin(-10,10,-10,-30)) #trbl
      
      #panel_height = unit(1,"npc") - sum(ggplotGrob(g)[["heights"]][-5]) - unit(2,"line")
      panel_height <- unit(panel_height_in ,"in")
      g <- g + guides(color = guide_colourbar(title.position = "left",
                                              ticks = T,
                                              frame.colour = 'black',
                                              frame.linewidth = 0.2,
                                              nbin = 1000, # high to set extreme ends
                                              direction='vertical',
                                              barwidth = .6,
                                              barheight = panel_height))
      g <- g + force_panelsizes(rows = panel_height,
                                cols = panel_height)
    }
    
  })
  
  p <- p[lengths(p) > 0]
  
  
  if(!missing(filename.prefix)){
    width <- ceiling(length(p)/nrow) * (panel_height_in + 1)
    height <- ifelse(nrow==1, (panel_height_in + 1), nrow * (panel_height_in + 1))
    if(as.pdf){
      f <- file.path(plot.dir,paste0(filename.prefix,".pdf"))
      message("Plotting to:",f)
      pdf(file = f, width = width, height = height, bg = 'transparent')
      do.call(grid.arrange, c(p, nrow=nrow))
      dev.off()
    } else{
      f <- file.path(plot.dir,paste0(filename.prefix,".png"))
      message("Plotting [",width,"x",height,"] to:",f)
      png(filename = f, width = width, height = height, bg = 'transparent', units = 'in', res=200)
      do.call(grid.arrange, c(p, nrow=nrow))
      dev.off()
    }
  }
  
  p
}

ViolinPlotEnhance <- function (object, features, type = "violin", idents = NULL, ncol = NULL,
                               sort = FALSE, assay = NULL, y.max = NULL, same.y.lims = FALSE,
                               adjust = 1, cols = NULL, pt.size = 0, group.by = NULL, split.by = NULL,
                               log = FALSE, slot = "data", stack = FALSE, combine = TRUE,
                               pal = gist_heat_r, 
                               y_axis_label="Expression Level (log normalized)",
                               title = NULL, 
                               x.axis.text.size=16, y.axis.text.size=12, strip.text.size=18,
                               x_axis_text_angle=0, y_text_hjust=0.5,
                               fill.by = NULL, flip = FALSE, raster = NULL)
{
  require(Seurat)
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  if (isTRUE(x = stack)) {
    if (!is.null(x = ncol)) {
      warning("'ncol' is ignored with 'stack' is TRUE",
              call. = FALSE, immediate. = TRUE)
    }
    if (!is.null(x = y.max)) {
      warning("'y.max' is ignored when 'stack' is TRUE",
              call. = FALSE, immediate. = TRUE)
    }
  }
  else {
    ncol <- ncol %||% ifelse(test = length(x = features) >
                               9, yes = 4, no = min(length(x = features), 3))
  }
  
  data <- NULL
  # features are genes or other slot
  if(!all(features %in% colnames(object@meta.data)) ){
    data <- FetchData(object = object, vars = features, slot = slot)
  }
  
  # check if features are in meta.data
  mdata <- NULL
  if( any(features %in% colnames(object@meta.data)) ){
    mfeatures <- features[which(features %in% colnames(object@meta.data))]
    mdata <- subset(object@meta.data, select = mfeatures)
    if(!is.null(data))
      data <- cbind(data, mdata)
    else data <- mdata
  }
  
  pt.size <- pt.size %||% AutoPointSize(data = object)
  features <- colnames(x = data)
  if (is.null(x = idents)) {
    cells <- colnames(x = object)
  }
  else {
    cells <- names(x = Idents(object = object)[Idents(object = object) %in%
                                                 idents])
  }
  data <- data[cells, , drop = FALSE]
  idents <- if (is.null(x = group.by)) {
    Idents(object = object)[cells]
  }
  else {
    object[[group.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = idents)) {
    idents <- factor(x = idents)
  }
  if (is.null(x = split.by)) {
    split <- NULL
  }
  else {
    split <- object[[split.by, drop = TRUE]][cells]
    if (!is.factor(x = split)) {
      split <- factor(x = split)
    }
    if (is.null(x = cols)) {
      cols <- hue_pal()(length(x = levels(x = idents)))
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    }
    else if (length(x = cols) == 1 && cols == "interaction") {
      split <- interaction(idents, split)
      cols <- hue_pal()(length(x = levels(x = idents)))
    }
    else {
      cols <- Col2Hex(cols)
    }
    if (length(x = cols) < length(x = levels(x = split))) {
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    }
    cols <- rep_len(x = cols, length.out = length(x = levels(x = split)))
    names(x = cols) <- levels(x = split)
    if ((length(x = cols) > 2) & (type == "splitViolin")) {
      warning("Split violin is only supported for <3 groups, using multi-violin.")
      type <- "violin"
    }
  }
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data)
  }
  if (isTRUE(x = stack)) {
    return(MultiViolinPlotEnhance(type = type, data = data, idents = idents,
                                  split = split, sort = sort, same.y.lims = same.y.lims,
                                  adjust = adjust, cols = cols, pt.size = pt.size,
                                  pal = pal, title = title,
                                  y_axis_label = y_axis_label,
                                  x.axis.text.size = x.axis.text.size,
                                  y.axis.text.size=y.axis.text.size,
                                  x_axis_text_angle=x_axis_text_angle,
                                  y_text_hjust=y_text_hjust,
                                  strip.text.size = strip.text.size,
                                  log = log, fill.by = fill.by, flip = flip))
  }
  plots <- lapply(X = features, FUN = function(x) {
    return(SingleExIPlot(type = type, data = data[, x, drop = FALSE],
                         idents = idents, split = split, sort = sort, y.max = y.max,
                         adjust = adjust, cols = cols, pt.size = pt.size,
                         log = log, raster = raster))
  })
  label.fxn <- switch(EXPR = type, violin = if (stack) {
    xlab
  } else {
    ylab
  }, splitViolin = if (stack) {
    xlab
  } else {
    ylab
  }, ridge = xlab, stop("Unknown ExIPlot type ", type, call. = FALSE))
  for (i in 1:length(x = plots)) {
    key <- paste0(unlist(x = strsplit(x = features[i], split = "_"))[1],
                  "_")
    obj <- names(x = which(x = Key(object = object) == key))
    if (length(x = obj) == 1) {
      if (inherits(x = object[[obj]], what = "DimReduc")) {
        plots[[i]] <- plots[[i]] + label.fxn(label = "Embeddings Value")
      }
      else if (inherits(x = object[[obj]], what = "Assay")) {
        next
      }
      else {
        warning("Unknown object type ", class(x = object),
                immediate. = TRUE, call. = FALSE)
        plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
      }
    }
    else if (!features[i] %in% rownames(x = object)) {
      plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
    }
  }
  if (combine) {
    require(patchwork)
    plots <- wrap_plots(plots, ncol = ncol)
    if (length(x = features) > 1) {
      plots <- plots & NoLegend()
    }
  }
  return(plots)
}

MultiViolinPlotEnhance <- function (data, idents, split = NULL, type = "violin", sort = FALSE,
                                    same.y.lims = same.y.lims, adjust = 1, pt.size = 0, cols = NULL,
                                    pal = gist_heat_r, title = NULL, y_axis_label=NULL,
                                    x.axis.text.size=9, y.axis.text.size=10,
                                    strip.text.size=24, x_axis_text_angle=45,
                                    y_text_hjust=0.5,
                                    seed.use = 42, log = FALSE, fill.by = NULL, flip = NULL)
{
  require(dplyr)
  require(cowplot)
  #theme_set(theme_cowplot())
  if (!(fill.by %in% c("feature", "ident"))) {
    stop("`fill.by` must be either `feature` or `ident`")
  }
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (!is.data.frame(x = data) || ncol(x = data) < 2) {
    stop("MultiExIPlot requires a data frame with >1 column")
  }
  #data <- Melt(x = data)
  data <- Seurat:::Melt(x = data)
  data <- data.frame(feature = data$cols, expression = data$vals,
                     ident = rep_len(x = idents, length.out = nrow(x = data)))
  if ((is.character(x = sort) && nchar(x = sort) > 0) || sort) {
    data$feature <- as.vector(x = data$feature)
    data$ident <- as.vector(x = data$ident)
    avgs.matrix <- sapply(X = split(x = data, f = data$ident),
                          FUN = function(df) {
                            return(tapply(X = df$expression, INDEX = df$feature,
                                          FUN = mean))
                          })
    idents.order <- hclust(d = dist(x = t(x = L2Norm(mat = avgs.matrix,
                                                     MARGIN = 2))))$order
    avgs.matrix <- avgs.matrix[, idents.order]
    avgs.matrix <- L2Norm(mat = avgs.matrix, MARGIN = 1)
    position <- apply(X = avgs.matrix, MARGIN = 1, FUN = which.max)
    mat <- hclust(d = dist(x = avgs.matrix))$merge
    orderings <- list()
    for (i in 1:nrow(mat)) {
      x <- if (mat[i, 1] < 0)
        -mat[i, 1]
      else orderings[[mat[i, 1]]]
      y <- if (mat[i, 2] < 0)
        -mat[i, 2]
      else orderings[[mat[i, 2]]]
      x.pos <- min(x = position[x])
      y.pos <- min(x = position[y])
      orderings[[i]] <- if (x.pos < y.pos) {
        c(x, y)
      }
      else {
        c(y, x)
      }
    }
    features.order <- orderings[[length(x = orderings)]]
    data$feature <- factor(x = data$feature, levels = unique(x = sort(x = data$feature))[features.order])
    data$ident <- factor(x = data$ident, levels = unique(x = sort(x = data$ident))[rev(x = idents.order)])
  }
  else {
    data$feature <- factor(x = data$feature, levels = unique(x = data$feature))
  }
  if (log) {
    noise <- rnorm(n = nrow(x = data))/200
    data$expression <- data$expression + 1
  }
  else {
    noise <- rnorm(n = nrow(x = data))/1e+05
  }
  for (f in unique(x = data$feature)) {
    if (all(data$expression[(data$feature == f)] == data$expression[(data$feature ==f)][1])) {
      warning("All cells have the same value of ", f, call. = FALSE,
              immediate. = TRUE)
    }
    else {
      data$expression[(data$feature == f)] <- data$expression[(data$feature ==f)] + noise[(data$feature == f)]
    }
  }
  if (type == "violin" && !is.null(x = split)) {
    data$split <- rep_len(x = split, length.out = nrow(data))
    vln.geom <- geom_violin
    fill.by <- "split"
  }
  else if (type == "splitViolin" && !is.null(x = split)) {
    data$split <- rep_len(x = split, length.out = nrow(data))
    vln.geom <- geom_split_violin
    fill.by <- "split"
    type <- "violin"
  }
  else {
    vln.geom <- geom_violin
  }
  switch(EXPR = type,
         violin = {
           geom <- list(vln.geom(scale = "width", adjust = adjust, trim = TRUE, linewidth = .1)) # size is border thickness of violin
         }, ridge = {
           geom <- list(geom_density_ridges(scale = 4), theme_ridges(),
                        scale_y_discrete(expand = c(0.01, 0)))
         }, stop("Unknown plot type: ", type))
  if (flip) {
    x <- "ident"
    x.label <- "Identity"
    y <- "expression"
    y.label <- y_axis_label #"Expression Level (log normalized)"
  }
  else {
    y <- "ident"
    y.label <- "Identity"
    x <- "expression"
    x.label <- y_axis_label #"Expression Level (log normalized)"
  }
  # plot <- ggplot(data = data, mapping = aes_string(x = x, y = y, fill = fill.by)[c(2, 3, 1)]) +
  #   labs(x = x.label, y = y.label, fill = NULL) + theme_cowplot()
  # plot <- do.call(what = "+", args = list(plot, geom))
  
  ### New here start
  setDT(data)
  data[,mean_expression:=mean(expression), by=list(feature, ident)]
  #View(data)
  plot <- ggplot(data = data, mapping = aes(x = data[[x]], y = as.numeric(data[[y]]), fill = mean_expression)[c(2, 3, 1)]) +
    labs(x = x.label, y = y.label, fill = NULL) + theme_cowplot()
  
  plot <- do.call(what = "+", args = list(plot, geom))
  ### New here end
  
  
  if (flip) {
    hjust <- 0.5
    if(x_axis_text_angle != 0)
      hjust <- 1
    plot <- plot + scale_y_continuous(expand = c(0, 0), labels = function(x) c(rep(x = "",times = length(x) - 2), x[length(x) - 1], "")) +
      facet_grid(feature ~ ., scales = (if (same.y.lims)
        "fixed"
        else "free")) + Seurat:::FacetTheme(panel.spacing = unit(0,"lines"),
                                            panel.background = element_rect(fill = NA,color = "black"),
                                            axis.text.y = element_text(size = y.axis.text.size, angle=0),
                                            axis.text.x = element_text(angle = x_axis_text_angle, vjust=0.5, hjust = hjust, size=x.axis.text.size),
                                            plot.title = element_text(hjust=0.5),
                                            strip.text.y.right = element_text(angle = 0, size=strip.text.size, hjust = y_text_hjust))
  }
  else {
    plot <- plot + scale_x_continuous(expand = c(0, 0), labels = function(x) c(rep(x = "",times = length(x) - 2), x[length(x) - 1], "")) +
      facet_grid(. ~ feature, scales = (if (same.y.lims)
        "fixed"
        else "free")) + Seurat:::FacetTheme(panel.spacing = unit(0,"lines"),
                                            panel.background = element_rect(fill = NA, color = "black"),
                                            axis.text.y = element_text(size = y.axis.text.size),
                                            axis.text.x = element_text(size = x.axis.text.size),
                                            plot.title = element_text(hjust=0.5),
                                            strip.text.x = element_text(angle = -90, size=strip.text.size))
  }
  if (log) {
    plot <- plot + scale_x_log10()
  }
  if (!is.null(x = cols)) {
    if (!is.null(x = split)) {
      idents <- unique(x = as.vector(x = data$ident))
      splits <- unique(x = as.vector(x = data$split))
      labels <- if (length(x = splits) == 2) {
        splits
      }
      else {
        unlist(x = lapply(X = idents, FUN = function(pattern,
                                                     x) {
          x.mod <- gsub(pattern = paste0(pattern, "."),
                        replacement = paste0(pattern, ": "), x = x,
                        fixed = TRUE)
          x.keep <- grep(pattern = ": ", x = x.mod, fixed = TRUE)
          x.return <- x.mod[x.keep]
          names(x = x.return) <- x[x.keep]
          return(x.return)
        }, x = unique(x = as.vector(x = data$split))))
      }
      if (is.null(x = names(x = labels))) {
        names(x = labels) <- labels
      }
    }
    else {
      #labels <- levels(x = droplevels(data$ident))
    }
    plot <- plot + scale_fill_manual(values = cols, labels = labels)
  } else {
    
    sum_func <- list(Minimum = min, Mean = mean, Maximum = max)
    
    ls_val <- unlist(lapply(sum_func, function(f) f(data$mean_expression)))
    
    # colors the violins across levels
    if(ls_val[1] > 0)
      ls_val[1] <- 0
    
    plot <- plot + scale_fill_gradientn(
      colours =  pal,
      limits=c(ls_val[1], ls_val[3]),
      breaks = c(ls_val[1],ls_val[3])
    )
    
  }
  
  plot <-  plot + ggtitle(title)
  
  #print(plot)
  
  #return(data)
  
  return(plot)
}

ViolinPlotQC <- function(seuratObj, features=c('nFeature_RNA','nCount_RNA','percent_mito','percent_ribo'), title="", ylim=NULL){
  gp <- lapply(features, function(i){
    ggplot(seuratObj@meta.data, aes(x=bioNames, y=get(i), fill=bioNames)) +
      geom_violin(trim = F, scale = 'width') +
      scale_fill_manual("", values = sample_palette) +
      geom_boxplot(fill='white', alpha=0.7, width=0.5, outlier.size = 1) +
      ggtitle(title) + xlab("") + ylab(i) + scale_y_continuous(limits = ylim) +
      theme_classic(14) + theme(axis.text.x = element_text(angle=90), legend.position = 'none')
  })
  #do.call(grid.arrange, c(gp, nrow=2))
}


PlotBarVariableInClusters <- function(seuratObj, 
                                      variable='orig.ident', 
                                      by.group='seurat_clusters', 
                                      legend.title="", 
                                      legend.column=1,
                                      palette=sample_palette, 
                                      addcounts=T, 
                                      orientation='horizontal', 
                                      width_multiplier=1,
                                      height_multiplier=1,
                                      text.angle=0,
                                      x.label,
                                      y.label,
                                      filename=""){
  dt <- as.data.table(seuratObj@meta.data[,c(eval(by.group),eval(variable))])
  setnames(dt, by.group, 'grouping')
  
  dt[,nCluster := .N, by=grouping]
  dt[,nClustervar := .N, by=list(grouping, get(variable))]
  dt <- unique(dt)
  dt <- dt[order(grouping, decreasing = F)]
  dt[,percClustervar := 100 * (nClustervar / nCluster) ]
  dtlab <- unique(dt[,list(grouping, nCluster)])
  
  fill <- scale_fill_manual(legend.title, values = palette)
  
  g <- NULL
  width <- 9
  height <- 4
  if(missing(y.label))
    y.label <- by.group
  if(missing(x.label))
    x.label <- "Percentage of cells"
  if(orientation=='horizontal'){
    g <- ggplot(dt, aes(x=grouping, y=percClustervar)) +
      geom_bar(stat='identity', aes(fill=get(variable)), linewidth = .95) +
      fill +
      xlab(x.label) +
      ylab(y.label) + 
      theme_bw(14) +
      coord_cartesian(ylim = c(0,100), expand = F, clip = 'off') +
      theme(panel.ontop = T,
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.grid.major.y = element_line(colour = 'black', size=.2, linetype = 'longdash'),
            legend.direction = "vertical",
            legend.position = "right",
            axis.text.y = element_text(size = 14),
            axis.text.x = element_text(size = 14, angle = text.angle, hjust=0.5, vjust=0.5),
            axis.title = element_text(size = 14),
            text = element_text(size = 12),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 14),
            legend.background = element_blank(),
            legend.margin = margin(0,0,0,1, unit = "cm"), #trbl
            plot.margin = unit(x = c(2,1,1,1), units = 'cm'),
            plot.background = element_blank())
    g <- g + guides(fill = guide_legend(title = legend.title, title.position = "top", ncol = legend.column))
    
    if(addcounts)
      g <- g + geom_text(data = dtlab, aes(label=paste0("n=",nCluster)), y=100, color='black',
                         vjust=0, hjust=0, angle=90, size=3.5)
  } else {
    width <- 4
    height <- 4
    if(missing(y.label))
      y.label <- "Percentage of cells"
    if(missing(x.label))
      x.label <- by.group
    dt <- dt[order(grouping, decreasing = T)]
    dt[,grouping:=factor(grouping, unique(grouping))]
    dtlab <- unique(dt[,list(grouping, nCluster)])
    g <- ggplot(dt, aes(y=grouping, x=percClustervar)) +
      geom_bar(stat='identity', aes(fill=get(variable)), linewidth = .95) +
      fill +
      xlab(x.label) + 
      ylab(y.label) +
      theme_bw(14) +
      coord_cartesian(xlim = c(0,100), expand = F, clip = 'off') +
      theme(panel.ontop = T,
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.grid.major.x = element_line(colour = 'black', linewidth=.2, linetype = 'longdash'),
            legend.direction = "horizontal",
            legend.position = "right",
            axis.text.x = element_text(size = 14),
            axis.title = element_text(size = 14),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 14),
            legend.margin = margin(0,0,0,1, unit = "cm"), #trbl
            legend.background = element_blank(),
            plot.margin = unit(x = c(1,2,1,1), units = 'cm'),
            plot.background = element_blank())
    g <- g + guides(fill = guide_legend(title = legend.title, title.position = "top", ncol = legend.column))
    
    if(addcounts)
      g <- g + geom_text(data = dtlab, aes(label=paste0("n=",nCluster)), x=101, color='black',
                         vjust=0.5, hjust=0, angle=0, size=3.5)
  }
  
  if(filename !=""){
    fout <- file.path(plot.dir,paste0(filename,'.pdf'))
    message("Plotting to:",fout)
    pdf(file = fout, width = width * width_multiplier, height = height * height_multiplier, bg = 'transparent')
    
    print(g)
    dev.off()
  }
  g
  
}

grid_arrange_list_shared_legend <- function(plots,
                                            ncol = length(plots),
                                            nrow = 1,
                                            legend.position = c("right","bottom","top")) {
  
  # plots <- list(...)
  position <- match.arg(legend.position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = legend.position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  
  if(ncol < length(plots))
    nrow <- ceiling(length(plots) / ncol)
  else if (nrow > 1)
    ncol <- ceiling(length(plots) / nrow)
  
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  library(grid)
  combined <- switch(position,
                     "top" = arrangeGrob(legend,
                                         do.call(arrangeGrob, gl),
                                         ncol = 1,
                                         heights = unit.c(lheight, unit(1, "npc") - lheight)),
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  # grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

