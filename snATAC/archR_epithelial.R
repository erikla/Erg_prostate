#devtools::install_github("GreenleafLab/ArchR", ref=" <- aster", repos = BiocManager::repositories())
#library(ArchR)
#ArchR::installExtraPackages()

library(ArchR)
library(ComplexHeatmap)
library(gridExtra)
library(ggrepel)
library(shadowtext)
library(magick)
library(ggplot2)
library(ggrastr)
# library(doMC)
# library(foreach)

library("BSgenome.Mmusculus.UCSC.mm10")

#registerDoMC(8)
addArchRThreads(threads = 8)
addArchRGenome("mm10")
# Used by geneRanges method (defined here) to extract gene ranges and symbols
#BiocManager::install("Mus.musculus")

sequencing.dir <- '/data/storage/sawyers/erg/multiome/processed/'
base.dir <- '/home/erik/projects/sawyers/erg/analysis/multiome_mouse_WF-2149_epithelial_no_sv/'
archr.dir <- file.path(base.dir,'archR')
data.dir <- file.path(archr.dir,'data')
plot.dir <- file.path(archr.dir,'Plots')
robj.dir <- file.path(archr.dir,'robj')
annotation.dir <- file.path(archr.dir,'Annotations')
spreadsheets.dir <- file.path(archr.dir,'spreadsheets')

## All murine specific
lum_mature_markers <- c("Cd24a","Foxa1","Gata3","Krt8","Krt18","Krt19","Nkx3-1")
basal_markers <- c("Itga6","Krt5","Krt14","Krt15","Krt17","Trp63")
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

chromatin_modifiers_list <- list(
  'H3K4 methyl' = c('Kmt2a','Kmt2b','Kmt2c','Kmt2d', 'Setd1a','Setd1b'),
  'H3K methyl' = c('Dot1l','Nsd1','Nsd2','Nsd3','Setd2'),
  #'Histone acetyl'
  'p300/CBP family' = c('Ep300', 'Crebbp'),
  'GCN5 family' = c('Nr6a1', 'Kat2b'),
  'MYST family' = c('Kat5', 'Kat6a', 'Kat6b', 'Kat7', 'Kat8')
)

## thresholded peaks
marker.cutOff = "FDR <= 0.05 & Log2FC >= 0.5"
abs.marker.cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.5"
gene.marker.cutOff = "FDR <= 0.05 & Log2FC >= 1"

# Sample Name	    Genotype	Time point	Lobe	Reporter allele	Note
# 5490-EPC-6m	    ERG+Pten-	6m	        DLP	  ERG-ires-GFP	  1 mouse
# 5503-PYC-6m	    Pten-	    6m	        DLP	  EYFP	          1 mouse
# 6045-EPC-3m	    ERG+Pten-	3m	        DLP	  ERG-ires-GFP	  1 mouse
# 6048-PhY-3m	    wt	      3m	        DLP	  n/a	            1 mouse
# 6051-PYC-3m	    Pten-	    3m	        DLP	  EYFP	          1 mouse
# 6390-c2-EPC-4w	ERG+Pten-	4w	        DLP	  ERG-ires-GFP	  pooled from two mice
# 6398-c2-PYC-4w	Pten-	    4w	        DLP	  EYFP	          pooled from two mice

INPUT_LABELS = list(
  "WF-2149_5490-EPC-6m_multiome" = "Erg+Pten- 6m",
  "WF-2149_5503-PYC-6m_multiome" = "Pten- 6m",
  "WF-2149_6045-EPC-3m_multiome" = "Erg+Pten- 3m",
  "WF-2149_6048-PhY-3m_multiome" = "WildType 3m",
  "WF-2149_6051-PYC-3m_multiome" = "Pten- 3m",
  "WF-2149_6390-c2-EPC-4w_multiome" = "Erg+Pten- 4w",
  "WF-2149_6398-c2-PYC-4w_multiome" = "Pten- 4w"
)

CONDITION_LABELS <- list(
  "WF-2149_5490-EPC-6m_multiome" = "Erg+Pten-",
  "WF-2149_5503-PYC-6m_multiome" = "Pten-",
  "WF-2149_6045-EPC-3m_multiome" = "Erg+Pten-",
  "WF-2149_6048-PhY-3m_multiome" = "WildType",
  "WF-2149_6051-PYC-3m_multiome" = "Pten-",
  "WF-2149_6390-c2-EPC-4w_multiome" = "Erg+Pten-",
  "WF-2149_6398-c2-PYC-4w_multiome" = "Pten-"
)

TIME_LABELS <- list(
  "WF-2149_5490-EPC-6m_multiome" = "6m",
  "WF-2149_5503-PYC-6m_multiome" = "3m",
  "WF-2149_6045-EPC-3m_multiome" = "3m",
  "WF-2149_6048-PhY-3m_multiome" = "3m",
  "WF-2149_6051-PYC-3m_multiome" = "3m",
  "WF-2149_6390-c2-EPC-4w_multiome" = "4w",
  "WF-2149_6398-c2-PYC-4w_multiome" = "4w"
)

samples_order <-c("WF-2149_5490-EPC-6m_multiome",
                  "WF-2149_5503-PYC-6m_multiome",
                  "WF-2149_6045-EPC-3m_multiome",
                  "WF-2149_6048-PhY-3m_multiome",
                  "WF-2149_6051-PYC-3m_multiome",
                  "WF-2149_6390-c2-EPC-4w_multiome",
                  "WF-2149_6398-c2-PYC-4w_multiome")

samples_palette_list = c('#E60012','#B27BB4','#F3AB19','#808000','#34BDEF','#0000ff','#63FFAC')
names(samples_palette_list) <- samples_order

bionames_order <- c("Erg+Pten- 6m",
                    "Erg+Pten- 3m",
                    "Erg+Pten- 4w",
                    "Pten- 6m",
                    "Pten- 3m",
                    "Pten- 4w",
                    "WildType 3m")
bionames_palette_list = c('#E60012','#F3AB19','#0000ff','#B27BB4','#34BDEF','#63FFAC','#808000') #'#4251A2',
names(bionames_palette_list) =  bionames_order

#cluster_palette <- c("#F0A0FF", "#0000ff", "#FFB500", "#4C005C", "#191919", "#94FFB5", "#FFFF00", "#8CD0FF", "#E60012", "#C0B9B2", "#005C31", "#FF8A9A")
cluster_palette <- c("#F0A0FF", "#0000ff", "#FFB500", "#4C005C", "#191919", "#94FFB5", "#FFFF00", "#8CD0FF", "#E60012", "#C0B9B2", "#005C31", "#FF8A9A")
condition_palette <- c('#D51F26','#272E6A','#CBA00D')
names(condition_palette) <- c("Erg+Pten-","Pten-","WildType")

epitype_palette <- c("dodgerblue","orangered1","chartreuse3","purple")
names(epitype_palette) <- c('Epi_Basal_1', 'Epi_Luminal_1', 'Epi_Luminal_2Psca', 'Epi_Luminal_3Foxi1')

condition_palette <- c('#CBA00D','#272E6A','#D51F26')
names(condition_palette) <- c("WildType","Pten-","Erg+Pten-")

cluster3_palette_list <- list("Basal" = "#1C86EE",
                              "EPC (IM)" = "#D51F26",
                              "L1" = "#FFC125",
                              "Tumor-Bas (IM)"= "#228B22",
                              "Tumor-L2" = "#FF90C9",
                              "L2" = "#4C005C")

cluster3_palette <- unlist(cluster3_palette_list)
names(cluster3_palette) <- names(cluster3_palette_list)

cluster_palette_list = cluster_palette
names(cluster_palette_list) =  paste0("C",seq_along(cluster_palette))

solarExtraWhite_palette <- c("#3361A5","#248AF3","#14B3FF","#88CEEF","white","#EAD397","#FDB31A","#E42A2A","#A31D1D")

rna_cluster_palette <- list(
  '1'= "#F0A0FF",
  '2'= "#0075DC",
  '3'= "#993F00",
  '4'= "#4C005C",
  '5'= "#191919",
  '6'= "#005C31",
  '7'= "#2BCE48",
  '8'= "#FFCC99",
  '9'= "#808080",
  '10'= "#94FFB5",
  '11'="#8F7C00",
  '12'="#9DCC00",
  '13'="#C20088",
  '14'="#003380",
  '15'="#FFA405",
  '16'="#FFA8BB",
  '17'="#426600",
  '18'="#FF0010",
  '19'="#5EF1F2",
  '20'="#00998F",
  '21'="#E0FF66",
  '22'="#740AFF",
  '23'="#990000",
  '24'="#FFFF80",
  '25'="#FFE100",
  '26'="#FF5005",
  '-1' = '#D9D9D9'
)
rna_cluster_palette <- unlist(rna_cluster_palette)


main <- function(){
  
  #BuildProj()
  
  #### Load object ####
  proj <- loadArchRProject(path = archr.dir)
  
  getAvailableMatrices(proj)
  
  # Look at data
  proj@cellColData
  
  table(proj@cellColData$Clusters)
  
  #### Plot UMAPS ####
  table(proj@cellColData$bioNames)
  #BioNames
  pclust <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", pal = cluster_palette_list, 
                          name = "Clusters", embedding = "UMAP")
  
  pBio <- plotEmbedding(ArchRProj = proj, 
                        colorBy = "cellColData", 
                        name = "bioNames", 
                        embedding = "UMAP", 
                        pal=bionames_palette_list,
                        baseSize = 14)
  pBio+pclust$layers
  
  pClust3 <- plotEmbedding(ArchRProj = proj,
                           colorBy = "cellColData",
                           pal = cluster3_palette,
                           name = "Cluster3",
                           embedding = "UMAP")
  pClust3 + pclust$layers
  
  #### Umaps by condition ####
  pCond <- plotEmbedding(ArchRProj = proj,
                         colorBy = "cellColData", 
                         name = "Condition",
                         embedding = "UMAP",
                         pal = condition_palette,
                         size=0.2, baseSize = 14)
  
  #### UMAP cluster3 ####
  pCelltype <- plotEmbedding(ArchRProj = proj,
                         colorBy = "cellColData",
                         name = "Cluster3",
                         embedding = "UMAP",
                         size=0.2, baseSize = 14)
  
  
  #### Fig. 6A: Facet condition ####
  pCond <- plotEmbedding(ArchRProj = proj,
                         colorBy = "cellColData",
                         name = "Cluster3",
                         embedding = "UMAP",
                         size=0.2, baseSize = 14)
  df <- pCond$data
  df$color <- gsub("[0-9]-","",df$color)
  df$color <- factor(x = df$color, levels = c("Tumor-Bas (IM)","Tumor-L2","L2","Basal","L1","EPC (IM)"), ordered = T)
  pCond$data <- df
  gp <- UmapPlotEnhance(ArchRProj = proj, p1 = pCond,
                        addLabels = F,
                        label.size = 4,
                        na.color = "gray50",
                        axis.label = 'Condition',
                        pt.size = 0.8,
                        pal=cluster3_palette_list,
                        title.size = 18,
                        #outline.color = 'gray80',
                        panel.border.size = 0)
  gp
  pCond

  
  #Cluster #
  plotEmbedding(ArchRProj = proj,
                colorBy = "cellColData",
                name = "Clusters",
                pal = cluster_palette_list,
                embedding = "UMAP")
  
  #------------------------------------------- Pseudobulk -------------------------------------------
  # Pseudo-Bulk for calling peaks per cluster #
  proj <- addGroupCoverages(ArchRProj = proj,
                            minReplicates = 2,
                            groupBy = "Clusters", 
                            force = T)
  
  # norm by ReadsInTSS is default; tilesk$reads <- tilesk$reads * 10^4/sum(normBy[cellGroupi, 1])
  
  #------------------------------------------- BigWigs -------------------------------------------  
  # Make BigWigs this method overwrites existing files.
  
  getGroupBW(ArchRProj = proj, groupBy = "Cluster3", ceiling = 10)
  
  getGroupBW(ArchRProj = proj, groupBy = "Cluster3", tileSize = 100, maxCells = 20000, normMethod = "nFrags")
  
  getGroupBW(ArchRProj = proj, groupBy = "Clusters", tileSize = 20, maxCells = 10000, normMethod = c("ReadsInTSS", "ReadsInPromoter"))
  
  getGroupBW(ArchRProj = proj, groupBy = "Clusters", tileSize = 100, maxCells = 20000) #, normMethod = 'none')
  
  getGroupBW(ArchRProj = proj, groupBy = "bioNames", tileSize = 100, maxCells = 20000) #, normMethod = 'none')
  
  #------------------------------------------- Call Peaks -------------------------------------------
  DT <- data.table(as.data.frame(proj@cellColData), keep.rownames = 'cell')
  table(DT[,list(bioNames, Clusters)])
  DT[,.N, by=Clusters][order(N)]
  
  # this is their implementation of Iterative Overlap
  # https://www.archrproject.com/bookdown/the-iterative-overlap-peak-merging-procedure.html
  pathToMacs2 <- findMacs2()
  proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    groupBy = "Clusters",
    peakMethod = 'Macs2',
    reproducibility = "2", # Requires only 1 sample to have peak call at this locus
    peaksPerCell = 500,
    maxPeaks = 150000,
    minCells = 25, # also used 250 here and same result
    shift = -75,
    extsize = 150,
    threads = 4, 
    force = T,
    pathToMacs2 = pathToMacs2
  )
  
  #------------------------------------------- Peak metrics -------------------------------------------
  ## Peaks were called using MACS2 above in each cluster
  proj <- addPeakMatrix(proj)
  PlotPeakAnnotations(ArchRProj = proj)
  
  #---- FRIPS ----#
  dt <- as.data.table(proj@cellColData[,c("bioNames","Clusters","FRIP","nFrags")])
  dt[,mean(FRIP), by=bioNames]
  cairo_pdf(filename = file.path(plot.dir, 'FRIPS_samples.pdf'), width = 6, height = 7, bg = 'transparent')
  ggplot(dt, aes(x=bioNames, y=FRIP, fill=bioNames)) + 
    geom_violin(draw_quantiles = T, trim = F, scale = 'width') + 
    geom_boxplot(fill='white', alpha=.6) +
    scale_fill_manual("", values = bionames_palette_list) +
    theme_ArchR(baseSize = 14) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1), plot.background = element_blank(), legend.position = 'none')
  dev.off()
  
  cairo_pdf(filename = file.path(plot.dir, 'FRIPS_cluster.pdf'), width = 6, height = 7, bg = 'transparent')
  ggplot(dt, aes(x=Clusters, y=FRIP, fill=Clusters)) + 
    geom_violin(draw_quantiles = T, trim = F, scale = 'width') + 
    geom_boxplot(fill='white', alpha=.6) +
    scale_fill_manual("", values = cluster_palette) +
    theme_ArchR(baseSize = 14) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1), plot.background = element_blank(), legend.position = 'none')
  dev.off()
  
  #---- ReadsInTSS ----#
  dt <- as.data.table(proj@cellColData[,c("bioNames","Clusters","ReadsInTSS")])
  dt[,mean(ReadsInTSS), by=bioNames]
  cairo_pdf(filename = file.path(plot.dir, 'ReadsInTSS_samples.pdf'), width = 6, height = 7, bg = 'transparent')
  ggplot(dt, aes(x=bioNames, y=ReadsInTSS, fill=bioNames)) + 
    geom_violin(draw_quantiles = T, trim = F, scale = 'width') + 
    geom_boxplot(fill='white', alpha=.6) +
    scale_fill_manual("", values = bionames_palette_list) +
    theme_ArchR(baseSize = 14) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1), plot.background = element_blank(), legend.position = 'none')
  dev.off()
  
  cairo_pdf(filename = file.path(plot.dir, 'ReadsInTSS_cluster.pdf'), width = 6, height = 7, bg = 'transparent')
  ggplot(dt, aes(x=Clusters, y=ReadsInTSS, fill=Clusters)) + 
    geom_violin(draw_quantiles = T, trim = F, scale = 'width') + 
    geom_boxplot(fill='white', alpha=.6) +
    scale_fill_manual("", values = cluster_palette) +
    theme_ArchR(baseSize = 14) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1), plot.background = element_blank(), legend.position = 'none')
  dev.off()
  
  saveArchRProject(ArchRProj = proj, load = TRUE)
  
  getAvailableMatrices(proj)
  
  #### Peak set (atlas) ####
  # This peak set contains an annotation for the group from which each peak originated. 
  # However, 
  # these annotations do not inherently mean that the given peak was only called in that group, 
  # rather that the annotated group had the highest normalized significance for that peak call.
  peakset <- getPeakSet(proj)
  
  table(peakset$peakType)
  
  which(duplicated(peakset))
  #### Write bed file of all peaks ####
  library(rtracklayer)
  export.bed(peakset, con = file.path(annotation.dir,'PeakAtlas.bed'))
  
  table(names(peakset))
  mclapply(unique(names(peakset)), function(n){
    idx <- which(names(peakset)==n)
    f <- file.path(annotation.dir,paste0(n,'.bed')) 
    message("Writing to ",f)
    export.bed(peakset[idx],con = f)
  }, mc.cores = 12)

  #### Co accessible peaks ####
  proj <- addCoAccessibility(
    ArchRProj = proj,
    reducedDims = "IterativeLSI"
  )
  
  #### Accessible Peaks####
  ##### One vs Others BioNames test #####
  # Identify marker peaks Per BioName
  # Performs a wilcoxon test of each group vs all others
  
  # may have to remove special characters (+/-) from names
  # https://github.com/GreenleafLab/ArchR/issues/1185
  sampleMarkersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "bioNames",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  ## Write bed files for all peaks
  library(rtracklayer)
  marker.cutOff
  no.sig.marker.cutOff <- "FDR <= 1"
  sampleMarkerList <- getMarkers(sampleMarkersPeaks, cutOff = no.sig.marker.cutOff) # used in plotbrowser track
  lengths(sampleMarkerList)
  (outDir <- file.path(getOutputDirectory(proj),'Annotations'))
  for(i in names(sampleMarkerList)){
    if(nrow(sampleMarkerList[[i]])==0)
      next
    gr <- GRanges(sampleMarkerList[[i]])
    short.name <- gsub(" ","_",i)
    gr@ranges@NAMES <- short.name
    sigLabel <- gsub("[<=&>=]","",no.sig.marker.cutOff) %>% gsub("[ ]+","_",.)
    fout <- file.path(outDir, paste0(short.name,'_',sigLabel,'.bed'))
    message("Writing file:",fout)
    bed <- data.table(data.frame(gr))[, list(seqnames, start, end)][order(seqnames, start, end)]
    write.table(bed, file = fout, sep = "\t", row.names = F, col.names = F, quote = F)
  }
  
  ##### Write bed files for significant peaks #####
  library(rtracklayer)
  marker.cutOff
  sampleMarkerList <- getMarkers(sampleMarkersPeaks, cutOff = marker.cutOff) # used in plotbrowser track
  lengths(sampleMarkerList)
  sum(lengths(sampleMarkerList))
  sampleMarkerList <- sampleMarkerList[names(sort(lengths(sampleMarkerList), decreasing = T))]
  sigLabel <- gsub("[<=&>=]","",marker.cutOff) %>% gsub("[ ]+","_",.)
  dt.bed <- lapply(names(sampleMarkerList), function(i){
    if(nrow(sampleMarkerList[[i]])==0)
      next
    gr <- GRanges(sampleMarkerList[[i]])
    short.name <- gsub(" ","_",i)
    gr@ranges@NAMES <- short.name
    dt <- data.table(data.frame(gr))
    
    ## Write Bed file
    fout <- file.path(outDir, paste0(short.name,'_',sigLabel,'.bed'))
    message("Writing file:",fout)
    bed <- dt[, list(seqnames, start, end)][order(seqnames, start, end)]
    write.table(bed, file = fout, sep = "\t", row.names = F, col.names = F, quote = F)
    
    ## Return
    bed <- dt[, list(chr=seqnames, start, end, name=".", 
                     score=MeanDiff, strand='.', signal=Log2FC, pvalue=FDR, 
                     qvalue=FDR, peak=round(width/2)) ][order(chr, start, end)]
  }) %>% set_names(names(sampleMarkerList))
  source("~/lib/scripts/R/GenomeHelpers.R")
  dt <- mergeMarkBedFiles.foverlaps(named.dt.list = dt.bed)
  dt <- reduceBed.GenomicRanges(bed = dt, bedKey = c('chr','start','end'))
  setnames(dt,"chr","#chr")
  write.table(dt, file = file.path(annotation.dir, paste0('bioNames_peaks_',sigLabel,'.bed')), 
              sep = "\t", row.names = F, col.names = T, quote = F)
  
  saveArchRProject(ArchRProj = proj, load = TRUE)
  
  
  ##### Plot bionames peak markers #####
  ## Heatmap
  heatmapSamplePeaksFilt <- plotMarkerHeatmap(
    seMarker = sampleMarkersPeaks,
    plotLog2FC = T,
    cutOff = marker.cutOff,
    returnMatrix = F,
    transpose = TRUE
  )
  heatmapSamplePeaksFilt
  
  # Write sample peaks and FDR to spreadsheets #
  dtSampleEnrichMarkerPeaks <- data.table(as.data.table(rowData(sampleMarkersPeaks)),
                                          "FDR"=do.call(cbind, assays(sampleMarkersPeaks)$FDR),
                                          "Log2FC"=do.call(cbind, assays(sampleMarkersPeaks)$Log2FC),
                                          "MeanDiff"=do.call(cbind, assays(sampleMarkersPeaks)$MeanDiff))
  writexl::write_xlsx(x = dtSampleEnrichMarkerPeaks, 
                      path = file.path(spreadsheets.dir,
                                       'archr_multiome_atac_marker_peak_FDR_by_sample.xlsx'), col_names = T, format_headers = T)
  dtSampleEnrichMarkerPeaks <- setDT( readxl::read_xlsx(file.path(spreadsheets.dir,
                                                                  'archr_multiome_atac_marker_peak_FDR_by_sample.xlsx')) )
  
  sigpeakcounts <- sapply(colnames(sampleMarkersPeaks), function(x){
    sig.str <- gsub("Log2FC",paste0("`Log2FC.",x,"`"),marker.cutOff) %>% 
      gsub("FDR",paste0("`FDR.",x,"`"),.) %>% 
      gsub("Mean",paste0("`Mean.",x,"`"), .) %>% 
      gsub("MeanDiff",paste0("`MeanDiff.",x,"`"), .)
    nrow(dtSampleEnrichMarkerPeaks[eval(parse(text = sig.str))])
  })
  sigpeakcounts
  #### Cluster peaks ####
  # Identify marker peaks Per Cluster
  # Performs a wilcoxon test of each cluster vs all others
  clusterMarkersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  
  heatmapClusterPeaksFilt <- plotMarkerHeatmap(
    seMarker = clusterMarkersPeaks,
    plotLog2FC = F,
    log2Norm = T,
    cutOff = marker.cutOff,
    returnMatrix = F,
    transpose = TRUE
  )
  heatmapClusterPeaksFilt
  
  # Shared peaks are present
  sigClusterMarkerCutoff <- getMarkers(clusterMarkersPeaks, cutOff =  marker.cutOff, returnGR = TRUE)
  
  ##### Plot all peaks by cluster #####
  htB <- DoHeatmapEnhance(marker.peaks = clusterMarkersPeaks, 
                          orientation.vertical = F,
                          marker.cutoff = "FDR <= 1 & MeanDiff > 0 & Mean > 0", 
                          sample.feature.ratio = .8, 
                          show.group.dendrogram = F, 
                          show.group.bar.count = F, 
                          group.bar = T,
                          feature.fontsize = 12, 
                          archRSort = T,
                          group.colors = cluster_palette_list,
                          heatmap.colors = ArchR::ArchRPalettes$solarExtra, 
                          #spreadsheets.path = file.path(spreadsheets.dir,'BioNames-Peak-Marker-Heatmap.xlsx'),
                          title = paste0("Peaks per cluster"))
  pdf(file = file.path(getOutputDirectory(proj),'Plots','Cluster-Peak-Marker-Heatmap_FDR1_MeanDiff.pdf'), width = 6, height = 4, bg = 'transparent')
  draw(htB, heatmap_legend_side = "bot")
  dev.off()
  
  #### Celltype peak markers ####
  # Identify marker peaks Per Cluster
  # Performs a wilcoxon test of each cluster vs all others
  celltypeMarkersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Cluster3",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  
  ##### Heatmaps of peaks by celltype #####
  heatmapCelltypePeaksFilt <- plotMarkerHeatmap(
    seMarker = celltypeMarkersPeaks,
    plotLog2FC = F,
    log2Norm = T,
    cutOff = sig.marker.cutOff,
    returnMatrix = F,
    transpose = TRUE
  )
  heatmapCelltypePeaksFilt
  
  # Shared peaks are present
  htB <- DoHeatmapEnhance(marker.peaks = celltypeMarkersPeaks,
                          orientation.vertical = F,
                          marker.cutoff = marker.cutOff,
                          sample.feature.ratio = .8, 
                          show.group.dendrogram = F,
                          feature.fontsize = 12, 
                          archRSort = T,
                          group.colors = cluster3_palette,
                          heatmap.colors = ArchR::ArchRPalettes$solarExtra, 
                          #spreadsheets.path = file.path(spreadsheets.dir,'Cluster-Peak-Marker-Heatmap.xlsx'),
                          title = paste0("Peaks per celltype"))
  pdf(file = file.path(getOutputDirectory(proj),'Plots','Celltype-Peak-Marker-Heatmap_All.pdf'), width = 6, height = 3, bg = 'transparent')
  draw(htB, heatmap_legend_side = "bot")
  dev.off()
  
  #------------------------------------------- Gene Scores -------------------------------------------
  # 1. Accessibility within the entire gene body contributes to the gene score.
  # 2. An exponential weighting function that accounts for the activity of putative distal regulatory elements in a distance-dependent fashion.
  # 3. Imposed gene boundaries that minimizes the contribution of unrelated regulatory elements to the gene score.
  # https://www.archrproject.com/bookdown/identification-of-marker-features.html
  
  
  # identify features (here gene scores) specific to groupby=bionames
  bionamesMarkersGS <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "bioNames",
    bias = c("TSSEnrichment", "log10(nFrags)"), # parameters for background generation
    testMethod = "wilcoxon"
  )
  
  ## returns log2 transformed data
  heatmapGS <- plotMarkerHeatmap(
    seMarker = bionamesMarkersGS,
    cutOff = gene.marker.cutOff,
    returnMatrix = T,
    transpose = F
  )
  
  # identify features (here gene scores) specific to groupby=clusters.
  clusterMarkersGS <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"), # parameters for background generation
    testMethod = "wilcoxon"
  )
  
  gene.marker.cutOff
  ## returns log2 transformed data
  clusterMarkersMAT <- plotMarkerHeatmap(
    seMarker = clusterMarkersGS,
    cutOff = gene.marker.cutOff,
    returnMatrix = T,
    #labelMarkers = unique(unlist(sapply(markerList, function(x){ x$name %>% head(10)} ))),
    transpose = F
  )
  dt <- as.data.table(clusterMarkersMAT, keep.rownames = "gene")
  dt <- dt[!gene %in% rem.genes]
  dt[order(C9, decreasing = T), head(gene, 100)]
  
  require(magick)
  topn <- 10
  htH <- plotEnrichGeneScoreHeatmap(clusterMarkersGS,
                                    topn = topn,
                                    cutOff = gene.marker.cutOff,
                                    title = paste0("Enriched gene scores in clusters\n","Top ",topn,
                                                   " per cluster (",gene.marker.cutOff,")"),
                                    remove.redundant = T,
                                    cluster.columns = T,
                                    labelFeatures = T, 
                                    add.score.to.name = F,
                                    labelFeaturesN = F,
                                    labelFeaturesFrac = F,
                                    outline = 0, 
                                    scaling = 'none',
                                    cell.outline = F,
                                    col = heatmap_palette, 
                                    na.col='white')
  fname <- file.path(plot.dir,"GeneScore-Enrichment-Cluster-Heatmap-Horizontal-Vierstra.pdf")
  message("Outputting file:",fname)
  
  pdf(file = fname, width = 20, height = 5.2, bg = 'transparent')
  draw(htH, heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 2), "mm")) #bltr
  dev.off()
  
  require(magick)
  topn <- 20
  htH <- plotEnrichGeneScoreHeatmap(clusterMarkersGS,
                                    topn = topn,
                                    cutOff = gene.marker.cutOff,
                                    title = paste0("Enriched gene scores in clusters\n","Top ",topn,
                                                   " per cluster (",gene.marker.cutOff,")"),
                                    remove.redundant = T,
                                    cluster.columns = T,
                                    labelFeatures = T, 
                                    add.score.to.name = F,
                                    labelFeaturesN = F,
                                    labelFeaturesFrac = F,
                                    outline = 0, 
                                    scaling = 'none',
                                    cell.outline = F,
                                    col = heatmap_palette, 
                                    na.col='white')
  fname <- file.path(plot.dir,"GeneScore-Enrichment-Cluster-Heatmap-Horizontal-Vierstra_top20.pdf")
  message("Outputting file:",fname)
  
  pdf(file = fname, width = 30, height = 5.2, bg = 'transparent')
  draw(htH, heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 2), "mm")) #bltr
  dev.off()
  
  # UMAP marker genes no imputation
  markerGenes2Label  <- unique(unlist(sapply(markerList, function(x){ head(x$name,20) } )))
  markerGenes2Label <- sort(unique(c(markerGenes2Label, unlist(markers_list))))
  ptop <- myplotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneScoreMatrix",
    name = markerGenes2Label,
    plotAs = 'points',
    embedding = "UMAP",
    size = 0.5,
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
  )
  
  #library(doParallel)
  dir.create(path = file.path(plot.dir,'marker_gene_scores_no_impute'), showWarnings = F)
  foreach(i = names(ptop)) %dopar%
    UmapPlotEnhance(ArchRProj = proj, 
                    p1 = ptop[[i]],
                    addLabels = F, 
                    pt.size = 1,
                    title = paste(i,"gene scores"),
                    legend.title = expression(Log[2]~('norm counts'+1)),
                    filename.prefix = paste0("/marker_gene_scores_no_impute/UMAP-ATAC-Marker-Genes-No-Imputation-",i))
  
  
  # UMAP Imputed marker genes
  proj <- addImputeWeights(proj)
  # UMAP marker genes no imputation
  addi <- c('Arid5a','Bach1','Foxa1','Foxb2','Foxi1','Foxl1','Hif3a','Nfkb1','Spdef','Trp63')
  markerGenes2Label  <- unique(unlist(sapply(markerList, function(x){ head(x$name,20) } )))
  markerGenes2Label <- c(markerGenes2Label, addi)
  markerGenes2Label <- sort(unique(c(markerGenes2Label, unlist(markers_list))))
  ptop <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneScoreMatrix",
    name = markerGenes2Label,
    plotAs = 'points',
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = getImputeWeights(proj)
  )
  
  
  dir.create(path = file.path(plot.dir,'marker_gene_scores_magic_impute'), showWarnings = F)
  foreach(i = names(ptop)) %dopar%
    UmapPlotEnhance(ArchRProj = proj, 
                    p1 = ptop[[i]],
                    addLabels = F, 
                    pt.size = 1,
                    title = paste("UMAP of InterativeLSI colored by\nGeneScoreMatrix:",i), 
                    legend.title = expression(Log[2]~('norm counts'+1)),
                    filename.prefix = paste0("/marker_gene_scores_magic_impute/UMAP-ATAC-Marker-Genes-Magic-Imputation-",i))
  
  
  #------------------------------------------- 3 month CCA scRNA & ATAC -------------------------------------------
  library(SingleCellExperiment)
  library(Seurat)
  # options(Seurat.object.assay.version = 'v5')
  # 
  multiRNA.file <- file.path(file.path(dirname(base.dir)),'scRNA/3months/seurat/robj/erg_epi_seurat.RDS')
  # 
  # stopped working with Seurat V5 otherwise just use this object
  seRNA <- readRDS(file = multiRNA.file)
  # 
  seATAC <- ArchRProjectToSeurat(build = T, ArchRProj = proj, matrix = "GeneScoreMatrix", output.file = NULL)
  # 
  # https://satijalab.org/seurat/articles/atacseq_integration_vignette.html
  seATAC <- ScaleData(seATAC, features = rownames(seATAC))
  
  # clear existing data
  if(length(grep("predict",colnames(proj@cellColData))))
    proj@cellColData[grep("predict",colnames(proj@cellColData))] <- NULL
  grep("predict",colnames(proj@cellColData))
  
  
  # Identify anchors
  transfer.anchors <- FindTransferAnchors(reference = seRNA,
                                          query = seATAC,
                                          features = VariableFeatures(object = seRNA),
                                          reference.assay = "RNA",
                                          query.assay = "ATAC",
                                          reduction = "cca")
  # Annotate cells
  epitype.predictions <- TransferData(anchorset = transfer.anchors,
                                      refdata = seRNA$wouter_class,
                                      weight.reduction = seATAC@reductions[['LSI']],
                                      dims = 1:29) # dim1 was removed as it is well correlated with depth
  colnames(epitype.predictions)[1] <- 'predicted.epitype.id'
  colnames(epitype.predictions)[6] <- 'predicted.score.max.epitype'
  proj <- addMetaDataArchR(ArchRProject = proj, metadata = epitype.predictions, add_prefix = "", force = T)
  
  celltype.predictions <- TransferData(anchorset = transfer.anchors,
                                       refdata = seRNA$seurat_clusters,
                                       weight.reduction = seATAC@reductions[['LSI']],
                                       dims = 1:29) # dim1 was removed as it is well correlated with depth
  
  # add predictions
  proj <- addMetaDataArchR(ArchRProject = proj, metadata = celltype.predictions, add_prefix = "", force = T)
  
  # Save annotation data
  prediction.table <- '/data/projects/sawyers/erg/analysis/multiome_mouse_WF-2149_epithelial_no_sv/archR/data/cellColData.txt'
  dt <- data.table( data.frame(proj@cellColData), keep.rownames = 'cellid')
  write.table(x = dt, file = prediction.table, col.names = T, row.names = F, sep = "\t")
  
  # non 3month epithelial cells not transferred
  proj@cellColData[is.na(proj$predicted.id),]$predicted.id <- ""
  proj@cellColData$predicted.label <- ""
  proj@cellColData[proj$predicted.id=='1',]$predicted.label <- "IM"
  proj@cellColData[proj$predicted.id=='2',]$predicted.label <- "Basal"
  proj@cellColData[proj$predicted.id=='3',]$predicted.label <- "L2 tumor"
  proj@cellColData[proj$predicted.id=='4',]$predicted.label <- "L2 tumor"
  proj@cellColData[proj$predicted.id=='5',]$predicted.label <- "L1 mixed"
  proj@cellColData[proj$predicted.id=='6',]$predicted.label <- "L1 mixed"
  proj@cellColData[proj$predicted.id=='7',]$predicted.label <- "L1 mixed"
  proj@cellColData[proj$predicted.id=='8',]$predicted.label <- "L3"
  proj@cellColData[proj$predicted.id=='9',]$predicted.label <- "L2 mixed"
  proj@cellColData[proj$predicted.id=='10',]$predicted.label <- "IM"
  proj@cellColData[proj$predicted.id=='11',]$predicted.label <- "L2 tumor"
  proj@cellColData[proj$predicted.id=='12',]$predicted.label <- "L1 mixed"
  
  proj@cellColData$Cluster2 <- ""
  proj@cellColData[proj@cellColData$Clusters %in% c("C1","C8"),]$Cluster2 <- "Basal mixed"
  proj@cellColData[proj@cellColData$Clusters %in% c("C3","C4"),]$Cluster2 <- "Basal tumor"
  proj@cellColData[proj@cellColData$Clusters %in% c("C2","C11"),]$Cluster2 <- "L2 tumor"
  proj@cellColData[proj@cellColData$Clusters %in% c("C5","C12"),]$Cluster2 <- "L2 mixed"
  proj@cellColData[proj@cellColData$Clusters %in% c("C7","C10"),]$Cluster2 <- "L1 mixed"
  proj@cellColData[proj@cellColData$Clusters %in% c("C6"),]$Cluster2 <- "IM Tp63(high)"
  proj@cellColData[proj@cellColData$Clusters %in% c("C9"),]$Cluster2 <- "IM Tp63(low)"
  
  proj@cellColData$Cluster3 <- ""
  proj@cellColData[proj@cellColData$Clusters %in% c("C1","C8"),]$Cluster3 <- "Basal"
  proj@cellColData[proj@cellColData$Clusters %in% c("C3","C4"),]$Cluster3 <- "Tumor-Bas (IM)"
  proj@cellColData[proj@cellColData$Clusters %in% c("C2","C11"),]$Cluster3 <- "Tumor-L2"
  proj@cellColData[proj@cellColData$Clusters %in% c("C5","C12"),]$Cluster3 <- "L2"
  proj@cellColData[proj@cellColData$Clusters %in% c("C7","C10"),]$Cluster3 <- "L1"
  proj@cellColData[proj@cellColData$Clusters %in% c("C6","C9"),]$Cluster3 <- "EPC (IM)"
  
  
  # save projects
  saveArchRProject(ArchRProj = proj)
  
  p1 <- plotEmbedding(
    proj, 
    colorBy = "cellColData",
    name = c("predicted.id","Clusters","Cluster2")
  )
  p1$data <- cbind(p1$predicted.id$data, label=p1$Cluster2$data$color)
  p1$data$color <- gsub("[0-9]+-","",p1$data$color)
  p1$data <- p1$data[order(p1$data$color, na.last = F),]
  p1$data$color[which(p1$data$color=="")] <- NA
  
  pal <- rev(as.vector(pals::alphabet(26)))
  names(pal) =  sort(unique(p1$data$color))
  pal <- sort(pal[unique(p1$data$color)])
  pal['8'] <- "#3A2465"
  pal['11'] <- "#F361A5"
  
  UmapPlotEnhance(proj, p1=p1, use.continuous = F, 
                  na.color = "gray60", outline.color = 'gray60',
                  pal = pal, pt.size = 0.6, feature = "predicted.id",
                  addLabels = T, label.size = 4, order.random = F,
                  legend.title = "scRNA integrated cluster",
                  filename.prefix = "RNA_constrained_integration_cluster")
  
  PlotBarVariableInClusters(proj, variable = 'predicted.id', 
                            palette = pal, orientation = 'vertical', 
                            filename = "scRNA_integrate_bar", legend.ncol = 4, 
                            height_multiplier = 0.8, width_multiplier = 0.9)
  
  
  PlotIntegrationHeatmap_ATAC_RNA(ArchRProj = proj, 
                                  rna_palette = pal,
                                  heatmap_palette = gist_heat_r,
                                  atac_cluster_name = "Clusters", 
                                  rna_cluster_name = "predicted.id")
  
  PlotIntegrationHeatmap_ATAC_RNA(ArchRProj = proj,
                                  rna_palette = pal,
                                  heatmap_palette = gist_heat_r,
                                  width_multiplier = 1.1,
                                  height_multiplier = 0.6,
                                  atac_cluster_name = "bioNames", 
                                  rna_cluster_name = "predicted.id")
  
  
  p1 <- plotEmbedding(
    proj, 
    colorBy = "cellColData",
    name = c("predicted.epitype.id","Clusters")
  )
  p1$data <- cbind(p1$predicted.epitype.id$data, label=p1$Clusters$data$color)
  p1$data$color <- gsub("[0-9]+-","",p1$data$color)
  p1$data <- p1$data[order(p1$data$color, na.last = F),]
  p1$data$color[which(p1$data$color=="")] <- NA
  
  UmapPlotEnhance(proj, p1=p1, use.continuous = F, 
                  na.color = "gray60", outline.color = 'gray60',
                  pal = epitype_palette, pt.size = 0.6, outline.size = 0.1,
                  feature = "predicted.epitype.id",
                  addLabels = T, label.size = 4, order.random = T,
                  legend.title = "Cluster epitype",
                  filename.prefix = "RNA_constrained_integration_epitype")
  
  p1 <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "cellColData",
    name = c("predicted.epitype.id","Clusters"),
    embedding = "FA",
    plotAs = 'points',
    size = 0.5,
    imputeWeights = NULL #getImputeWeights(proj)
  )
  p1$data <- cbind(p1$predicted.epitype.id$data, label=p1$Clusters$data$color)
  p1$data$color <- gsub("[0-9]+-","",p1$data$color)
  p1$data <- p1$data[order(p1$data$color, na.last = F),]
  p1$data$color[which(p1$data$color=="")] <- NA
  
  PlotBarVariableInClusters(proj, variable = 'predicted.epitype.id', 
                            palette = pal, orientation = 'vertical', 
                            filename = "scRNA_integrate_epitype_bar", legend.ncol = 2, 
                            height_multiplier = 0.8, width_multiplier = 0.9)
  
  #------------------------------------------- Define "Cluster Groups"  -------------------------------------------
  colnames(proj@cellColData)
  proj@cellColData$ClusterGroups <- ""
  proj@cellColData[proj@cellColData$Clusters %in% c("C7","C10"),]$ClusterGroups <- "L1_mixed"
  proj@cellColData[proj@cellColData$Clusters %in% c("C2","C11"),]$ClusterGroups <- "L2_tumor"
  proj@cellColData[proj@cellColData$Clusters %in% c("C9"),]$ClusterGroups <- "IM_Tp63-"
  proj@cellColData[proj@cellColData$Clusters %in% c("C6"),]$ClusterGroups <- "IM_Tp63+"
  proj@cellColData[proj@cellColData$Clusters %in% c("C3","C4"),]$ClusterGroups <- "Basal_tumor"
  proj@cellColData[proj@cellColData$Clusters %in% c("C1","C8"),]$ClusterGroups <- "Basal_mixed"
  proj@cellColData[proj@cellColData$Clusters %in% c("C5","C12"),]$ClusterGroups <- "L2_mixed"
  
  saveArchRProject(ArchRProj = proj)
  
  proj <- loadArchRProject(path = archr.dir)
  
  #------------------------------------------- 1. Motif calling -------------------------------------------
  getAvailableMatrices(proj)
  
  library(TFBSTools)
  
  system(paste("wget https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/Vierstra-Human-Motifs.rds -O", file.path(robj.dir,"Vierstra-Human-Motifs.rds")))
  mPWMs <- readRDS( file.path(robj.dir,"Vierstra-Human-Motifs.rds") )
  names(mPWMs)
  # update names
  names(mPWMs) <- gsub(":","-",names(mPWMs))
  motifPWMs <- lapply(mPWMs, function(i){
    i@name <- gsub(":","-",i@name)
    i
  })
  motifPWMs <- do.call(PWMatrixList, motifPWMs)
  
  vierstramotifs <- gsub("_.*","",names(motifPWMs)) %>% unique %>% toupper()
  cisbpmotifs <- gsub(".*_LINE[0-9]+_","",names(chromVARmotifs::mouse_pwms_v2)) %>% gsub("_.*","",.) %>% unique %>% toupper()
  missing <- setdiff(cisbpmotifs, vierstramotifs) %>% sort
  
  saveRDS(object = motifPWMs, file = file.path(robj.dir,"Vierstra-Human-Motifs-Update.rds"))
  
  proj <- addMotifAnnotations(ArchRProj = proj,
                              motifSet = "cisbp",
                              cutOff = 5e-05, # p-value cutoff against bg set
                              width = 7, # motif width
                              name = "cisbp",
                              force = T)
  
  # Scan for all motifs (human+mouse)
  motifPWMs <- readRDS(file = file.path('/home/erik/annotation/vierstra/Vierstra-Human-Motifs-Update.rds'))
  names(motifPWMs)
  proj <- addMotifAnnotations(proj, 
                              motifPWMs = motifPWMs, 
                              annoName = "Vierstra",
                              cutOff = 5e-05, # p-value cutoff against bg set
                              width = 7, # motif width
                              force = T)
  saveArchRProject(ArchRProj = proj, overwrite = T)
  
  getAvailableMatrices(proj)
  
  # Filter for motifs that are not well defined or expression is missing
  if(!'VierstraMotifToGene' %in% names(attributes(proj))){
    attr(proj, 'VierstraMotifToGene') <- Motif2GeneLookupDT(proj = proj, motifpwm = proj@peakAnnotation@listData$Vierstra$motifs)
  }
  (dtM2Glookup <- proj@VierstraMotifToGene)
  
  #-------------------------------------------2. Chromvar Deviations -------------------------------------------
  # chromVar is designed for predicting enrichment of TF activity on a per-cell basis from sparse chromatin accessibility data. The two primary outputs of chromVAR are:
  # 1. “deviations” - A deviation is a bias-corrected measurement of how far the per-cell accessibility of a given feature (i.e motif) 
  # deviates from the expected accessibility based on the average of all cells or samples.
  # 2. “z-score” - The z-score, also known as a “deviation score” is the z-score for each bias-corrected deviation across all cells.
  #.    The absolute value of the deviation score is correlated with the per-cell read depth. 
  #     This is because, with more reads, you have higher confidence that the difference in per-cell accessibility of the given feature (i.e. motif) from the expectation is greater than would occur by chance.
  ## stores both deviations and z-scores (deviation score)
  getAvailableMatrices(proj)
  proj <- addDeviationsMatrix(proj, peakAnnotation = "cisbp", force = T, threads = 8)
  
  getPeakAnnotation(proj)$Name
  
  proj <- addDeviationsMatrix(ArchRProj = proj, 
                              peakAnnotation = "Vierstra",
                              force = TRUE, 
                              threads = 6)
  
  saveArchRProject(ArchRProj = proj, overwrite = T)
  
  
  archrMotifNames <- data.table(rownames(proj@peakAnnotation$cisbp$motifSummary)[
    match(names(chromVARmotifs::mouse_pwms_v2), proj@peakAnnotation$cisbp$motifSummary$ID)
  ], proj@peakAnnotation$cisbp$motifSummary$ID)
  
  #------------------------------------------- 3. Motif Significant Enrichment -------------------------------------------
  
  getAvailableMatrices(proj)
  #-----sample enrichment summary ----#
  #### Per Sample - get enriched peaks overlapped with motifs ####
  sampleMarkersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "bioNames",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  
  #-----cluster enrichment summary ----#
  #### Per Cluster - get enriched peaks overlapped with motifs ####
  clusterMarkersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  # Per Cluster - get motifs within enriched peaks
  clusterEnrichMotifVierstra <- peakAnnoEnrichment(
    seMarker = clusterMarkersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Vierstra",
    cutOff = marker.cutOff # cutoff applies to seMarker
  )
  
  allClusterEnrichMotifVierstra <- peakAnnoEnrichment(
    seMarker = clusterMarkersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Vierstra",
    cutOff = "FDR <= 1" # cutoff applies to seMarker
  )
  
  ## Motif differences ##
  clusterMarkersMotifs <- getMarkerFeatures(
    ArchRProj = proj, 
    testMethod = "wilcoxon",
    binarize = FALSE,
    useMatrix = "VierstraMatrix",
    groupBy = "Clusters",
    useSeqnames="z"
  )
  
  #-----cluster groups enrichment summary ----#
  #### Per cluster groups class - get enriched peaks ####
  ClusterGroupsMarkersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "ClusterGroups",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
  )
  
  ClusterGroupsEnrichMotifVierstra <- peakAnnoEnrichment(
    seMarker = ClusterGroupsMarkersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Vierstra",
    cutOff = marker.cutOff # cutoff applies to seMarker
  )
  
  ## Motif differences ##
  ClusterGroupsMarkersMotifs <- getMarkerFeatures(
    ArchRProj = proj, 
    testMethod = "wilcoxon",
    binarize = FALSE,
    useMatrix = "VierstraMatrix",
    groupBy = "ClusterGroups",
    useSeqnames="z"
  )
  ClusterGroupsMarkersMotifs
  
  
  #-----cluster3 enrichment summary ----#
  #### Per cluster groups class - get enriched peaks ####
  Cluster3MarkersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Cluster3",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  
  marker.cutOff
  Cluster3EnrichMotifVierstra <- peakAnnoEnrichment(
    seMarker = Cluster3MarkersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Vierstra",
    cutOff = marker.cutOff # cutoff applies to seMarker
  )
  
  # Write enrichment to spreadsheets #
  dtCluster3EnrichMotifVierstra <- data.table("motifID"=names(Cluster3EnrichMotifVierstra),
                                              "mlog10Padj"=do.call(cbind, assays(Cluster3EnrichMotifVierstra)$mlog10Padj),
                                              "enrichment"=do.call(cbind, assays(Cluster3EnrichMotifVierstra)$Enrichment))
  writexl::write_xlsx(x = dtCluster3EnrichMotifVierstra,
                      path = file.path(spreadsheets.dir, 'archr_multiome_atac_motifs_EP_Cluster3EnrichMotifVierstra.xlsx'), 
                      col_names = T, format_headers = T)
  
  
  #### Plot Cluster3 heatmap ####
  require(magick)
  pval.cutoff <- 5e-5
  min.cutoff <- -log10(pval.cutoff)

  require(magick)
  pval.cutoff <- 5e-5
  min.cutoff <- -log10(pval.cutoff)
  htH <- plotEnrichMotifHeatmap(Cluster3EnrichMotifVierstra,
                                n = Inf,
                                title = paste0("Enriched T.F. motifs in cluster group\n(padj <",pval.cutoff,")"),
                                remove.redundant = T, 
                                cluster.columns = T,
                                mlog10Padj.cutOff = min.cutoff,
                                labelMotifs = F, 
                                outline = 0, 
                                scaling = 'none',
                                cell.outline = F,
                                col = ArchRPalettes$comet, 
                                na.col='white')
  fname <- file.path(plot.dir,"Motif-Enrichment-Cluster3-Heatmap-Horizontal-Vierstra.png")
  message("Outputting file:",fname)
  res <- 200
  w <- res*8 # res * inches
  h <- res*2.5
  fig <- image_graph(width = w, height = h, res = res)
  draw(htH, heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 5), "mm")) #bltr
  dev.off()
  image_write(image = fig, path = fname, format = 'png')
  
  
  #------------------------------------------- Variance Rank Filter Grid by CellTypes -------------------------------------------
  dtVierstraVarDevFilt <- data.table(read_xlsx(file.path(spreadsheets.dir,
                                                         'archr_multiome_atac_motifs_vierstra.xlsx')))
  dtCluster3EnrichMotifVierstra <- readxl::read_xlsx(path = file.path(spreadsheets.dir,
                                                                      'archr_multiome_atac_motifs_EP_Cluster3EnrichMotifVierstra.xlsx')) %>% setDT()
  dtMotifGeneCorCluster3NoImpute <- correlateAccessibilityExpression(build=F,
                                                                     ArchRproj = proj, 
                                                                     dtInputChromVar = dtVierstraVarDevFilt,
                                                                     vassay = vassay, 
                                                                     expassay = expassay, 
                                                                     geneMatrix = "GeneScoreMatrix",
                                                                     #matGE = matGE,
                                                                     min.expression = 0.0,
                                                                     groupby = 'Cluster3',
                                                                     ncores = 16, 
                                                                     cor.method = 'spearman',
                                                                     xlsx='dtMotifGeneCorCluster3NoImputeSpearman.xlsx')
  dtCor <- dtMotifGeneCorCluster3NoImpute
  
  # This simply combines correlation and chromvar statistics with sample enrichment
  dt <- dtCluster3EnrichMotifVierstra[,c("motifID", grep("mlog10Padj", colnames(dtCluster3EnrichMotifVierstra), value = T)), with=F]
  setnames(dt, colnames(dt), gsub("mlog10Padj.","",colnames(dt)))
  dt <- melt(dt, id.vars = 'motifID', variable.name = 'Cluster3', value.name = 'mlog10Padj')
  ms <- merge(dtCor, dt, by=c('motifID','Cluster3'))
  
  ms[,fisher_p := fisherCpvalLog10(c(cor_mlog10Padj, mlog10Padj))$p, by=1:nrow(ms) ]
  ms[,m10log_fisher_p := -log10(max(c(fisher_p, .Machine$double.xmin))), by=1:nrow(ms)]
  
  ms[,score:= mlog10Padj * sign(avg_accessibility) * frac_exp]
  ms <- ms[order(abs(score), m10log_fisher_p, mlog10Padj, cor_mlog10Padj, decreasing = T)]
  ms[,ridx := 1:nrow(ms)]
  ms[,short_motifID := gsub("_.*","", toupper( motifID) )] # remove some redundancy
  ms[,short_familyID := gsub("/.*","",family)]
  
  dtClusterData <- fread('/data/projects/sawyers/erg/analysis/scRNA/3months/seurat//data/average_expression_singlecell_rna_seq_cluster_data.csv.gz')
  setnames(dtClusterData, old = colnames(dtClusterData), new = gsub("seurat_cluster_","",colnames(dtClusterData)))
  
  expData <- RelateExpressionData(ArchRProj = proj, RnaAverageExpression = dtClusterData, groupby = 'Cluster3')
  expData <- melt(expData, id.vars='gene', value.name = "AvgExpression")
  
  ## ms contains enrichment per Cluster3 vs others in mlog10Padj from dtCluster3EnrichMotifVierstra
  ## and it contains geneCorrelated, coefficient, cor, cor_pval, etc from fn: correlateAccessibilityExpression
  ## we rank score calculated by mlog10Padj * sign(avg_accessibility) * fraction of cells expressing gene
  
  # 186660
  nrow(ms)
  ms <- merge(ms, expData, by.x=c("Cluster3", "motifGene"), by.y=c("variable", "gene"))
  nrow(ms)
  
  top_familyN <- 1
  message("Getting top ",top_familyN," motifs per family")
  # set pval threshold
  #pval.cutoff <- 1e-25
  pval.cutoff <- 1e-20
  
  min.cutoff <- round(-log10(pval.cutoff))
  ms[,plot := F]
  mark <- ms[cor_padj < 0.01 & 
               score >= min.cutoff &
               frac_both > 0.1 &
               ((sign(avg_accessibility) == sign(avg_expression) & AvgExpression >=1 & coefficient > 0) |
                  (sign(avg_accessibility) != sign(avg_expression) & AvgExpression >= 1 & coefficient < 0))]
  mark <- mark[order(score, decreasing = T)]
  #mark <- mark[mark[, .I[head(1)] , by = list(short_motifID, family)]$V1][order(score, decreasing=T),head(.SD, top_familyN),by=list(family)][,ridx]
  mark <- mark[mark[, .I[head(1)] , by = list(short_motifID, family)]$V1][order(score, decreasing=T),head(.SD, top_familyN),by=list(family, Cluster3)][,ridx]
  mark <- ms[ridx %in% mark]
  setkey(ms, motifID, family, motifGene, geneCorrelated)
  mark <- lapply(mark$ridx, function(i){
    r <- mark[ridx==i, list(motifID, family, motifGene, geneCorrelated)]
    ms[r, unlist(ridx)]
  }) %>% unlist
  
  ms[ridx %in% mark, plot := T]
  nrow(ms[plot==T])
  ms[plot==T]
  
  writexl::write_xlsx(x = ms, path = file.path(spreadsheets.dir,'Cluster3_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx'), format_headers = T, col_names = T)
  ms <- setDT(readxl::read_xlsx(path = file.path(spreadsheets.dir,'Cluster3_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx')))
  
  msmat <- as.matrix(dcast.data.table(ms[plot==T], motifID ~ Cluster3, value.var = 'score', fun.aggregate = absmax), rownames="motifID")
  msmat[is.na(msmat)] <- 0
  msmat[is.infinite(msmat)] <- 0
  
  msgroup <- ms[ms[plot==T, .I[which.max(mlog10Padj)], by=list(motifID, Cluster3)]$V1]
  msgroup <- msgroup[,.(same=ifelse(motifGene==geneCorrelated,'*','')), by=list(motifID, Cluster3)]
  msgroup <- as.matrix(dcast.data.table(msgroup, motifID ~ Cluster3, value.var = 'same'), rownames="motifID")
  msgroup <- msgroup[unique(ms[plot==T, motifID]),]
  msgroup[is.na(msgroup)] <- ""
  
  identical(rownames(msmat), rownames(msgroup))
  identical(colnames(msmat), colnames(msgroup))
  
  #### Fig 6D ####
  nrow(msmat)
  htV <- plotEnrichMotifHeatmap(msmat,
                                geneMotifSameMat=msgroup,
                                n = Inf,
                                orientation_horizontal = F,
                                remove.redundant = T,
                                cluster.columns = T, #F
                                cluster.rows = T,
                                labelMotifClusters= '',
                                label.cluster.family=T,
                                labelMotifN = T, 
                                mlog10Padj.cutOff = min.cutoff,
                                labelMotifs = T,
                                title = paste0("Enriched TF motifs per Group\npadj < ",pval.cutoff), 
                                title.legend = "Enrichment score",
                                scaling = 'motif',
                                row_name_fontsize = 12,
                                column_name_fontsize = 16,
                                na.col = 'white',
                                legend.direction = 'horizontal',
                                col = solarExtraWhite_palette,
                                label.col = godsnot_102) #godsnot_102
  
  
  
  cairo_pdf(filename = file.path(plot.dir,paste0("Motif-Heatmap-Cluster3-vertical_Vierstra_top",top_familyN,"_family.pdf")),
            width = 6, height = 14, bg = 'transparent', onefile = T, family = 'Helvetica')
  draw(htV, heatmap_legend_side = "bottom", padding = unit(c(5, 2, 2, 8), "mm")) #bltr
  dev.off()
  
  
  #------------------------------------------- Plot Reductions -------------------------------------------
  (g <- grep( "ETS1_400" , names(proj@peakAnnotation$Vierstra$motifs), value=T, ignore.case = T))
  
  motifId <- g[1]
  smotifId <- gsub("#.*","",motifId)
  gp <- UmapPlotEnhance(ArchRProj = proj,
                        p1 = NULL, 
                        feature = paste0('z:',motifId),
                        colorBy = 'VierstraMatrix',
                        pt.size = 1,
                        addLabels = F, 
                        embedding = 'UMAP',
                        title = smotifId, 
                        label.size = 14, 
                        title.size = 18, 
                        break.unit = NA, 
                        panel_width_in = 3,
                        scale_color_bar = T,
                        pal = solarExtraWhite_palette,
                        legend.title = 'Motif deviation zscore',
                        legend.size = 18, 
                        legend.height.factor = .6,
                        filename.prefix = paste0('Umap_zDeviation_',smotifId))
  gp

  
  #------------------------------------------- Motif co-occurrences -------------------------------------------
  getAvailableMatrices(proj)
  
  peaks <- getPeakSet(proj)
  peakRS <- getMatrixFromProject(proj, useMatrix='PeakMatrix')
  summary(peakRS)
  
  # rows are peaks and columns are cells
  # dim(peakRS)
  #peakmat <- assays(peakRS)[[1]]
  
  # Genes linked to significantly different peaks #
  # Rows are peaks and columns are motifs
  motif_matches <- getMatches(proj, "Vierstra")
  ## easier but same above 
  ## motif_matches <- readRDS(proj@peakAnnotation$Vierstra$Matches)
  summary(motif_matches)
  gr <- rowRanges(motif_matches)
  
  #### Cell type marker peaks ####
  # These are all peaks in the peak matrix grouped by celltype
  # The getMarkers function does not return disjoint peaks
  groupBy <- "Cluster3"
  celltypeMarkersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = groupBy,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  # genomic locations
  rowData(celltypeMarkersPeaks)
  
  ## Cell type aka Cluster3 cutoff
  allcelltypeMarkerCutoff <- getMarkers(celltypeMarkersPeaks, cutOff =  "FDR <= 1", returnGR = TRUE)
  allcelltypeMarkerCutoff <- lapply(allcelltypeMarkerCutoff, sort)
  lengths(allcelltypeMarkerCutoff)
  # add grouping to row
  # not disjoint!
  allcelltypeMarkerCutoff <- lapply( names(allcelltypeMarkerCutoff), function(x){
    mcols(allcelltypeMarkerCutoff[[x]])[eval(groupBy)] <- x
    sort(allcelltypeMarkerCutoff[[x]])
  }) %>% set_names( names(allcelltypeMarkerCutoff) )
  ## unlist
  allcelltypeMarkerCutoff <- unlist(as(allcelltypeMarkerCutoff, "GRangesList"), use.names = F)
  ## Make disjoint
  ## sort by location and FDR, then get unique
  allcelltypeMarkerCutoff <- sort(allcelltypeMarkerCutoff, by=~seqnames+start+end+FDR)
  alldisjointcelltypeMarkerCutoff <- unique(allcelltypeMarkerCutoff)
  
  
  sig.marker.cutOff <- "FDR <= 0.05 & Log2FC >= 1"
  sigLabel <- gsub("[<=&>=]","",sig.marker.cutOff) %>% gsub("[ ]+","_",.) %>% paste0("disjoint_",.)
  sigcelltypeMarkerCutoff <- getMarkers(celltypeMarkersPeaks, cutOff = sig.marker.cutOff, returnGR = TRUE)
  
  # add cluster name to row
  # not disjoint!
  sigcelltypeMarkerCutoff <- lapply( names(sigcelltypeMarkerCutoff), function(x){
    mcols(sigcelltypeMarkerCutoff[[x]])[eval(groupBy)] <- x
    sort(sigcelltypeMarkerCutoff[[x]])
  }) %>% set_names( names(sigcelltypeMarkerCutoff) )
  lapply(sigcelltypeMarkerCutoff, function(x){ x[seqnames(x)=="chr1" & start(x)==13197067] })
  
  # make disjoint
  sigcelltypeMarkerCutoff <- unlist(as(sigcelltypeMarkerCutoff, "GRangesList"), use.names = F)
  ## sort by location and FDR, then get unique
  sigcelltypeMarkerCutoff <- sort(sigcelltypeMarkerCutoff, by=~seqnames+start+end+FDR)
  sigcelltypeMarkerCutoff <- unique(sigcelltypeMarkerCutoff)
  ##
  table(mcols(sigcelltypeMarkerCutoff)[eval(groupBy)]) # disjoint
  
  motif_cooccur_peak_count_matrix <- function(motif_matches, markers, save_obj_path="", overwrite=F, no.cores=12){
    
    if(file.exists(save_obj_path) & overwrite != T)
      stop("File exists..exiting. Please change overwrite to True or save to different file")
    
    gr <- rowRanges(motif_matches)
    ov <- findOverlaps(query = markers, subject = gr)
    qidx <- queryHits(ov)
    sidx <- subjectHits(ov)
    
    # hould be same as markers
    marker_hits <- markers[qidx]
    
    # subset query
    gr <- unname(gr[sidx])
    
    identical(ranges(gr), ranges(marker_hits))
    
    M <- (assays(motif_matches)$matches)[sidx, ]
    
    motifs <- colnames(motif_matches)
    H <- matrix(data = 0, nrow = length(motifs), ncol = length(motifs))
    counts <- NULL
    # A sparse matrix: [peaks x motif], indicating if motif is found
    #m <- lapply(motifs, function(motif_of_interest1){
    for(i in seq_along(motifs)){
      motif_of_interest1 <- motifs[i]
      message(i,") hits with ",motif_of_interest1)
      M1 <- M[, motif_of_interest1]
      motifs2 <- motifs[seq(i,length(motifs), by=1)]
      rn <- mclapply(motifs2, function(motif_of_interest2){
        M2 <- M[, motif_of_interest2]
        d <- unname(which( (M1 & M2) > 0 ))
        length(d)
      }, mc.cores = no.cores)
      counts <- c(counts, do.call(what = c, rn))
    }
    H[lower.tri(H, diag=TRUE)] <- as.integer(counts)
    H <- t(H)
    rownames(H) <- colnames(H) <- motifs
    if(save_obj_path != ""){
      message("Saving counts to ",save_obj_path)
      saveRDS(H, save_obj_path)
    }
    return(H)
  }
  
  
  motif_cooccur_peak_count_matrix <- function(motif_matches, markers, motif_of_interest, save_obj_path="", overwrite=F, no.cores=12){
    
    if(file.exists(save_obj_path) & overwrite != T)
      stop("File exists..exiting. Please change overwrite to True or save to different file")
    
    gr <- rowRanges(motif_matches)
    ov <- findOverlaps(query = markers, subject = gr)
    qidx <- queryHits(ov)
    sidx <- subjectHits(ov)
    
    # hould be same as markers
    marker_hits <- markers[qidx]
    
    # subset query
    gr <- unname(gr[sidx])
    
    identical(ranges(gr), ranges(marker_hits))
    
    M <- (assays(motif_matches)$matches)[sidx, ]
    
    motifs <- colnames(motif_matches)
    counts <- NULL
    # A sparse matrix: [peaks x motif], indicating if motif is found
    #m <- lapply(motifs, function(motif_of_interest1){
    
    message("Finding hits with ",motif_of_interest1)
    M1 <- M[, motif_of_interest1]
    counts <- mclapply(motifs, function(motif_of_interest2){
      M2 <- M[, motif_of_interest2]
      d <- unname(which( (M1 & M2) > 0 ))
      length(d)
    }, mc.cores = no.cores)
    counts <- unlist(counts)
    names(counts) <- motifs
    counts
    
  }
  
  motif_peak_hits <- function(motif_matches, markers, peakRS=NULL, motif_of_interest1, motif_of_interest2=NULL){
    
    gr <- rowRanges(motif_matches)
    ov <- findOverlaps(query = markers, subject = gr)
    qidx <- queryHits(ov)
    sidx <- subjectHits(ov)
    
    # should be same as markers
    marker_hits <- markers[qidx]
    
    # subset query
    gr <- sort(unname(gr[sidx]))
    marker_hits <- sort(unname(marker_hits))
    
    if( !identical(ranges(gr), ranges(marker_hits)) )
      stop("Stopping..motif ranges and hits do not contain same sets")
    
    # peak x motif
    M <- (assays(motif_matches)$matches)[sidx, ]
    M1 <- NULL
    # A sparse matrix: [peaks x motif], indicating if motif is found
    if(length(motif_of_interest1)>1){
      m <- lapply(motif_of_interest1, function(i){
        M[, i]
      })
      M1 <- Reduce('+', m)
    } else {
      M1 <- M[, motif_of_interest1] 
    }
    if(!is.null(motif_of_interest2)){
      M2 <- M[, motif_of_interest2]
      d <- unname(which( (M1 & M2) > 0 ))
      marker_hits <- marker_hits[d,]
      message(length(marker_hits)," hits for ", motif_of_interest1," and ",motif_of_interest2)
    } else {
      d <- unname(which( M1 > 0 ))
      marker_hits <- marker_hits[d,]
      message(length(marker_hits)," hits for ", motif_of_interest1)
    }
    if(!is.null(peakRS)){
      ov <- findOverlaps(query = marker_hits, subject = peakRS)
      peakSS <- peakRS[subjectHits(ov),]
      if(length(peakSS) != length(marker_hits))
        stop("Peak matrix and marker hits are different")
      group.name <- names(mcols(marker_hits)[ncol(mcols(marker_hits))])
      mcols(marker_hits)['ncells'] <- 0
      mcols(marker_hits)['fraccells'] <- 0
      mcols(marker_hits)['nhits'] <- 0
      for(i in unique(mcols(marker_hits)[[group.name]])){
        message("Getting per cell peak counts for ",i)
        # get peak locations
        pidx <- which(mcols(marker_hits)[[group.name]]==i)
        # get cells
        total.cells.group <- length(colnames(peakRS)[which(peakRS[[group.name]]==i)])
        cells <- colnames(peakSS)[which(peakSS[[group.name]]==i)]
        ncells <- rowSums2(sign(assays(peakSS)[[1]][pidx, cells]))
        mcols(marker_hits)[pidx,'ncells'] <- ncells
        mcols(marker_hits)[pidx,'fraccells'] <- ncells/total.cells.group
        # num hits
        nhits <- rowSums2(assays(peakSS)[[1]][pidx, cells])
        mcols(marker_hits)[pidx,'nhits'] <- nhits
        
      }
    }
    
    marker_hits
  }
  
  # Overlap GR list of significant peaks with motif positions
  # retrieve motif hits that match an interest
  
  ## Significant peaks in cell types
  sig.marker.cutOff
  
  motif.dir <- file.path(annotation.dir,'motif_hits')
 
  #### Co-occur motif lists ####
  peaks <- getPeakSet(proj)
  peaks <- unname(peaks) # most similar to disjoint peaks
  
  ## Motifs identified with chromvar in Cluster3 that correlate with gene score
  ms <- setDT(readxl::read_xlsx(path = file.path(spreadsheets.dir,'Cluster3_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx')))
  sigmotifs <- unique(ms[plot==T, motifID])
  grep("\\+",sigmotifs, value=T)
  sigmotifs <- gsub(".+\\+","",sigmotifs)
  
  ## This is pairwise
  list_cooccupy_motif_list <- lapply(sigmotifs, function(i){
    (m <- grep(i, colnames(motif_matches), ignore.case = T, value = T))
    ys <- sapply(sigmotifs, function(j) grep(j, colnames(motif_matches), ignore.case = T, value = T) )
    cooccupy_motif_list <- mclapply(ys, function(i){
      message("Evaluating ",m," with ",i)
      motif_peak_hits(motif_matches = motif_matches, markers = peaks, 
                      motif_of_interest1 = m, motif_of_interest2 = i)
    }, mc.cores = 4)
    names(cooccupy_motif_list) <- paste( gsub("#.*","",m), gsub("#.*","",ys), sep = "_and_")
    cooccupy_motif_list
  })
  names(list_cooccupy_motif_list) <- gsub("#.*","",sigmotifs)
  
  names(list_cooccupy_motif_list$Gata2_383)
  
  list_cooccupy_motif_list <- unlist(unname(list_cooccupy_motif_list), recursive = FALSE)
  
  ## motif hits not overlapping ETS1_400
  grep("NFKB2_1188", colnames(motif_matches), value = T)
  region_interest <- motif_peak_hits(motif_matches = motif_matches,
                                     markers = peaks, 
                                     motif_of_interest1 = "P63_1342#P53-like/1-P53L-257-0",
                                     motif_of_interest2 = NULL)
  
  #### ETS1_400 motif ####
  (ets1.motif <- grep("ETS1_400", colnames(motif_matches), value = T))
  ets1_regions <- motif_peak_hits(motif_matches = motif_matches,
                                  markers = peaks, 
                                  motif_of_interest1 = ets1.motif,
                                  motif_of_interest2 = NULL)
  no_ets1_motifs <- setdiff(sigmotifs, ets1.motif)
  # No overlap with ETS1
  list_motif_not_with_ETS1 <- lapply(no_ets1_motifs, function(i){
    ## need the exact motif name
    (m <- grep(i, colnames(motif_matches), ignore.case = T, value = T))
    message("Evaluating ",m," without ETS")
    x <- motif_peak_hits(motif_matches = motif_matches, 
                         markers = peaks, 
                         motif_of_interest1 = m, motif_of_interest2 = NULL)
    ## Not in ETS regions
    x[!x %over% ets1_regions,]
  }) %>% set_names(., paste0(gsub("#.*","",no_ets1_motifs),"_no_ETS1"))
  # Cooccur with ETS1
  list_cooccupy_motif_ETS1 <- lapply(no_ets1_motifs, function(i){
    (m <- grep(i, colnames(motif_matches), ignore.case = T, value = T))
    message("Evaluating ",m," with ETS")
    x <- motif_peak_hits(motif_matches = motif_matches, markers = peaks, 
                         motif_of_interest1 = m, motif_of_interest2 = NULL)
    # Overlap with Ets
    x[x %over% ets1_regions,]
  }) %>% set_names(., paste( gsub("#.*","",no_ets1_motifs), "ETS1", sep = "_and_"))
  
  #### Several ETS factors ####
  (ets.motifs <- grep("^ETS|^ETV|^ERG|^FLI", colnames(motif_matches), value = T))
  ets_regions <- lapply(ets.motifs, function(i){
    motif_peak_hits(motif_matches = motif_matches,
                    markers = peaks, 
                    motif_of_interest1 = i,
                    motif_of_interest2 = NULL)
  })
  ets_regions <- unlist(as(ets_regions, "GRangesList"), use.names = F)
  ## sort by location and FDR, then get unique
  ets_regions <- unique(sort(ets_regions, by=~seqnames+start+end))
  
  no_ets_motifs <- setdiff(sigmotifs, ets.motifs)
  list_motif_not_with_ets <- lapply(no_ets_motifs, function(i){
    ## need the exact motif name
    (m <- grep(i, colnames(motif_matches), ignore.case = T, value = T))
    message("Evaluating ",m," without ETS")
    x <- motif_peak_hits(motif_matches = motif_matches, 
                         markers = peaks, 
                         motif_of_interest1 = m, motif_of_interest2 = NULL)
    ## Not in ETS regions
    x[!x %over% ets_regions,]
  }) %>% set_names(., paste0(gsub("#.*","",no_ets_motifs),"_no_ETS"))
  
  ## Overlap significant with ETS
  list_cooccupy_motif_ets <- lapply(no_ets_motifs, function(i){
    (m <- grep(i, colnames(motif_matches), ignore.case = T, value = T))
    message("Evaluating ",m," with ETS")
    x <- motif_peak_hits(motif_matches = motif_matches, markers = peaks, 
                         motif_of_interest1 = m, motif_of_interest2 = NULL)
    # Overlap with Ets
    x[x %over% ets_regions,]
  }) %>% set_names(., paste( gsub("#.*","",no_ets_motifs), "ETS", sep = "_and_"))
  
  ## ETS motifs without candidates
  list_ets_without_motif <- lapply(no_ets_motifs, function(i){
    (m <- grep(i, colnames(motif_matches), ignore.case = T, value = T))
    message("Evaluating ",m," with ETS")
    x <- motif_peak_hits(motif_matches = motif_matches, markers = peaks, 
                         motif_of_interest1 = m, motif_of_interest2 = NULL)
    # Ets without motif
    ets_regions[!ets_regions %over% x]
  }) %>% set_names(., paste( "ETS", gsub("#.*","",no_ets_motifs), sep = "_no_"))
  ## Add all the ETS factors
  list_ets_without_motif[['ETS_motifs']] <- ets_regions
  
  
  
  ## Plot these chromvar
  gp <- lapply(ets.motifs, function(motifId){
    smotifId <- gsub("#.*","",motifId)
    g <- UmapPlotEnhance(ArchRProj = proj,
                         p1 = NULL, 
                         feature = paste0('z:',motifId),
                         colorBy = 'VierstraMatrix',
                         pt.size = 1,
                         addLabels = F, 
                         embedding = 'UMAP',
                         title = smotifId, 
                         label.size = 14, 
                         title.size = 18, 
                         break.unit = NA, 
                         panel_width_in = 3,
                         scale_color_bar = T,
                         pal = solarExtraWhite_palette,
                         legend.title = 'Motif deviation zscore',
                         legend.size = 18, 
                         legend.height.factor = .6)
    #filename.prefix = paste0('Umap_zDeviation_',smotifId))
  })
  length(gp)
  f <- file.path(plot.dir,"ETS-like_motif_chromvar_umaps.png")
  require(magick)
  message("Plotting to:",f)
  res <- 300
  fig <- image_graph(width = 30*res, height = 30*res, res = res, bg = 'transparent')
  do.call(grid.arrange, c(gp, nrow=8))
  dev.off()
  image_write(image = fig, path = f, format = 'png')
  
  
  ## ECDF
  dt <- as.data.table(ets_regions)
  short.name <- gsub("#.*","",motif_of_interest1[1])
  short.name <- 'ETS regions'
  #cairo_pdf(filename = file.path(plot.dir,paste0(short.name,"_ecdf.pdf")), width = 8.5, height = 6)
  png(filename = file.path(plot.dir,paste0(short.name,"_ecdf.png")), width = 7, height = 4, units = 'in', res = 300)
  ggplot(dt, aes(x=Log2FC, color=get(groupBy))) +
    stat_ecdf(geom = "step", linewidth=1.2) +
    xlab(bquote(Log[2]*Fc)) + ylab("Fn(x)") +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_x_continuous(limits = c(-3,3), expand = c(0,0)) +
    #facet_wrap(~condition) +
    ggtitle(paste(short.name,"peak accessibility")) +
    geom_hline(yintercept = 0.5, linetype='dotted') +
    geom_vline(xintercept = 0, linetype='dotted') +
    scale_color_manual("", values = cluster3_palette) +
    theme_bw(18) + theme(text = element_text(family='Arial'),
                         axis.text = element_text(size=18),
                         legend.text = element_text(size=18),
                         panel.grid = element_blank(),
                         plot.background = element_blank(),
                         panel.background = element_blank(),
                         legend.background = element_blank(),
                         strip.background = element_blank(),
                         strip.text = element_text(size=12, face='bold'))
  dev.off()
  
  #### Co-occur ChromVar  ####
  
  # 1. create object
  proj <- addPeakAnnotations(ArchRProj = proj, regions = list_cooccupy_motif_list, name = "vierstra_cooccur_motif", force = T)
  
  proj <- addPeakAnnotations(ArchRProj = proj, regions = list_motif_not_with_ETS1, name = "vierstra_motif_nooverlap_ETS1", force = T)
  
  proj <- addPeakAnnotations(ArchRProj = proj, regions = list_cooccupy_motif_ETS1, name = "vierstra_cooccur_motif_with_ETS1", force = T)
  
  proj <- addPeakAnnotations(ArchRProj = proj, regions = list_motif_not_with_ets, name = "vierstra_motif_nooverlap_ets", force = T)
  
  proj <- addPeakAnnotations(ArchRProj = proj, regions = list_cooccupy_motif_ets, name = "vierstra_cooccur_motif_with_ets", force = T)
  
  proj <- addPeakAnnotations(ArchRProj = proj, regions = list_ets_without_motif, name = "vierstra_motif_ets_nooverlap_candidate", force = T)
  
  # 2. Run chromvar on these and store to project
  proj <- addDeviationsMatrix(proj, "vierstra_cooccur_motif", force = T)
  
  proj <- addDeviationsMatrix(proj, "vierstra_motif_nooverlap_ETS1", force = T)
  
  proj <- addDeviationsMatrix(proj, "vierstra_cooccur_motif_with_ETS1", force = T)
  
  proj <- addDeviationsMatrix(proj, "vierstra_motif_nooverlap_ets", force = T)
  
  proj <- addDeviationsMatrix(proj, "vierstra_cooccur_motif_with_ets", force = T)
  
  proj <- addDeviationsMatrix(proj, "vierstra_motif_ets_nooverlap_candidate", force = T)
  
  saveArchRProject(ArchRProj = proj, overwrite = T)
  
  
  # 3. plot deviations
  # this combines chromvar deviations
  plotCooccurVarDev <- getVarDeviations(ArchRProj = proj, name = "vierstra_cooccur_motifMatrix", plot = TRUE, n = 100)
  # this gets all motifs
  dtplotCooccurVarDev <- data.table(plotCooccurVarDev$data) %>% setnames(., 'name','motifID')
  png(filename = file.path(plot.dir,"default_chromvar_top_deviation_coccur.png"), width = 6, height = 10, units = 'in', res = 150)
  (PlotDeviations(dt = dtplotCooccurVarDev, n=50, text.adjust = T)) # Plots Top deviations
  dev.off()
  
  plotCoETS1VarDev <- getVarDeviations(ArchRProj = proj, name = "vierstra_cooccur_motif_with_ETS1Matrix", plot = TRUE, n = 100)
  # this gets all motifs
  dtplotCoETS1VarDev <- data.table(plotCoETS1VarDev$data) %>% setnames(., 'name','motifID')
  (PlotDeviations(dt = dtplotCoETS1VarDev, n=50, text.adjust = T)) # Plots Top deviations  
  
  plotCoetsVarDev <- getVarDeviations(ArchRProj = proj, name = "vierstra_cooccur_motif_with_etsMatrix", plot = TRUE, n = 100)
  # this gets all motifs
  dtplotCoetsVarDev <- data.table(plotCoetsVarDev$data) %>% setnames(., 'name','motifID')
  (PlotDeviations(dt = dtplotCoetsVarDev, n=50, text.adjust = T)) # Plots Top deviations
  
  ## what's available
  getAvailableMatrices(proj)
  
  
  ## Motifs identified with chromvar in Cluster3 that correlate with gene score
  ms <- setDT(readxl::read_xlsx(path = file.path(spreadsheets.dir,'Cluster3_Motif_Heatmap_Vierstra_Filt_Spearman.xlsx')))
  sigmotifs <- unique(ms[plot==T, motifID])
  grep("\\+",sigmotifs, value=T)
  sigmotifs <- gsub(".+\\+","",sigmotifs)
  
  #### Load Co-occurrence matrix of Chromvar scores and perform tests ####
  assay <- getMatrixFromProject(ArchRProj = proj, useMatrix = 'VierstraMatrix', binarize = F, threads = 8)
  
  coetsassay <- getMatrixFromProject(ArchRProj = proj, useMatrix = 'vierstra_cooccur_motif_with_etsMatrix', binarize = F, threads = 8)
  
  nonetsassay <- getMatrixFromProject(ArchRProj = proj, useMatrix = 'vierstra_motif_nooverlap_etsMatrix', binarize = F, threads = 8)
  
  etsnocandidateassay <- getMatrixFromProject(ArchRProj = proj, useMatrix = 'vierstra_motif_ets_nooverlap_candidateMatrix', binarize = F, threads = 8)
  
  z <- assays(assay[rownames(assay),])$z
  
  zn <- assays(nonetsassay[rownames(nonetsassay),])$z
  
  znc <- assays(etsnocandidateassay[rownames(etsnocandidateassay),])$z
  
  ze <- assays(coetsassay[rownames(coetsassay),])$z
  
  znames <- (grep(paste(sigmotifs, collapse = "|"), rownames(z), ignore.case = T, value = T))
  
  ze <- rbind(z[znames,], zn, znc, ze) # individual, with ETS, without ETS
  short.names <- strsplit(rownames(ze),"_and_") %>% lapply(.,function(i) gsub("#.*","",i)) %>% lapply(., function(i) paste(i,collapse = ' & ')) %>% unlist
  (rownames(ze) <- short.names)
  
  #### Fig. 6B ####
  smotifId <- "ETS_motifs"
  g0 <- UmapPlotEnhance(ArchRProj = proj,
                        p1 = NULL, 
                        feature = paste0('z:',smotifId),
                        colorBy = 'vierstra_motif_ets_nooverlap_candidateMatrix',
                        pt.size = 1,
                        outline.size=.1,
                        addLabels = F, 
                        embedding = 'UMAP',
                        title = "Erg/Ets motif", 
                        label.size = 14, 
                        title.size = 18, 
                        break.unit = NA, 
                        panel_width_in = 3,
                        scale_color_bar = T,
                        pal = solarExtraWhite_palette,
                        legend.title = 'Motif deviation zscore', 
                        legend.size = 14, 
                        legend.height.factor = .6,
                        filename.prefix = paste0('Umap_zDeviation_',smotifId))
 
  # use sigmotifs
  x <- c("ETS_motifs",sigmotifs)
  idx <- sapply(gsub("#.*","", x), function(i){
    grep( i, rownames(ze) )
  }) %>% unlist
  rownames(ze[idx,])
  
  dt <- as.data.table(data.frame(proj@cellColData), keep.rownames = "cells")  
  dtplot <- as.data.table(as.matrix(ze[idx,]), keep.rownames = 'motifs') %>% 
    melt(., id.vars='motifs', value.name = "zscores", variable.name = "cells")
  dtplot <- merge(dtplot, dt[,.(cells, Cluster3)], by='cells')
  dtplot[,motifs := factor(x = motifs, levels = unique(motifs), ordered = T)]
  dtplot[,quantile(x = zscores, c(.5,.95)), by=list(motifs, Cluster3)]
  dtplot[,motif_base := gsub(" .*|_no.*","", motifs)]
  dtplot[grep("^ETS_no_", motifs), motif_base := gsub("^ETS_no_","", motifs)]
  dtplot[,z_mean:=mean(zscores), by=list(Cluster3,motifs)]
  
  #### Wilcoxon single sample signed rank test ####
  # Wilcoxon 1 sample signed rank test of the null that distribution of x is symmetric about mu (0 default)
  the.motifs <- as.character(unique(dtplot$motifs))
  dtwilcox <- dtplot[, .(wilcox.pval = wilcox.test(zscores, alternative='greater')$p.value), by=list(Cluster3, motifs)]
  dtwilcox[,wilcox.pval.adjust := p.adjust(wilcox.pval)]
  dtwilcox[,type:="one-sample"]
  
  # t-test
  the.motifs <- as.character(unique(dtplot$motifs))
  dttest <- dtplot[, .(ttest.pval = t.test(zscores, alternative='greater')$p.value), by=list(Cluster3, motifs)]
  dttest[,ttest.pval.pval.adjust := p.adjust(ttest.pval)]
  dttest[,type:="one-sample"]
  
  ## combine tests
  dtstats <- na.omit(Reduce(function(x, y) x[y, on = .(Cluster3,motifs,type)], list(dtwilcox, dttest) ))
  unique(dtstats$motifs)
  dtstats[,motif_base := gsub(" .*|_no.*","", motifs)]
  dtstats[grep("^ETS_no_", motifs), motif_base := gsub("^ETS_no_","", motifs)]

  #### Wilcoxon 2 sample rank sum (Mann-Whitney) test ####
  dtwilcox_pw <- lapply(as.character(unique(dtplot$motif_base)), function(i){
    lapply(unique(dtplot$Cluster3), function(j){
      dt_test_base <- dtplot[Cluster3==j & motif_base==i][grep("&|_no_",motifs)]
      x <- dt_test_base[,unique(grep("&",motifs, value=T))]
      
      dtx <- dt_test_base[unique(grep("&",motifs))]
      ys <- dt_test_base[,unique(grep("_no_",motifs, value=T)) ]
      
      dt_res <- lapply(ys, function(k){
        dty <- dt_test_base[motifs==k]
        wilcox.pval <- wilcox.test(x = dtx$zscores, y = dty$zscores, alternative='greater')$p.value
        data.table(motif_base=i, Cluster3=j, group1 = x, group2 = k, wilcox.pval)
      }) %>% rbindlist()
    }) %>% rbindlist()
  }) %>% rbindlist()
  dtwilcox_pw[,wilcox.pval.adjust := p.adjust(wilcox.pval)]
  dtwilcox_pw[motif_base=='NFAC1_1145']
  
  #### Fig. 6F ####
  ## significant candidates above baseline (i.e. mean 0)
  k <- "EPC (IM)"
  candidates <- dtstats[grep("&|_no_", motifs, invert = T)][Cluster3==k][wilcox.pval.adjust < 0.05, motif_base]
  candidates <- intersect(candidates , c("NFKB2_1188","NFAC1_1145","STAT3_1819"))

  library(ggpubr) # add statistics within plot
  pmotif_base <- as.character(unique(dtplot[motif_base %in% candidates]$motif_base))
  gp <- lapply(pmotif_base, function(m){
    stat_table <- dtwilcox_pw[Cluster3==k & (motif_base==m | motif_base=="ETS_motifs")]
    selection <- unique(c(stat_table$group1, stat_table$group2))
    dt_plot_table <- dtplot[Cluster3==k & (motif_base==m | motif_base=="ETS_motifs"),list(motifs, Cluster3, zscores, motif_base)]
    dt_plot_table <- dt_plot_table[motifs %in% selection]
    dtylim <- dt_plot_table[,.(ymin=min(zscores), ymax=max(zscores)), by=motifs]
    stat_table[,wilcox.pval.adjust.label := formatC(wilcox.pval.adjust, format="G", digits=2)]
    p_meds <- dt_plot_table[,.(med=round(median(zscores), digits=2)), by=list(motifs)]
    p_meds <- p_meds[order(med, decreasing = T)]
    dt_plot_table[,motifs:=factor(motifs, levels = p_meds$motifs, ordered = T)]
    dt_plot_table[,cfill:="indianred2"]
    dt_plot_table[grep("&",motifs), cfill:="gray80"]
    
    sp <- ggplot(dt_plot_table, aes(x=motifs, y=zscores)) +
      geom_hline(yintercept = 0, linewidth=.5, linetype='dashed', color='gray50') +
      geom_violin(linewidth=1, scale = "width", aes(fill=cfill)) +
      geom_boxplot(color='black', linewidth=.5, fill='white', width=.5, outlier.shape = NA) +
      geom_text(data = p_meds, aes(x = motifs, y = med, label = med), size = 4, vjust = -.2) + # label medians
      ggtitle(paste0(m," in ",k)) + 
      ylab("") + xlab("") +
      scale_y_continuous(limits = c(floor(min(dtylim$ymin))-1, ceiling(max(dtylim$ymax)) + c(1))) +
      scale_color_manual("celltypes",values = "black") + 
      scale_fill_identity("") + 
      theme_bw(14) + 
      theme(strip.background = element_blank(), 
            plot.title = element_text(size=12),
            plot.margin = unit(x = c(.1,1,0,0), units = "cm"), #trbl
            panel.grid.minor = element_blank(), 
            panel.grid.major = element_line(linewidth = 0.2), 
            axis.text.x = element_text(size=12, angle = -45, hjust=0))
    sp + stat_pvalue_manual(data = stat_table, 
                            x = "group2", y.position = dtylim[motifs %in% stat_table$group2, ymax],
                            label = "wilcox.pval.adjust.label", 
                            size = 3,
                            color = "black")
  }) %>% set_names(pmotif_base)
  
  base.size <- c(3,4) # individual plot size
  ncol <- 3
  nrow <- 1
  do.call(grid.arrange, c(gp, ncol=ncol))
  
}

#### Helpers ####
absmax <- function(x) { x[which(abs(x)==max(abs(x)))][1] }

RebuildProj <- function(proj){
  
  message("After subset project must again add Iterative LSI and perform clustering")
  message("Consider the use of scaleDims")
  proj <- addIterativeLSI(ArchRProj = proj,
                          useMatrix = "TileMatrix",
                          force = T,
                          corCutOff = 0.75, # removes dimensions having more that this correlation with sequencing depth
                          iterations = 3,
                          clusterParams = list(resolution = c(0.2), sampleCells = 10000, n.start=10), # See Seurate:FindClusters, uses KNN->SNN
                          varFeatures = 25000,
                          scaleDims = F,  # Set to False so that we remove highly correlated dimensions
                          dimsToUse = 1:30)
  
  scaleDimsFalseClustering <- EvaluateUMAPClustering(ArchRProj = proj,
                                                     ncores = 6,
                                                     file.prefix = 'scaleDimsFalse',
                                                     scaleClusterDims = F,
                                                     resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5,
                                                                     0.6, 0.7, 0.8, 0.9, 1.0,
                                                                     1.1, 1.2, 1.3, 1.4, 1.5,
                                                                     1.6, 1.7, 1.8, 1.9, 2))
  res <- scaleDimsFalseClustering$rlist
  
  ## Speed up silhouette calculation by calculating distance matrix once
  sample.ratio <- 1
  ncells <- nrow(res[[1]]@cellColData)
  num_sample <- ceiling(nrow(res[[1]]@cellColData) * sample.ratio)
  sample_cells <- sample(x = res[[1]]$cellNames, size = num_sample)
  message("Sampled ",100*sample.ratio,"% cells: ",num_sample," out of ",ncells," for silhouette calculation")
  
  library(cluster)
  library(distances)
  
  d <- distances(res[[1]]@reducedDims$IterativeLSI$matSVD[sample_cells,])
  dist.matrix <- distances::distance_matrix(d)
  
  sils <- lapply(res, function(r){
    sil <- silhouetteScore(ArchrProj = r, dist.matrix = dist.matrix)
    #summary(sil)$avg.width
  }) %>% set_names( names(res) )
  
  um <- lapply(names(scaleDimsFalseClustering$rplots), function(r){
    score <- summary(sils[[r]])$avg.width
    g <- scaleDimsFalseClustering$rplots[[r]] + 
      ggtitle(paste("Cluster res.:", r,"\nmean silhouette ", formatC(score, digits = 3)))
  }) %>% set_names(names(scaleDimsFalseClustering$rplots))
  
  png(filename = file.path(plot.dir, 'silhouette_scores.png'), width = 30, height = 20, units = 'in', res=150, bg = 'transparent')
  do.call(grid.arrange, c(um, nrow=4))
  dev.off()
  
  pdf(file = file.path(plot.dir, 'silhouette_scores.pdf'), width = 10, height = 8, onefile=T)
  for(r in names(res)){
    g <- um[[r]]
    p1 <- ~plot(sils[[r]], main=paste("resolution:",r), cex.names = par("cex.axis")) # RStudio sometimes does not display silhouette plots correctly
    grid.arrange(cowplot::as_grob(p1), g, nrow=1, widths=c(.35,.65))
  }
  dev.off()
  
  rm(dist.matrix)
  rm(d)
  gc()
  
  #### Preferred resolution ####
  r <- scaleDimsFalseClustering$rlist[['0.7']]
  
  # combine C1 and C10
  r@cellColData[r@cellColData$Clusters=='C10',]$Clusters <- 'C1'
  r@cellColData[r@cellColData$Clusters=='C11',]$Clusters <- 'C8'
  r@cellColData[r@cellColData$Clusters=='C15',]$Clusters <- 'C12'
  
  ndesc <- sort(table(r$Clusters), decreasing = T)
  # Re-Number by #cells
  for(i in seq_along(ndesc)){
    # relabel clusters
    r@cellColData[r@cellColData$Clusters==names(ndesc)[i],]$Clusters <- paste0("X",i)
  }
  r@cellColData$Clusters <- gsub("X","C",r@cellColData$Clusters)
  table(r$Clusters)
  
  #### Save project ####
  proj <- saveArchRProject(ArchRProj = r, load = T)
  
  #### Plot final UMAPS ####
  table(proj@cellColData$bioNames)
  #BioNames
  plabels <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", 
                           name = "Clusters", embedding = "UMAP") %>% 
    UmapPlotEnhance(ArchRProj = proj, 
                    p1 = ., 
                    addLabels = T, 
                    title.size = 16,
                    pal = NA, 
                    label.size = 6)
  
  pBio <- plotEmbedding(ArchRProj = proj, 
                        colorBy = "cellColData", 
                        name = "bioNames", 
                        embedding = "UMAP", 
                        size=0.2, baseSize = 14) %>% 
    UmapPlotEnhance(ArchRProj = proj, p1 = ., addLabels = F, 
                    label.size = 4, 
                    pal=bionames_palette_list, 
                    title.size = 18, 
                    panel.border.size = 0,
                    filename.prefix = "UMAP-ATAC-BioNames")
  
  png(filename = file.path(plot.dir,'UMAP-ATAC-BioNamesClusters.png'), width = 5.5, height = 3.5, units = 'in', bg = 'transparent', res = 250)
  pBio+plabels$layers
  dev.off()
  
  #Cluster #
  plotEmbedding(ArchRProj = proj, colorBy = "cellColData", 
                name = "Clusters", embedding = "UMAP") %>% 
    UmapPlotEnhance(ArchRProj = proj, 
                    p1 = ., 
                    addLabels = T, 
                    title.size = 16,
                    pal = cluster_palette_list, 
                    label.size = 6, 
                    order.random = T,
                    panel.border.size = 0,
                    filename.prefix = "UMAP-ATAC-ClusterNumber")
  
  #nFrag
  dtl <- list(data=data.table(proj@embeddings$UMAP$df, log10(proj@cellColData$nFrags), keep.rownames = 'cells') %>% setnames(.,c('cells',"x","y","color")))
  p <- UmapPlotEnhance(ArchRProj = proj,
                       feature = "nFrags", 
                       break.unit = 1,
                       p1 = dtl, 
                       addLabels = F, 
                       panel.border.size = 0,
                       title.size = 16, 
                       scale_color_bar = F, 
                       legend.title = "log10(nFrags)", 
                       filename.prefix = "UMAP-ATAC-log10nFRAG")
  
  p <- UmapPlotEnhance(ArchRProj = proj,
                       feature = "nFrags", 
                       break.unit = 10000,
                       #p1 = dtl, 
                       addLabels = F, 
                       panel.border.size = 0,
                       title.size = 16, 
                       scale_color_bar = F, 
                       legend.title = "nFrags", 
                       filename.prefix = "UMAP-ATAC-nFRAG")
  
  
  #### Plot bar abundance ####
  PlotBarVariableInClusters(proj, variable = 'bioNames', 
                            palette = bionames_palette_list, orientation = 'vertical', 
                            filename = "bionames_bar", legend.ncol = 2, 
                            height_multiplier = 0.8, width_multiplier = 0.9)
  
  PlotBarVariableInClusters(proj, variable = 'Condition', palette = condition_palette, 
                            orientation = 'vertical', filename = "condition_bar", legend.ncol = 3, 
                            height_multiplier = 0.8, width_multiplier = 0.9)
  
  PlotBarVariableInClusters(proj, variable = 'Cluster3', yvariable = "Condition", palette = cluster3_palette, 
                            orientation = 'vertical', filename = "cluster3_bar", legend.ncol = 3, xlabel = "Condition",
                            height_multiplier = 0.5, width_multiplier = 1.5)
  
  PlotBarVariableInClusters(proj, variable = 'Cluster3', yvariable = "Clusters", palette = cluster3_palette, 
                            orientation = 'vertical', filename = "cluster3_bar", legend.ncol = 3, xlabel = "Condition",
                            height_multiplier = 0.5, width_multiplier = 1.5)
  
  pdf(file = file.path(plot.dir,'DepthDimensionCor.pdf'), width = 8, height = 4)
  PlotCorDimensionDepth(proj)
  dev.off()
  
  
  #### Confusion Matrix Sample and Cluster ####
  table(proj$Clusters)
  cM <- confusionMatrix(paste0(proj$Clusters), paste0(proj$bioNames))
  cM
  library(pheatmap)
  # apply(cM, 1, function(x){ x/sum(x)})
  cM <- t(cM) / rowSums(t(cM))
  cM
  p <- pheatmap::pheatmap(
    mat = as.matrix(cM),
    color = ArchRPalettes$solarExtra, #
    clustering_method = 'ward.D',
    treeheight_row = 20,
    treeheight_col = 20,
    cellwidth = 20,
    cellheight =20,
    #color = RColorBrewer::brewer.pal(n=8, name = "YlOrRd"), # paletteContinuous(set = "Set1"),
    border_color = "black",
    number_color = 'white',
    display_numbers = T,
    fontsize_row = 14,
    fontsize_col = 14,
    number_format="%.2f"
  )
  
  plotPDF(p, name = "confusion-matrix-wardmethod-cluster-sample.pdf", ArchRProj = proj, addDOC = FALSE, width = 8, height = 4)
  
  proj@cellColData$bioNames <- factor(x = proj$bioNames, levels = bionames_order, ordered = T)
  #### Save object ####
  proj <- saveArchRProject(ArchRProj = proj, load = T)
  
}

addClustersIgraph <- function (input = NULL, reducedDims = "IterativeLSI", name = "Clusters", 
                               sampleCells = NULL, seed = 1, method = "Seurat", dimsToUse = NULL, 
                               scaleDims = NULL, corCutOff = 0.75, knnAssign = 10, nOutlier = 5, 
                               maxClusters = 25, testBias = TRUE, filterBias = FALSE, biasClusters = 0.01, 
                               biasCol = "nFrags", biasVals = NULL, biasQuantiles = c(0.05, 
                                                                                      0.95), biasEnrich = 10, biasProportion = 0.5, biasPval = 0.05, 
                               nPerm = 500, prefix = "C", ArchRProj = NULL, verbose = TRUE, 
                               tstart = NULL, force = FALSE, logFile = createLogFile("addClusters"), 
                               ...) 
{
  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj", "null"))
  if (is(ArchRProj, "ArchRProject")) {
    message("When running addClusters 'input' param should be used for 'ArchRProj'. Replacing 'input' param with user 'ArchRPRoj'...")
    input <- ArchRProj
    rm(ArchRProj)
    gc()
  }
  ArchR:::.validInput(input = input, name = "input", valid = c("ArchRProj", "matrix"))
  ArchR:::.validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  ArchR:::.validInput(input = name, name = "name", valid = c("character"))
  ArchR:::.validInput(input = sampleCells, name = "sampleCells", valid = c("integer", "null"))
  ArchR:::.validInput(input = seed, name = "seed", valid = c("integer"))
  ArchR:::.validInput(input = method, name = "method", valid = c("character"))
  ArchR:::.validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  ArchR:::.validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  ArchR:::.validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  ArchR:::.validInput(input = knnAssign, name = "knnAssign", valid = c("integer"))
  ArchR:::.validInput(input = nOutlier, name = "nOutlier", valid = c("integer"))
  ArchR:::.validInput(input = testBias, name = "testBias", valid = c("boolean"))
  ArchR:::.validInput(input = filterBias, name = "filterBias", valid = c("boolean"))
  ArchR:::.validInput(input = biasClusters, name = "biasClusters", valid = c("numeric"))
  ArchR:::.validInput(input = biasCol, name = "biasCol", valid = c("character"))
  ArchR:::.validInput(input = biasQuantiles, name = "biasQuantiles", valid = c("numeric"))
  ArchR:::.validInput(input = biasEnrich, name = "biasEnrich", valid = c("numeric"))
  ArchR:::.validInput(input = biasProportion, name = "biasProportion", valid = c("numeric"))
  ArchR:::.validInput(input = biasPval, name = "biasPval", valid = c("numeric"))
  ArchR:::.validInput(input = nPerm, name = "nPerm", valid = c("integer"))
  ArchR:::.validInput(input = prefix, name = "prefix", valid = c("character"))
  ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
  ArchR:::.validInput(input = tstart, name = "tstart", valid = c("timestamp", "null"))
  ArchR:::.validInput(input = force, name = "force", valid = c("boolean"))
  ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))
  ArchR:::.startLogging(logFile = logFile)
  ArchR:::.logThis(append(args, mget(names(formals()), sys.frame(sys.nframe()))), 
                   "addClusters Input-Parameters", logFile = logFile)
  if (is.null(tstart)) {
    tstart <- Sys.time()
  }
  if (inherits(input, "ArchRProject")) {
    input <- addCellColData(ArchRProj = input, data = rep(NA, 
                                                          nCells(input)), name = name, cells = getCellNames(input), 
                            force = force)
    if (reducedDims %ni% names(input@reducedDims)) {
      stop("Error reducedDims not available!")
    }
    matDR <- getReducedDims(ArchRProj = input, reducedDims = reducedDims, 
                            dimsToUse = dimsToUse, corCutOff = corCutOff, scaleDims = scaleDims)
  }
  else if (inherits(input, "matrix")) {
    matDR <- input
  }
  else {
    stop("Input an ArchRProject or Cell by Reduced Dims Matrix!")
  }
  set.seed(seed)
  nr <- nrow(matDR)
  if (!is.null(sampleCells)) {
    if (sampleCells < nrow(matDR)) {
      ArchR:::.logDiffTime("Estimating Clusters by Sampling", tstart, 
                           verbose = verbose, logFile = logFile)
      estimatingClusters <- 1
      idx <- sample(seq_len(nrow(matDR)), sampleCells)
      matDRAll <- matDR
      matDR <- matDR[idx, , drop = FALSE]
    }
    else {
      estimatingClusters <- 0
    }
  }
  else {
    estimatingClusters <- 0
  }
  if (grepl("seurat", tolower(method))) {
  }
  else if (grepl("scran", tolower(method))) {
  }
  else {
    stop("Clustering Method Not Recognized!")
  }
  clust <- tryCatch({
    if (grepl("seurat", tolower(method))) {
      clustParams <- list(...)
      clustParams$verbose <- verbose
      clustParams$tstart <- tstart
      clustParams$method <- 'igraph'
      message("...Using igraph...")
      clust <- ArchR:::.clustSeurat(mat = matDR, clustParams = clustParams, 
                                    logFile = logFile)
    }
    else if (grepl("scran", tolower(method))) {
      clustParams <- list(...)
      clustParams$verbose <- verbose
      clustParams$tstart <- tstart
      clustParams$x <- t(matDR)
      clustParams$d <- ncol(matDR)
      clustParams$k <- ifelse(exists("...$k"), ...$k, 25)
      clust <- .clustScran(clustParams = clustParams, logFile = logFile)
    }
  }, error = function(e) {
    errorList <- clustParams
    ArchR:::.logError(e, fn = "runClusters", info = "", errorList = errorList, 
                      logFile = logFile)
  })
  if (estimatingClusters == 1) {
    ArchR:::.logDiffTime("Finding Nearest Clusters", tstart, verbose = verbose, 
                         logFile = logFile)
    knnAssigni <- as.matrix(.computeKNN(matDR, matDRAll[-idx, , drop = FALSE], knnAssign))
    clustUnique <- unique(clust)
    clustMatch <- match(clust, clustUnique)
    knnAssigni <- matrix(apply(knnAssigni, 2, function(x) clustMatch[x]), ncol = knnAssign)
    ArchR:::.logDiffTime("Assigning Nearest Clusters", tstart, verbose = verbose, logFile = logFile)
    clustAssign <- lapply(seq_along(clustUnique), function(x) {
      rowSums(knnAssigni == x)
    }) %>% Reduce("cbind", .) %>% apply(., 1, which.max)
    clustOld <- clust
    clust <- rep(NA, nr)
    clust[idx] <- clustOld
    clust[-idx] <- clustUnique[clustAssign]
    matDR <- matDRAll
    remove(matDRAll)
    gc()
  }
  if (testBias) {
    if (inherits(input, "ArchRProject")) {
      if (is.null(biasVals)) {
        biasDF <- getCellColData(input, select = biasCol)
      }
      else {
        biasDF <- DataFrame(row.names = rownames(matDR), 
                            bias = biasVals)
      }
    }
    else {
      if (!is.null(biasVals)) {
        biasDF <- DataFrame(row.names = rownames(matDR), 
                            bias = biasVals)
      }
      else {
        message("No biasVals for testing bias continuing without bias detection")
        testBias <- FALSE
      }
    }
  }
  if (testBias) {
    clust <- tryCatch({
      biasDF$Q <- ArchR:::.getQuantiles(biasDF[, 1])
      tabClust <- table(clust)
      tabClustP <- tabClust/sum(tabClust)
      idxTest <- which(tabClustP < biasClusters)
      names(clust) <- rownames(matDR)
      if (length(idxTest) > 0) {
        ArchR:::.logDiffTime("Testing Biased Clusters", tstart, 
                             verbose = verbose, logFile = logFile)
        testDF <- lapply(seq_along(idxTest), function(i) {
          clustTesti <- names(tabClustP)[idxTest[i]]
          biasQ <- biasDF[names(clust)[which(clust == 
                                               clustTesti)], 2]
          biasBgd <- matrix(sample(x = biasDF[names(clust)[which(clust != 
                                                                   clustTesti)], 2], size = nPerm * length(biasQ), 
                                   replace = if (nPerm * length(biasQ) > nrow(biasDF[names(clust)[which(clust != 
                                                                                                        clustTesti)], ])) 
                                     TRUE
                                   else FALSE), nrow = length(biasQ), ncol = nPerm)
          n1 <- colSums(biasBgd >= max(biasQuantiles))
          n2 <- colSums(biasBgd <= min(biasQuantiles))
          pval1 <- max(sum(sum(biasQ >= max(biasQuantiles)) < 
                             n1) * 2, 1)/length(n1)
          pval2 <- max(sum(sum(biasQ <= min(biasQuantiles)) < 
                             n2) * 2, 1)/length(n2)
          enrich1 <- sum(biasQ >= max(biasQuantiles))/max(median(n1), 
                                                          1)
          enrich2 <- sum(biasQ <= min(biasQuantiles))/max(median(n2), 
                                                          1)
          per1 <- sum(biasQ >= max(biasQuantiles))/length(biasQ)
          per2 <- sum(biasQ <= min(biasQuantiles))/length(biasQ)
          if (enrich1 > enrich2) {
            enrichClust <- enrich1
            enrichPval <- min(pval1, 1)
            enrichPer <- per1
          }
          else {
            enrichClust <- enrich2
            enrichPval <- min(pval2, 1)
            enrichPer <- per2
          }
          DataFrame(Cluster = clustTesti, enrichClust = enrichClust, 
                    enrichPval = enrichPval, enrichProportion = enrichPer)
        }) %>% Reduce("rbind", .)
        clustAssign <- testDF[which(testDF$enrichClust > 
                                      biasEnrich & testDF$enrichProportion > biasProportion & 
                                      testDF$enrichPval <= biasPval), 1]
        if (length(clustAssign) > 0) {
          if (filterBias) {
            ArchR:::.logDiffTime(sprintf("Assigning Biased Clusters (n = %s) to Neighbors", 
                                         length(clustAssign)), tstart, verbose = verbose, 
                                 logFile = logFile)
            for (i in seq_along(clustAssign)) {
              clusti <- clustAssign[i]
              idxi <- which(clust == clusti)
              knni <- ArchR:::.computeKNN(matDR[-idxi, , drop = FALSE], 
                                          matDR[idxi, , drop = FALSE], knnAssign)
              clustf <- unlist(lapply(seq_len(nrow(knni)), 
                                      function(x) names(sort(table(clust[-idxi][knni[x, 
                                      ]]), decreasing = TRUE)[1])))
              clust[idxi] <- clustf
            }
          }
          else {
            ArchR:::.logDiffTime(sprintf("Identified Biased Clusters (n = %s), set filterBias = TRUE to re-assign these cells: ", 
                                         length(clustAssign)), tstart, verbose = verbose, 
                                 logFile = logFile)
            message("Biased Clusters : ", appendLF = FALSE)
            for (i in seq_along(clustAssign)) {
              message(clustAssign[i], " ", appendLF = FALSE)
            }
            message("")
          }
        }
      }
      clust
    }, error = function(e) {
      errorList <- list(idxTest = if (exists("testDF", 
                                             inherits = FALSE)) fragx else "Error with idxTest!", 
                        biasDF = if (exists("testDF", inherits = FALSE)) fragx else "Error with biasDF!", 
                        testDF = if (exists("testDF", inherits = FALSE)) fragx else "Error with testDF!", 
                        clustAssign = if (exists("idf", inherits = FALSE)) fragx else "Error with clustAssign!")
      ArchR:::.logError(e, fn = "testBias", info = "", errorList = errorList, 
                        logFile = logFile)
    })
  }
  ArchR:::.logDiffTime("Testing Outlier Clusters", tstart, verbose = verbose, 
                       logFile = logFile)
  tabClust <- table(clust)
  clustAssign <- which(tabClust < nOutlier)
  if (length(clustAssign) > 0) {
    ArchR:::.logDiffTime(sprintf("Assigning Outlier Clusters (n = %s, nOutlier < %s cells) to Neighbors", 
                                 length(clustAssign), nOutlier), tstart, verbose = verbose, 
                         logFile = logFile)
    for (i in seq_along(clustAssign)) {
      clusti <- names(clustAssign[i])
      idxi <- which(clust == clusti)
      knni <- ArchR:::.computeKNN(matDR[-idxi, ], matDR[idxi, ], 
                                  knnAssign)
      clustf <- unlist(lapply(seq_len(nrow(knni)), function(x) names(sort(table(clust[-idxi][knni[x, 
      ]]), decreasing = TRUE)[1])))
      clust[idxi] <- clustf
    }
  }
  if (!is.null(maxClusters)) {
    if (length(unique(clust)) > maxClusters) {
      ArchR:::.logDiffTime(sprintf("Identified more clusters than maxClusters allowed (n = %s). Merging clusters to maxClusters (n = %s).\nIf this is not desired set maxClusters = NULL!", 
                                   length(clustAssign), maxClusters), tstart, verbose = verbose, 
                           logFile = logFile)
      meanDR <- t(ArchR:::.groupMeans(t(matDR), clust))
      hc <- hclust(dist(as.matrix(meanDR)))
      ct <- cutree(hc, maxClusters)
      clust <- mapLabels(labels = clust, oldLabels = names(ct), 
                         newLabels = paste0(prefix, ct))
    }
  }
  ArchR:::.logDiffTime(sprintf("Assigning Cluster Names to %s Clusters", 
                               length(unique(clust))), tstart, verbose = verbose, logFile = logFile)
  if (length(unique(clust)) > 1) {
    meanDR <- t(ArchR:::.groupMeans(t(matDR), clust))
    hc <- hclust(dist(as.matrix(meanDR)))
    out <- mapLabels(labels = clust, oldLabels = hc$labels[hc$order], 
                     newLabels = paste0(prefix, seq_along(hc$labels)))
  }
  else {
    out <- rep(paste0(prefix, "1"), length(clust))
  }
  if (inherits(input, "ArchRProject")) {
    input <- ArchR:::.suppressAll(addCellColData(input, data = out, 
                                                 name = name, cells = rownames(matDR), force = TRUE))
    ArchR:::.logDiffTime("Finished addClusters", t1 = tstart, verbose = verbose, 
                         logFile = logFile)
    return(input)
  }
  else if (!inherits(input, "ArchRProject")) {
    ArchR:::.logDiffTime("Finished addClusters", t1 = tstart, verbose = verbose, 
                         logFile = logFile)
    return(out)
  }
}


EvaluateUMAPClustering <- function(ArchRProj, file.prefix="clusters", minDist=0.5, spread=1, scaleClusterDims=F, ncores=8, 
                                   resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5, 
                                                   0.6, 0.7, 0.8, 0.9, 1.0,
                                                   1.1, 1.2, 1.3, 1.4, 1.5)) {
  library(cluster)
  
  dt <- data.table(ArchRProj@reducedDims$IterativeLSI$matSVD[,c(1,2,3)])
  dt[,ReadsInPromoter := ArchRProj@cellColData$ReadsInPromoter]
  dt[,labels := as.character(ArchRProj@cellColData$bioNames)]
  colors <-  sapply(dt$labels, function(x) { bionames_palette_list[[x]] })
  
  png(filename = file.path(plot.dir,paste0(file.prefix,'_LSI_Correlation.png')), 
      width = 6, height = 6, units = 'in', res = 200)
  plot(dt[,list(LSI1, LSI2, LSI3, ReadsInPromoter)], col=colors, pch=1, lwd=1, cex=.5)
  dev.off()
  
  png(filename = file.path(plot.dir,paste0(file.prefix,'_clusters_pca.png')),
      width = 6, height = 6, units = 'in', res = 200)
  plot(dt$LSI2, dt$LSI3, col=colors, pch=1, lwd=1, cex=.5)
  # Add a legend
  legend("bottomleft", inset=.02, title="Sample", bg = NULL,
         names(bionames_palette_list), col=unlist(bionames_palette_list), horiz=FALSE, cex=.8, pch=1)
  dev.off()
  # 
  # #### Batch correction ####
  # ArchRProj <- addHarmony(
  #   ArchRProj = ArchRProj,
  #   reducedDims = "IterativeLSI",
  #   name = "Harmony",
  #   groupBy = "library_batch"
  # )
  
  #### Clustering ####
  # quick clustering approach
  # Uses Seurat::FindClusters(), RAM hungry because "Seurat" casts data as dense 'matrix' rather than keeping it sparse
  # must use my own method for clustering, it uses Seurat sparse matrix, but it is slow
  rclust <- mclapply(resolutions, function(r){
    message("computing for resolution: ",r)
    addClustersIgraph(input = ArchRProj,
                      reducedDims = "IterativeLSI",
                      method = "Seurat", 
                      name = "Clusters", 
                      algorithm = 4, #leiden not louvain
                      resolution = r, 
                      ncores = 12,
                      scaleDims = scaleClusterDims, # Uses value from LSI ArchRProj@reducedDims$IterativeLSI$scaleDims
                      force = T)
  }, mc.cores = ncores ) %>% set_names(resolutions)
  
  if(!dir.exists(robj.dir))
    dir.create(robj.dir)
  saveRDS(object = rclust, file = file.path(robj.dir,'archrR_clustering_resolutions.RDS'))
  
  #rclust <- readRDS(file = file.path(robj.dir,'archrR_clustering_resolutions.RDS'))
  
  message("Adding umaps for resolutions")
  rlist <- mclapply(rclust, function(r){
    addUMAP(ArchRProj = r,
            force = T,
            reducedDims = "IterativeLSI",
            name = "UMAP",
            nNeighbors = 30, 
            seed = 42,
            minDist = minDist,
            spread = spread,
            metric = "cosine")
  }, mc.cores = ncores) %>% set_names( names(rclust) ) # efficiently operates on sparse matrix
  
  rplots <- lapply(names(rlist), function(r){
    p <- plotEmbedding(ArchRProj = rlist[[r]], 
                       colorBy = "cellColData", 
                       name = "Clusters", 
                       embedding = "UMAP", 
                       size=0.2, 
                       baseSize = 14)
    p <- UmapPlotEnhance(ArchRProj = rlist[[r]], 
                         p1 = p, 
                         addLabels = T, 
                         title.size = 16,
                         pal =  cluster_palette_list, 
                         legend.title = paste("resolution",r))
    
  })
  names(rplots) <- names(rlist)
  
  # spread changed from default(=1) to 0.6
  # minDist default(=0.5)
  # spread can therefore be used to control the inter-cluster distances to some extent,
  # where as min_dist controls the size of the clusters.
  # https://jlmelville.github.io/uwot/abparams.html
  psample <- plotEmbedding(ArchRProj = rlist[[1]], 
                           colorBy = "cellColData", 
                           name = "bioNames", 
                           embedding = "UMAP")
  psample <- UmapPlotEnhance(ArchRProj = rlist[[1]], 
                             p1 = psample, 
                             addLabels = F, 
                             title.size = 16,
                             label.size = 3, 
                             pal =  bionames_palette_list, 
                             legend.title = "", 
                             filename.prefix = "Umap_bioNames")
  psample
  
  png(filename = file.path(plot.dir,paste0(file.prefix,'_clusters_assess_spread_',spread,'.png')), 
      width = 30, height = 20, units = 'in', res = 200)
  do.call(grid.arrange, c(rplots, nrow=4))
  dev.off()
  
  ## Decide best resolution
  return(list('rlist'=rlist, 'rplots'=rplots))
  
}

#### Bug in filter doublets
filterDoubletsFix <- function(ArchRProj = NULL, cutEnrich = 1, cutScore = -Inf, filterRatio = 1, dryrun=FALSE){
  
  fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
  for (i in seq_along(fn)) {
    tryCatch({
      eval(parse(text = paste0(fn[i], "<-ArchR:::", fn[i])))
    }, error = function(x) {
    })
  }
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = cutEnrich, name = "cutEnrich", valid = c("numeric"))
  .validInput(input = cutScore, name = "cutScore", valid = c("numeric"))
  .validInput(input = filterRatio, name = "filterRatio", valid = c("numeric"))
  
  if(any(grepl("filterDoublets", names(ArchRProj@projectSummary)))){
    stop("Already ran filterDoublets on ArchRProject! Cannot be re-ran on an ArchRProject!")
  }
  
  df <- getCellColData(ArchRProj, c("Sample", "DoubletEnrichment", "DoubletScore"))
  splitDF <- split(seq_len(nrow(df)), as.character(df$Sample))
  
  cellsFilter <- lapply(splitDF, function(y){
    
    x <- df[y, ,drop = FALSE]
    
    n <- nrow(x)
    
    x <- x[order(x$DoubletEnrichment, decreasing = TRUE), ]
    
    if(!is.null(cutEnrich)){
      x <- x[which(x$DoubletEnrichment >= cutEnrich), ]
    }
    
    if(!is.null(cutScore)){
      x <- x[which(x$DoubletScore >= cutScore), ]
    }
    
    if(nrow(x) > 0){
      head(rownames(x), filterRatio * n * (n / 100000))
    }else{
      NULL
    }
    
  }) %>% unlist(use.names=FALSE)
  
  message("Total ", length(cellsFilter), " cells to be filtered from ArchRProject!")
  tabRemove <- table(df[cellsFilter,]$Sample)
  tabAll <- table(df$Sample)
  samples <- unique(df$Sample)
  for(i in seq_along(samples)){
    if(!is.na(tabRemove[samples[i]])){
      message("\t", samples[i], " : ", tabRemove[samples[i]], " of ", tabAll[samples[i]], " (", round(100 * tabRemove[samples[i]] / tabAll[samples[i]], 1),"%)")
    }else{
      message("\t", samples[i], " : ", 0, " of ", tabAll[samples[i]], " (0%)")
    }
  }
  
  if(length(cellsFilter) > 0){
    
    
    if(dryrun){
      message("Dry run, not filtering. Returning cells to filter")
      return(df[cellsFilter,])
    }
    
    
    else
      ArchRProj@cellColData <- ArchRProj@cellColData[rownames(ArchRProj@cellColData) %ni% cellsFilter,,drop=FALSE]
    
  }
  
  ArchRProj
  
}

## This code is in progress
silhouetteScore <- function(ArchrProj, dist.matrix=NULL, sample.ratio=0.5, reduction='IterativeLSI', clusters='Clusters'){
  library(cluster)
  library(distances)
  
  if(is.null(dist.matrix)){
    message("Calculating distance matrix")
    ncells <- nrow(ArchrProj@cellColData)
    num_sample <- ceiling(nrow(ArchrProj@cellColData) * sample.ratio) # subsample cells
    sample_cells <- sample(x = ArchrProj$cellNames, size = num_sample)
    message("Sampled ",100*sample.ratio,"% cells: ",num_sample," out of ",ncells," for distance matrix")
    
    d <- distances(ArchrProj@reducedDims[[reduction]]$matSVD[sample_cells,])
    dist.matrix <- distances::distance_matrix(d)
  } else{
    message("Using input distance matrix")
  }
  
  clusters <- gsub("^C","",ArchrProj@cellColData[sample_cells,clusters]) %>% as.numeric()
  clusters <- as.factor(clusters)
  message("Calculating silhouette information for ",length(levels(clusters))," clusters")
  sil <- silhouette(x = as.numeric(x = clusters), dist = dist.matrix)
}


geneRanges <- function(db, column="ENTREZID")
{
  g <- genes(db, columns=column)
  col <- mcols(g)[[column]]
  genes <- granges(g)[rep(seq_along(g), lengths(col))]
  mcols(genes)[[column]] <- as.character(unlist(col))
  genes
}

correlateAccessibilityExpression <- function(ArchRproj,
                                             dtInputChromVar,
                                             vassay,
                                             expassay,
                                             matGE,
                                             geneMatrix='GeneExpressionMatrix',
                                             min.dev.zscore=-Inf,
                                             min.expression=0,
                                             min.accessibililty=-Inf,
                                             groupby=NULL,
                                             imputeValues=F,
                                             ncores=16,
                                             build=T,
                                             cor.method='pearson',
                                             xicor.symmetric=T, # If xicor should we perform symmetric test?
                                             xlsx=NULL){
  stopifnot("cor.method must be one of spearman, pearson, xicor"=cor.method %in% c('spearman','pearson','xicor'))
  
  require(readxl)
  
  if(!build)
    return( setDT(read_xlsx(path = file.path(spreadsheets.dir, xlsx))) )
  
  message("#---------------------------------- Correlatating TF and Gene -----------------------------------------")
  
  require(XICOR)
  
  if(missing(vassay))
    vassay <- getMatrixFromProject(ArchRProj = ArchRproj, useMatrix = 'VierstraMatrix', binarize = F, threads = ncores)
  
  if(missing(expassay))
    expassay <- getMatrixFromProject(ArchRProj = ArchRproj, useMatrix = geneMatrix, binarize = F, threads = ncores)
  
  ## make function
  DT <- copy(dtInputChromVar)
  DT[,gene:=gsub(",.*","",gene)] # pick one from list , separated
  setkey(DT, idx)
  
  DTcellColData <- data.table(data.frame(ArchRproj@cellColData), keep.rownames = 'cell')
  
  dt <- NULL
  
  if(missing(matGE)){
    if(imputeValues){
      # Imputed values
      message("Getting imputed values")
      matGE <- imputeMatrix(assay(expassay), getImputeWeights(ArchRproj))
      rownames(matGE) <- rowData(expassay)$name
    } else {
      message("Using Unimputed expression matrix values")
      matGE <- expassay@assays@data[[geneMatrix]]
      rownames(matGE) <- rowData(expassay)$name
    }
  } else {
    message("Using provided gene expression matrix values")
  }
  
  
  #1. “deviations” - A deviation is a bias-corrected measurement of how far the per-cell accessibility of a given feature (i.e motif) 
  # deviates from the expected accessibility based on the average of all cells or samples.
  
  #2. “z-score” - The z-score, also known as a “deviation score” is the z-score for each bias-corrected 
  # deviation across all cells. The absolute value of the deviation score is correlated with the per-cell read depth. 
  # This is because, with more reads, you have higher confidence that the difference in per-cell accessibility of the given feature (i.e. motif) from the expectation is greater than would occur by chance.
  
  # deviations
  X <- vassay@assays@data$deviations[DT$motifID,]
  # z-score
  Z <- vassay@assays@data$z[DT$motifID,]
  # expression
  Y <- as.matrix(matGE[which(rownames(matGE) %in% DT$gene), ,drop=T])
  
  DT <- DT[DT$motifID %in% rownames(X)]
  DT <- DT[DT$gene %in% rownames(Y)]
  
  if( !is.null(groupby) ){
    if(!groupby %in% colnames(DTcellColData))
      stop("Error ",groupby," is not a valid value for grouping")
    message("Grouping correlations by ",groupby)
    
    # we will scan a motif/gene combination based on the motif 'family'
    # i.e. Tcf7l1, Tcfl2, and Tcf7 are all found within the TCF/LEF-HMG-29-0 family.
    # correlate each of these genes with each motif member in the family
    #dt <- foreach(i=DT$idx) %dopar% {
    families <- unique(DT$family)
    dt <- mclapply(families, function(f){
      #dt <- lapply(families, function(f){
      dtf <- DT[family == f]
      # for each motif
      dtm <- lapply( unique(dtf$motifID), function(m){
        # for each gene in Family
        dg <- lapply(unique( dtf$gene ), function(g){
          message("Correlating deviations of ",m," with gene ",g)
          motif.assigned.gene <- dtf[motifID==m, unique(gene)]
          
          x <- X[m,]
          x <- x[x > min.accessibililty]   # insist that cell accessibility score is above some threshold
          y <- Y[g,]
          y <- y[y > min.expression]       # insist that cell expresses gene above the min
          z <- Z[m,]
          
          
          # for each groupby
          result <- lapply( unique(DTcellColData[,get(groupby)]), function(j){
            sc <- DTcellColData[get(groupby) == j, cell]
            gn <- length(sc) # number of cells in group
            yn <- length(intersect(sc, names(y)))
            xn <- length(intersect(sc, names(x)))
            # cells with access. and express gene
            common <- intersect(sc, intersect(names(y), names(x)))
            sy <- y[which(names(y) %in% common)] # cell expression
            sx <- x[which(names(x) %in% common)] # cell deviations or accessibility
            sz <- z[which(names(z) %in% common)] # cell z-score deviations
            expr <- tryCatch({
              if(cor.method=='xicor'){
                # is y ~ f(x), i.e. is expression a function of accessibility?
                t1 <- xicor(x = sx, y = sy, pvalue = T, ties = T)
                names(t1) <- c("estimate","sd","p.value")
                out <- t1
                if(xicor.symmetric){
                  t2 <- xicor(x = sy, y = sx, pvalue = T, ties = T)
                  names(t2) <- c("estimate","sd","p.value")
                  if(t2$estimate > out$estimate)
                    out <- t2
                }
                out
              }
              else # correlating cell accessibility with expression
                cor.test(sx, sy, method = cor.method, exact=F)
            }, warning = function(warning_condition) {
              # message(warning_condition)
              list(estimate=expr$estimate, p.value=ifelse(is.na(expr$p.value), 1, expr$p.value))
            }, error = function(error_condition) {
              # message(error_condition)
              list(estimate=NA, p.value=1)
              #stop(error_condition)
            }, finally = {
              # cleanup-code
            })
            d <- data.table(j, 
                            family=f, 
                            motifID=m, 
                            motifGene=motif.assigned.gene, 
                            geneCorrelated=g,
                            coefficient=expr$estimate,
                            cor_pval=expr$p.value,
                            avg_accessibility=mean(sx),
                            avg_expression=mean(sy),
                            avg_zscore_deviation=mean(sz),
                            ncell_group=gn,
                            ncell_both=length(common),
                            ncell_motif=xn,
                            ncell_exp=yn,
                            frac_both=length(common)/gn,
                            frac_motif=xn/gn,
                            frac_exp=yn/gn
            )
            setnames(d, 'j', eval(groupby))
            d
          }) %>% rbindlist #groupby
        }) %>% rbindlist #gene
        dg[,cor_padj := p.adjust(cor_pval), by=geneCorrelated]
      }) %>% rbindlist #motif
    }, mc.cores = ncores) # family                
    #})
    dt <- rbindlist(dt)
    dt[,cor_mlog10Padj := -log(cor_padj+.Machine$double.xmin)/log(10)]
  } else {
    message("Calculating correlations for all entries")
    dt <- mclapply( DT$idx, function(i){
      f <- DT[i, family] %>% unique
      m <- DT[i, motifID] %>% unique
      g <- DT[i, gene] %>% unique
      gn <- length(DT$idx)
      # Imputed values
      # seGS <- getMatrixFromProject(ArchRProj)
      # matGS <- imputeMatrix(assay(seGS), getImputeWeights(ArchRProj))
      # Raw values
      #x <- vassay@assays@data$deviations[m,]
      #x <- vassay@assays@data$z[m,]
      #x <- x[which(x > min.dev.zscore)] # z-scores cutoff, otherwise may correlate with many 0s
      #x <- vassay@assays@data$z[m,]
      #y <- expassay@assays@data$GeneExpressionMatrix[which(rowData(expassay)$name==g),]
      #y <- matGE[g,]
      x <- X[m,]
      z <- Z[m,]
      y <- Y[g,]
      y <- y[y > min.expression]
      common <- intersect(names(y), names(x))
      sy <- y[names(y) %in% common]
      sx <- x[names(x) %in% common]
      sz <- x[names(z) %in% common]
      expr <- tryCatch({
        if(cor.method=='xicor'){
          out <- xicor(sx, sy, pvalue = T)
          names(out) <- c("estimate","sd","p.value")
          out
        }
        else # correlating cell accessibility with expression
          cor.test(sx, sy, method = cor.method, exact=F)
        
      }, warning = function(warning_condition) {
        # message(warning_condition)
        list(estimate=expr$estimate, p.value=ifelse(is.na(expr$p.value), 1, expr$p.value))
      }, error = function(error_condition) {
        # message(warning_condition)
        list(estimate=NA, p.value=1)
      }, finally={
        # cleanup-code
      })
      d <- data.table(family=f, 
                      motifID=m, 
                      geneCorrelated=g, 
                      coefficient=expr$estimate, 
                      cor_pval=expr$p.value,
                      avg_accessibility=mean(sx),
                      avg_expression=mean(sy),
                      avg_zscore_deviation=mean(sz),
                      ncell_both=length(common),
                      ncell_motif=length(x),
                      ncell_exp=length(y),
                      frac_both=length(common)/gn,
                      frac_motif=length(x)/gn,
                      frac_exp=length(y)/gn)
    }, mc.cores = ncores)
    dt <- rbindlist(dt)
    dt[,cor_padj := p.adjust(cor_pval), by=geneCorrelated]
    dt[,cor_mlog10Padj := -log(cor_padj+.Machine$double.xmin)/log(10)]
  }
  dt <- dt[order(coefficient, decreasing = T)]
  
  message("#---------------------------------- Finished correlations -----------------------------------------")
  
  if(!is.null(xlsx)){
    out.file <- file.path(spreadsheets.dir, xlsx)
    message("Writing to : ",out.file)
    writexl::write_xlsx(x = dt, path = out.file)
  }
  
  dt
}


Motif2GeneLookupDT <- function(proj, motifpwm, geneMatrix="GeneExpressionMatrix"){
  require(biomaRt)
  ### Gene relateing to these motifs NOT expressed ##
  dtmotif <- lapply(names(motifpwm), function(i){
    motif.names <- gsub("_.*","", toupper(i))
    x <- unlist(strsplit(motif.names, '\\+'), use.names = F)
    data.table(motifID=i, motif=x)
  }) %>% rbindlist()
  
  dtmotif <- dtmotif[order(motifID)]
  ## Update a few known ones
  dtmotif[grep("^ZN[0-9]+", motif, value = F),
          motif := gsub("ZN","ZNF", motif)]
  
  # replace NKX21 with NKX2-1
  dtmotif[grep("^NKX[0-9]+$",motif, value = F),
          motif := gsub('^(NKX[0-9]+)([0-9]+)$', '\\1-\\2', motif)]
  
  dtmotif[grep("^HX[A-C][0-9]+$", motif, value = F),
          motif := gsub('HX', 'HOX', motif)]
  
  dtmotif[grep("^PRD[0-9]+$", motif, value = F),
          motif := gsub('PRD', 'PRDM', motif)]
  
  dtmotif[grep('ZSCAN4', motif, value=F), motif := 'Zscan4d']
  
  # First assume motifs are human and translate human to mouse homologs
  source("~/lib/scripts/R/GenomeHelpers.R")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  h2m <- convertHumanGeneList2Mouse(x = unique(dtmotif$motif), human = human, mouse = mouse)
  
  # Create a look up table from motif to gene
  lookup <- data.table(motif=h2m$HGNC.symbol, gene=h2m$Gene.name)
  
  # Several genes were not found
  missing <- data.table(motif=setdiff(dtmotif$motif, h2m$HGNC.symbol), gene="")
  
  ## make manual adjustments to names
  # Try uniprot for lookup
  # https://www.uniprot.org/uniprot/?query=suh&sort=score
  missing[motif=='ANDR', gene:='Ar']
  missing[motif=='AP2A', gene:='Ap2a1']
  missing[motif=='AP2B', gene:='Ap2b1']
  missing[motif=='AP2C', gene:='Ap2c1']
  missing[motif=='TCFAP2A', gene:='Ap2a1']
  missing[motif=='TCFAP2B', gene:='Ap2b1']
  missing[motif=='TCFAP2C', gene:='Ap2c1']
  missing[motif=='ARI5B', gene:='Arid5b']
  missing[motif=='ATF6A', gene:='Atf6']
  missing[motif=='ATF7', gene:='Atf7']
  
  missing[motif=='BC11A', gene:='Bcl11a']
  missing[motif=='BHA15', gene:='Bhlha15']
  missing[motif=='BHLHB2', gene:='Bhlhe40']
  missing[motif=='BHLHB3', gene:='Bhlhe41']
  missing[motif=='BMAL1', gene:='Arntl']
  missing[motif=='BRAC', gene:='T']
  
  missing[motif=='CART1', gene:='Alx1']
  missing[motif=='COT1', gene:='Map3k8']
  missing[motif=='COT2', gene:='Nr2f2']
  
  missing[motif=='ERR1', gene:='Esrra']
  missing[motif=='ERR2', gene:='Esrrb']
  missing[motif=='ERR3', gene:='Esrrg']
  missing[motif=='EVI1', gene:='Mecom']
  missing[motif=='EWSR1-FLI1', gene:='Ewsr1']
  missing[motif=='FOXP1', gene:='Foxp1']
  
  missing[motif=='GCR', gene:='Nr3c1']
  missing[motif=='Hen1', gene:='Henmt1']
  missing[motif=='HINFP1', gene:='Hinfp']
  missing[motif=='HNF6', gene:='Onecut1']
  missing[motif=='HTF4', gene:='Tcf12']
  missing[motif=='ITF2', gene:='Tcf4']
  
  missing[motif=='MYBA', gene:='Mybl1']
  missing[motif=='NR1A4', gene:='Nr4a1']
  missing[motif=='NDF1', gene:='Neurod1']
  missing[motif=='NDF2', gene:='Neurod2']
  
  missing[motif=='NF2L1', gene:='Nfe2l1']
  missing[motif=='NF2L2', gene:='Nfe2l2']
  
  missing[motif=='NFAC1', gene:='Nfatc1']
  missing[motif=='NFAC2', gene:='Nfatc2']
  missing[motif=='NFAC3', gene:='Nfatc3']
  missing[motif=='NFAC4', gene:='Nfatc4']
  missing[motif=='NGN2', gene:='Neurog2']
  
  missing[motif=='P53', gene:='Trp53']
  missing[motif=='P63', gene:='Trp63']
  missing[motif=='P73', gene:='Trp73']
  missing[motif=='PO2F1', gene:='Pou2f1']
  missing[motif=='PO2F1', gene:='Pou2f1']
  missing[motif=='PO2F2', gene:='Pou2f2']
  missing[motif=='PO3F1', gene:='Pou3f1']
  missing[motif=='PO3F2', gene:='Pou3f2']
  missing[motif=='PO5F1', gene:='Pou5f1']
  missing[motif=='POU5F1P1', gene:='Pou5f1']
  
  missing[motif=='PRGR', gene:='Pgr']
  
  missing[motif=='RHOX11', gene:='Rhox11']
  missing[motif=='SMCA1', gene:='Smc1a']
  missing[motif=='SMCA5', gene:='Smc5a']
  missing[motif=='SRBP1', gene:='Srebf1']
  missing[motif=='SRBP2', gene:='Srebf2']
  missing[motif=='STA5A', gene:='Stat5a']
  missing[motif=='STA5B', gene:='Stat5b']
  missing[motif=='STF1', gene:='Pdx1']
  
  missing[motif=='T', gene:='T']
  missing[motif=='TAF1', gene:='Taf1']
  missing[motif=='TF2L1', gene:='Tfcp2l1']
  missing[motif=='TF7L1', gene:='Tcf7l1']
  missing[motif=='TF7L2', gene:='Tcf7l2']
  missing[motif=='TFE2', gene:='E2f1']
  missing[motif=='TRP63', gene:='Trp63']
  missing[motif=='TRP73', gene:='Trp73']
  missing[motif=='TWST1', gene:='Twist1']
  missing[motif=='TYY1', gene:='Yy1']
  missing[motif=='YY2', gene:='Mbtps2']
  missing[motif=='ZKSC1', gene:='Zkscan1']
  missing[motif=='ZSCAN4', gene:='Zscan4d']
  
  ##gmex <- getMatrixFromProject(ArchRProj = proj, useMatrix = "geneExpressionMatrix")
  gmex <- ArchR:::.getFeatureDF(ArrowFiles = getArrowFiles(proj), subGroup = geneMatrix)
  y <- sapply(missing[gene=="",motif], function(x){
    grep(paste0("^",x), gmex$name, ignore.case = T, value = T) %>% paste(., collapse = ',')
  })
  missing[gene==""]$gene <- y
  
  lookup <- rbindlist(list(lookup, missing))[order(motif)]
  lookup <- unique(lookup)
  
  expsed <- sapply(lookup$gene, function(i){
    x <- NULL
    if(grepl(",", i, fixed = T))
      x <- sapply( unlist(strsplit(i,',')), function(j){
        grep(paste0("^",j,"$"), gmex$name, ignore.case = T, value = T) %>% paste(., collapse = ',')
      })
    else
      x <- grep(paste0("^",i,"$"), gmex$name, ignore.case = T, value = T) %>% paste(., collapse = ',')
    any(nchar(x)) > 0
  })
  
  lookup[,expressed := unlist(expsed)]
  DT <- merge(dtmotif, lookup, on='motif', all.x=T, all.y=T)
}


# Given -log10pvalues will return combined statistic
fisherCpvalLog10 <- function(x){
  chi <- 2*sum(x*log(10))
  # (2*2.811076*log(10)) + (2*255.9578*log(10))
  df <- 2*length(x) # df = 2k
  p <- pchisq(q = chi, df = df, lower.tail = F)
  list(chisq=chi, df=df, p=p)
}


ArchRProjectToSeurat <- function(build=T,
                                 ArchRProj,
                                 matrix = 'PeakMatrix',
                                 output.file=file.path(robj.dir,'Erg_epi_no_sv_archR_atac_peaks_multiome.rds')){
  
  if(!build){
    return(readRDS(output.file))
  }
  
  library(SingleCellExperiment)
  library(Seurat)
  #getAvailableMatrices(ArchRProj)
  message("Using ",matrix)
  
  rse <- getMatrixFromProject(ArchRProj, matrix, threads = 8)
  
  if(matrix=='GeneScoreMatrix'){
    names(assays(rse)) <- "counts"
    assays(rse) <- list(counts = assays(rse)$counts,
                        logcounts = log1p(assays(rse)$counts) )
    rownames(rse) <- rowData(rse)$name
  }
  else if(matrix=='PeakMatrix'){
    names(assays(rse)) <- "counts"
    assays(rse) <- list(counts = assays(rse)$counts,
                        logcounts = log1p(assays(rse)$counts) )
    rownames(rse) <- paste0("Peak",seq_along(rse))
  }
  else if(matrix=='VierstraMatrix'){
    i <- names(assays(rse))
    i <- sapply(i, function(j){
      ifelse(j=="deviations","logcounts",ifelse(j=="z","counts",j)) })
    names(assays(rse)) <- i
    assays(rse) <- list(counts = assays(rse)$counts,
                        logcounts = assays(rse)$logcounts )
  }
  
  ## create single cell experiment
  sce <- as(rse, "SingleCellExperiment")
  
  
  ## create seurat object
  #obj <- as.Seurat(x = sce, data='data')
  obj <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))
  obj <- SetAssayData(object = obj, slot = "data", new.data = logcounts(sce))
  obj@assays
  obj <- RenameAssays(obj, RNA = 'ATAC')
  
  lsi <- ArchRProj@reducedDims$IterativeLSI$matSVD[,-1]
  umap <- as.matrix(ArchRProj@embeddings$UMAP$df)
  
  obj[['LSI']] <- Seurat::CreateDimReducObject(embeddings = lsi,
                                               key = "LSI_", 
                                               assay = "ATAC")
  obj[['UMAP']] <- Seurat::CreateDimReducObject(embeddings = umap,
                                                key = "UMAP_", 
                                                assay = "ATAC")
  
  f <- FindNeighbors(object = obj, k.param = 30, reduction = 'LSI', annoy.metric = 'euclidean')
  unique(rowSums(f@graphs$ATAC_nn))
  names(f@graphs)
  f@reductions
  
  if(length(output.file) > 0)
    if(!is.na(output.file)){
      message("Saving output to ",output.file)
      saveRDS(object = f, file=output.file) 
    }
  
  f
  
}


RelateExpressionData <- function(ArchRProj, RnaAverageExpression, groupby="Cluster3"){
  ## relate the groupby and predicted id
  
  mat <- table(ArchRProj@cellColData[[groupby]], ArchRProj$predicted.id)
  # divide each element by sum of column (i.e. this is fraction of cells in RNA cluster mapped that to this ATAC cluster)
  rnafrac <- apply(mat, 2, function(x){ x/sum(x) })
  # % of atac cluster that maps to this rna cluster
  atacweight <- t(apply(mat, 1, function(x){ 100 * x/sum(x) }))
  
  MClusterData <- as.matrix(RnaAverageExpression[,-1])[,colnames(mat)]
  rownames(MClusterData) <- RnaAverageExpression$gene
  
  groups <- unique(ArchRProj@cellColData[[groupby]])
  
  ## The average expression per group of all genes
  avgExpGroup <- sapply(groups, function(x){
    #mean(rowSums(MClusterData[] * atacweight[x,]))
    #rowMeans(MClusterData[] * atacweight[x,])
    mean(rowMeans(MClusterData[] * atacweight[x,]))
  }) 
  avgExpGroupByGene <- lapply(rownames(MClusterData), function(g){
    sapply(groups, function(x){ mean(MClusterData[g,] * atacweight[x,]) })
  })
  names(avgExpGroupByGene) <- rownames(MClusterData)
  
  dt <- data.table(do.call(rbind, avgExpGroupByGene), keep.rownames = 'gene')
  
  #dt[gene=='Stat1']
}

#### Plots ####
UmapPlotEnhance <- function(ArchRProj, 
                            p1=NULL,
                            useMatrix="VierstraMatrix",
                            colorBy = "cellColData",
                            feature,
                            embedding='UMAP',
                            pal=as.vector(ArchRPalettes$solarExtra),
                            panel.border.size=0,
                            na.color="gray80",
                            axis.label=NA,
                            addLabels=F, 
                            outline.size=.1,
                            outline.color='gray20', #'grey80',
                            label.size=3,
                            panel_width_in=3.5,
                            order.random=FALSE,
                            pt.size=1.2,
                            pt.alpha=1,
                            break.unit=1,
                            vmin=NULL,
                            title, 
                            title.size=16, 
                            show.legend=T, 
                            legend.title="", 
                            legend.size=11,
                            legend.title.size=10,
                            legend.title.rotation=90,
                            legend.height.factor=.4,
                            legend.position="right",
                            legend.justification="top",
                            legend.direction="vertical",
                            use.continuous=T,
                            scale_color_bar=T,
                            fixed_coord=T,
                            small_axis_title=T,
                            filename.prefix, 
                            nrow=1,
                            as.pdf=F){
  # Plot the integration UMAP
  require(shadowtext)
  require(ggrastr)
  require(ggh4x)
  library(gtable)
  
  if(panel_width_in <= 1)
    stop("panel_width_in must be larger than 1")
  
  panel_width <- unit(panel_width_in, "in")
  
  if(is.null(p1)){
    p1 <- plotEmbedding(ArchRProj = ArchRProj,
                        colorBy = colorBy,
                        name = feature,
                        embedding = embedding,
                        plotAs = 'points',
                        size = 0.5,
                        imputeWeights = NULL)
    if(use.continuous)
      p1$data$color <- as.numeric(p1$data$color)
    else
      p1$data$color <- as.character(p1$data$color)
  }
  if(is.list(p1)){
    if(!"data" %in% names(p1)){
      thiscall <- match.call(expand.dots = F)
      p <- lapply(names(p1), function(n){
        thiscall[["p1"]] <- p1[[n]]
        thiscall[["title"]] <- n
        thiscall[["filename.prefix"]] <- n
        eval.parent(thiscall)
      })
      return(p)
    }
  }
  if(is.data.frame(p1)){
    p1 <- list("data"=p1)
  }
  
  dt <- data.table(copy(p1$data), keep.rownames = 'rn')
  dt[,rn:=as.integer(rn)]
  dt <- dt[order(rn, decreasing = F)]
  
  ylabel <- gsub("Dimension ","",p1$labels$y)
  xlabel <- gsub("Dimension ","",p1$labels$x)
  
  ## This creates 0 valued points not in faceting window
  if(axis.label %in% colnames(ArchRProj@cellColData)){
    dt[,axis.label := as.character(ArchRProj@cellColData[, eval(axis.label)])]
    #dt[,axis.label := ArchRProj@cellColData[, eval(axis.label)]]
    dt <- lapply(unique(dt$axis.label), function(i){
      ds <- copy(dt[axis.label!=i])
      dt_zero <- copy(ds)
      dt_zero[,color:=as.integer(as.character(color))]
      dt_zero$color <- NA
      dt_zero$axis.label <- i
      rbindlist(list(dt_zero, ds))
    }) %>% rbindlist
  }
  
  if(order.random){
    dt$color <- gsub("[0-9]+-","", dt$color)
    dt <- dt[order(sample(color), na.last = F)] 
  }
  else if(!is.numeric(dt$color) & length(levels(dt$color)) > 0){
    lvls <- levels(dt$color)
    lvls <- gsub("[0-9]+-","",lvls)
    dt$color <- gsub("[0-9]+-","", dt$color)
    dt[,color:=factor(color, levels = lvls, ordered = T)]
    dt <- dt[order(color, na.last = F)]
    cat(lvls,"\n")
  }
  else {
    message("Sorting colors")
    dt <- dt[order(color, decreasing = F, na.last = F)]
  }
  
  if(!"label" %in% colnames(dt)){
    dt$label <- dt$color
  }
  
  dt$label <- gsub("[0-9]+-","", dt$label)
  median.labels <- dt
  median.labels <- median.labels[,.(color, x=median(x),y=median(y)), by=label]
  median.labels <- median.labels[order(label)]
  
  p1 <- NULL
  
  if(all(is.na(pal))){
    p1 <- ggplot(dt, aes(x=x, y=y)) + theme_ArchR()
    outline.size <- 0
  }
  else{
    p1 <- ggplot(dt, aes(x=x, y=y, fill=color)) + theme_ArchR()
  }
  
  if(outline.size > 0){
    stroke <- outline.size # deintensified stroke
    p1 <- p1 + geom_point_rast(size=pt.size, shape=1, stroke=outline.size, 
                               col=outline.color, alpha=pt.alpha) # standard stroke #'grey40'
  }
  
  if(small_axis_title){
    xmax <- dt[,max(x)]
    ymax <- dt[,max(y)]
    xmin <- dt[,min(x)]*1.05
    ymin <- dt[,min(y)]*1.05
    expansion <- .10
    xd <- (xmax-xmin) * expansion
    yd <- (ymax-ymin) * expansion
    
    d <- max(xd, yd)
    
    p1 <- p1 +
      geom_segment(x=xmin, y=ymin, xend=xmin+d, yend=ymin, linewidth=.5, arrow = arrow(length = unit(.2, "cm"), type = 'closed')) +
      geom_segment(x=xmin, y=ymin, xend=xmin, yend=ymin+d, linewidth=.5, arrow = arrow(length = unit(.2, "cm"), type = 'closed')) +
      theme(axis.title.x =  element_text(hjust = 0, vjust=1, size = 12), 
            axis.title.y = element_text(hjust = 0, vjust=1, size = 12))
    
  }
  
  p1 <- p1 + geom_point_rast(size=pt.size, shape=21, stroke=0, alpha=pt.alpha) +
    xlab(xlabel) + ylab(ylabel) +
    theme(text = element_text(size = 15),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(size=title.size, hjust = 0.5),
          legend.text = element_text(size = 12),
          panel.border = element_rect(linewidth = panel.border.size, color='black'),
          panel.spacing = unit(0, "lines"), # in case facet
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
  
  if('axis.label' %in% colnames(dt)){
    message("Facetting by axis.label")
    if(!is.null(nrow))
      p1 <- p1 + facet_wrap(~ axis.label, nrow = nrow)
    else
      p1 <- p1 + facet_wrap(~ axis.label)
  }
  
  if(addLabels)
    p1 <- p1 + geom_shadowtext(data = median.labels, aes(x=x, y=y, label=label), 
                               color='black', bg.color='white', bg.r = 0.05, size=label.size)
  
  ## Set dimensions of plot
  gt <- ggplotGrob(p1)
  panel <- gtable_filter(gt, "panel")
  # height is scaled compared to width
  panel_height_in <- panel_width_in*as.numeric(gsub("null","",panel$heights))
  # set the panel height based on width
  panel_height <- unit(panel_height_in, "in") # reduce for title
  
  
  if(is.numeric(dt$color)){
    message("continuous: limits")
    
    lim <- range(dt$color)
    message("Data limits: [",lim[1],":", lim[2],']')
    ## round the limit up to nearest 10th
    lim <- c(floor(1000*lim[1])/1000, ceiling(1000*lim[2])/1000)
    
    my.breaks <- NULL
    
    maxval <- ceiling(max(lim))
    numbreaks <- 4
    if((maxval %% 2) != 0) #odd number
      numbreaks <- 3  
    
    if(scale_color_bar){
      # lim <- c(-max(abs(lim)), max(abs(lim)))
      lim <- c(- ceiling(max(lim)), ceiling(max(lim))) # we cap colors from the high end
      break.unit <- ceiling(lim[2]/numbreaks)
      my.breaks <- seq(from=-break.unit*numbreaks, 
                       to=break.unit*numbreaks, 
                       by = break.unit)
    }
    
    if(!missing(pal) & !is.null(pal)){
      #names(pal) <- levels(dt$color)
      names(pal) <- levels(lim)
    } else {
      pal <- ArchRPalettes$solarExtra
    }
    
    message("Limits: [",lim[1],":", lim[2],']')
    if(is.null(break.unit) | is.na(break.unit)){
      break.unit <- round((lim[2] - lim[1]) / numbreaks, digits = 0)
    }
    if(lim[2] < 1 | break.unit < 1){
      break.unit <- 0.2
    }
    
    message("break unit:",break.unit)
    
    # if limit is smaller than break
    if(break.unit > lim[2]){ # & scale_color_bar){
      # center at 0
      my.breaks <- seq(from=floor(lim[1]/break.unit)*break.unit, 
                       to=ceiling(lim[2]/break.unit)*break.unit, 
                       by=break.unit)
    } else if(is.null(my.breaks)) {
      my.breaks <- seq(from=lim[1], 
                       to=lim[2], 
                       by = break.unit)
    }
    
    if(!is.null(vmin)){
      my.breaks[1] <- lim[1] <- vmin
    }
    
    message("Breaks: [",paste(my.breaks,collapse = ","),']')
    
    p1 <- p1 + scale_fill_gradientn(legend.title,
                                    guide = "colourbar", 
                                    colours = pal, 
                                    na.value = na.color,
                                    limits = lim, 
                                    breaks = my.breaks)
    
    p1 <- p1 + theme(legend.title = element_text(angle = 90, size=legend.title.size))
    
    p1 <- p1 + guides(fill = guide_colourbar(title.position = "left",
                                             ticks = T,
                                             ticks.linewidth = .5,
                                             ticks.colour = 'black',
                                             frame.colour = 'black',
                                             frame.linewidth = 0.4,
                                             nbin = 1000, # set high to put ticks at extreme ends
                                             direction='vertical',
                                             barwidth = .5,
                                             barheight = panel_height*legend.height.factor))
    panel_height_in <- panel_height_in+.5
  }
  #discrete
  else{
    message("discrete")
    v <- na.omit(unique(dt$color))
    n <- length(v)
    if(length(pal) != n)
      pal <- unname(pal[v])
    message("discrete color : ", n)
    if(missing(pal))
      pal <- pals::brewer.accent(n = n)
    else if(is.null(pal)){
      message("No colors defined")
    }
    else{
      #p1 <- p1 + scale_fill_manual(legend.title, values = unname(pal[v]), na.value = "gray80")
      p1 <- p1 + scale_fill_manual(legend.title, values = pal, na.value = "gray80")
      ncol.legend <- ifelse(n>10,ifelse(n>20,3,2),1)
      p1 <- p1 + guides(fill = guide_legend(override.aes = list(size=5), 
                                            ncol=ncol.legend, 
                                            title.position = ifelse(legend.title.rotation!=0,"left", "center")))
      p1 <- p1 + theme(legend.box.margin=margin(-10,10,-10,-10)) #trbl
      
      if(ncol.legend>1)
        panel_width_in <- panel_width_in+.5
    }
  }
  
  if(!show.legend){
    p1 <- p1 + theme(legend.position = 'none')
    panel_width_in <- panel_width_in+1
    panel_height_in <- panel_height_in+.5
  }
  else{
    panel_width_in <- panel_width_in+2.5
    panel_height_in <- panel_height_in+.5
  }
  if(legend.title!=""){
    p1 <- p1 + theme(legend.box.margin=margin(-10,10,-10,-10), legend.title = element_text(angle = legend.title.rotation)) #trbl
  }
  if(!missing(title)){
    p1 <- p1 + ggtitle(label = title) + theme(plot.title = element_text(size = title.size)) 
    panel_height_in <- panel_height_in+.5 # bit more
  }
  
  message("Dimensions: [",panel_width,",",panel_height,"]")
  
  p1 <- p1 + force_panelsizes(cols = panel_width, rows = panel_height)
  
  if(!missing(filename.prefix)){
    if(as.pdf){
      pdf(file = file.path(getOutputDirectory(ArchRProj),'Plots',paste0(filename.prefix,".pdf")), width = panel_width_in, height = panel_height_in, )
      print(p1)
      dev.off()
    } else{
      png(file = file.path(getOutputDirectory(ArchRProj),'Plots',paste0(filename.prefix,".png")), 
          width = panel_width_in, height = panel_height_in, units = 'in', bg = 'transparent', res = 300)
      print(p1)
      dev.off()
    }
  }
  p1
}

UmapPlotDensity <- function(ArchrProj,
                            feature='Phase',
                            title.size=16,
                            pt.size=0.8,
                            panel_width_in=3.5,
                            pal=pals::kovesi.linear_kryw_5_100_c64(10), #pal=pals::viridis(10), #pal=pals::parula(10)
                            na.color="gray85",
                            contour.color = 'blue',
                            filename.prefix,
                            fixed_coord=T,
                            embedding='UMAP',
                            show.axis=T,
                            coord1='IterativeLSI#UMAP_Dimension_1',
                            coord2='IterativeLSI#UMAP_Dimension_2',
                            legend.height.factor=1,
                            nrow=1,
                            as.pdf=F){
  
  
  .get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  
  dt <- merge( data.table( ArchrProj@embeddings[[embedding]][['df']], keep.rownames = 'cellid' ),
               data.table( as.data.frame(ArchrProj@cellColData), keep.rownames = 'cellid' ), by='cellid')
  dt[,temp_features := get(feature)]
  
  dt <- dt[!is.na(get(coord1))]
  
  require(ggh4x)
  require(gtable)
  
  panel_height_in <<- 0
  panel_height <<- 0
  panel_width <<- panel_width_in
  
  xmin <<- Inf
  xmax <<- -Inf
  ymin <<- Inf
  ymax <<- -Inf
  
  #fsort <- sort(unique(dt[,temp_features]))
  fsort <- rev(bionames_order)
  p <- lapply(fsort, function(s){
    # select data of interest
    toplot <- dt[temp_features == s]
    if(nrow(toplot) > 1){
      
      # increases space between plots but also allows for contours not being cutoff
      xrange <- range(dt[,get(coord1)])
      xrange[1] <- xrange[1]*1.1
      xrange[2] <- xrange[2]*1.4
      yrange <- range(dt[,get(coord2)])
      yrange[1] <- yrange[1]*1.1
      yrange[2] <- yrange[2]*1.18
      
      toplot[,density := .get_density(x = get(coord1), y= get(coord2), n=500)]
      toplot[,density := density/max(density)]
      
      g <- ggplot(toplot, aes(x=get(coord1), y=get(coord2), color=density)) +
        geom_point_rast(data = dt, size = pt.size, fill = na.color, shape=21, color='black', stroke=0.5) + # include all points and color gray
        geom_point_rast(data = dt, size = pt.size, color = na.color, shape=1) + # include all points and color gray
        geom_point_rast(size = pt.size, stroke=0.5) +
        scale_x_continuous(limits=xrange) +
        scale_y_continuous(limits=yrange) +
        scale_color_gradientn('', guide = "colourbar", colours = pal, 
                              limits=c(0,1),
                              breaks=seq(from=0, to=1, by=.2)) +
        #ggtitle(paste("Density", gsub("_"," ",s) )) + 
        ggtitle(gsub("_"," ",s)) + 
        xlab( gsub(".*#|_Dimension","",coord1) ) + ylab( gsub(".*#|_Dimension","",coord2) ) +
        theme_bw(12) +
        theme(panel.grid = element_blank(),
              plot.title = element_text(size=title.size, hjust = 0.5, vjust = -1),
              plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
              plot.margin = margin(0,0,0,0),
              strip.background = element_blank(),
              panel.border = element_blank(),
              panel.background = element_rect(fill = "transparent"), # bg of the panel
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              legend.background = element_blank(),
              legend.key.size = unit(16, units='points'),
              legend.title = element_text(angle = 90, hjust = 0.5),
              legend.text = element_text(size = 12),
              legend.position = 'right',
              legend.justification = "top",
              legend.margin=margin(0,0,0,0),
              legend.box.margin=margin(-10,10,-10,-10)) #trbl
      
      if(!is.null(contour.color))
        if(!is.na(contour.color))
          g <- g + geom_density_2d(data=toplot, linewidth=.2, color=contour.color) #pal[length(pal)])
      
      if(fixed_coord)
        g <- g + coord_fixed(clip = 'off')
      
      #  get dimensions of plot
      gt <- ggplotGrob(g)
      panel <- gtable_filter(gt, "panel")
      # height is scaled compared to width
      panel_height_in <<- panel_width_in*as.numeric(gsub("null","",panel$heights))
      # set the panel height based on width
      panel_height <<- max(panel_height, unit(panel_height_in, "in"))
      panel_width <<- max(panel_width, unique(panel_width_in, "in"))
      
      xmin <<- min(xmin, toplot[,min(get(coord1))])
      xmax <<- max(xmax, toplot[,max(get(coord1))])
      ymin <<- min(ymin, toplot[,min(get(coord2))])
      ymax <<- max(ymax, toplot[,max(get(coord2))])
      
      g <- g + guides(color = guide_colourbar(title.position = "left",
                                              ticks = T,
                                              frame.colour = 'black',
                                              frame.linewidth = 0.2,
                                              nbin = 1000, # high to set extreme ends
                                              direction='vertical',
                                              barwidth = .6,
                                              barheight = panel_height*legend.height.factor))
    }
  })
  
  p <- p[lengths(p) > 0]
  
  p <- lapply(p, function(g){
    
    if(show.axis){
      xd <- xmax-xmin
      yd <- ymax-ymin
      
      xminx <- xmin -(abs(xmin)*.05)
      yminy <- ymin -(abs(ymin)*.05)
      
      len <- (xminx+xd*.1) - xminx
      yend <- yminy + len
      
      g <- g +
        geom_segment(x=xminx, y=yminy, xend=xminx+xd*.1, yend=yminy, linewidth=.25,
                     arrow = arrow(length = unit(.15, "cm"), type = 'closed'), color='black') +
        geom_segment(x=xminx, y=yminy, xend=xminx, yend=yend, linewidth=.25,
                     arrow = arrow(length = unit(.15, "cm"), type = 'closed'), color='black') +
        theme(axis.title.x =  element_text(hjust = 0, size = 12),
              axis.title.y = element_text(hjust = 0, size = 12))
      
    } else {
      g <- g + xlab("") + ylab("")
    }
    
    g <- g + force_panelsizes(cols = panel_width, rows = panel_height)
  })
  
  
  panel_width_in <- (panel_width_in+.05) * (length(p)/nrow)
  panel_height_in <- (panel_height_in+.05) * nrow
  
  #panel_width_in <- (panel_width_in+1) * (length(p)/nrow)
  #panel_height_in <- (panel_height_in+.5) * nrow
  
  if(!missing(filename.prefix)){
    if(as.pdf){
      f <- file.path(plot.dir,paste0(filename.prefix,".pdf"))
      message("Plotting to:",f)
      pdf(file = f, width = panel_width_in, height = panel_height_in, bg = 'transparent')
      grid_arrange_list_shared_legend(p, nrow = nrow, legend.position = "topright")
      dev.off()
    } else{
      f <- file.path(plot.dir,paste0(filename.prefix,".png"))
      message("Plotting [",panel_width_in,"x",panel_height_in,"] to:",f)
      png(filename = f, width = panel_width_in, height = panel_height_in, bg = 'transparent', units = 'in', res=200)
      grid_arrange_list_shared_legend(p, nrow = nrow, legend.position = "topright")
      dev.off()
    }
  }
  
  p
}

DoHeatmapEnhance <- function(marker.peaks, 
                             marker.cutoff="!is.na(FDR)", 
                             filt.heatmap.peak = NULL,
                             label.features,
                             subset.features,
                             orientation.vertical = TRUE,
                             heatmap.colors = heatmap_palette,
                             sample.feature.ratio = 1,
                             plot.log2fc = FALSE,
                             title="Single cell Heatmap",
                             legend.title="Peak accessibility",
                             spreadsheets.path=NULL,
                             group.bar = TRUE, 
                             group.colors = NULL, 
                             show.group.bar.count = TRUE,
                             disp.min = -2.5,
                             show.feature.names = FALSE, 
                             show_column_names = TRUE,
                             show.group.dendrogram = TRUE, 
                             show.cell.dendrogram = FALSE,
                             show.feature.dendrogram = FALSE,
                             show.group.legend = FALSE,
                             group.fontsize = 10, 
                             feature.fontsize = 10, 
                             legend.fontsize = 9,
                             raster = TRUE, 
                             draw.lines = TRUE,
                             lines.width = NULL, 
                             disp.identity.legend = TRUE,
                             clustering.within = T, # this only matters if plotting individual single cells
                             archRSort = TRUE,
                             sample.order = NULL,
                             group.bar.height = 1)
{
  
  
  if(!missing(subset.features)){
    message("Subsetting peaks for features")
    idx <- which(overlapsAny(GRanges(rowData(marker.peaks)), subset.features))
    marker.peaks <- marker.peaks[idx]
    #idx <- which(rowData(marker.peaks)$name %in% unlist(subset.features))
    #marker.peaks <- marker.peaks[idx,]
    message("Recommend setting both archRsort and clustering.within to FALSE")
  }
  
  ## returns a normalized filtered list ##
  if(is.null(filt.heatmap.peak)){
    message("Generating matrix of log2 normalized peaks")
    filt.heatmap.peak <- plotMarkerHeatmap(
      seMarker = marker.peaks,
      log2Norm = T,
      plotLog2FC = plot.log2fc,
      returnMatrix = T,
      cutOff = marker.cutoff,
      transpose = TRUE
    )
    legend.title <- bquote(.(legend.title)~log[2]~"FC")
  }
  
  if(plot.log2fc){
    legend.title <- paste(legend.title,"logFC")
  }
  
  # Write peaks and FDR to spreadsheets #
  if(class(marker.peaks)=='SummarizedExperiment'){
    dtheatmap_peaks <- data.table(as.data.table(rowData(marker.peaks)),
                                  "FDR"=do.call(cbind, assays(marker.peaks)$FDR),
                                  "Log2FC"=do.call(cbind, assays(marker.peaks)$Log2FC),
                                  "Mean"=do.call(cbind, assays(marker.peaks)$Mean),
                                  "MeanDiff"=do.call(cbind, assays(marker.peaks)$MeanDiff))
    # Write sample peaks and FDR to spreadsheets #
    if(!is.null(spreadsheets.path)){
      message("Writing to spreadsheet ",spreadsheets.path)
      writexl::write_xlsx(x = dtheatmap_peaks, path = file.path(spreadsheets.path), col_names = T, format_headers = T)
    }
    
    message("Filtering peaks with ",marker.cutoff)
    
    sig.marker.peaks <- getMarkers(marker.peaks, cutOff =  marker.cutoff, returnGR = TRUE)
    # add cluster name to row
    sig.marker.peaks <- lapply( names(sig.marker.peaks), function(x){
      sig.marker.peaks[[x]]$Cluster <- x
      sort(sig.marker.peaks[[x]])
    })
    sig.marker.peaks <- unlist(as(sig.marker.peaks, "GRangesList"), use.names = F)
    ## sort by location and FDR, then get unique
    sig.marker.peaks <- sort(sig.marker.peaks, by=~seqnames+start+end+FDR)
    sig.marker.peaks <- unique(sig.marker.peaks)
    sigpeakcounts <- table(sig.marker.peaks$Cluster)
    
    # message("Filtering peaks with ",marker.cutoff)
    # sigpeaks <- lapply(colnames(marker.peaks), function(x){
    #   sig.str <- gsub("Log2FC",paste0("`Log2FC.",x,"`"),marker.cutoff) %>% 
    #     gsub("FDR",paste0("`FDR.",x,"`"),.) %>% 
    #     gsub("Mean",paste0("`Mean.",x,"`"), .)
    #   dtheatmap_peaks[eval(parse(text = sig.str))]
    #   sigpeakcounts <- sapply(sigpeaks, nrow)
    #}) %>% set_names(colnames(marker.peaks))
  } else if(class(marker.peaks)=='list'){
    sigpeaks <- marker.peaks
    sigpeakcounts <- sapply(sigpeaks, nrow)
  }
  
  sigpeakcounts <- sigpeakcounts[order(sigpeakcounts, decreasing = T)]
  print("sig peak counts:")
  print(sigpeakcounts)
  
  # order largest to smallest
  nlabel <- sapply(names(sigpeakcounts), function(i){ paste0(i," (n=",sigpeakcounts[i],")") })
  print("Label counts:")
  print(nlabel)
  
  # peaks are rows and cells/clusters/samples on columns
  mat <- t(filt.heatmap.peak[names(sigpeakcounts),])
  
  peaks <- rownames(mat)
  print(paste("number of peaks:",length(peaks)))
  
  group_anno <- factor(rep(colnames(mat), nrow(mat)))
  levels(group_anno) <- names(sigpeakcounts)
  print("group anno:")
  print(levels(group_anno))
  ngroup <- length(unique(group_anno))
  
  if(sample.feature.ratio<1){
    message("Sampling ",100*sample.feature.ratio,"% of features in each group")
    d <- dim(mat)
    peaks <- sample(x = rownames(mat), size = d[1]*sample.feature.ratio)
    
    #r <- sapply( levels(group_anno), function(i){
    #  idx <- sample(x = seq_len(length.out = nrow(sigpeaks[[i]])), size = nrow(sigpeaks[[i]])*sample.feature.ratio)
    #  sigpeaks[[i]][idx][,paste0(seqnames,":",start,"-",end)]
    #})
    #print(lengths(r))
    #peaks <- unlist(r, use.names = F)
    
    # subset matrix to peaks
    mat <- mat[peaks,]
    # subset cluster annotation to subsampled peaks
    group_anno <- factor(rep(colnames(mat), nrow(mat)))
    # re assign group size
    ngroup <- length(unique(group_anno))
  } else
    message("Using all peaks")
  
  # reorder matrix columns to match group_anno
  mat <- mat[,levels(group_anno)]
  
  if(archRSort){
    bS <- ArchR:::.binarySort(m = mat, clusterCols = TRUE)
    mat <- bS[[1]][, colnames(mat)]
    clusterCols <- bS[[2]][['order']]
    mat <- mat[, clusterCols, drop = FALSE]
    group_anno <- factor(colnames(mat), levels = colnames(mat))
  } 
  if(!is.null(sample.order)){
    bS <- ArchR:::.binarySort(m = mat, clusterCols = FALSE)
    mat <- bS[[1]][, colnames(mat)]
    clusterCols <- bS[[2]][['order']]
    mat <- mat[, sample.order, drop = FALSE]
    group_anno <- factor(colnames(mat), levels = colnames(mat))
  }
  
  print("Group Anno:")
  print(levels(group_anno))
  #Don't scale data as we get log2FC and a certain threshold see: plotMarkerHeatmap
  #message("Scaling data")
  #mat <- scale(mat)
  # what's the value range in the matrix
  q <- quantile(mat, c(0.1, 0.95))
  q <- c(-max(abs(q)), max(abs(q)))
  
  print("Matrix ranges:")
  print(paste(q, collapse = ":"))
  
  if(is.null(heatmap.colors)){
    heatmap.colors <- colorRampPalette(c("royalblue4","royalblue2","deepskyblue",
                                         "grey90","orangered3","red3","darkred"), space="rgb")(1024);
  }
  col_fun = circlize::colorRamp2(seq(floor(q[1]), ceiling(q[2]), length = length(heatmap.colors)), heatmap.colors)
  library(grid)
  library(ComplexHeatmap)
  top_anno <- NULL
  
  if(group.bar){
    
    lvls <- levels(group_anno)
    
    if(is.null(group.colors)){
      group.colors <- pals::alphabet(length(lvls))
      names(group.colors) <- lvls
    }
    else if(length(names(group.colors))==length(group.colors)){
      group.colors <- group.colors[lvls]
      cat(group.colors,"\n")
    }
    
    df <- data.frame(sample = lvls)
    anno.counts <- sapply(df$sample, function(x){ sigpeakcounts[x] })
    df$sample <- sapply(df$sample, function(x){ nlabel[x] })
    names(group.colors) <- sapply(names(group.colors), function(x){ nlabel[x] })
    
    group.bar.dist <- NULL
    if(show.group.bar.count){
      if(orientation.vertical){
        group.bar.dist <- anno_barplot(
          x = anno.counts,
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = "gray40"), #group.colors
          width = unit(x = group.bar.height, units = 'cm'),
        )
      } else {
        group.bar.dist <- anno_barplot(
          x = anno.counts,
          which = 'row',
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = "gray40"), #group.colors
          height = unit(x = group.bar.height, units = 'cm'),
        )
      }
    }
    
    group_anno <- NULL
    
    if(orientation.vertical){
      group_anno = HeatmapAnnotation(
        show_annotation_name = F, 
        dist = group.bar.dist,
        df = df,
        col = list(sample = c(group.colors)),
        border = TRUE,
        show_legend = show.group.legend)
      colnames(mat) <- sapply(colnames(mat), function(x){ nlabel[x] })
      cluster_anno <- factor(rep(colnames(mat), nrow(mat))) 
      
    } else {
      group_anno = rowAnnotation(
        df = df,
        col = list(sample = c(group.colors)),
        #height = unit(x = group.bar.height, units = 'cm'),
        border = TRUE,
        show_legend = show.group.legend,
        show_annotation_name = F,
        dist = group.bar.dist)
      colnames(mat) <- sapply(colnames(mat), function(x){ nlabel[x] })
      cluster_anno <- factor(rep(colnames(mat), nrow(mat)))
    }
    
  }
  
  # prefer group
  dendc <- F
  column_title <- NULL
  column_split <- NULL
  show_column_dend <- FALSE
  cluster_rows <- FALSE
  cluster_column_slices <- FALSE
  ## We don't input single cell so clustering within is same as group clustering
  if(archRSort | nrow(mat) >= 65000){
    msg <- ifelse(nrow(mat) >= 65000,
                  "Complexheatmap cannot handle dendrogram for large matrices, turning off dendrogram.",
                  "Using archR Sort")
    message(msg)
    #dendc <- as.dendrogram(hclust(dist(t(mat))))
    #show_column_dend <- T
    #column_split <- ngroup
  }
  else{
    if(clustering.within){
      message("Performing within & between group clustering")
      dendc <- cluster_within_group(mat, levels(cluster_anno)) # both groups and within groups are reordered, better for discerning clusters, but looks rougher 
      # For large heatmaps clustering rows with dendrogram works best
      cluster_rows <- as.dendrogram(hclust(dist(mat, method = 'binary')))
    }
    else{
      message("Performing between group clustering")
      dendc <- cluster_between_groups(mat, levels(cluster_anno)) # only reordered between groups, looks more smooth
      # For smaller heatmaps (with labels) clustering rows based on dendc works better
      cluster_rows <- TRUE
      #cluster_rows <- as.dendrogram(hclust(dist(mat, method = 'binary')))
      show_column_names <- FALSE # dendrogram shows these
    }
  }
  
  if(show.group.dendrogram){
    message("Show group dendrogram")
    if(!dendc){ # archR sort
      dendc <- T #as.dendrogram(hclust(dist(t(mat))))
      #dendc <- as.dendrogram(hclust(dist(t(mat))))
      column_title <- NULL
    }
    cluster_column_slices <- TRUE
    show_column_dend <- T
    column_split <- ngroup
  } else if (show.cell.dendrogram){
    message("Show cell dendrogram")
    cluster_column_slices <- TRUE
    show_column_dend <- T
    column_split <- cluster_anno
  }
  
  
  heatmap_legend_param <- NULL
  if(orientation.vertical){
    heatmap_legend_param <- list(title = legend.title, 
                                 title_position = 'topleft',
                                 direction = 'vertical',
                                 legend_position="right",
                                 legend_gp = gpar(fontsize = legend.fontsize))
  } else {
    heatmap_legend_param = list(#at = range(mat),
      title = legend.title,
      title_position = 'topcenter',
      legend_gp = gpar(fontsize = legend.fontsize),
      legend_width = unit(3, "cm"),
      grid_height = unit(.3, "cm"),
      legend_position="bottom",
      legend_direction = 'horizontal')
  }
  
  feature_anno <- NULL
  if(!missing(label.features)){
    if(orientation.vertical)
      feature_anno <- rowAnnotation(foo = anno_mark(at = match(x = label.features, table = rownames(mat)), labels = label.features))
    else
      feature_anno <- HeatmapAnnotation(foo = anno_mark(at = match(x = label.features, table = rownames(mat)), labels = label.features))
  }
  
  message("Generating heatmap orientation ",ifelse(orientation.vertical,'vertical','horizontal'))
  cat(dim(mat))
  ht <- NULL
  if(orientation.vertical)
    ht <- Heatmap(mat,
                  name = title,
                  row_title = marker.cutoff, 
                  column_title = column_title,
                  cluster_columns = dendc,
                  column_split = column_split,
                  show_column_dend = show_column_dend,
                  cluster_column_slices = cluster_column_slices,
                  column_title_gp = gpar(fontsize = group.fontsize),
                  column_gap = unit(0.5, "mm"),
                  cluster_rows = cluster_rows,
                  show_row_dend = show.feature.dendrogram,
                  col = col_fun,
                  row_names_gp = gpar(fontsize = feature.fontsize),
                  column_title_rot = 90,
                  top_annotation = group_anno,
                  right_annotation = feature_anno,
                  heatmap_legend_param = heatmap_legend_param,
                  show_column_names = show_column_names,
                  show_row_names = show.feature.names,
                  use_raster = raster,
                  border = 'black',
                  raster_quality = 5)
  else
    ht <- Heatmap(t(mat),
                  name = title,
                  column_title = title,
                  row_title = marker.cutoff,
                  cluster_rows = dendc,
                  row_split = column_split,
                  show_row_dend = show_column_dend,
                  cluster_row_slices = cluster_column_slices,
                  row_title_gp = gpar(fontsize = group.fontsize),
                  row_gap = unit(0.5, "mm"),
                  cluster_columns = cluster_rows,
                  show_column_dend = show.feature.dendrogram,
                  col = col_fun,
                  column_names_gp = gpar(fontsize = feature.fontsize),
                  row_title_rot = 90,
                  top_annotation = feature_anno,
                  right_annotation = group_anno,
                  heatmap_legend_param = heatmap_legend_param,
                  show_row_names = show_column_names,
                  show_column_names = show.feature.names,
                  use_raster = raster,
                  border = 'black',
                  raster_quality = 5)
  ht
  
}

PlotBarVariableInClusters <- function(ArchRProj, 
                                      variable='Condition', 
                                      yvariable='Clusters',
                                      palette=condition_palette, 
                                      addcounts=T, 
                                      filename="", 
                                      orientation='horizontal', 
                                      xlabel = "Cluster number",
                                      width_multiplier=1,
                                      height_multiplier=1,
                                      na.color='gray60',
                                      legend.direction="horizontal",
                                      legend.position="bottom",
                                      legend.fontsize=11,
                                      legend.ncol=1){
  
  dt <- as.data.table(ArchRProj@cellColData)
  dt <- dt[,c("xvar", "variable") := list(get(yvariable), get(variable))]
  dt <- dt[,list(xvar,variable)]
  
  dt[,nCluster := .N, by=xvar]
  dt[,nClusterExperiment := .N, by=list(xvar, variable)]
  dt <- unique(dt)
  
  g <- NULL
  if(orientation=='horizontal'){
    dt <- dt[order(nCluster, decreasing = T),]
    dt[,Clusters:=factor(x = xvar, levels = unique(xvar), ordered = T)]
    dt[,percClusterExperiment := 100 * (nClusterExperiment / nCluster) ]
    dtlab <- unique(dt[,list(Clusters, nCluster)])
    g <- ggplot(dt, aes(x=Clusters, y=percClusterExperiment)) +
      geom_bar(stat='identity', aes(fill=variable), width = .95) +
      scale_fill_manual(eval(variable),values = palette) +
      ylab("Percentage of cells") + xlab(xlabel) +
      theme_ArchR() +
      coord_cartesian(ylim = c(0,100), expand = F, clip = 'off') +
      theme(panel.ontop = T,
            panel.grid = element_blank(),
            panel.grid.major.y = element_line(colour = 'black', linewidth = .2, linetype = 'longdash'),
            legend.direction = legend.direction, 
            legend.position = legend.position,
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 12),
            text = element_text(size = 14),
            legend.text = element_text(size = legend.fontsize),
            legend.background = element_blank(),
            legend.margin = margin(0,1,0,1, unit = "cm"), #trbl,
            plot.background = element_blank())
    if(addcounts)
      g <- g + geom_text(data = dtlab, aes(label=paste0("n=",nCluster)), y=100, color='black', 
                         vjust=0, hjust=0, angle=45, size=3.5)
    
    g <- g + guides(fill=guide_legend(ncol = legend.ncol, title.position = "top", title.hjust = 0.5))
    
    if(filename !=""){
      fout <- file.path(getOutputDirectory(ArchRProj),'Plots',paste0(filename,'.pdf'))
      message("Plotting to:",fout)
      pdf(file = fout, width = 7, height = 3.5, bg = 'transparent')
      print(g)
      dev.off()
    }
  } else { #vertical orientation
    dt <- dt[order(nCluster, decreasing = F),]
    dt[,Clusters:=factor(x = xvar, levels = unique(xvar), ordered = T)]
    dt[,percClusterExperiment := 100 * (nClusterExperiment / nCluster) ]
    dtlab <- unique(dt[,list(Clusters, nCluster)])
    
    message("Plotting ",sum(dtlab$nCluster), " cells")
    
    g <- ggplot(dt, aes(x=percClusterExperiment, y=Clusters)) +
      geom_bar(stat='identity', aes(fill=variable), width = .95) +
      scale_fill_manual(variable, values = palette, na.value = na.color) +
      xlab("Percentage of cells") +  ylab(xlabel) +
      theme_ArchR() +
      coord_cartesian(xlim = c(0,100), expand = F, clip = 'off') +
      theme(panel.ontop = T,
            panel.grid = element_blank(),
            panel.grid.major.x = element_line(colour = 'black', linewidth = .2, linetype = 'longdash'),
            legend.direction = legend.direction, 
            legend.position = legend.position,
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title = element_text(size = 14),
            text = element_text(size = 14),
            legend.text = element_text(size = legend.fontsize),
            legend.background = element_blank(), 
            plot.margin = unit(c(1,2,1,1), "cm"), #trbl
            legend.margin = margin(0,1,0,1, unit = "cm"), #trbl,
            plot.background = element_blank())
    g <- g + guides(fill=guide_legend(ncol = legend.ncol, title.position = "top", title.hjust = 0.5))
    if(addcounts)
      g <- g + geom_text(data = dtlab, aes(label=paste0("n=",nCluster)), 
                         x=100, color='black', 
                         vjust=0.5, hjust=0, angle=0, size=4)
    if(filename !=""){
      width <- 4 * width_multiplier
      height <- 6 * height_multiplier
      fout <- file.path(getOutputDirectory(ArchRProj),'Plots',paste0(filename,'_vertical.pdf'))
      message("Plotting to:",fout)
      pdf(file = fout, width = width, height = height, bg = 'transparent')
      print(g)
      dev.off()
    }
  }
  g
}



PlotCorDimensionDepth <- function(ArchRProj){
  dt <- data.table(do.call(cbind, ArchRProj@reducedDims$IterativeLSI$corToDepth), keep.rownames = 'dimension')
  dt <- melt(dt, id.vars = 'dimension')
  dt[,dimension:= gsub("LSI","",dimension)]
  dt[,dimension := factor(x = dimension, levels = unique(dimension), ordered = T)]
  ggplot(dt, aes(x=dimension, y=value, color=variable)) + geom_point(size=3) + ggtitle(label = 'PCC Depth & dimension') + 
    scale_color_brewer("LSI Dimension", palette = 'Set1') + xlab("LSI dimension") + ylab("PCC") + theme_bw()
}


PlotIntegrationHeatmap_ATAC_RNA <- function(ArchRProj, 
                                            atac_palette,
                                            rna_palette,
                                            atac_cluster_name = 'Clusters',
                                            rna_cluster_name = 'gex_seurat_clusters',
                                            heatmap_palette=gist_heat_r,
                                            width_multiplier = 1,
                                            height_multiplier = 1,
                                            show.percentages=FALSE){
  
  dt <- data.table(table(data.frame(rna_group=as.character(ArchRProj@cellColData[[rna_cluster_name]]), # ), 
                                    atac_cluster=ArchRProj@cellColData[[atac_cluster_name]])))
  
  dt[,rna_group:=as.character(rna_group)]
  dt <- dt[rna_group!=-1]
  dt <- dt[order(atac_cluster)]
  dt[rna_group=="", rna_group:=NA]
  dt[,rna_percent:= 100*N/sum(N), by=rna_group]
  dtmat <- dcast.data.table(dt, atac_cluster ~ rna_group, value.var = 'rna_percent')
  mat <- as.matrix(dtmat[,-c('atac_cluster'), with=F])
  rownames(mat) <- dtmat$atac_cluster
  
  cm = cor(mat)
  
  od =  hclust(dist(mat))$order
  
  mat = mat[od, ]
  
  bS <- ArchR:::.binarySort(m = mat, clusterCols = TRUE)
  mat <- bS[[1]][, colnames(mat)]
  clusterCols <- 1:ncol(mat)
  clusterCols <- bS[[2]][['order']]
  mat <- mat[, clusterCols, drop = FALSE]
  
  default.cell.outline <- gpar(col = "gray80", lwd = 1)
  
  rna_palette <- rna_palette[!is.na(names(rna_palette))]
  ha_column = HeatmapAnnotation(df = data.frame(rna_group=colnames(mat)),
                                show_annotation_name = F,
                                show_legend = F,
                                col = list(rna_group = rna_palette[colnames(mat)]),
                                annotation_legend_param = list(title="", labels_gp = gpar(fontsize = 14)))  
  heatmap <- ComplexHeatmap::Heatmap(name = "RNA % shared in ATAC",
                                     matrix = mat,
                                     row_title = "snATAC clusters", 
                                     row_title_side = "left", 
                                     cluster_rows = T, 
                                     cluster_columns = F,
                                     column_title = "scRNA clusters", 
                                     column_title_side = "bottom",
                                     rect_gp = default.cell.outline, # creates "space" between cells
                                     #rect_gp = gpar(col = "black", lwd = .5),
                                     #col = unname(ArchR::ArchRPalettes$whitePurple),
                                     col = heatmap_palette,
                                     bottom_annotation = ha_column,
                                     show_column_dend = F, 
                                     show_row_dend = F,
                                     cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                                       if(show.percentages==TRUE){
                                         ifelse(mat[i, j] > 40,
                                                shadowtext::grid.shadowtext(paste0(round(mat[i, j], digits = 1),"%"), x, y, 
                                                                            bg.colour = 'grey60', gp = gpar(col='white'), bg.r = 0.05),
                                                ifelse(mat[i, j] <= 0,
                                                       grid.text("-", x, y,  gp=gpar(col='grey40')),
                                                       grid.text(paste0(round(mat[i, j]),"%"), x, y,  gp=gpar(col='grey10'))
                                                )
                                         )
                                       }
                                       grid.rect(x, y, w, h, gp = gpar(fill = "transparent", col = "white", lwd = .5)) # outlines cells
                                     },
                                     heatmap_legend_param = list(at = c(0, 100),
                                                                 labels = c("0%", "100%"),
                                                                 title_position = 'leftcenter-rot',
                                                                 title_gp = gpar(fontsize = 12),
                                                                 grid_width = unit(.5, "cm"),
                                                                 border = 'black',
                                                                 legend_height = unit(4, "cm")),
                                     row_names_side = "left")
  pdf(file = file.path(getOutputDirectory(ArchRProj),'Plots',
                       paste0(atac_cluster_name,'_',atac_cluster_name,'_Integration_Heatmap.pdf')), 
      width = 4.1*width_multiplier, height = 3.6*height_multiplier, bg = 'transparent')
  draw(heatmap)
  dev.off()
  
  heatmap
}

plotEnrichMotifHeatmap <- function(enrichMotifs,
                                   geneMotifSameMat=NULL,
                                   group.by=NULL,
                                   n=10,
                                   mlog10Padj.cutOff=0,
                                   labelMotifs=T,
                                   labelMotifN=T,
                                   labelMotifFrac=F,
                                   labelSpecificMotifs,
                                   labelMotifClusters,
                                   remove.redundant=T,
                                   remove.redundant.family=F,
                                   cluster.columns=T,
                                   cluster.rows=F,
                                   scaling = 'None', #column
                                   add.score.to.name=T,
                                   label.cluster.family=T,
                                   col=unname(ArchRPalettes$whitePurple),
                                   na.col='gray90',
                                   cell.outline=T,
                                   label.col=pals::polychrome(36),
                                   orientation_horizontal=T,
                                   title.legend ="Average Accessibility",
                                   title = "Enriched transcription factor motifs per family",
                                   row_title = "",
                                   row_name_fontsize=12,
                                   column_name_fontsize=16,
                                   outline=T,
                                   legend.direction="horizontal",
                                   quantile.min=.01,
                                   quantile.max=.99,
                                   label.outline=0.5){
  
  require(circlize)
  require(ComplexHeatmap)
  
  default.cell.outline <- gpar(col = "gray30", lwd = 1)
  
  if(!cell.outline)
    default.cell.outline <- gpar(col = NA)
  
  
  mat <- NULL
  if ('matrix' %in% class(enrichMotifs))
    mat <- enrichMotifs
  else if('SummarizedExperiment' %in% class(enrichMotifs))
    mat <- as.matrix(assays(enrichMotifs)[["mlog10Padj"]])
  else if ('data.table' %in% class(enrichMotifs)){
    mat <- as.matrix(enrichMotifs, rownames='motifID')
    mat <- mat[,grep("mlog10Padj", colnames(mat))]
    colnames(mat) <- gsub('mlog10Padj','',colnames(mat))
  }
  
  # remove empty
  ridx <- which(colSums(mat)==0)
  if(length(ridx)>0)
    mat <- mat[,-ridx]
  
  # remove.redundant
  # Keeps n items with padj value > mlog10Padj.cutOff
  keep <- lapply(seq_len(ncol(mat)), function(x) {
    sorted <- mat[order(mat[, x], decreasing = T),]
    if(remove.redundant){
      mnames <- gsub("#.*","",rownames(sorted))
      dups <- which(duplicated(mnames))
      if(length(dups) > 0){
        sorted <- sorted[-dups,]
      }
    }
    if(remove.redundant.family){
      fnames <- gsub(".*#","",rownames(sorted))
      dups <- which(duplicated(fnames))
      if(length(dups) > 0){
        sorted <- sorted[-dups,]
      }
    }
    pass <- which(sorted[, abs(x)] >= abs(mlog10Padj.cutOff))
    head(names(pass), n)
  }) %>% unlist %>% unique
  if(length(keep)==0)
    stop("All entries filtered out, consider reducing cutoff")
  mat <- mat[keep, , drop = FALSE]
  
  mat[is.infinite(mat)] <- NA
  
  samples <- colnames(mat)
  # mat[features, samples]
  if(cluster.columns){
    bS <- ArchR:::.binarySort(m = mat, clusterCols = cluster.columns)
    mat <- bS[[1]][, colnames(mat)]
    clusterCols <- 1:ncol(mat)
    clusterCols <- bS[[2]][['order']]
    mat <- mat[, clusterCols, drop = FALSE]
  } else{ # order by score
    cnames <- melt(data.table(mat, keep.rownames = 'motifid'), id.vars = 'motifid')[order(value, decreasing = T), motifid] %>% unique()
    mat <- mat[cnames,]
  }
  
  if(orientation_horizontal)
    mat <- mat[,samples]
  
  if(!is.null(geneMotifSameMat)){
    geneMotifSameMat <- geneMotifSameMat[rownames(mat), colnames(mat)]
  } else {
    geneMotifSameMat <- matrix(data = NA, nrow = nrow(mat), ncol = ncol(mat), dimnames = list(rownames(mat), colnames(mat)))
  }
  
  # get families
  
  cluster.names <- gsub(".*#","", rownames(mat)) %>% gsub("-[0-9]+-[01].*","",.)
  #cluster.names <- sapply(strsplit(x = rownames(mat), split='-'), function(x){ x[[ (length(x)-2) ]] })
  #cat(cluster.names)
  
  # clean up motif names
  rownames(mat) <- gsub("#.*","",rownames(mat))
  
  # clean up row names
  if(add.score.to.name){
    rownames(mat) <- round(apply(mat, 1, max)) %>% paste0(names(.)," (",.,")")
  }
  
  rowScale <- function(mat=NULL, min=NULL, max=NULL){
    if (!is.null(min)) {
      rMin <- min
    }
    else {
      rMin <- matrixStats::rowMins(mat, na.rm = T)
    }
    if (!is.null(max)) {
      rMax <- max
    }
    else {
      rMax <- matrixStats::rowMaxs(mat, na.rm = T)
    }
    # rScale <- rMax - rMin
    # matDiff <- mat - rMin
    # matScale <- matDiff/rScale
    #cat("rMax",rMax,"\n")
    cat("length rMin",length(rMin),"\n")
    rScale <- pmax(rMax, abs(rMin)) # centering point
    cat("length rScale",length(rScale),"\n")
    matDiff <- mat # - rMin
    matScale <- matDiff/rScale # centering
    out <- list(mat = matScale, min = rMin, max = rMax)
    return(out)
  }
  
  # scaling
  if(scaling=='sample'){
    #mat <- ArchR:::.rowScale(t(mat), min = 0)[[1]] %>% t(.) # scales across clusters
    mat <- rowScale(t(mat))[[1]] %>% t(.) # scales across clusters
    scaling <- "\nscaled by sample"
  }
  else if(scaling=='motif'){
    ##mat <- ArchR:::.rowScale(mat, min = 0)[[1]] # scales TF across motifs
    mat <- rowScale(mat)[[1]] # scales TF across motifs
    scaling <- "\nscaled by motif"
  }
  else{
    scaling <- "\nNo scaling"
  }
  message(scaling)
  val.range <- c(quantile(mat, quantile.min, na.rm=T), quantile(mat, quantile.max, na.rm=T))
  #val.range <- range(mat, na.rm = T)
  cat("scale.range = ",val.range,"\n")
  heatmap.legend.ticks <- c(val.range[1],val.range[2])
  heatmap.legend.labels <- c("min", "max")
  
  if(min(val.range, na.rm = T) < 0){
    val.range <- c(-max(abs(val.range), na.rm = T),max(abs(val.range), na.rm = T))
    heatmap.legend.ticks <- c(val.range[1],0,val.range[2])
    heatmap.legend.labels <- c("min",0,"max")  
  }
  
  cat("new scale.range = ",val.range,"\n")
  
  heatcolors <- colorRamp2(seq(val.range[1], val.range[2], length = length(col)), col)
  
  #ncl <- apply(mat, 2, function(x) { length(which(!is.na(x))) } )
  nclup <- apply(mat, 2, function(x) { length(which(x>0)) } )
  nclempty <- apply(mat, 2, function(x) { length(which(is.na(x))) } )
  ncldown <- apply(mat, 2, function(x) { length(which(x<0)) } )
  fracup <- nclup / (nclup + ncldown + nclempty)
  fracdown <- ncldown / (nclup + ncldown + nclempty)
  cat(paste(colnames(mat),nclup+ncldown,collapse = ", "),"\n")
  clusterRows <- F
  if(cluster.rows){
    cmat <- mat
    cmat[is.na(cmat)] <- 0
    cmat[is.nan(cmat)] <- 0
    cmat[is.infinite(cmat)] <- 0
    clusterRows <- as.dendrogram(hclust(dist(cmat))) 
  }
  
  clusterCols <- 1:ncol(mat)
  
  cluster_anno <- clusterCols
  
  if(!is.null(group.by))
    cluster_anno<- factor(proj@cellColData[colnames(mat),group.by])
  cat(levels(cluster_anno))
  
  # get families
  #cluster.names <- gsub(".*#","", rownames(mat)) %>% gsub("-[0-9]+-[01].*","",.)
  #cat(cluster.names)
  
  # clean up motif names
  #rownames(mat) <- gsub("#.*","",rownames(mat))
  
  ht <- NULL
  
  #### horizontal orientation ####
  if(orientation_horizontal){
    message("\nCreating horizontal heatmap with ",nrow(mat)," motifs")
    # Barplot annotation
    ha <- NULL
    # these are number of events
    if(labelMotifN){
      ha = rowAnnotation(
        dist0 = anno_barplot(
          which = 'row',
          x = nclup[clusterCols]+ncldown[clusterCols],
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = 'black'),
          width = unit(2, "cm")
        )
        # dist1 = anno_barplot(
        #   which = 'row',
        #   x = nclup[clusterCols],
        #   bar_width = 1,
        #   border = FALSE,
        #   gp = gpar(col = "white", fill = tail(col, n=1)),
        #   width = unit(2, "cm")
        # ),
        # dist2 = anno_barplot(
        #   which = 'row',
        #   x = ncldown[clusterCols],
        #   bar_width = 1,
        #   border = FALSE,
        #   gp = gpar(col = "white", fill = col[1]),
        #   width = unit(2, "cm")
        # )
        ,show_annotation_name = FALSE )
    } else if(labelMotifFrac){
      # This is fraction of events
      ha = rowAnnotation(
        dist3 = anno_barplot(
          which = 'row',
          x = cbind(fracup[clusterCols], fracdown[clusterCols]),
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = c(tail(col, n=1), head(col, n=1)) ), #up color, down color
          #gp = gpar(col = "white", fill = tail(col, n=1)),
          width = unit(2, "cm")
        ), 
        show_annotation_name = FALSE)
    }
    
    # label motif family
    colAnn <- NULL
    if(!missing(labelMotifClusters)){
      # Bottom side motif family labels
      rloc <- grep(paste(labelMotifClusters, collapse = "-"), cluster.names, ignore.case = T)
      label.names <- unique(cluster.names[rloc])
      labels <- label.col[seq_along(label.names)]
      names(labels) <- label.names
      alllabels <- cluster.names
      alllabels[setdiff(seq_along(alllabels), rloc)] <- NA
      nc <- ifelse(length(labels) < 20, 2, 3)
      colAnn <- HeatmapAnnotation(df=data.frame('Clustered'=alllabels),
                                  col = list('Clustered'=labels),
                                  which='col',
                                  show_annotation_name = T, 
                                  annotation_label = "Motif family", 
                                  annotation_name_side = "left", 
                                  annotation_name_gp = gpar(fontface='italic'),
                                  na_col = na.col,
                                  gp = gpar(col = "black", lwd=label.outline),
                                  show_legend = F,
                                  # legend for cluster annotation
                                  annotation_legend_param = list(legend_height = unit(5, "in"),
                                                                 labels_gp = gpar(fontsize = 12),
                                                                 border = TRUE,
                                                                 ncol = nc),
                                  text=anno_text(alllabels, just = "left", rot = 45, location = 0,
                                                 gp = gpar(fontsize = 11, col = 'black'))
      )
    }
    if(outline)
      outline <- gpar(col = "grey10", lwd = 0.1)
    else outline <- gpar(col = NA)
    
    legend.height <- unit(1.5, "cm")
    if(legend.direction!='vertical')
      legend.height <- unit(.5, "cm")
    ht <- ComplexHeatmap::Heatmap(name=paste0(title.legend, scaling),
                                  matrix = t(mat),
                                  column_title = title,
                                  #col = unname(ArchRPalettes$comet),
                                  col = heatcolors,
                                  na_col = na.col,
                                  rect_gp = default.cell.outline, # creates "space" between cells
                                  border = T,
                                  border_gp = gpar(col='black', lwd=.2),
                                  show_column_dend = F,
                                  cluster_columns = cluster.columns,
                                  cluster_rows = cluster.rows,
                                  use_raster = F,
                                  show_column_names = labelMotifs,
                                  column_names_gp = gpar(fontsize=column_name_fontsize),
                                  column_names_side = 'bottom',
                                  row_names_side = 'left',
                                  heatmap_legend_param = list(at = heatmap.legend.ticks,
                                                              labels = heatmap.legend.labels,
                                                              title_position = 'topcenter',
                                                              title_gp = gpar(fontsize = 12),
                                                              grid_width = unit(2, "cm"),
                                                              grid_height = legend.height,
                                                              grid_border = 'black',
                                                              border='black',
                                                              gp = gpar(col = "black", lwd=label.outline),
                                                              legend_position = 'bottom',
                                                              legend_direction = legend.direction),
                                  right_annotation = ha,
                                  top_annotation = colAnn,
                                  ## add text to each grid
                                  layer_fun = function(j, i, x, y, w, h, col) {
                                    #grid.points(x, y, pch=ifelse(geneMotifSameMat[j, i]=="","","*"), size = unit(4, "mm"), gp=gpar(col='white'))
                                    if(cell.outline)
                                      grid.rect(x, y, w, h, gp = gpar(fill = "transparent", col = "white", lwd = .5)) # outlines cells
                                  }
    )
    if(labelMotifN){
      # just text
      cluster_anno = anno_mark(at = seq_len(ncol(mat)),
                               labels_gp = gpar(fontsize = row_name_fontsize),
                               labels = paste("n=",nclup[clusterCols]+ncldown[clusterCols]),
                               which = "row",
                               link_width = unit(3, "mm"))
      ht <- ht + rowAnnotation(mark = cluster_anno)
    }
  }
  #### vertical orientation ####
  else {
    message("\nCreating vertical heatmap with ",nrow(mat)," motifs")
    # Barplot indicating number of motifs top annotation
    ha <- NULL
    if(labelMotifN){
      ha = columnAnnotation(
        dist0 = anno_barplot(
          which = 'column',
          x = nclup[clusterCols]+ncldown[clusterCols],
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = 'black', fontsize = column_name_fontsize),
          width = unit(2, "cm")
        )
        # dist1 = anno_barplot(
        #   which = 'row',
        #   x = nclup[clusterCols],
        #   bar_width = 1,
        #   border = FALSE,
        #   gp = gpar(col = "white", fill = tail(col, n=1)),
        #   width = unit(2, "cm")
        # ),
        # dist2 = anno_barplot(
        #   which = 'row',
        #   x = ncldown[clusterCols],
        #   bar_width = 1,
        #   border = FALSE,
        #   gp = gpar(col = "white", fill = col[1]),
        #   width = unit(2, "cm")
        # )
        ,show_annotation_name = FALSE )
    } else if(labelMotifFrac){
      # This is fraction of events
      ha = columnAnnotation(
        dist3 = anno_barplot(
          which = 'column',
          x = cbind(fracup[clusterCols], fracdown[clusterCols]), #up, down
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = c(tail(col, n=1), head(col, n=1)) ), #up color, down color
          width = unit(2, "cm")
        ), 
        show_annotation_name = FALSE)
    }
    
    legend.height <- unit(1.5, "cm")
    if(legend.direction!='vertical')
      legend.height <- unit(.5, "cm")
    ht <- ComplexHeatmap::Heatmap(name=paste0(title.legend, scaling),
                                  matrix = mat,
                                  column_title = title,
                                  row_title = row_title,
                                  top_annotation = ha,
                                  col = heatcolors,
                                  na_col = na.col,
                                  border = T,
                                  border_gp = gpar(col='black', lwd=.2),
                                  rect_gp = default.cell.outline, # creates "space" between cells
                                  show_column_dend = T,
                                  cluster_columns = F,
                                  cluster_rows = clusterRows,
                                  show_row_dend = F,
                                  use_raster = F,
                                  show_row_names = labelMotifs,
                                  row_names_side = 'left',
                                  row_names_gp = gpar(fontsize = row_name_fontsize),
                                  column_names_gp = gpar(fontsize = column_name_fontsize),
                                  column_names_side = 'top',
                                  column_names_centered = F,
                                  heatmap_legend_param = list(at = heatmap.legend.ticks,
                                                              labels = heatmap.legend.labels,
                                                              title_position = 'topleft',
                                                              title_gp = gpar(fontsize = 12),
                                                              labels_gp = gpar(fontsize = 12),
                                                              grid_width = unit(.5, "cm"),
                                                              grid_height = legend.height,
                                                              grid_border = 'black',
                                                              border='black',
                                                              legend_side='bottom',
                                                              legend_direction = legend.direction),
                                  # add text to each grid
                                  layer_fun = function(j, i, x, y, w, h, col) {
                                    #grid.points(x, y, pch=ifelse(geneMotifSameMat[i, j]=="","","*"), size = unit(4, "mm"), gp=gpar(col='white'))
                                    if(cell.outline)
                                      grid.rect(x, y, w, h, gp = gpar(fill = "transparent", col = "white", lwd = .5)) # outlines cells
                                  }
    )
    
    # Label
    if(!missing(labelSpecificMotifs)){
      rloc <- grep(paste(labelSpecificMotifs, collapse = "-"), rownames(mat), ignore.case = T)
      labels <- rownames(mat)[rloc]
      motif_anno = anno_mark(at = rloc,
                             labels_gp = gpar(fontsize = 12),
                             labels = labels,
                             which = "row",
                             link_width = unit(8, "mm"))
      ht <- ht + rowAnnotation(mark = motif_anno)
    }
    
    # Right side motif family labels
    if(!missing(labelMotifClusters)){
      rloc <- grep(paste(labelMotifClusters, collapse = "-"), cluster.names, ignore.case = T)
      label.names <- unique(cluster.names[rloc])
      labels <- label.col[seq_along(label.names)]
      names(labels) <- label.names
      alllabels <- cluster.names
      alllabels[setdiff(seq_along(alllabels), rloc)] <- NA
      ht <- ht + HeatmapAnnotation(#df=data.frame('Clustered'=alllabels),
        bar=alllabels,
        col = list('bar'=labels), # controls bar color
        text=anno_text(alllabels, 
                       gp = gpar(fontsize = 12, just = "left",
                                 location = 0.5,
                                 col = 'black')),
        which='row',
        show_annotation_name = FALSE,
        na_col = na.col,
        gp = gpar(col = "black", lwd=label.outline), 
        show_legend = F,
        annotation_legend_param = list(legend_height = unit(5, "in"),
                                       labels_gp = gpar(fontsize = 13)))
    }
  }
  
  ht
}

plotEnrichGeneScoreHeatmap <- function(featureGS,
                                       geneFeaturesSameMat=NULL,
                                       group.by=NULL,
                                       topn=10,
                                       cutOff=NULL,
                                       labelFeatures=T,
                                       labelFeaturesN=T,
                                       labelFeaturesFrac=F,
                                       labelSpecificFeatures,
                                       labelFeaturesClusters,
                                       remove.redundant=T,
                                       cluster.columns=T,
                                       cluster.rows=F,
                                       scaling = 'None', #column
                                       add.score.to.name=T,
                                       col=unname(ArchRPalettes$whitePurple),
                                       na.col='gray90',
                                       cell.outline=T,
                                       label.col=pals::polychrome(36),
                                       orientation_horizontal=T,
                                       title.legend ="Average Accessibility",
                                       title = "Enriched Gene scores",
                                       row_title = "",
                                       row_name_fontsize=12,
                                       column_name_fontsize=11,
                                       outline=T,
                                       quantile.min=.01,
                                       quantile.max=.99,
                                       label.outline=0.5){
  
  require(circlize)
  require(ComplexHeatmap)
  
  default.cell.outline <- gpar(col = "gray30", lwd = 1)
  
  if(!cell.outline)
    default.cell.outline <- gpar(col = NA)
  
  
  mat <- NULL
  if ('matrix' %in% class(featureGS))
    mat <- featureGS
  else if('SummarizedExperiment' %in% class(featureGS)){
    #gene.marker.cutOff <- "FDR <= 0.01 & Log2FC >= 1"
    #topn <- 20
    cat("using cutoff:",cutOff)
    # returns matrix containing mean log2 'expression' in features
    featureMatGS <- plotMarkerHeatmap(
      seMarker = featureGS,
      cutOff = cutOff, # apply cutoff
      returnMatrix = T,
      transpose = F
    )
    title.legend <- "log2(average accessibility)"
    # large matrix is returned so we filter for topn, apply cutoff anyway
    featureListGS <- getMarkers(featureGS, cutOff = cutOff,  n = topn) # top number per group
    topfeatures <- lapply(featureListGS, function(x) { x$name }) %>% unlist %>% unique
    mat <- featureMatGS[topfeatures,]
  }
  else{
    return("Unknown type of feature set")
  }
  
  # remove empty
  ridx <- which(colSums(mat)==0)
  if(length(ridx)>0)
    mat <- mat[,-ridx]
  
  mat[is.infinite(mat)] <- NA
  
  # subsetting
  if(cluster.columns){
    bS <- ArchR:::.binarySort(m = mat, clusterCols = cluster.columns)
    mat <- bS[[1]][, colnames(mat)]
    clusterCols <- 1:ncol(mat)
    clusterCols <- bS[[2]][['order']]
    mat <- mat[, clusterCols, drop = FALSE]
  } else{ # order by score
    cnames <- melt(data.table(mat, keep.rownames = 'geneID'), id.vars = 'geneID')[order(value, decreasing = T), geneID] %>% unique()
    mat <- mat[cnames,]
  }
  
  
  if(!is.null(geneFeaturesSameMat)){
    geneFeaturesSameMat <- geneFeaturesSameMat[rownames(mat), colnames(mat)]
  } else {
    geneFeaturesSameMat <- matrix(data = NA, nrow = nrow(mat), ncol = ncol(mat), dimnames = list(rownames(mat), colnames(mat)))
  }
  
  # get families
  #cluster.names <- sapply(strsplit(x = rownames(mat), split='-'), function(x){ x[[ (length(x)-2) ]] })
  cluster.names <- sapply(strsplit(x = gsub(".*#","",rownames(mat)), split='/'), function(x){ x[[1]] })
  message("Clustered families:",cluster.names)
  
  # clean up Features names
  rownames(mat) <- gsub("#.*","",rownames(mat))
  
  # clean up row names
  if(add.score.to.name){
    rownames(mat) <- round(apply(mat, 1, max)) %>% paste0(names(.)," (",.,")")
  }
  
  rowScale <- function(mat=NULL, min=NULL, max=NULL){
    if (!is.null(min)) {
      rMin <- min
    }
    else {
      rMin <- matrixStats::rowMins(mat, na.rm = T)
    }
    if (!is.null(max)) {
      rMax <- max
    }
    else {
      rMax <- matrixStats::rowMaxs(mat, na.rm = T)
    }
    # rScale <- rMax - rMin
    # matDiff <- mat - rMin
    # matScale <- matDiff/rScale
    #cat("rMax",rMax,"\n")
    cat("length rMin",length(rMin),"\n")
    rScale <- pmax(rMax, abs(rMin)) # centering point
    cat("length rScale",length(rScale),"\n")
    matDiff <- mat # - rMin
    matScale <- matDiff/rScale # centering
    out <- list(mat = matScale, min = rMin, max = rMax)
    return(out)
  }
  
  # scaling
  if(scaling=='sample'){
    #mat <- ArchR:::.rowScale(t(mat), min = 0)[[1]] %>% t(.) # scales across clusters
    mat <- rowScale(t(mat))[[1]] %>% t(.) # scales across clusters
    scaling <- "\nscaled by sample"
  }
  else if(scaling=='gene'){
    ##mat <- ArchR:::.rowScale(mat, min = 0)[[1]] # scales TF across Featuress
    mat <- rowScale(mat)[[1]] # scales TF across Featuress
    scaling <- "\nscaled by gene"
  }
  else{
    message("No scaling")
    scaling <- ""
  }
  val.range <- c(quantile(mat, quantile.min, na.rm=T), quantile(mat, quantile.max, na.rm=T))
  #val.range <- range(mat, na.rm = T)
  cat("scale.range = ",val.range,"\n")
  heatmap.legend.ticks <- c(val.range[1],val.range[2])
  heatmap.legend.labels <- c("min", "max")
  
  if(min(val.range, na.rm = T) < 0){
    val.range <- c(-max(abs(val.range), na.rm = T),max(abs(val.range), na.rm = T))
    heatmap.legend.ticks <- c(val.range[1],0,val.range[2])
    heatmap.legend.labels <- c("min",0,"max")  
  }
  
  cat("new scale.range = ",val.range,"\n")
  
  heatcolors <- colorRamp2(seq(val.range[1], val.range[2], length = length(col)), col)
  
  #ncl <- apply(mat, 2, function(x) { length(which(!is.na(x))) } )
  nclup <- apply(mat, 2, function(x) { length(which(x>0)) } )
  nclempty <- apply(mat, 2, function(x) { length(which(is.na(x))) } )
  ncldown <- apply(mat, 2, function(x) { length(which(x<0)) } )
  fracup <- nclup / (nclup + ncldown + nclempty)
  fracdown <- ncldown / (nclup + ncldown + nclempty)
  cat(paste(colnames(mat),nclup+ncldown,collapse = ", "),"\n")
  clusterRows <- F
  if(cluster.rows){
    cmat <- mat
    cmat[is.na(cmat)] <- 0
    cmat[is.nan(cmat)] <- 0
    cmat[is.infinite(cmat)] <- 0
    clusterRows <- as.dendrogram(hclust(dist(cmat))) 
  }
  
  clusterCols <- 1:ncol(mat)
  
  cluster_anno <- clusterCols
  
  if(!is.null(group.by))
    cluster_anno<- factor(proj@cellColData[colnames(mat),group.by])
  cat(levels(cluster_anno))
  
  # get families
  #cluster.names <- gsub(".*#","", rownames(mat)) %>% gsub("-[0-9]+-[01].*","",.)
  #cat(cluster.names)
  
  # clean up Features names
  #rownames(mat) <- gsub("#.*","",rownames(mat))
  
  ht <- NULL
  
  # horizontal
  if(orientation_horizontal){
    message("\nCreating horizontal heatmap with ",ncol(mat)," entries and ",nrow(mat)," groups")
    # Barplot annotation
    ha <- NULL
    # these are number of events
    if(labelFeaturesN){
      ha = rowAnnotation(
        dist0 = anno_barplot(
          which = 'row',
          x = nclup[clusterCols]+ncldown[clusterCols],
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = 'black'),
          width = unit(2, "cm")
        )
        # dist1 = anno_barplot(
        #   which = 'row',
        #   x = nclup[clusterCols],
        #   bar_width = 1,
        #   border = FALSE,
        #   gp = gpar(col = "white", fill = tail(col, n=1)),
        #   width = unit(2, "cm")
        # ),
        # dist2 = anno_barplot(
        #   which = 'row',
        #   x = ncldown[clusterCols],
        #   bar_width = 1,
        #   border = FALSE,
        #   gp = gpar(col = "white", fill = col[1]),
        #   width = unit(2, "cm")
        # )
        ,show_annotation_name = FALSE )
    } else if(labelFeaturesFrac){
      # This is fraction of events
      ha = rowAnnotation(
        dist3 = anno_barplot(
          which = 'row',
          x = cbind(fracup[clusterCols], fracdown[clusterCols]),
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = c(tail(col, n=1), head(col, n=1)) ), #up color, down color
          #gp = gpar(col = "white", fill = tail(col, n=1)),
          width = unit(2, "cm")
        ), 
        show_annotation_name = FALSE)
    }
    
    # label clusters
    colAnn <- NULL
    if(!missing(labelFeaturesClusters)){
      # cluster annotation
      rloc <- grep(paste(labelFeaturesClusters, collapse = "-"), cluster.names, ignore.case = T)
      label.names <- unique(cluster.names[rloc])
      labels <- label.col[seq_along(label.names)]
      names(labels) <- label.names
      alllabels <- cluster.names
      alllabels[setdiff(seq_along(alllabels), rloc)] <- NA
      nc <- ifelse(length(labels) < 20, 2, 3)
      colAnn <- HeatmapAnnotation(df=data.frame('Clustered'=alllabels),
                                  col = list('Clustered'=labels),
                                  which='col',
                                  show_annotation_name = FALSE,
                                  na_col = na.col,
                                  gp = gpar(col = "black", lwd=label.outline),
                                  # legend for cluster annotation
                                  annotation_legend_param = list(legend_height = unit(5, "in"), 
                                                                 labels_gp = gpar(fontsize = 12),
                                                                 border = TRUE,
                                                                 ncol = nc))
      
    }
    if(outline)
      outline <- gpar(col = "grey10", lwd = 0.1)
    else outline <- gpar(col = NA)
    
    ht <- ComplexHeatmap::Heatmap(name=paste0(title.legend, scaling),
                                  matrix = t(mat),
                                  column_title = title,
                                  #col = unname(ArchRPalettes$comet),
                                  col = heatcolors,
                                  na_col = na.col,
                                  rect_gp = default.cell.outline, # creates "space" between cells
                                  border = T,
                                  border_gp = gpar(col='black', lwd=.2),
                                  show_column_dend = F,
                                  cluster_columns = F,
                                  cluster_rows = F,
                                  use_raster = F,
                                  show_column_names = labelFeatures,
                                  column_names_gp = gpar(fontsize=column_name_fontsize),
                                  column_names_side = 'top',
                                  row_names_side = 'left',
                                  heatmap_legend_param = list(at = heatmap.legend.ticks,
                                                              labels = heatmap.legend.labels,
                                                              title_position = 'topcenter',
                                                              title_gp = gpar(fontsize = 12),
                                                              grid_width = unit(2, "cm"),
                                                              grid_height = unit(.5, "cm"),
                                                              grid_border = 'black',
                                                              border='black',
                                                              gp = gpar(col = "black", lwd=label.outline),
                                                              legend_position = 'right',
                                                              legend_direction = 'horizontal'),
                                  right_annotation = ha,
                                  top_annotation=colAnn,
                                  ## add text to each grid
                                  layer_fun = function(j, i, x, y, w, h, col) {
                                    #grid.points(x, y, pch=ifelse(geneFeaturesSameMat[j, i]=="","","*"), size = unit(4, "mm"), gp=gpar(col='white'))
                                    if(cell.outline)
                                      grid.rect(x, y, w, h, gp = gpar(fill = "transparent", col = "white", lwd = .5)) # outlines cells
                                  }
    )
    
    if(labelFeaturesN){
      # just text
      cluster_anno = anno_mark(at = seq_len(ncol(mat)),
                               labels_gp = gpar(fontsize = row_name_fontsize),
                               labels = paste("n=",nclup[clusterCols]+ncldown[clusterCols]),
                               which = "row",
                               link_width = unit(3, "mm"))
      ht <- ht + rowAnnotation(mark = cluster_anno)
    }
    ht <- draw(ht, heatmap_legend_side = "bottom", padding = unit(c(3, 2, 5, 8), "mm")) #bltr
    
    
    # vertical orientation
  } else {
    message("\nCreating vertical heatmap with ",nrow(mat)," entries and ",ncol(mat)," groups")
    # Barplot indicating number of Featuress top annotation
    ha <- NULL
    if(labelFeaturesN){
      ha = columnAnnotation(
        dist0 = anno_barplot(
          which = 'column',
          x = nclup[clusterCols]+ncldown[clusterCols],
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = 'black'),
          width = unit(2, "cm")
        )
        # dist1 = anno_barplot(
        #   which = 'row',
        #   x = nclup[clusterCols],
        #   bar_width = 1,
        #   border = FALSE,
        #   gp = gpar(col = "white", fill = tail(col, n=1)),
        #   width = unit(2, "cm")
        # ),
        # dist2 = anno_barplot(
        #   which = 'row',
        #   x = ncldown[clusterCols],
        #   bar_width = 1,
        #   border = FALSE,
        #   gp = gpar(col = "white", fill = col[1]),
        #   width = unit(2, "cm")
        # )
        ,show_annotation_name = FALSE )
    } else if(labelFeaturesFrac){
      # This is fraction of events
      ha = columnAnnotation(
        dist3 = anno_barplot(
          which = 'column',
          x = cbind(fracup[clusterCols], fracdown[clusterCols]), #up, down
          bar_width = 1,
          border = FALSE,
          gp = gpar(col = "white", fill = c(tail(col, n=1), head(col, n=1)) ), #up color, down color
          width = unit(2, "cm")
        ), 
        show_annotation_name = FALSE)
    }
    
    ht <- ComplexHeatmap::Heatmap(name=paste0(title.legend, scaling),
                                  matrix = mat,
                                  column_title = title,
                                  row_title = row_title,
                                  top_annotation = ha,
                                  col = heatcolors,
                                  na_col = na.col,
                                  border = T, 
                                  border_gp = gpar(col='black', lwd=.2),
                                  rect_gp = default.cell.outline, # creates "space" between cells
                                  show_column_dend = T,
                                  cluster_columns = F,
                                  cluster_rows = clusterRows,
                                  show_row_dend = F,
                                  use_raster = F,
                                  show_row_names = labelFeatures,
                                  row_names_gp = gpar(fontsize = row_name_fontsize),
                                  column_names_side = 'top',
                                  column_names_centered = F,
                                  row_names_side = 'left',
                                  heatmap_legend_param = list(at = heatmap.legend.ticks,
                                                              labels = heatmap.legend.labels,
                                                              title_position = 'topleft',
                                                              #title_position = "leftcenter-rot",
                                                              title_gp = gpar(fontsize = 12),
                                                              labels_gp = gpar(fontsize = 12),
                                                              grid_width = unit(.5, "cm"),
                                                              grid_height = unit(1.5, "cm"),
                                                              grid_border = 'black',
                                                              border='black',
                                                              legend_side='bottom',
                                                              legend_direction = 'vertical'),
                                  # add text to each grid
                                  # add text to each grid
                                  layer_fun = function(j, i, x, y, w, h, col) {
                                    #grid.points(x, y, pch=ifelse(geneFeaturesSameMat[i, j]=="","","*"), size = unit(4, "mm"), gp=gpar(col='white'))
                                    if(cell.outline)
                                      grid.rect(x, y, w, h, gp = gpar(fill = "transparent", col = "white", lwd = .5)) # outlines cells
                                  }
    )
    
    # Label
    if(!missing(labelSpecificFeaturess)){
      rloc <- grep(paste(labelSpecificFeaturess, collapse = "-"), rownames(mat), ignore.case = T)
      labels <- rownames(mat)[rloc]
      Features_anno = anno_mark(at = rloc,
                                labels_gp = gpar(fontsize = 12),
                                labels = labels,
                                which = "row",
                                link_width = unit(8, "mm"))
      ht <- ht + rowAnnotation(mark = Features_anno)
    }
    
    # Cluster colored tiles
    if(!missing(labelFeaturesClusters)){
      rloc <- grep(paste(labelFeaturesClusters, collapse = "-"), cluster.names, ignore.case = T)
      label.names <- unique(cluster.names[rloc])
      labels <- label.col[seq_along(label.names)]
      names(labels) <- label.names
      alllabels <- cluster.names
      alllabels[setdiff(seq_along(alllabels), rloc)] <- NA
      ht <- ht + HeatmapAnnotation(df=data.frame('Clustered'=alllabels),
                                   col = list('Clustered'=labels),
                                   which='row',
                                   show_annotation_name = FALSE,
                                   na_col = na.col,
                                   gp = gpar(col = "black", lwd=label.outline),
                                   annotation_legend_param = list(legend_height = unit(5, "in"),
                                                                  labels_gp = gpar(fontsize = 13)))
    }
    ht <- draw(ht, heatmap_legend_side = "right", padding = unit(c(3, 2, 5, 8), "mm")) #bltr
  }
  
  ht
}

PlotPeakAnnotations <- function(ArchRProj){
  
  lookup <- unique(data.table(data.frame(ArchRProj@cellColData[,c("Clusters","Cluster3")])))
  
  peakDT <- data.table(metadata(ArchRProj@peakSet)$PeakCallSummary)
  peakDT[,cluster := gsub("\\(.*","",Group)]
  peakDT$cluster_name <- lookup[match(x = peakDT$cluster, table = lookup$Clusters), Clusters]
  peakDT[is.na(cluster_name), cluster_name := ""]
  peakDT[,cluster_name := paste(cluster_name, gsub(".*\\(","(",Group))]
  peakDT[,clust := as.integer(gsub("^C","",gsub("-.*","",cluster)))]
  peakDT <- peakDT[order(clust),]
  peakDT[,cluster_name := factor(x = cluster_name, levels = unique(cluster_name), ordered = T)]
  
  pal <- paletteDiscrete(values = peakDT$Var1)
  pal <- c(Distal = "#60BA64", Exonic = "#73C6FF", Intronic = "#620FA3", Promoter = "#FFC554", pal)
  
  maxp <- peakDT[,.(x=sum(Freq)), by=cluster][,max(x)]
  p1 <- ggplot(peakDT, aes(x=cluster_name, y=Freq, fill=Var1)) +
    ylab(expression("Number of Peaks "~(x10^3))) + xlab("") +
    geom_col(col='black') +
    scale_y_continuous(limits = c(0, maxp), expand = c(0,0.5)) + # expansion(0,0)) +
    scale_fill_manual("Annotation",values = pal) +
    theme_ArchR(baseSize = 12) + theme(text = element_text(family = 'Arial'),
                                       panel.border = element_blank(),
                                       axis.line = element_line(linewidth = .4, colour = 'black'),
                                       plot.background = element_blank(),
                                       plot.margin = margin(.5, 0, .5, 2, "cm"),
                                       axis.text.x = element_text(angle = 45, hjust=1),
                                       legend.position = 'right',
                                       legend.text = element_text(size=11),
                                       legend.box.background = element_blank(),
                                       legend.background = element_blank())
  
  cairo_pdf(filename = file.path(getOutputDirectory(ArchRProj),'Plots','Peak_annotation_region.pdf'),
            width = 8, height = 4.5,
            bg = 'transparent')
  print(p1)
  dev.off()
}

PlotBarDEPeaks <- function(ArchRProj, seDETests, absLog2FC=1, fdr=0.05, width=15, height=6, outfile.name="Peaks_in_pairwise_de_clusters_Barplot.pdf"){
  dtSummaryPairPlots <- lapply(names(seDETests), function(i){
    dtsummary <- data.table(FDR=seDETests[[i]]@assays@data$FDR$x,
                            Log2FC=seDETests[[i]]@assays@data$Log2FC$x)
    dtsummary <- dtsummary[FDR <= fdr & abs(Log2FC) >= absLog2FC, list(direction=sign(Log2FC))][, .N, by=list(direction)]
    dtsummary[,name:=i]
    dtsummary[,numerator:=gsub("_vs_.*","",name)]
    dtsummary[,denominator:=gsub(".*_vs_","",name)]
  }) %>% rbindlist()
  
  dtSummaryPairPlotsRev <- copy(dtSummaryPairPlots)
  dtSummaryPairPlotsRev[,direction := direction*-1]
  dtSummaryPairPlotsRev[,temp:=numerator]
  dtSummaryPairPlotsRev[,numerator:=denominator]
  dtSummaryPairPlotsRev[,denominator:=temp]
  dtSummaryPairPlotsRev[,temp:=NULL]
  
  dtSummaryPair <- rbindlist(list(dtSummaryPairPlots, dtSummaryPairPlotsRev))
  dtSummaryPair <- dtSummaryPair[order(numerator, denominator)]
  dtSummaryPair[,clust := as.integer(gsub("^C","",gsub("-.*","",denominator)))]
  dtSummaryPair <- dtSummaryPair[order(clust, direction)]
  dtSummaryPair[,denominator := factor(x = denominator, levels = unique(denominator), ordered = T)]
  
  dtSummaryPair[,clust := as.integer(gsub("^C","",gsub("-.*","",numerator)))]
  dtSummaryPair <- dtSummaryPair[order(clust, direction)]
  dtSummaryPair[,numerator := factor(x = numerator, levels = unique(numerator), ordered = T)]
  
  ## Barplot
  ylim <- 110000 #round(max(abs(dtSummaryPair$N)))
  g <- ggplot(dtSummaryPair, aes(x=denominator, y=direction)) +
    geom_col(data = dtSummaryPair[direction > 0], aes(y=direction*N), fill='indianred') +
    geom_col(data = dtSummaryPair[direction < 0], aes(y=direction*N), fill='dodgerblue') +
    ggtitle("Significantly Differential Peak Accessibility", subtitle = bquote(FDR<= ~ .(fdr)~","~Log2FC>= ~ .(absLog2FC))) +
    theme_ArchR() +
    geom_hline(yintercept = 0, linetype='solid') +
    scale_y_continuous(breaks = seq(from=-ylim, to=ylim, by=10000), labels = abs(seq(from=-ylim/1000, to=ylim/1000, by=10))) +
    ylab(expression("Number of Peaks (x"~10^3~")")) + xlab("") +
    #facet_wrap( ~group, nrow = 1, scales = "free_x") +
    facet_grid(~numerator, scales = "free_x", space = 'free_x') +
    theme(strip.placement = "inside",
          plot.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, angle = 90),
          axis.text.x = element_text(angle=90, hjust = 1, vjust=.5),
          panel.spacing.y = unit(.1, "lines"),
          panel.spacing.x = unit(.1, "lines"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank())
  #g <- g + ggforce::facet_row(vars(numerator), scales = 'free_x', space = 'free', strip.position = 'right')
  pdf(file = file.path(getOutputDirectory(ArchRProj),'Plots',outfile.name), width = width, height = height, bg = 'transparent')
  print(g)
  dev.off()
  dtSummaryPair
}

PlotDeviations <- function(dt, n=30, text.adjust=T, splitstr="_and_", nudge.x.text=300, title){
  
  dt <- setDT(dt)
  dt <- dt[order(rank, decreasing = F)]
  if(length(dt[,grepl(splitstr,motifID)])>0){
    dt[, c("label1","label2") :=tstrsplit(motifID,splitstr)]
    dt[,label := paste(
      paste(gsub("#.*","", label1),gsub("-.*","",gsub('.*#','',label1)), sep = ':'),
      paste(gsub("#.*","", label2),gsub("-.*","",gsub('.*#','',label2)), sep = ':'), sep=splitstr)]
  } else{
    dt[,label := paste(gsub("#.*","", motifID),
                       gsub("-.*","",gsub('.*#','',motifID)), sep = ':')]
  }
  
  if(missing(title))
    title <- paste("Top",n,"deviations")
  
  g <- ggplot(dt, aes(x = rank, y = combinedVars, color = combinedVars)) + 
    geom_point(size = 1) + scale_color_gradientn(colors = paletteContinuous(set = "comet"))
  
  if(text.adjust)
    g <- g + ggrepel::geom_label_repel(data = dt[seq_len(n) ], 
                                       aes(x = rank, y = combinedVars, label = label), max.overlaps = 100,
                                       size = 3, box.padding = 0.1, label.padding = 0.1, color = "black", 
                                       segment.size = .1,
                                       nudge_x = nudge.x.text, direction = "y", hjust = "left")
  else
    g <- g + ggrepel::geom_label_repel(data = dt[seq_len(n)], 
                                       aes(x = rank, y = combinedVars, label = label), max.overlaps = 100,
                                       nudge_x = nudge.x.text, min.segment.length = 2, direction = "both",
                                       size = 3, box.padding = 0.1, label.padding = 0.1, color = "black", segment.size = .1)
  
  g <- g + ylab("Variability") + xlab("Rank Sorted Annotations") +
    guides(color = guide_colourbar(title = "Scaled score",direction='horizontal', barwidth = 10, barheight = 1)) + 
    ggtitle(title) +
    theme_ArchR(plotMarginCm = .1)
}

addDeviationsMatrix <- function (ArchRProj = NULL, peakAnnotation = NULL, matches = NULL, 
                                 bgdPeaks = getBgdPeaks(ArchRProj, method = "chromVAR"), matrixName = NULL, 
                                 out = c("z", "deviations"), binarize = FALSE, threads = getArchRThreads(), 
                                 verbose = TRUE, parallelParam = NULL, force = FALSE, logFile = createLogFile("addDeviationsMatrix")) 
{
  message("Using my own implementation to fix following error:\n")
  message('as(<lgCMatrix>, "dgCMatrix") is deprecated since Matrix 1.5-0; do as(., "dMatrix") instead')
  
  
  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  ArchR:::.validInput(input = peakAnnotation, name = "peakAnnotation", 
                      valid = c("character", "null"))
  ArchR:::.validInput(input = matches, name = "matches", valid = c("SummarizedExperiment", 
                                                                   "null"))
  ArchR:::.validInput(input = bgdPeaks, name = "bgdPeaks", valid = c("SummarizedExperiment"))
  ArchR:::.validInput(input = matrixName, name = "matrixName", valid = c("character", 
                                                                         "null"))
  ArchR:::.validInput(input = out, name = "out", valid = c("character"))
  ArchR:::.validInput(input = binarize, name = "binarize", valid = c("boolean"))
  ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))
  ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
  ArchR:::.validInput(input = parallelParam, name = "parallelParam", 
                      valid = c("parallelparam", "null"))
  ArchR:::.validInput(input = force, name = "force", valid = c("boolean"))
  ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))
  tstart <- Sys.time()
  ArchR:::.startLogging(logFile = logFile)
  if (!inherits(ArchRProj, "ArchRProject")) {
    stop("Error Needs to be ArchR Project for Input!")
  }
  ArrowFiles <- getSampleColData(ArchRProj)$ArrowFiles
  threads <- min(length(ArrowFiles), threads)
  allCells <- rownames(getCellColData(ArchRProj))
  outDir <- getOutputDirectory(ArchRProj)
  if (!all(file.exists(ArrowFiles))) {
    stop("Error Input Arrow Files do not all exist!")
  }
  print(matches)
  if (is.null(matches)) {
    anno <- getPeakAnnotation(ArchRProj, peakAnnotation)
    matches <- readRDS(anno$Matches)
    if (is.null(matrixName)) {
      matrixName <- paste0(anno$Name, "Matrix")
    }
  }
  else {
    if (is.null(matrixName)) {
      matrixName <- paste0("MotifMatrix")
    }
  }
  annotationsMatrix <- SummarizedExperiment::assay(matches)
  rownames(annotationsMatrix) <- paste0(seqnames(matches), 
                                        "_", start(matches), "_", end(matches))
  #annotationsMatrix <- as(annotationsMatrix, "dgCMatrix")
  annotationsMatrix <- as(annotationsMatrix, "dMatrix")
  gc()
  ArrowFiles <- getArrowFiles(ArchRProj)
  useMatrix <- "PeakMatrix"
  availableChr <- ArchR:::.availableSeqnames(ArrowFiles, useMatrix)
  rS <- suppressMessages(ArchR:::.getRowSums(ArrowFiles = ArrowFiles, 
                                             seqnames = availableChr, useMatrix = useMatrix, filter0 = FALSE))
  rownames(rS) <- paste0(rS$seqnames, "_", rS$idx)
  rS <- rS[paste0(seqnames(ArchRProj@peakSet), "_", mcols(ArchRProj@peakSet)$idx), 
  ]
  rS$start <- start(ArchRProj@peakSet)
  rS$end <- end(ArchRProj@peakSet)
  rS$GC <- ArchRProj@peakSet$GC
  rownames(rS) <- paste0(rS$seqnames, "_", rS$start, "_", rS$end)
  annotationsMatrix <- annotationsMatrix[rownames(rS), ]
  args <- mget(names(formals()), sys.frame(sys.nframe()))
  args$peakAnnotation <- NULL
  rm(peakAnnotation)
  args$annotationsMatrix <- annotationsMatrix
  args$featureDF <- rS
  args$useMatrix <- useMatrix
  args$ArrowFiles <- ArrowFiles
  args$allCells <- allCells
  args$matrixName <- matrixName
  args$X <- seq_along(ArrowFiles)
  args$FUN <- ArchR:::.addDeviationsMatrix
  args$logFile <- logFile
  args$registryDir <- file.path(getOutputDirectory(ArchRProj), 
                                paste0(matrixName, "DeviationsRegistry"))
  args$ArchRProj <- NULL
  args$matches <- NULL
  ArchR:::.logThis(args, "Deviations-Args", logFile = logFile)
  outList <- ArchR:::.batchlapply(args)
  ArchR:::.logDiffTime("Completed Computing Deviations!", tstart, addHeader = TRUE, 
                       logFile = logFile)
  ArchR:::.endLogging(logFile = logFile)
  gc()
  return(ArchRProj)
}


addMetaDataArchR <- function(ArchRProject, metadata, add_prefix="", force=T){
  
  input <- ArchRProject
  
  dtmetadata <- data.table(metadata, keep.rownames = "cellid")
  
  #order metadata to match project
  dtmetadata <- dtmetadata[order(match(cellid, rownames(input@cellColData)))]
  
  #if( !identical(dtmetadata$cellid, rownames(input@cellColData)) )
  #  stop("CellIds between ArchRProject and metadata are different")
  
  
  transfer_cols <- setdiff(colnames(dtmetadata), "cellid")
  
  setnames(x = dtmetadata, old = transfer_cols, new = paste0(add_prefix,transfer_cols))
  
  if(!force)
    transfer_cols <- setdiff(transfer_cols, colnames(input@cellColData))
  else
    transfer_cols <- setdiff(colnames(dtmetadata), "cellid")
  
  for (i in transfer_cols){
    message("Adding column ",i)
    input <- addCellColData(ArchRProj = input, 
                            data = as.vector(unlist(dtmetadata[,i,with=F], use.names = F)), 
                            name = i, 
                            cells = dtmetadata$cellid,
                            force = force)
  }
  
  input
}

grid_arrange_list_shared_legend <- function(plots, ncol = length(plots), nrow = 1,
                                            legend.position = c("bottom", "right","topright")) {
  
  # plots <- list(...)
  position <- match.arg(legend.position)
  position.justification <- "center"
  if(position=="topright"){
    position <- "right"
    position.justification <- c(1,.9)
  }
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position, legend.justification = position.justification))$grobs
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
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight),
                                            padding = unit.c(unit(0,"lines"))),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(.8, "npc") - lwidth, lwidth*1.2),
                                           padding = unit.c(unit(0,"lines"))))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}
