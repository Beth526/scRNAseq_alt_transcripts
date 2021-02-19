
library(Seurat)
library(sctransform)
library(scater)
library(cluster, quietly = TRUE)
library(purrr)

setwd('/Users/beth/Desktop')

orig <- Read10X_h5('Parent_NGSC3_DI_PBMC_filtered_feature_bc_matrix.h5',use.names=T,unique.features=T)
alt<- Read10X_h5('Alternate_transcripts_added_PBMC.h5',use.names=T,unique.features = T)

orig <- CreateSeuratObject(counts = orig, project = "orig")
alt <- CreateSeuratObject(counts = alt, project = "alt")

pipeline <- function(SeuratObj){
  
  #Filter out cells by high percent.mt and counts using 3 median deviations with scater
  MT.genes <- grep("^MT-", rownames(SeuratObj), value=TRUE)
  SeuratObj[["percent.mt"]] <- PercentageFeatureSet(SeuratObj, features = MT.genes)
  temp1 <- scater::isOutlier(SeuratObj$percent.mt, nmads = 3, type = "higher")
  mt.threshold <- min(SeuratObj$percent.mt[temp1])
  temp2 <- scater::isOutlier(SeuratObj$nCount_RNA, nmads = 3, type = "higher")
  Count.threshold <- min(SeuratObj$nCount_RNA[temp2])
  SeuratObj <- subset(SeuratObj, subset = percent.mt < mt.threshold & 
                        nCount_RNA < Count.threshold)
  
  #SCTransform to get normalized gene counts corrected for cell sequencing depth and mt.percent
  SeuratObj <- SCTransform(SeuratObj, return.only.var.genes = FALSE, vars.to.regress = c('percent.mt'))
  
  return(SeuratObj)
}

orig <- pipeline(orig)

#plot PCA orig
orig <- RunPCA(orig, verbose = TRUE)

stdevs <- orig@reductions$pca@stdev
orig_stdev <- DataFrame(stdevs, seq(1:50))
orig_stdev <- as.data.frame(orig_stdev)

p = ggplot(orig_stdev, aes(x=seq.1.50., y=stdevs)) + geom_point(stat="identity",color="blue") + theme_classic() + xlab('PC') + ylab('Standard Deviation') + ylim(0,35) + theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12))
ggsave(filename = "PC_orig.jpg", plot=p, width = 4, height = 3, dpi = 300)

alt <- pipeline(alt)

#plot PCA alt
alt <- RunPCA(alt, verbose = TRUE)
#ElbowPlot(alt, ndims=50)

stdevs <- alt@reductions$pca@stdev
alt_stdev <- DataFrame(stdevs, seq(1:50))
alt_stdev <- as.data.frame(alt_stdev)

p= ggplot(alt_stdev, aes(x=seq.1.50., y=stdevs)) + geom_point(stat="identity",color="blue") + theme_classic() + xlab('PC') + ylab('Standard Deviation') + ylim(0,35) + theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12))
ggsave(filename = "PC_alt.jpg", plot=p, width = 4, height = 3, dpi = 300)

#plot loadings
p = VizDimLoadings(alt, dims = 1:2, reduction = "pca",nfeatures=20, balanced=T)
ggsave(filename = "loadings_alt.jpg", plot=p, width=5, height=6, dpi=300)

p = VizDimLoadings(orig, dims = 1:2, reduction = "pca",nfeatures=20, balanced=T)
ggsave(filename = "loadings_orig.jpg", plot=p, width=5, height=6, dpi=300)

alt <- FindNeighbors(alt, dims = 1:10, verbose = TRUE)
orig <- FindNeighbors(orig, dims = 1:10, verbose = TRUE)

#because alt transcripts bias the PCs, probably should not use alt genes for clustering 
#unless maybe all highly correlated alt transcripts are removed 

#mean silhoutte score function that takes r, resolution for Louvain algorithm in FindClusters()
mean_sil <- function(SeuratObj, r) {
  SeuratObj <-  FindClusters(SeuratObj, verbose = FALSE, resolution = r)  
  dist.matrix <- dist(x = Embeddings(object = SeuratObj[['pca']])[, 1:25])
  clusters <- SeuratObj$seurat_clusters
  sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
  return(mean(sil[,3]))
}


#number of clusters function that takes r, resolution for Louvain algorithm in FindClusters()
number_clusters <- function(SeuratObj, r) {
  SeuratObj <-  FindClusters(SeuratObj, verbose = FALSE, resolution = r)  
  return(length(levels(SeuratObj@active.ident)))
}

#plot mean sil score and number clusters for orig
num_clusters <-  map(seq(0.1,2,0.1), number_clusters, SeuratObj = orig)
num_clusters <- unlist(num_clusters)

mean_sil_score <- map(seq(0.1,2,0.1), mean_sil, SeuratObj = orig)
mean_sil_score <- unlist(mean_sil_score)

orig_df <- DataFrame(seq(0.1,2,0.1),num_clusters,mean_sil_score)
orig_df <- as.data.frame(orig_df)
colnames(orig_df) <- c('resolution','num_clusters', 'mean_sil_score')

p = ggplot(orig_df, aes(x=resolution)) +
  
  geom_line( aes(y=num_clusters), color = 'blue') + 
  geom_line( aes(y=mean_sil_score*100), color = 'red') +
  
  scale_y_continuous(
    limits = c(5,45),
    # Features of the first axis
    name = "Number of Clusters",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./100, name="Mean Sil Score") ) +
  
    theme_classic() +
    
    theme(axis.title.y = element_text(color = 'blue', size=12),
          axis.title.y.right = element_text(color = 'red', size=12)
      ) 

ggsave(filename = "resolution_orig.jpg", plot=p, width=5, height=3, dpi=300)


#plot mean sil score and number clusters for alt 
num_clusters <-  map(seq(0.1,2,0.1), number_clusters, SeuratObj = alt)
num_clusters <- unlist(num_clusters)

mean_sil_score <- map(seq(0.1,2,0.1), mean_sil, SeuratObj = alt)
mean_sil_score <- unlist(mean_sil_score)

alt_df <- DataFrame(seq(0.1,2,0.1),num_clusters,mean_sil_score)
alt_df <- as.data.frame(alt_df)
colnames(alt_df) <- c('resolution','num_clusters', 'mean_sil_score')

p = ggplot(alt_df, aes(x=resolution)) +
  
  geom_line( aes(y=num_clusters), color = 'blue') + 
  geom_line( aes(y=mean_sil_score*100), color = 'red') +
  scale_y_continuous(
    
    # Features of the first axis
    name = "Number of Clusters",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./100, name="Mean Sil Score") ) +
  
  theme_classic() +
  
  theme(axis.title.y = element_text(color = 'blue', size=12),
        axis.title.y.right = element_text(color = 'red', size=12)
  ) 

ggsave(filename = "resolution_alt.jpg", plot=p, width=5, height=3, dpi=300)

#find clusters and run umap  
orig <- FindClusters(orig, verbose = TRUE, resolution = 0.3)
alt <- FindClusters(alt, verbose = TRUE, resolution = 0.3)

alt <- RunUMAP(alt, dims = 1:15, verbose = TRUE)
DimPlot(alt, label = TRUE, reduction = "umap")
ggsave(filename = "alt_0.3_umap.jpg", plot=p, width = 4, height = 4, dpi = 300)

orig <- RunUMAP(orig, dims = 1:7, verbose = TRUE)
p=DimPlot(orig, label = TRUE, reduction = "umap")
ggsave(filename = "orig_0.3_umap.jpg", plot=p, width = 4, height = 4, dpi = 300)


#find the number of alternate transcript markers that are unique cluster markers
#the parent gene and other alternate transcipts are not

all_markers <- FindAllMarkers(alt, verbose = FALSE,
                              assay = "SCT", slot = "data",
                              min.pct = 0.25, logfc.threshold = 0.25)

test <- all_markers[(all_markers$cluster==2) & (grepl("[\'exons\',\'others\']", all_markers$gene) == T),]
test = test[test$p_val_adj<0.05,]
test[order(test$gene),]

how_many_alt <- function(all_markers, cluster){
  test <- all_markers[(all_markers$cluster==cluster),]
  test = test[test$p_val_adj<0.05,]
  test$gene <- sub('-others-.*','',test$gene)
  test$gene <- sub('-exons-.*','',test$gene)
  test = test[order(test$gene),]
  test = test[!(duplicated(test$gene) | duplicated(test$gene, fromLast = TRUE)),]
  test = test[grepl("(exons|others)", rownames(test)) == T,]
  return(nrow(test))
  #return(test)
}

unique_alt_markers <- map(seq(0,10,1),how_many_alt,all_markers = all_markers)

#there was an average of 19 unique alt transcript markers per cluster!
