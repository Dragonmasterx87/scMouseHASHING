

pbmc.umis <- Read10X(data.dir = "C:/Users/mqadir/Box/Lab 2301/sCell RefrenceDatasets/scISLET mice data 04142020/Raw data/Islet_2_mRNA/umi count/", strip.suffix = TRUE)
pbmc.htos <- Read10X(data.dir = "C:/Users/mqadir/Box/Lab 2301/sCell RefrenceDatasets/scISLET mice data 04142020/Raw data/Islet_2_HTO/umi_count/", gene.column = 1)

pbmc.umis
pbmc.htos

# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

# Confirm that the HTO have the correct names
rownames(pbmc.htos)

# Setup Seurat object
pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis)

# Normalize RNA data with log normalization
pbmc.hashtag <- NormalizeData(pbmc.hashtag)

# Find and scale variable features
pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = "mean.var.plot")
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))

# Add HTO data as a new assay independent from RNA
pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")

# If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
# clustering function for large applications You can also play with additional parameters (see
# documentation for HTODemux()) to adjust the threshold for classification Here we are using the
# default settings
pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)

# Global classification results
table(pbmc.hashtag$HTO_classification.global)

# Group cells based on the max HTO signal
Idents(pbmc.hashtag) <- "HTO_maxID"
RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]])[1:2], ncol = 2)
FeatureScatter(pbmc.hashtag, feature1 = "HTO-3-Male-CTTGCCGCATGTCAT", feature2 = "HTO-4-Female-AAAGCATTCTTCACG")

Idents(pbmc.hashtag) <- "HTO_classification.global"
VlnPlot(pbmc.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# First, we will remove negative cells from the object
pbmc.hashtag.subset <- subset(pbmc.hashtag, idents = "Negative", invert = TRUE)

# Calculate a distance matrix using HTO
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = pbmc.hashtag.subset, assay = "HTO"))))

# Calculate tSNE embeddings with a distance matrix
pbmc.hashtag.subset <- RunTSNE(pbmc.hashtag.subset, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot(pbmc.hashtag.subset)

# To increase the efficiency of plotting, you can subsample cells using the num.cells argument
HTOHeatmap(pbmc.hashtag, assay = "HTO", ncells = 5000)

# Extract the singlets
pbmc.singlet <- subset(pbmc.hashtag, idents = "Singlet")
pbmc.singlet1 <- subset(pbmc.hashtag, subset = HTO_maxID == c("HTO-3-Male-CTTGCCGCATGTCAT", "HTO-4-Female-AAAGCATTCTTCACG"))
alk3n3.integrated <- pbmc.singlet1

# PCA analysis data will be stored in the "reductions' slot
alk3n3.integrated <- RunPCA(object = alk3n3.integrated, features = VariableFeatures(object = alk3n3.integrated))

# Examine and visualize PCA results a few different ways
print(x = alk3n3.integrated[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object = alk3n3.integrated, dims = 1:2, reduction = "pca")
DimPlot(object = alk3n3.integrated, reduction = "pca")
DimHeatmap(object = alk3n3.integrated, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(object = alk3n3.integrated, dims = 1:15, cells = 500, balanced = TRUE)

# Using JacStraw plots to visualize LDR measurements
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
# alk3n3.integrated <- JackStraw(object = alk3n3.integrated, num.replicate = 100)
# alk3n3.integrated <- ScoreJackStraw(object = alk3n3.integrated, dims = 1:10)
# JackStrawPlot(object = alk3n3.integrated, dims = 1:10)

# Elbow plot
ElbowPlot(object = alk3n3.integrated)

# Clustering
alk3n3.integrated <- FindNeighbors(object = alk3n3.integrated, dims = 1:20)
alk3n3.integrated <- FindClusters(object = alk3n3.integrated, resolution = 0.4)

# NON-LINEAR DIMENSIONALITY REDUCTION ####
# RunUMAP
# reticulate::py_install(packages = 'umap-learn')

# We ran RunUMAP on the parameters: metric = cosine and umap.method = UWOT
# you are welcome to try UMAp-learn, the interpretation of the outcome doesnt
# adversally affect the cell clustering on superficial observation
alk3n3.integrated <- RunUMAP(alk3n3.integrated, dims = 1:20, metric = 'cosine', umap.method = 'uwot')

# Visualization
# Note: you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
# UMAP Visualization
p20 <- DimPlot(alk3n3.integrated, group.by = c("HTO_maxID", "seurat_clusters"), combine = FALSE, pt.size = 1)
p20 <- lapply(X = p20, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3))))
CombinePlots(p20)

# DIFFERENTIALLY EXPRESSED-GENE ANALYSIS ####
# Find markers for every cluster compared to all remaining cells, report only the positive ones
# Here we define a DE gene as a gene which has:
# Fold Change of >1.5
# Atleast 10% of cells express that gene
# p value < 0.001
alk3n3.integrated.markers <- FindAllMarkers(object = alk3n3.integrated, only.pos = TRUE, logfc.threshold = 0.41, slot = 'data', test.use = 'wilcox')
alk3n3.integrated.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(alk3n3.integrated.markers, 'C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/scISLET mice data 04142020/output data/mouseislets.csv')

# Look at your default assay
DefaultAssay(object = alk3n3.integrated)

# Change default assay to RNA, to visualize data on graph
# You can toggle between integrated, SCT and RNA to see different expression profiles/different normalizations
DefaultAssay(object = alk3n3.integrated) <- "RNA"

# You can look at the expression of a particular gene across an entire data set here
# As an example we need to remove all mesenchymal cells from the analysis
# We use COL1A1 and THY1 as mesenchymal identifiers
# Violin plot
VlnPlot(object = alk3n3.integrated, features = c("Cyp19a2"), group.by = "seurat_clusters", ncol = 1)

plots <- VlnPlot(alk3n3.integrated, features = c("Srd5a1"), split.by = "HTO_maxID", group.by = "seurat_clusters", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

# UMAP expression plot
FeaturePlot(object = alk3n3.integrated, 
            features = c("Pgr", "Ins1", "Gcg", "Sst", "Ghrl", "Ppy"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 100,
            order = TRUE)


# Look at a heatmap of top 10 most differentially expressed genes
# Change default assay to RNA, save information in the "RNA" assay
DefaultAssay(object = alk3n3.integrated) <- "integrated"

# Create heatmap using doheatmap
top100.all <- alk3n3.integrated.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
DoHeatmap(object = alk3n3.integrated, features = top100.all$gene) + NoLegend()
