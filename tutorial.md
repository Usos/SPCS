# Tutorial for SPCS
_Yusong Liu_


Spatial and Pattern Combined Smoothing (SPCS) is a novel two-factor smoothing technique, that employs k-nearest neighbor technique to utilize associations from transcriptome and Euclidean space from the Spatial Transcriptomic (ST) data. Here we present an R implementation of this method and provide a step-by-step example using a PDAC slide (PDACA1) in (Moncada R, et.al 2020).

## Before the smoothing

Before performing the smoothing, we have to introduce the dependent packages and example data we will use.

### Dependencies

All the codes are built with R 4.1.2. Here are the packages we used in our implementation. Note that different versions of the packages may also work well with our codes.




Table: Dependent Packages

|package    |version |
|:----------|:-------|
|Matrix     |1.3.4   |
|rsvd       |1.0.5   |
|factoextra |1.0.7   |
|foreach    |1.5.1   |
|dplyr      |1.0.7   |
|doParallel |1.0.16  |
|ggplot2    |3.3.5   |

### Example Data

Data we used is PDACA1 slide in (Moncada R, et.al 2020). There are two tables, data table and coordinate table, for each slide. Data table is a matrix that each column is representing a spot and each row representing a gene, while coordinate table including spatial coordinate for each spot. We also provided histopathlogical labels as a reference.


```r
data<-read.csv(file = "PDACA1.txt", sep=",")
coord<-read.csv(file = "PDACA1_coord.txt", sep="," )

data[1:10,1:7]
```

```
##         ST1 ST2 ST3 ST4 ST5 ST6 ST7
## A1BG      0   0   0   0   0   0   0
## A1CF      0   0   0   0   0   0   0
## A2M      13   0   4   0   0   0   9
## A2ML1     0   0   0   0   0   0   0
## A3GALT2   0   0   0   0   0   0   0
## A4GALT    1   0   0   0   0   0   0
## A4GNT     0   0   1   0   0   0   0
## AAAS      0   0   0   0   0   0   0
## AACS      0   0   0   0   0   0   0
## AADAC     0   0   0   0   0   0   0
```

```r
coord[1:7,]
```

```
##     coord1 coord2
## ST1     10     10
## ST2     10     13
## ST3     10     14
## ST4     10     15
## ST5     10     16
## ST6     10     17
## ST7     10     19
```

Since ST expression is extremely sparse, we have to do gene filtering on expression data. In this project, we only kept genes with less than 70% zero expressed spots in our analysis. We also perform Count Per Million (CPM) normalization and logarithm transformation as all ST RNA-seq based analysis will do.


```r
gene.zero.cutoff = 0.7 
gene.var.cutoff = 0.0

filtered.genes <- selectGenes(data, gene.zero.cutoff, gene.var.cutoff)

data <- filtered.genes$data
coord <- coord[colnames(data),]

data <- NormalizeLibrarySize(data)
data <- log2(data+1)

data[1:10,1:7]
```

```
##                 ST1      ST2      ST3      ST4      ST5      ST6      ST7
## ABCC3      5.240812 4.519153 4.730248 0.000000 0.000000 0.000000 2.551839
## AC004556.1 0.000000 0.000000 0.000000 0.000000 0.000000 3.883147 2.551839
## ACADVL     3.822477 5.082989 3.783602 6.348425 0.000000 3.883147 5.318697
## ACTB       6.800112 6.280852 6.688896 8.036555 6.805612 4.833412 7.198832
## ACTG1      6.169857 0.000000 0.000000 4.400587 4.534370 6.794952 5.318697
## ACTN1      4.068414 0.000000 3.783602 5.366023 4.534370 4.833412 3.423259
## ACTN4      3.822477 5.082989 6.019190 7.148172 6.805612 6.794952 3.962682
## ADAR       0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
## AEBP1      5.942012 0.000000 5.296978 0.000000 0.000000 0.000000 4.915659
## AES        4.068414 4.519153 3.783602 5.366023 0.000000 0.000000 5.318697
```

## Performing SPCS Smoothing

SPCS smooths the ST expression date using both spatial position and expression pattern information. Let $X_i$ be a vector of gene expression values for spot $i$, smoothed expression $X_i^\prime$ can be calculated by:
$$X_i^\prime=(1-\alpha)X_i+\alpha(\beta\frac{\sum_{j\in N_s(i)}{c_{ji}^sX_j}}{\sum_{k\in N_s(i)}{c_{ki}^s}}+(1-\beta)\frac{\sum_{j\in N_p(i)}{c_{ji}^pX_j}}{\sum_{k\in N_p(i)}{c_{ki}^p}})$$
There are 4 hyper-parameters in the equation above, we determine them as followings:


```r
tau.p=16
tau.s=2
alpha=0.6
beta=0.4
```

If the slide going to smoothed is a Visium slides, please set following parameters:

```r
is.hexa=T
```
If missing spot padding is needed, please set:

```r
is.padding=T
```

Hence, there are two steps before get the smoothed expression: finding neighbors and calculating smoothed expression.

### Finding neighbors

For pattern neighbors, we first transform the expression of spots into a 10-dimensional principal component (PC) space (i.e. pattern space)and then select $\tau_p$ nearest neighbor spots as pattern neighbors for each spot. Last, for each neighbor of object spot, calculating its pattern contribution to the object spot. 

For spatial neighbors, we can calculate their contributions according to its spatial coordinates.

To speed up the neighborhood detection step, we employed doParallel package as parallel calculating engine. The number of threads used to find neighbors is set to 8 in default and can be set according to the number of CPU cores.


```r
nThreads=8

data.mat<-Matrix(as.matrix(data))

rpca.res <- calcRandPCA(t(data.mat)) 
pdist.mat <- calcPDist(rpca.res)  
rmat <- calcPContribution(pdist.mat)

positions <- coord
colnames(positions) <- c("st.x","st.y")
positions$barcodes <- rownames(positions)

p.neighbors <- findPNeighbors(rmat, positions, tau.p, nThreads)
s.neighbors <- findSNeighbors(positions, tau.s, nThreads)
```

### Calculating smoothed expression

After obatined all the pattern and spatial neighbors for each spot, we can calculate smoothed expression using the equation we provided.


```r
data.combinedsmooth <- getCSmoothedExp(data.mat, s.neighbors, p.neighbors, alpha, beta, nThreads)
data.combinedsmooth <- as.data.frame(as.matrix(data.combinedsmooth))
```

### One-step SPCS smoothing

For ease of use, we further provided an one-step smoothing function:


```r
nThreads=8
data.combinedsmooth <- SPCS(data, coord, tau.p, tau.s, alpha, beta, nThreads)

data.combinedsmooth[1:10,1:7]
```

```
##                  ST1       ST2       ST3       ST4       ST5       ST6      ST7
## ABCC3      2.5779190 3.9695155 3.8595167 1.3628916 1.2578928 1.9193818 2.981912
## AC004556.1 0.9148087 0.7979809 0.6368954 0.2624510 0.3813050 2.3590790 2.056393
## ACADVL     3.0153712 4.0926595 3.5918232 4.0227667 1.4830451 3.5395401 4.188449
## ACTB       6.9431678 6.4961383 6.6606753 6.8408532 6.1897022 6.0293188 6.776505
## ACTG1      4.6120628 2.2392488 2.1856797 3.4641784 3.2047171 4.9142529 4.982735
## ACTN1      3.3222626 1.1208215 3.0313268 3.2411354 2.8409670 3.7509756 3.156533
## ACTN4      2.9270816 4.2051710 5.1258645 5.6093018 4.7290141 4.5223957 4.393630
## ADAR       0.9933134 1.3720080 0.9272210 0.3628613 0.3888164 0.9838539 1.536320
## AEBP1      4.7695336 1.4467837 3.4597273 2.0850280 1.2248232 1.0053940 3.093675
## AES        3.8656010 3.4119207 3.1890870 3.5477942 1.0710624 1.5261867 3.681295
```

## Draw Expression heatmaps 

To visualize the result of smoothing, we can draw heatmap of specific gene expression. We show the heatmap of an important oncogene _TM4SF1_ for example. To make heatmap of different genes comparable, expression of genes first linearly transformed into [0,1] by dividing each value by the maximum expression value. 

```r
gene<-"TM4SF1"
barcodes<-colnames(data.combinedsmooth)
data.norm<-as.data.frame(t(data.combinedsmooth)) %>% select(gene)
data.norm<-sweep(data.norm,2,apply(data.norm, 2, max),"/")
data.norm<-data.norm %>% mutate(barcode=barcodes)
coord.heatmap<-coord %>% mutate(barcode=rownames(coord))
data.using<-inner_join(x=coord.heatmap, y=data.norm, copy=T)

ggplot(data.using, aes(coord1, coord2, fill=!!sym(gene))) + ggtitle(gene) + 
      geom_tile() + coord_equal() +
      scale_fill_viridis_c() +
      theme_bw()+
      theme(legend.position="bottom", plot.title=element_text(hjust = 0.5), 
            panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())
```
![heatmap-1](https://user-images.githubusercontent.com/5370174/139857309-89cdba15-b9b8-4a8c-bd1d-a5d7355fd5b0.png)

## Citation
If any code in this reposition is used in any publishable works, please citing:
- **Liu Y, Wang T, Duggan B _et al._**, "SPCS: A Spatial and Pattern Combined Smoothing Method for Spatial Transcriptomic Expression", _Briefings in Bioinformatics_ (2022), bbac116, doi: https://doi.org/10.1093/bib/bbac116.  
      
If the test data in this reposition is used in any publishable works, please citing:
  - **Moncada R, Barkley D, Wagner F _et al._**, "Integrating microarray-based spatial transcriptomics and single-cell RNA-seq reveals tissue architecture in pancreatic ductal adenocarcinomas", _Nature Biotechnology_ (2020), 38:333-342, doi: https://doi.org/10.1038/s41587-019-0392-8.
