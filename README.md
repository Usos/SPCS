# SPCS: A Spatial and Pattern Combined Smoothing Method for Spatial Transcriptomic Expression
_Yusong Liu, Tongxin Wang, Ben Duggan, Kun Huang, Jie Zhang, Xiufen Ye, Travis S. Johnson_  
College of Intelligent Systems Science and Engineering, Harbin Engineering University  
Indiana University School of Medicine

If there's any question about this project, please send emails to johnstr@iu.edu or yexiufen@hrbeu.edu.cn. 

## The method
Spatial and Pattern Combined Smoothing (SPCS) is a novel two-factor smoothing technique, that employs k-nearest neighbor technique to utilize associations from transcriptome and Euclidean space from the Spatial Transcriptomic (ST) data. This reposition is an R implementation of SPCS method, including tutorial and a test slide. The test slide is one of pancreatic ductal adenocarcinoma (PDAC) slide provided in (Moncada R, et.al 2020).

![workflow](https://user-images.githubusercontent.com/5370174/150074625-862c5bca-947f-4440-aa93-6656a12abee3.png)


## The dependencies
This is an R implementation based on following packages:

- Matrix 1.3.4
- rsvd 1.0.5
- factoextra 1.0.7
- foreach 1.5.1
- dplyr 1.0.7
- doParallel 1.0.16

To draw the heatmap of gene expression, please also include:

- ggplot2 3.3.5

P.S: Other versions of these packages may also work well with our implementation.

## The tutorial
We have provided a step-by-step tutorial for our SPCS R implementation. Please see "tutorial.md" for details.

## The citation
If any code in this reposition is used in any publishable works, please citing:
  - **Liu Y, Wang T, Duggan B _et al._**, "SPCS: A Spatial and Pattern Combined Smoothing Method for Spatial Transcriptomic Expression", _Briefings in Bioinformatics_ (2022), bbac116, doi: https://doi.org/10.1093/bib/bbac116.       
If the test data in this reposition is used in any publishable works, please citing:
  - **Moncada R, Barkley D, Wagner F _et al._**, "Integrating microarray-based spatial transcriptomics and single-cell RNA-seq reveals tissue architecture in pancreatic ductal adenocarcinomas", _Nature Biotechnology_ (2020), 38:333-342, doi: https://doi.org/10.1038/s41587-019-0392-8.

