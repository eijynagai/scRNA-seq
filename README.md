# Single Cell Trancriptomics 

A personal collection focused on scRNA-seq analysis.
<br>

## Significance

scRNA-seq has emerged as a standard technology to study and enhance the 
understanding of transcription. The space here is dedicated to collect
papers, methods, tutorials, and stuff related to the technology. It has 
been manually curated and content checked since 2020 April.

## Welcome to single cell world

### Getting started: scRNA-seq analysis 101

* [Andrews et al., Nature Protocols, 2020. Tutorial: guidelines for the computational analysis of single-cell RNA sequencing data](https://www.nature.com/articles/s41596-020-00409-w)
* [BioTuring's Blog: Single-cell RNA-seq tutorials](https://blog.bioturing.com/category/single-cell-rna-seq-tutorials/)
* [Full introductory course of scRNA-seq data analysis from Hemberg lab](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html)
* [Seurat tutorial from Satija lab](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
* [Orchestrating Single-Cell Analysis with Bioconductor](http://bioconductor.org/books/release/OSCA/)
* [Single cell RNA-seq tutorial from Theis' lab](https://github.com/theislab/single-cell-tutorial)
* [Complete tutorial from Broad Institute](https://broadinstitute.github.io/2019_scWorkshop/)
* [Luecken and Theis, MolSysBio., 2019. Current best practices in single‐cell RNA‐seq analysis: a tutorial](https://www.embopress.org/doi/full/10.15252/msb.20188746)
* [2017/2018 Single Cell RNA Sequecing Analysis Workshop at UCD,UCB,UCSF](https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/)
* [HBC Single-cell RNA-seq: raw sequencing data to counts](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_pre-QC.html)
* [2016 YeoLab single cell bioinformatics](https://github.com/YeoLab/single-cell-bioinformatics)
* [Lun et al., F1000research, A step-by-step workflow for low-level analysis of single-cell RNA-seq data](https://f1000research.com/articles/5-2122/v1)

### Advanced analysis

* [Seurat Wrappers: seurat integrated functions of popular tools such Monocle3, Harmony, Velocito, etc](https://github.com/satijalab/seurat-wrappers)
* 
* [RNA Velocyto video](https://www.youtube.com/watch?v=EPTgF4EA2zY)
* [Workshop on Machine learning for Single Cell Analysis from Krishnaswamy lab 2021, available for free on YouTube s2](https://www.youtube.com/playlist?list=PLW5ccXvlqgHjcCWuxMxdp6rtVs8fop8tM)
<br>

### Review papers
* [Andrew E. Teschendorff & Andrew P. Feinberg, Nature Reviews Genetics, 2021. Statistical mechanics meets single-cell biology](https://www.nature.com/articles/s41576-021-00341-z)
* [Vilas Menon, Briefings in Functional Genomics, 2018. Clustering single cells: a review of approaches on
high-and low-depth single-cell RNA-seq data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6063268/pdf/elx044.pdf)
* [Hwang et al, Nature EMM, 2018. Single-cell RNA sequencing technologies and bioinformatics pipelines](https://www.nature.com/articles/s12276-018-0071-8)
* [Vieth et al, Nature Comm, 2019. A systematic evaluation of single cell RNA-seq analysis pipelines](https://www.nature.com/articles/s41467-019-12266-7)
* [Luecken and Theis, Mol Sys Biol, 2019. Current best practices in single‐cell RNA‐seq analysis: a tutorial](https://www.embopress.org/doi/10.15252/msb.20188746) 
* [Cheng et al., Frontier in Genetics, 2019. Single-Cell RNA-Seq Technologies and Related Computational Data Analysis](https://www.frontiersin.org/articles/10.3389/fgene.2019.00317/full)
* [Haque et al., Genome Medicine, 2017. A practical guide to single-cell RNA-sequencing for biomedical research and clinical applications](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0467-4)
* [Lafzi et al., Nature Protocols, 2018. Tutorial: guidelines for the experimental design of single-cell RNA sequencing studies](https://www.nature.com/articles/s41596-018-0073-y)
* [Soneson and Robinson, Nature Methods, 2018. Bias, robustness and scalability in single-cell differential expression analysis](https://www.nature.com/articles/nmeth.4612)

### Biclustering

* [Pontes et al., Journal of Biological Informatics, 2015. Biclustering on expression data: A review](https://reader.elsevier.com/reader/sd/pii/S1532046415001380?token=48BFDC8133FBF0BA66678534D47DD8A14AADA2D663504FDBCFFA72450181161CBF4D5638C0DB37919F4CAE75AF0A638F)


### Challenges and open questions
* [Lahnemann et al., Genome Biology, 2020. Eleven grand challenges in single-cell data science](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1926-6)

### Multiomics integration
* [Gorin et al., Genome Biology, 2020. Protein velocity and acceleration from single-cell multiomics experiments](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1945-3)

<br>

Let's divide the scRNA-seq workflow in two parts:

1. Generation of cells and the expression matrix

2. Computational analysis of the expression matrix

<br><br><br>

## 1- Generation of cells and the expression matrix


### Experimental design

* [Domenico et al., STAR Protocols, 2020. Optimized workflow for single-cell transcriptomics on infectious diseases including COVID-19](https://www.sciencedirect.com/science/article/pii/S2666166720302203)
* [Molin and Camillo, Briefings in Bioinformatics, 2019. How to design a single-cell RNA-sequencing experiment: pitfalls, challenges and perspectives](https://academic.oup.com/bib/article/20/4/1384/4831233)
* [Zhang et al., Nature Communications, 2020. Determining sequencing depth in a single-cell RNA-seq experiment](https://www.nature.com/articles/s41467-020-14482-y)
* [Lafzi et al., Nature Protocols, 2018. Tutorial: guidelines for the experimental design of single-cell RNA sequencing studies](https://www.nature.com/articles/s41596-018-0073-y)


### 2- Computational analysis of the expression matrix

### Droplet identification
* [Lun et al., Genome Biology, 2019. EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y)


### Doublet inference
* [Wolock et al., Cell Systems, 2019. Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data](https://www.sciencedirect.com/science/article/pii/S2405471218304745)
* [McGinnis et al., Cell Systems, 2019. DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors](https://www.sciencedirect.com/science/article/pii/S2405471219300730)


### Normalization
* [SCTransform package repository](https://github.com/ChristophH/sctransform)
* [Seurat SCTransform current vignette](https://satijalab.org/seurat/articles/sctransform_vignette.html)
* [Hafemeister and Satija, Genome Biology, 2019. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1)
* [(scran) Lun et al., Genome Biology, 2016. Pooling across cells to normalize single-cell RNA sequencing data with many zero counts](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0947-7)
* [Bacher et al., Nature Methods, 2017. SCnorm: robust normalization of single-cell RNA-seq data](http://pages.cs.wisc.edu/~newton/papers/publications/nmeth.4263.pdf)
* [Tang et al., Bioinformatics, 2020. bayNorm: Bayesian gene expression recovery, imputation and normalization for single-cell RNA-sequencing data](https://academic.oup.com/bioinformatics/article/36/4/1174/5581401)


### Batch effect correction
Batch effect are technical noise such as the time the experiment was done, the person carrying out the experiment or differences in reagents, etc. To correclty apply batch effect correction to a dataset, the experiment cannot be confounded (i.e., each bach must contatin at least two biological conditions).
* [Ilicic et al., Genome Biolgy, 2016. Classification of low quality cells from single-cell RNA-seq data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4758103/)
* [Brennecke et al., Nature Methods, 2013. Accounting for technical noise in single-cell RNA-seq experiments](https://www.nature.com/articles/nmeth.2645)
* [Hildreth Robert Frost, NAR, 2020. Variance-adjusted Mahalanobis (VAM): a fast and accurate method for cell-specific gene set scoring](https://academic.oup.com/nar/article/48/16/e94/5868339)

Recommended tools as in [Tutorial: guidelines for the computational analysis of single-cell RNA sequencing data. Nature Protocols, 2021](https://pubmed.ncbi.nlm.nih.gov/33288955/). One important point is that even applying mutual-nearest-neighbors (mnn) tools between cells, one may incorrectly remove real biological signals if applied to a confounded experiment.

* [mmnCorrect. Haghverdi et atl., Nature Biotech, 2018. Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors](https://www.nature.com/articles/nbt.4091)
* [Seurat. Butler et al., Nature Biotech., 2018. Integrating single-cell transcriptomic data across different conditions, technologies, and species](https://www.nature.com/articles/nbt.4096)


### Data imputation


* [Talwar et al., Scientific Reports, 2018. AutoImpute: Autoencoder based imputation of single-cell RNA-seq data](https://www.nature.com/articles/s41598-018-34688-x)
* [Andrews and Hemberg, F1000Research, 2019. False signals induced by single-cell imputation](https://f1000research.com/articles/7-1740)
* [Huang et al., Nat. Methods, 2018. SAVER: gene expression recovery for single-cell RNA sequencing](https://www.nature.com/articles/s41592-018-0033-z)
* [Peng et al., Genome Biology, 2019. SCRABBLE: single-cell RNA-seq imputation constrained by bulk RNA-seq data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1681-8)
* [Patruno et at., Briefings in Bioinformatics, 2020. A review of computational strategies for denoising and imputation of single-cell transcriptomic data](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbaa222/5916940)
* [Tang et al., Bioinformatics, 2020. bayNorm: Bayesian gene expression recovery, imputation and normalization for single-cell RNA-sequencing data](https://academic.oup.com/bioinformatics/article/36/4/1174/5581401)




### 


### Co-expression
* [Farahbod M. and Pavlidis P., Genome Research, 2020. Untangling the effects of cellular composition
on coexpression analysis](https://genome.cshlp.org/content/30/6/849.full.pdf+html)



### Deconvolution

* [Dong et. al, Briefings in Bioinformatics, 2020. SCDC: bulk gene expression deconvolution by multiple single-cell RNA sequencing references.](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbz166/5699815)


### List of collections

* [Packages collected by Sean Davis](https://github.com/seandavi/awesome-single-cell)
* [Collection from CrazyHotTommy](https://github.com/crazyhottommy/RNA-seq-analysis#single-cell-rna-seq)
* [scRNA-tools](https://www.scrna-tools.org/)
* [A collection of 55 trajectory inference methods](https://github.com/dynverse/dynmethods#list-of-included-methods)
* [Anthony Gitter. Single-cell RNA-seq pseudotime estimation algorithms. 2018. doi:10.5281/zenodo.1297422](https://github.com/agitter/single-cell-pseudotime)


### Data resources
* [The Human Cell Atlas](https://data.humancellatlas.org/analyze/portals/single-cell-expression-atlas)
* [Datasets from 10XGenomics](https://support.10xgenomics.com/single-cell-gene-expression/datasets)
* [Manually curated database from Svensson's lab](http://www.nxn.se/single-cell-studies/)
* [Human and mouse data from Hemberg's lab](https://hemberg-lab.github.io/scRNA.seq.datasets/)

### List of packages
* [scRNA-tools: A collection of scRNA-seq tools](https://www.scrna-tools.org/)
* [Seurat3](https://satijalab.org/seurat/)
* [Monocle3](https://github.com/cole-trapnell-lab/monocle3)
* [Velocity](http://velocyto.org/)
* [future: Unified Parallel and Distributed Processing in R for Everyone](https://github.com/HenrikBengtsson/future/tree/master)

### Marker genes identification
* [Lorenzo et al., NAR, 2020. Combining single-cell rna-sequencing with a molecular atlas unveils new markers for caenorhabditis elegans neurons classes](https://academic.oup.com/nar/article/48/13/7119/5857708)
* [Descartes: human gene expression during development. Database for cell identification](https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/dataset/liver)
* [Azimuth: App for reference-based single cell analysis](https://azimuth.hubmapconsortium.org/)
* [SingleR: bioconductor package to annotate single-cell RNA-seq data](https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html)

### Big projects, Atlas
* [The Human Cell Atlas](https://data.humancellatlas.org/analyze/portals/single-cell-expression-atlas)
* [Saunders et al, Cell, 2018. Molecular Diversity and Specializations among the Cells of the Adult Mouse Brain](https://www.sciencedirect.com/science/article/pii/S0092867418309553?via%3Dihub#sec4)


### Comparison with bulk RNA-seq
* [Bacher and Kendziorski, Genome Biology, 2016. Design and computational analysis of
single-cell RNA-sequencing experiments](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-016-0927-y)



