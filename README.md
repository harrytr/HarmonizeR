# HARMONIZER

## Information
CCLE ML &amp; Network Optimization

HARMONIZER is an Shiny app developed in R and combines Python scripts to facilitate in depth analysis of the Cancel Cell Line Encyclopedia [CCLE](https://depmap.org/portal/).

It mainly enables the following avenues of analysis:

1) Exploratory analysis using clustering, PCA, and statistics by subsetting on a specific disease and other clinical features (gene, cell lines, mutations etc).
2) Machine Learning Classification based on the [RENOIR](https://github.com/alebarberis/renoir) package - it enables the use of a multitude of learning methods for binomial and multinomial classification using CCLE datasets across different settings
3) Network Optimization and Reconstruction using the [CARNIVAL](https://github.com/saezlab/carnival) package, and additionally performing network analytics (such as community detection, centralities etc.)

## Getting Started

Download the three required  input ".csv" files from the link [here](https://depmap.org/portal/data_page/?tab=currentRelease).

You will need : 

- OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv
- OmicsSomaticMutations.csv
- OmicsCNGene.csv

Put the renamed files into a folder inside /inputs/, with the name of the releaase of these files. For example "/inputs/24Q2").

### Installing

Clone this Repo locally. Load the Shiny App by running the app.R from Rstudio. You will need R.4.5.0 and Python 3.11.7 (tested and working ; anything more recent might not work). 

## Running HARMONIZER

```R
runApp('Harmonizer.R')
```

It should look like this:

<img width="904" alt="Screenshot 2024-11-22 at 14 45 02" src="https://github.com/user-attachments/assets/4a49ddf1-ec56-489c-89ef-ccc36d08d17e">

In the UI from the bottom right check the box "Read new CCLE version" and use the "24Q2" in the box to let the platform know where to find the files (it will look now then in /inputs/24Q2".

Once the files have been read, a "24Q2.RData" will be created with all the matrices included, so you don't have to re-read the .CSV files each time you run the platform but just load the .Rdata file (much faster) in memory.

### Prerequisites

CARNIVAL requires the interactive version of IBM Cplex The IBM ILOG Cplex is freely available through Academic Initiative [here](https://www.ibm.com/products/ilog-cplex-optimization-studio?S_PKG=CoG&cm_mmc=Search_Google-_-Data+Science_Data+Science-_-WW_IDA-_-+IBM++CPLEX_Broad_CoG&cm_mmca1=000000RE&cm_mmca2=10000668&cm_mmca7=9041989&cm_mmca8=kwd-412296208719&cm_mmca9=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_&cm_mmca10=267798126431&cm_mmca11=b&mkwid=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_|470|135655&cvosrc=ppc.google.%2Bibm%20%2Bcplex&cvo_campaign=000000RE&cvo_crid=267798126431&Matchtype=b&gclid=Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB). The [CBC](https://projects.coin-or.org/Cbc) solver is open source and freely available for any user. 


## License

Distributed under the GNU GPLv2 License. See accompanying file.

## References

[Triantafyllidis et al.](https://www.biorxiv.org/content/10.1101/2024.12.28.630070v3.abstract):
> Causality-aware graph neural networks for functional stratification and phenotype prediction at scale, Charalampos P. Triantafyllidis, Ricardo Aguas, bioRxiv 2024.12.28.630070; doi: https://doi.org/10.1101/2024.12.28.630070

Harmonizer is highly updated version of the initial release of a similar platform from:

[Triantafyllidis et al.]([https://pubs.rsc.org/en/content/articlehtml/2015/ib/c4ib00294f](https://www.cell.com/iscience/fulltext/S2589-0042(23)02368-4?uuid=uuid%3A7b7fb9c3-4515-46e2-8e6e-fe22489b11b9)):

> Triantafyllidis, C.P. et al. (2023) ‘A machine learning and directed network optimization approach to uncover TP53 regulatory patterns’, iScience, 26(12). Available at: https://doi.org/10.1016/j.isci.2023.108291.

Credits also go to [Enio Gjerga](https://scholar.google.com/citations?user=IzQPpf0AAAAJ&hl=en), for the functions implemented (as indicated within each function).

The main differences are:

- improved data compatibility (now working with latest versions of CCLE (2024+)
- Removed TCGA support
- Improved community detection on optimized networks (directed supported)
- Improved gene set signature extraction from the networks as now included state of gene (up/down regulated)
- Improved reading of .DOT figure optimized CARNIVAL network, now being read using Python (PyDoT)
- Multiple user options on both optimizing, or classifying
- Data are being drawn now directly from the dataset and populated real-time on the UI
- Stratification of analyses now available across user selected features (for both classification or optimization)
- Updated generated figures (heatmaps, corrplots etc.)
- Exported optimized networks now also in .graphml format and visNetwork format
- Deep Learning integration


