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

Clone this Repo locally. Load the Shiny App by running the app.R from Rstudio. You will need R.4.2.2 (tested and working ; anything more recent might not work). 

## Running HARMONIZER

```R
runApp('Harmonizer.R')
```

It should look like this:

![harmonizer_UI](https://github.com/user-attachments/assets/30459019-cdc7-4378-a8fe-5850c83b7c76)

In the UI from the bottom right check the box "Read new CCLE version" and use the "24Q2" in the box to let the platform know where to find the files (it will look now then in /inputs/24Q2".

Once the files have been read, a "24Q2.RData" will be created with all the matrices included, so you don't have to re-read the .CSV files each time you run the platform but just load the .Rdata file (much faster) in memory.

### Prerequisites

CARNIVAL requires the interactive version of IBM Cplex The IBM ILOG Cplex is freely available through Academic Initiative [here](https://www.ibm.com/products/ilog-cplex-optimization-studio?S_PKG=CoG&cm_mmc=Search_Google-_-Data+Science_Data+Science-_-WW_IDA-_-+IBM++CPLEX_Broad_CoG&cm_mmca1=000000RE&cm_mmca2=10000668&cm_mmca7=9041989&cm_mmca8=kwd-412296208719&cm_mmca9=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_&cm_mmca10=267798126431&cm_mmca11=b&mkwid=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_|470|135655&cvosrc=ppc.google.%2Bibm%20%2Bcplex&cvo_campaign=000000RE&cvo_crid=267798126431&Matchtype=b&gclid=Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB). The [CBC](https://projects.coin-or.org/Cbc) solver is open source and freely available for any user. 


## License

Distributed under the GNU GPLv2 License. See accompanying file.

## References

Harmonizer is highly updated version of the initial release of a similar platform from:

[Triantafyllidis et al.]([https://pubs.rsc.org/en/content/articlehtml/2015/ib/c4ib00294f](https://www.cell.com/iscience/fulltext/S2589-0042(23)02368-4?uuid=uuid%3A7b7fb9c3-4515-46e2-8e6e-fe22489b11b9)):

> Triantafyllidis, C.P. et al. (2023) ‘A machine learning and directed network optimization approach to uncover TP53 regulatory patterns’, iScience, 26(12). Available at: https://doi.org/10.1016/j.isci.2023.108291.

Credits also go to Enio Gjerga, for the functions implemented (as indicated within each function).

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

## Examples of Results

![Plots_cancers_TNBC_24Q2 pdf-image-002](https://github.com/user-attachments/assets/a5a9e624-91b0-4837-ae09-7edec9fc50c7)
![Plots_cancers_TNBC_24Q2 pdf-image-003](https://github.com/user-attachments/assets/b2c61a23-79fa-46bc-bf3f-7fa3e6e0bec9)
![Plots_cancers_TNBC_24Q2 pdf-image-008](https://github.com/user-attachments/assets/cb3a57ec-9cfb-4c1e-af4d-4703557c2fe9)
![Plots_cancers_TNBC_24Q2 pdf-image-010](https://github.com/user-attachments/assets/afbf00da-002c-4b72-895a-3f65029fc28f)
![Plots_cancers_TNBC_24Q2 pdf-image-012](https://github.com/user-attachments/assets/5be0d6ad-7cb7-41c3-bc4b-dde8d6ab10e9)
![Plots_cancers_TNBC_24Q2 pdf-image-017](https://github.com/user-attachments/assets/63b7e2ef-ea01-416b-9d22-17a21e2cff82)
![Plots_cancers_TNBC_24Q2 pdf-image-016](https://github.com/user-attachments/assets/f9925b51-0344-4383-b283-915ebb421378)
![net](https://github.com/user-attachments/assets/d2422c61-2984-4030-850d-aae43e381960)
![net2](https://github.com/user-attachments/assets/88117405-3259-4c21-a021-2472cf3ad37c)


