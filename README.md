## Human Microbiome Research Manual

> Materials of the repository were prepared on the basis of materials of the educational intensive from the Blastim<br>
> It can be used as a helpful repository with cheat-sheets for 16s rRNA amplicon based metagenome studies.

## Introduction
The gut microbiome plays a crucial role in human health, influencing immune responses, metabolism, and disease susceptibility. Crohn’s Disease, a chronic inflammatory bowel disease, has been associated with microbial dysbiosis. Understanding the microbial alterations in Crohn’s Disease can provide insights into disease mechanisms and potential therapeutic targets.

## Abstract
In this study, the gut microbiome composition in individuals with Crohn’s Disease (CD) and Healthy Controls (HC) was investigated. The goal was to identify taxonomic shifts, potential functional implications, and associations with disease. High-throughput sequencing data and various bioinformatics tools were employed to analyze microbial diversity and abundance.

## Materials and Methods

Data obtained from the article [«_A microbial signature for Crohn's disease_»](https://pubmed.ncbi.nlm.nih.gov/28179361/) was used. These data include information on microbial samples from healthy individuals and patients with Crohn's disease.

### Data Preprocessing
- Raw 16S rRNA gene sequencing data were obtained from a cohort of CD patients and HC.
- Metadata, including sample information and diagnosis labels, were organized and linked to the corresponding samples.

### Filtering Rare Taxa
- Samples with zero counts or low coverage were filtered out.
- The final dataset included microbial counts for each sample.
- Taxa that were present in at least 30% of the samples were retaines.

### Coverage Analysis
- The sequencing depth and coverage for each sample was evaluated.
- Samples with sufficient coverage were considered for subsequent analyses.

### Alpha Diversity Metrics
- Computed alpha diversity indices (`Shannon index`) with the R package `MicrobiomeR`. The index was calculated after thinning 5 times to 19000 reads and making a mean, then a boxplot was constructed.
- Assessed statistical significance using non-parametric tests (`Wilcoxon rank-sum test`).

### Beta Diversity Analysis
- Conducted Principal Coordinate Analysis (PCoA) based on Bray-Curtis dissimilarity.
- Tested group differences using permutational multivariate analysis of variance (`adonis`).

### Software and Packages
- R version 4.2.1
- Utilized `data.table`, `openxlsx`, `MicrobeR`, `ggplot2`, `zCompositions`, `NearestBalance`, `GUniFrac`, `vegan`, `ape` and `selbal` packages.

## Results

The abundance of various bacterial taxa in different samples is illustrated in the Figure 1. The samples are divided into two categories: CD and HC. Each vertical bar represents a specific sample, and the colors within each bar correspond to different types of bacteria. The height of each color segment indicates the relative abundance of that bacterial group within the sample. CD samples show distinct microbial composition compared to HC samples.

<div style='justify-content: center'>
<img src="https://github.com/iliapopov17/Human-Microbiome-Research-User-Manual/blob/main/img/microbiome_barplot.jpg" align='center', width="50%">
</div>

_Figure 1. Microbiome Barplot._

The coverage quality of the samples is represented in Figure 2. The x-axis shows the number of samples, and the y-axis represents the count of samples falling into each coverage range.  The quality of sequencing can be regarded as sufficient.

<div style='justify-content: center'>
<img src="https://github.com/iliapopov17/Human-Microbiome-Research-User-Manual/blob/main/img/coverage_quality.jpg" align='center', width="50%">
</div>

_Figure 2. Coverage quality._

The post-filtration coverage of the samples is represented in Figure 3 and in Table 1.

<div style='justify-content: center'>
<img src="https://github.com/iliapopov17/Human-Microbiome-Research-User-Manual/blob/main/img/post-filtration%20coverage.jpg" align='center', width="50%">
</div>

_Figure 3. Post-filtration coverage._

|   |N microbes|Minimum coverage|
|---|----------|----------------|
|Before|210|19412|
|After|89|12981|

_Table 1. Mocrobial samples and coverage before and after filtration._

The proportion of microbes remaining in the assay across different samples is illustrated in Figure 4.

<div style='justify-content: center'>
<img src="https://github.com/iliapopov17/Human-Microbiome-Research-User-Manual/blob/main/img/remaining%20proportion.jpg" align='center', width="50%">
</div>

_Figure 4. Remaining proportion._

The relative abundance of various microbial organisms in different samples is visualized in Figure 5. The x-axis represents the sample categories (CD and HC), while the y-axis lists the microbial taxa. The color intensity within each cell indicates the abundance percentage, with darker colors representing higher relative abundance.

<div style='justify-content: center'>
<img src="https://github.com/iliapopov17/Human-Microbiome-Research-User-Manual/blob/main/img/relative%20abundance%20of%20the%20main%20microbial%20organisms.jpg" align='center', width="50%">
</div>

_Figure 5. Relative abundance of the main microbial organisms._

The boxplot in Figure 6. compares the alpha diversity (measured using the `Shannon index`) between two groups: CD and HC. The x-axis represents the group labels, while the y-axis shows the `Shannon index` values. The red box corresponds to CD samples, and the blue box corresponds to HC samples.

The _p_-value of `1.795217e-09` indicates a significant difference in alpha diversity between the two groups. CD samples exhibit lower alpha diversity compared to HC samples.

<div style='justify-content: center'>
<img src="https://github.com/iliapopov17/Human-Microbiome-Research-User-Manual/blob/main/img/alpha%20diversity.jpg" align='center', width="50%">
</div>

_Figure 6. Alpha diversity._

The scatter plot in Figure 7 visualizes the beta diversity using the first two axes of a principal coordinate analysis (`PCoA`). Each point represents a sample, and the colors differentiate between the CD and HC groups. The percentage of variance explained by each axis is indicated in the axis labels.

The statistical analysis using the `adonis2` test shows that the difference in beta diversity between the two groups is statistically significant (p-value = 0.001).

<div style='justify-content: center'>
<img src="https://github.com/iliapopov17/Human-Microbiome-Research-User-Manual/blob/main/img/beta%20diversity.jpg" align='center', width="50%">
</div>

_Figure 7. Beta diversity._

The heatmap in Figure 8 visualizes the relative abundance of microbial taxa in different samples, split into two categories: “health-related” and “disease-related.” The color intensity represents the abundance percentage, with darker colors indicating higher abundance.

<div style='justify-content: center'>
<img src="https://github.com/iliapopov17/Human-Microbiome-Research-User-Manual/blob/main/img/heatmap%20with%20split.jpg" align='center', width="50%">
</div>

_Figure 8. Heatmap with split._

The bar plot in Figure 9 illustrates the mean difference between CD and HC individuals for various microbial taxa. Each bar represents a specific taxon, and the x-axis shows the CLR(v) values. The dark blue bars indicate the direction and magnitude of the difference in abundance between the two groups.

<div style='justify-content: center'>
<img src="https://github.com/iliapopov17/Human-Microbiome-Research-User-Manual/blob/main/img/mean%20difference.jpg" align='center', width="50%">
</div>

_Figure 9. Mean difference._

The bar plot in Figure 10 illustrates the approximate, simplified difference between CD and HC individuals for various microbial taxa. Each bar represents a specific taxon, and the x-axis shows the CLR(b) values. The dark blue bars indicate the direction and magnitude of the difference in abundance between the two groups.

<div style='justify-content: center'>
<img src="https://github.com/iliapopov17/Human-Microbiome-Research-User-Manual/blob/main/img/approximate%20difference.jpg" align='center', width="50%">
</div>

_Figure 10. Approximate difference._

The boxplot in Figure 11 illustrates the balance value in each sample for two different categories: CD and HC. The x-axis represents the group labels, while the y-axis shows the balance values. The p-value of `2.621864e-17` indicates a significant difference in balance between the two groups.

<div style='justify-content: center'>
<img src="https://github.com/iliapopov17/Human-Microbiome-Research-User-Manual/blob/main/img/balance%20value.jpg" align='center', width="50%">
</div>

_Figure 11. Balance value._

## Discussion

### Taxonomic Shifts
The heatmap (Figure 8) revealed distinct patterns of microbial abundance between CD and HC samples. Taxa clustering within each group suggests shared ecological niches or functional roles.
Notably, certain taxa may be enriched in CD (e.g., Firmicutes) or depleted (e.g., Bacteroidetes). These shifts could impact host health and immune responses.
The “Approximate Difference” plot (Figure 10) highlights specific taxa with significant differences. For instance, some taxa exhibit higher CLR(b) values in CD, while others are more abundant in HC.

### Functional Implications
Taxa with higher abundance in CD might contribute to inflammation or immune dysregulation. For example, increased Proteobacteria may correlate with disease severity.
Conversely, taxa enriched in HC (e.g., Bacteroides) could play protective roles by promoting gut barrier integrity or modulating immune responses.
Functional metagenomic analyses are warranted to explore specific pathways associated with these taxa.

### Balance and Dysbiosis
The “Balance Value” boxplot (Figure 11) indicates differences in microbial balance between CD and HC. Dysbiosis (imbalance) may disrupt homeostasis and contribute to disease pathogenesis.
CD samples exhibit lower balance values, suggesting altered microbial interactions and potential loss of beneficial cross-feeding relationships.

## Data availability

In this repository you can find:

1) `Human Microbiome Research.Rmd` - contains the whole script of the study
2) `Human Microbiome Research.html` - contains the report on the study
3) `data_Crohns_disease` - contains the data (`counts.csv`, `metadata.csv`)
