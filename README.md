# MicrobiomeDiversityKit

R functions for microbiome beta diversity analysis, including data preprocessing, normalization, ADONIS/PERMANOVA testing, and visualization.

## Requirements
```r
# CRAN packages
install.packages(c("vegan", "ade4", "cluster", "psych", "ggplot2", 
                   "ggrepel", "ggpubr", "gridExtra", "viridis",
                   "RColorBrewer", "wesanderson", "corrplot", 
                   "dunn.test", "imputeTS"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("phyloseq", "DESeq2", "microbiome"))

# GitHub packages
devtools::install_github("ggloor/CoDaSeq")
devtools::install_github("umerijaz/microbiomeSeq")
devtools::install_github("jfq3/ggordiplots")
devtools::install_github("jchen1981/GMPR")
```

## Usage

Load the functions:
```r
# Set seed for reproducibility
set.seed(12345)

# Source the functions
source("./Functions/beta_diver_functions.R")
source("./Functions/plot_beta_diver_functions.R")
```

---

## Data Input Functions

### filter_otu_table()

Filter a phyloseq object by sample subset and prevalence threshold.
```r
taxa_matrix <- filter_otu_table(
  in_phylo = phyloseq_object,
  Samples = sample_names_vector,
  prev = 0.2,              # Minimum prevalence (20%)
  Tax_level = "Genus"      # Options: "Genus", "Species", "SGB", "none"
)
```

### read.infile.data()

Read abundance data from multiple formats.
```r
data_list <- read.infile.data(
  infileGenus = "genus_table.tsv",   # For table format
  format = "table",                   # Options: "table", "phyloseq", "DADA2"
  min.rar = TRUE,                     # Use minimum sample depth for rarefaction
  min.num.seqs = 10000                # Minimum number of sequences
)

# Returns a list with:
# - data_list$Genus: Raw counts
# - data_list$Genus.rar: Rarefied counts
```

### dada2_to_phyloseq()

Convert DADA2 output files to a phyloseq object.
```r
physeq <- dada2_to_phyloseq(
  abundance.file = "seqtab_nochim.txt",
  tax.file = "tax_silva.txt",
  sep = " "
)
```

---

## Data Transformation

### CLR.transformation()

Centered Log-Ratio transformation for compositional data. Supports two methods depending on data sparsity.
```r
# Standard CLR (for low sparsity data)
clr_matrix <- CLR.transformation(
  in.table = abundance_table,
  min.reads = 0,
  min.prop = 0.001,
  GMPR_aproximation = FALSE
)

# GMPR-based normalization (for high sparsity data)
gmpr_matrix <- CLR.transformation(
  in.table = abundance_table,
  GMPR_aproximation = TRUE,
  GMPR_aproximation_log = FALSE,
  GMPR_aproximation_square_root = TRUE,
  min_ct = 2,
  intersect_no = 4
)
```

---

## Filtering Functions

### filter.table()

Filter taxa by minimum mean relative abundance.
```r
filtered_table <- filter.table(
  in.table = abundance_table,
  min.percentage = 0.01    # Minimum 0.01% mean abundance
)
```

### filter.prevalence()

Filter taxa by maximum prevalence of zeros.
```r
filtered_table <- filter.prevalence(
  in.table = abundance_table,
  max.prev = 0.8           # Maximum 80% zeros allowed
)
```

### filter_low_prevalence()

Filter taxa by zero percentage, optionally within groups.
```r
# Global filtering
filtered_table <- filter_low_prevalence(
  in.table = abundance_table,
  max_percentage_0 = 80
)

# Group-specific filtering
filtered_table <- filter_low_prevalence(
  in.table = abundance_table,
  Categories_groups = metadata["Group"],
  max_percentage_0 = 80
)
```

---

## Beta Diversity Analysis

### ADONIS_func()

Perform PERMANOVA (ADONIS) test for multiple variables.
```r
adonis_results <- ADONIS_func(
  in.matrix = taxa_matrix,
  Distance = "bray",          # Distance metric: "bray", "euclidean", etc.
  in.Metadata = metadata,
  prefix = "output",
  permutations = 1000
)

# Returns data.frame with columns:
# - Variable, Fmodel, R2, p.value, BH.adj.p.value, N
```

### capscale_cum_variance()

Calculate cumulative variance explained by metadata variables using db-RDA.
```r
variance_results <- capscale_cum_variance(
  in.Metadata = metadata,
  in.matrix = taxa_matrix,
  Distance = "bray",
  prefix = "output",
  permutations = 1000,
  adj.pval.cutof = 0.1
)
```

### non_redundant_variance()

Calculate non-redundant variance for each variable.
```r
nr_variance <- non_redundant_variance(
  in.Metadata = metadata,
  in.matrix = taxa_matrix,
  Distance = "bray"
)

# Returns list with:
# - nr: Non-redundant variance table
# - nr_sig: Formatted for plotting
```

---

## Ordination Functions

### vegan_PCoA_envfit()

Perform PCoA with environmental fitting.
```r
pcoa_results <- vegan_PCoA_envfit(
  in.matrix = taxa_matrix,
  Metadata2enfit = metadata,
  choices = c(1, 2),
  distance = "euclidean",
  perm = 999
)

# Returns list with:
# - PCoA: Ordination object
# - xlab, ylab: Axis labels with variance explained
# - fit: envfit results
# - df_ord: Ordination scores
```

---

## Plotting Functions

### PCoA_grapper()

Create PCoA plot with ADONIS statistics.
```r
pcoa_plot <- PCoA_grapper(
  Var = "Treatment",
  tempMetadata = metadata,
  table.ADONIS = adonis_results,
  taxa_matrix = taxa_matrix,
  distance = "bray",
  Ellipses = TRUE,
  Label = TRUE,
  colors = c("Control" = "blue", "Treatment" = "red")
)
```

### total_variance_donut_plot()

Create donut plot showing total variance explained.
```r
donut_plot <- total_variance_donut_plot(capscale_results)
```

---

## Complete Workflow Example
```r
# 1. Load functions
source("./Functions/beta_diver_functions.R")
source("./Functions/plot_beta_diver_functions.R")

# 2. Load data
physeq <- readRDS("phyloseq_object.rds")
metadata <- read.table("metadata.tsv", header = TRUE, row.names = 1, sep = "\t")

# 3. Filter and prepare data
taxa_matrix <- filter_otu_table(
  in_phylo = physeq,
  Samples = rownames(metadata),
  prev = 0.1,
  Tax_level = "Genus"
)

# 4. Transform data (CLR)
taxa_clr <- CLR.transformation(in.table = t(taxa_matrix))

# 5. Run ADONIS
adonis_results <- ADONIS_func(
  in.matrix = taxa_clr,
  Distance = "euclidean",
  in.Metadata = metadata,
  prefix = "beta_diversity"
)

# 6. Create PCoA plots for significant variables
significant_vars <- adonis_results$Variable[adonis_results$BH.adj.p.value < 0.05]

pdf("PCoA_plots.pdf", width = 10, height = 8)
for (var in significant_vars) {
  p <- PCoA_grapper(
    Var = var,
    tempMetadata = metadata,
    table.ADONIS = adonis_results,
    taxa_matrix = taxa_clr,
    distance = "euclidean",
    Ellipses = TRUE
  )
  print(p)
}
dev.off()
```

---

## Functions Reference

| Function | Description |
|----------|-------------|
| `filter_otu_table()` | Filter phyloseq object by samples and prevalence |
| `dada2_to_phyloseq()` | Convert DADA2 output to phyloseq |
| `dada2genusTable()` | Create genus table from DADA2 with rarefaction |
| `read.infile.data()` | Read abundance data from multiple formats |
| `CLR.transformation()` | CLR or GMPR-based normalization |
| `filter.table()` | Filter by minimum mean abundance |
| `filter.prevalence()` | Filter by maximum zero prevalence |
| `filter_low_prevalence()` | Filter by zero percentage (global or by group) |
| `estimate0.min()` | Estimate minimum values for zero replacement |
| `adonis.func()` | Core ADONIS function |
| `ADONIS_func()` | ADONIS wrapper with output |
| `capscale_cum_variance()` | Cumulative variance with db-RDA |
| `non_redundant_variance()` | Non-redundant variance calculation |
| `vegan_PCoA_envfit()` | PCoA with environmental fitting |
| `ggplot_envfit()` | Create ggplot from envfit results |
| `ord_labels()` | Generate ordination axis labels |
| `PCoA_grapper()` | PCoA plot with ADONIS stats |
| `total_variance_donut_plot()` | Donut plot for variance explained |
| `median_per_condition()` | Calculate median abundance per condition |

---

## Citation

If you use this package, please cite:

> Vázquez-Castellanos JF. MicrobiomeDiversityKit: R functions for microbiome beta diversity analysis. GitHub repository: https://github.com/jorgevazcast/microbiomediversitykit
