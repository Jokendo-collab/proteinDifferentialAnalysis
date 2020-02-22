
#Install the required library
#BiocManager::install("DEP")

#Set working directory
setwd("/home/javan/Desktop/dlovuData")
# Loading a package required for data handling
library("dplyr")
library("DEP")

# For LFQ analysis
run_app("LFQ")

# Load the data

data <- UbiLength

# We filter for contaminant proteins and decoy database hits, which are indicated by "+" in the columns "Potential.contaminants" and "Reverse", respectively. 
data <- filter(data, Reverse != "+", Potential.contaminant != "+")

#The dataframe has the following dimensions

dim(data)

colnames(data)

#Data preprocessing

# Are there any duplicated gene names?
data$Gene.names %>% duplicated() %>% any()

# Make a table of duplicated gene names
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)

# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

# Are there any duplicated names?
data$name %>% duplicated() %>% any()


#Generate a SummarizedExperiment object

# Generate a SummarizedExperiment object using an experimental design
LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
experimental_design <- UbiLength_ExpDesign
data_se <- make_se(data_unique, LFQ_columns, experimental_design)

# Generate a SummarizedExperiment object by parsing condition information from the column names
LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
data_se_parsed <- make_se_parse(data_unique, LFQ_columns)

# Let's have a look at the SummarizedExperiment object
data_se


#Filter on missing values
# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

# Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 0)

# Less stringent filtering:
# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt2 <- filter_missval(data_se, thr = 1)

# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

#Normalization
#The data is background corrected and normalized by variance stabilizing transformation (vsn).

# Normalize the data
data_norm <- normalize_vsn(data_filt)

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)


#Impute the data for the missing valuee

# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

# To check whether missing values are biased to lower intense proteins, the densities and 
# cumulative fractions are plotted for proteins with and without missing values.

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)

# All possible imputation methods are printed in an error, if an invalid function name is given.
impute(data_norm, fun = "")

# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)

# Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)

#Differential enrichment analysis

# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "Ctrl")

# Test all possible comparisons of samples
data_diff_all_contrasts <- test_diff(data_imp, type = "all")


# Test manually defined comparisons
data_diff_manual <- test_diff(data_imp, type = "manual", 
                              test = c("Ubi4_vs_Ctrl", "Ubi6_vs_Ctrl"))

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

#Visualization of the previuos results

#PCA plot
# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

#Correlation matrix
#A correlation matrix can be plotted as a heatmap, to visualize the Pearson correlations between the different samples.

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Spectral")

#Heatmap of all significant proteins
# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))

# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)

#Volcano plots of specific contrasts
# Plot a volcano plot for the contrast "Ubi6 vs Ctrl""
plot_volcano(dep, contrast = "Ubi6_vs_Ctrl", label_size = 2, add_names = TRUE)

#Barplot of protein of interest
# Plot a barplot for USP15 and IKBKG
plot_single(dep, proteins = c("USP15", "IKBKG"))

# Plot a barplot for the protein USP15 with the data centered
plot_single(dep, proteins = "USP15", type = "centered")

#Frequency plot of significant proteins and overlap of conditions

# Plot a frequency plot of significant proteins for the different conditions
plot_cond(dep)



# The data is provided with the package 
data <- UbiLength
experimental_design <- UbiLength_ExpDesign

# The wrapper function performs the full analysis
data_results <- LFQ(data, experimental_design, fun = "MinProb", 
                    type = "control", control = "Ctrl", alpha = 0.05, lfc = 1)

#Prints the entire results
report(data_results)






















