library(limma)
library(ggrepel)
library(dplyr)
library(gplots)
library(openxlsx) # save excel file
library(ggplot2)
library(devtools)
library(ggvenn)
library(reshape2)
library(FactoMineR) # til PCA
library(factoextra) # visualisering PCA (fviz)
# library(promor) # prediction model
library(tidyr) # unite function
library(ggpmisc)
library(plotly)

# Project structure:
# FALL_study/
# ├── data/      # Input data
# ├── scripts/   # Analysis scripts
# └── output/    # Figures and output files


########################## Plasma ##############################################################################

###### PCA ######

#### from http://www.sthda.com/english/articles/22-principal-component-methods-videos/65-pca-in-r-using-factominer-quick-scripts-and-videos/ ###

# import baseline characteristics from "merged" document created with "p_baseline_characteristics.R"

pca_p <- read.csv("../data/plasma_characteristics.csv")

pca_p$group <- factor(pca_p$group)
table(pca_p$group)
pca_p$group <- factor(pca_p$group, labels=c("Healthy","MASLD","PBC")) # Rename NAFLD to MASLD
table(pca_p$group)

# Split PBC into PBC with MASLD and PBC without MASLD
pca_p$group_1 <- pca_p$group 

pca_p$group_1 <- ifelse(pca_p$group == "PBC" & pca_p$con_masld == "Yes", "PBC with MASLD", 
                        ifelse(pca_p$group == "PBC" & pca_p$con_masld == "No", "PBC without MASLD", 
                               ifelse(pca_p$group == "Healthy","Healthy",
                                      ifelse(pca_p$group=="MASLD","MASLD",
                                             pca_p$group))))

pca_p$group_1 <- factor(pca_p$group_1)
table(pca_p$group_1)


### Split PBC by Ursochol treatment (urso = Yes / No) ####

# 1) Create new grouping variable in pca_p
pca_p$group_urso <- pca_p$group

pca_p$group_urso <- ifelse(
  pca_p$group == "PBC" & pca_p$urso == "yes", "PBC_urso_yes",
  ifelse(
    pca_p$group == "PBC" & pca_p$urso == "no",  "PBC_urso_no",
    as.character(pca_p$group)
  )
)

pca_p$group_urso <- factor(pca_p$group_urso)
table(pca_p$group_urso, useNA = "ifany")


# Subset dataframe
dim(pca_p)
pca_p[1:3,45:55] # proteins from column 52:633
pca_p[1:3,630:635] # proteins from column 52:633
pca_p <- pca_p %>% 
  select(all_of(c("subject", "age", "sex", "bmi", "cirrose", "con_masld", "group", "group_1", "group_urso")), 52:633)

# Correct rest of structure in PCA 
str(pca_p)
pca_p$sex <- factor(pca_p$sex)
pca_p$cirrose <- factor(pca_p$cirrose)
pca_p$con_masld <- factor(pca_p$con_masld)
pca_p$group_1   <- factor(pca_p$group_1)
pca_p$group_urso<- factor(pca_p$group_urso)


# Now to PCA
res.pca.p <- PCA(pca_p, 
               graph = FALSE, 
               quali.sup = c("subject","sex","cirrose","con_masld","group","group_1","group_urso"),
               quanti.sup = c("bmi","age"))

# Show the percentage of variances explained by each principal component
eig.val.p <- res.pca.p$eig
#
#n.pc <- nrow(eig.val.p)
#
#barplot(eig.val.p[1:n.pc, 2], 
#        names.arg = 1:n.pc,
#        main = "Variances Explained by PCs (%)",
#        xlab = "Principal Components",
#        ylab = "Percentage of variance",
#        col = "steelblue")

# Visualize the graph of individuals. Individuals with a similar profile are grouped together.

# By diagnosis

# save plot
# jpeg(file = "../output/PCA/pca_plasma_diagnosis_pbc_masld_healthy.jpeg")

fviz_pca_ind(res.pca.p, 
             geom.ind = "point", 
             col.ind = pca_p$group, 
             palette = c("#58aadb","grey","#AA0233"),
             pointsize = 3, 
             pointshape = 19, 
             invisible = "quali", 
             legend.title = "Diagnosis") +
  ggtitle("Individuals - PCA plasma") +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),  # Increase legend text size
    axis.text = element_text(size = 14),    # Increase x- and y-axis text size
    axis.title = element_text(size = 16),   # Increase x- and y-axis title size
    plot.title = element_text(size = 18),   # Increase plot title size
    legend.position = "bottom"              # Move legend to the bottom
  )

dev.off()


###### 3D PCA ######


# Extract PCA coordinates
pca_3d_p <- as.data.frame(res.pca.p$ind$coord)

# Add metadata
pca_3d_p$subject <- pca_p$subject
pca_3d_p$group <- pca_p$group

# Create 3D PCA plot
pca_3d_plot_p <- plot_ly(
  data = pca_3d_p,
  x = ~Dim.1,
  y = ~Dim.2,
  z = ~Dim.3,
  color = ~group,
  colors = c("#58aadb", "grey", "#AA0233"),
  text = ~subject,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 5)
) %>%
  layout(
    title = "3D PCA plasma",
    scene = list(
      xaxis = list(title = paste0("PC1 (", round(eig.val.p[1, 2], 1), "%)")),
      yaxis = list(title = paste0("PC2 (", round(eig.val.p[2, 2], 1), "%)")),
      zaxis = list(title = paste0("PC3 (", round(eig.val.p[3, 2], 1), "%)"))
    ),
    legend = list(title = list(text = "Diagnosis"))
  )

# View plot
pca_3d_plot_p


# save plot
# htmlwidgets::saveWidget(
#   pca_3d_plot_p,
#   "../output/PCA/pca_3d_plasma_diagnosis.html",
#   selfcontained = TRUE
# )


###### 3D PCA without subject IDs ######

pca_3d_plot_p_no_ids <- plot_ly(
  data = pca_3d_p,
  x = ~Dim.1,
  y = ~Dim.2,
  z = ~Dim.3,
  color = ~group,
  colors = c("#58aadb", "grey", "#AA0233"),
  text = ~paste("Group:", group),
  hoverinfo = "text",
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 5)
) %>%
  layout(
    title = "3D PCA plasma",
    scene = list(
      xaxis = list(title = paste0("PC1 (", round(eig.val.p[1, 2], 1), "%)")),
      yaxis = list(title = paste0("PC2 (", round(eig.val.p[2, 2], 1), "%)")),
      zaxis = list(title = paste0("PC3 (", round(eig.val.p[3, 2], 1), "%)"))
    ),
    legend = list(title = list(text = "Diagnosis"))
  )

pca_3d_plot_p_no_ids

# save plot
# htmlwidgets::saveWidget(
#   pca_3d_plot_p_no_ids,
#   "../output/PCA/pca_3d_plasma_diagnosis_no_IDs.html",
#   selfcontained = TRUE
# )


############ Identify potential outliers ###############


# PCA plot with sample IDs
# Save PCA plot with sample IDs
# jpeg(
#   file = "../output/PCA/pca_plasma_diagnosis_with_IDs.jpeg",
#   width = 3000,
#   height = 2500,
#   res = 300
# )


fviz_pca_ind(
  res.pca.p,
  geom.ind = "point",
  col.ind = pca_p$group,
  palette = c("#58aadb", "grey", "#AA0233"),
  pointsize = 3,
  pointshape = 19,
  invisible = "quali",
  legend.title = "Diagnosis",
  repel = TRUE
) +
  geom_text_repel(
    aes(label = pca_p$subject),
    size = 3,
    max.overlaps = Inf
  ) +
  ggtitle("Individuals - PCA plasma with sample IDs") +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.position = "bottom"
  )

dev.off()


# PCA coordinates
pca_coord <- as.data.frame(res.pca.p$ind$coord[,1:5])

# Mahalanobis distance
md <- mahalanobis(
  pca_coord,
  colMeans(pca_coord),
  cov(pca_coord)
)

# Add to dataframe
pca_coord$Mahalanobis <- md
pca_coord$subject <- pca_p$subject

# Threshold using chi-square
threshold <- qchisq(0.99, df = 5)

# Potential outliers
outliers <- pca_coord[pca_coord$Mahalanobis > threshold, ]


# Define plasma outliers 
plasma_outlier_ids <- c("FALL_1", "FALL_36")



###### PCA without outliers ######

# Remove outliers
pca_p_no_outliers <- pca_p %>%
  filter(!subject %in% plasma_outlier_ids)

# Check dimensions
dim(pca_p)
dim(pca_p_no_outliers)

# Check that outliers are removed
subset(pca_p_no_outliers, subject %in% plasma_outlier_ids)


# Now to PCA
res.pca.p.no_outliers <- PCA(
  pca_p_no_outliers,
  graph = FALSE,
  quali.sup = c("subject","sex","cirrose","con_masld",
                "group","group_1","group_urso"),
  quanti.sup = c("bmi","age")
)

# Show the percentage of variances explained by each principal component
eig.val.p.no_outliers <- res.pca.p.no_outliers$eig

#n.pc.p.no_outliers <- min(10, nrow(eig.val.p.no_outliers))
#
#barplot(
#  eig.val.p.no_outliers[1:n.pc.p.no_outliers, 2],
#  names.arg = 1:n.pc.p.no_outliers,
#  main = "Variances Explained by PCs (%)",
#  xlab = "Principal Components",
#  ylab = "Percentage of variance",
#  col = "steelblue"
#)

# Visualize the graph of individuals. Individuals with a similar profile are grouped together.

# By diagnosis

fviz_pca_ind(
  res.pca.p.no_outliers,
  geom.ind = "point",
  col.ind = pca_p_no_outliers$group,
  palette = c("#58aadb","grey","#AA0233"),
  pointsize = 3,
  pointshape = 19,
  invisible = "quali",
  legend.title = "Diagnosis"
) +
  ggtitle("Individuals - PCA plasma without outliers") +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.position = "bottom"
  )

# Save plot
# ggsave(
#   "../output/PCA/pca_plasma_diagnosis_without_outliers.jpeg",
#   width = 10,
#   height = 8,
#   device = "jpeg"
# )


###### PCA with and without MASLD 

# save plot
# jpeg(file = "../output/PCA/pca_plasma_diagnosis_pbc_con.masld_masld_healthy.jpeg")

fviz_pca_ind(res.pca.p, 
             geom.ind = "point", # default is point + text and text is 1:n samples
             col.ind = pca_p$group_1, # color by group
             palette = c("#dba458","#58aadb","grey","#AA0233"),             
             pointsize = 3, # increased size
             pointshape = 19, # makes circles
             invisible = "quali", # removes centroids (group mean dots/points)
             legend.title = "Diagnosis") +
  ggtitle("Individuals - PCA plasma")

dev.off()

# By sex

# save plot
# jpeg(file = "../output/PCA/pca_plasma_sex_pbc_masld_healthy.jpeg")

fviz_pca_ind(res.pca.p, 
             geom.ind = "point", # default er c("point","text") 
             col.ind = pca_p$sex,
             pointsize = 3,
             pointshape = 19,
             invisible = "quali",
             legend.title = "Sex") +
  ggtitle("Individuals - PCA plasma")

# dev.off()

# By cirrhosis

# save plot
# jpeg(file = "../output/PCA/pca_plasma_cirrosis_pbc_masld_healthy.jpeg")
 
fviz_pca_ind(res.pca.p, 
             geom.ind = "point", # default er c("point","text") 
             col.ind = pca_p$cirrose,
             pointsize = 3,
             pointshape = 19,
             invisible = "quali",
             legend.title = "Cirrhosis") +
  ggtitle("Individuals - PCA plasma")+
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),  # Increase legend text size
    axis.text = element_text(size = 14),    # Increase x- and y-axis text size
    axis.title = element_text(size = 16),   # Increase x- and y-axis title size
    plot.title = element_text(size = 18),   # Increase plot title size
    legend.position = "bottom"              # Move legend to the bottom
  )

# dev.off()

# Cirrhosis cluster together

# subset only pbc and healthy
pca_p_pbc <- subset(pca_p, group %in% c("PBC","Healthy"))
levels(pca_p_pbc$group)
pca_p_pbc$group <- droplevels(pca_p_pbc$group)
table(pca_p_pbc$group)

# Now to PCA

res.pca.p.pbc <- PCA(pca_p_pbc, 
               graph = FALSE, 
               quali.sup = c("subject","sex","cirrose","con_masld","group","group_1","group_urso"),
               quanti.sup = c("bmi","age"))

# Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component
eig.val.p.pbc <- res.pca.p.pbc$eig
#
#barplot(eig.val.p.pbc[1:10, 2], 
#        names.arg = 1:10, 
#        main = "Variances Explained by PCs (%)",
#        xlab = "Principal Components",
#        ylab = "Percentage of variances",
#        col ="steelblue")

# Visualize the graph of individuals. Individuals with a similar profile are grouped together.

# By diagnosis
# jpeg(file = "../output/PCA/pca_plasma_diagnosis_pbc_healthy.jpeg")

fviz_pca_ind(res.pca.p.pbc, 
             geom.ind = "point", # default is point + text and text is 1:n samples
             col.ind = pca_p_pbc$group, # color by group
             palette = c("#58aadb","grey"),
             pointsize = 3, # increased size
             pointshape = 19, # makes circles
             invisible = "quali",# removes centroids (group mean dots/points)
             legend.title = "Diagnosis") +
  ggtitle("Individuals - PCA plasma")

# dev.off()

# By sex
# jpeg(file = "../output/PCA/pca_plasma_sex_pbc_healthy.jpeg")
 

fviz_pca_ind(res.pca.p.pbc, 
             geom.ind = "point", # default er c("point","text") 
             col.ind = pca_p_pbc$sex,
             pointsize = 3,
             pointshape = 19,
             invisible = "quali",
             legend.title = "Sex") +
  ggtitle("Individuals - PCA plasma")

# dev.off()

# By cirrhosis
# jpeg(file = "../output/PCA/pca_plasma_cirrhosis_pbc_healthy.jpeg")
 
fviz_pca_ind(res.pca.p.pbc, 
             geom.ind = "point", # default er c("point","text") 
             col.ind = pca_p_pbc$cirrose,
             pointsize = 3,
             pointshape = 19,
             invisible = "quali",
             legend.title = "Cirrhosis") +
  ggtitle("Individuals - PCA plasma")

# dev.off()

# Cirrhosis cluster together


################## Differential expression analysis plasma ########################################################

p <- read.csv("../data/plasma_minprob_data_for_proteomics.tsv", sep = "\t")

p[1:5, 1:5] #instead of head because the df is too big
dim(p)
str(p[1:5])
p$group <- as.factor(p$group)
table(p$group)

# Check if there are any missing values
sum(is.na(p))

# Check if any neg values
any(p[,3:584] < 0)

# Preparing for other group comparisons
p$group_1 <- p$group

p <- p %>%
  mutate(group_1 = case_when(
    group_1 %in% c("Healthy", "NAFLD") ~ "Control",    # Healthy and MASLD as a single Control group
    TRUE ~ group_1
  ))

# Check the new group counts
table(p$group_1)

p$group_2 <- p$group

p <- p %>%
  mutate(group_2 = case_when(
    group_2 %in% c("Healthy", "NAFLD","AIH","PSC") ~ "Others",   # All groups except PBC as a single group
    TRUE ~ group_2
  ))

# Check the new group counts
table(p$group_2)


# Preparing for Ursochol group comparisons
p$group_urso <- p$group

# Subject vectors for Ursochol groups (from metadata)
subjects_pbc_urso_yes <- pca_p$subject[pca_p$group_urso == "PBC_urso_yes"]
subjects_pbc_urso_no  <- pca_p$subject[pca_p$group_urso == "PBC_urso_no"]


p <- p %>%
  mutate(group_urso = case_when(
    subject %in% subjects_pbc_urso_yes ~ "PBC_urso_yes",
    subject %in% subjects_pbc_urso_no  ~ "PBC_urso_no",
    TRUE ~ as.character(group)
  ))

p$group_urso <- factor(p$group_urso)


# Check the new group counts
table(p$group_urso, useNA = "ifany")

# Dataset without outliers 
p_no_outliers <- subset(p, !subject %in% plasma_outlier_ids)

# Check dimensions
dim(p)
dim(p_no_outliers)

# Dataset without cirrhosis

pca_p$subject[pca_p$cirrose=="Yes"]

p_no_c <- subset(p, ! subject %in% c("FALL_6", "FALL_75"))
dim(p)
dim(p_no_c)

# Dataset PBC without MASLD

p_no_masld <- p %>%
  filter(!(subject %in% pca_p$subject[pca_p$group_1 == "PBC with MASLD"]))

dim(p)
dim(p_no_masld)



############### PBC vs. healthy - all individuals ################### 

p_pbc <- subset(p, group %in% c("PBC","Healthy"))
levels(p_pbc$group)
table(p_pbc$group)
p_pbc$group <- droplevels(p_pbc$group)
p_pbc[1:3,1:5]
dim(p_pbc)

design_p_pbc <- model.matrix(~ 0 + group, p_pbc)
design_p_pbc[1:20,]

dim(p_pbc)
p_pbc[1:3,1:5]
p_pbc[1:3,584:586]
p_pbc <- p_pbc[,3:584] # create only numeric values for p_pbc, used for model underneath
dim(p_pbc)

# lmFit expects input array to have structure: protein x sample
# lmFit fits a linear model using weighted least squares for each protein:
# cite: Phipson, B, Lee, S, Majewski, IJ, Alexander, WS, and Smyth, GK (2016). Robust
# hyperparameter estimation protects against hypervariable genes and improves power to
#detect differential expression. Annals of Applied Statistics 10(2), 946–963.

fit_p_pbc <- lmFit(t(p_pbc), design_p_pbc)
head(coef(fit_p_pbc))

# Comparisons between groups (log fold-changes) are obtained as contrasts of
# these fitted linear models:
# Samples are grouped based on experimental condition
# The variability of protein expression is compared between these groups
contr_p_pbc <- makeContrasts(groupPBC-groupHealthy, levels=design_p_pbc)

contr_p_pbc

# Estimate contrast for each protein
tmp_p_pbc <- contrasts.fit(fit_p_pbc, contr_p_pbc)

# Empirical Bayes smoothing of standard errors (shrinks standard errors
# that are much larger or smaller than those from other proteins towards the average standard error)
tmp_p_pbc <- eBayes(tmp_p_pbc)

tmp_p_pbc
# could also set robust = T and trend = T in eBayes (almost same results). 
# Check: https://support.bioconductor.org/p/56560/

### Diagnostic plots:
# Scatterplot of residual-variances vs average log-expression
# save plot
# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy/All/plotSA_pbc_vs_healthy_plasma.jpeg")

plotSA(tmp_p_pbc, main = "Residual variances vs. average log-expression, PBC vs. Healthy, Plasma")

# dev.off()


# Mean Difference plot. Log-intensity ratios (differences, y) versus log-intensity averages (means, x):
# save plot
# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy/All/MD-plot_pbc_vs_healthy_plasma.jpeg")

plotMD(tmp_p_pbc, main = "Mean-difference plot, PBC vs. Healthy, Plasma") 

# dev.off()


# Extract results
top.table.p.pbc <- topTable(tmp_p_pbc, coef = 1, sort.by = "P", n = Inf) #default is "BH" which is alias of "fdr"

results.p.pbc <-  as.data.frame(top.table.p.pbc)
results.p.pbc$protein <- row.names(results.p.pbc)

names(results.p.pbc)

results.p.pbc[1:20,1:6]

# Function to extract UniProt ID
extract_uniprot_id <- function(gene_name_uniprot_id) {
  # Extract the part after the last dot
  uniprot_id <- sub(".*\\.", "", gene_name_uniprot_id)
  
  # Check if the extracted part contains letters
  if (grepl("[A-Za-z]", uniprot_id)) {
    return(uniprot_id)
  } else {
    # Extract the part before the last dot
    return(sub(".*\\.(.*)\\..*", "\\1", gene_name_uniprot_id))
  }
}


# Function to extract GeneName before UniProt ID
extract_gene_name <- function(gene_name_uniprot_id) {
  # Split the string by dots
  parts <- unlist(strsplit(gene_name_uniprot_id, "\\."))
  
  # Initialize the gene name
  gene_name <- parts[1]  # Default to the first part
  
  # Check if there's a UniProt ID pattern (letter followed by numbers/letters) after any dot
  for (i in 2:(length(parts) - 1)) {
    if (grepl("^[A-Za-z][A-Za-z0-9]+$", parts[i])) {
      gene_name <- paste(parts[1:(i-1)], collapse = ".")
      break
    }
  }
  
  # Check the last part for UniProt ID pattern if not found earlier
  if (grepl("^[A-Za-z][A-Za-z0-9]+$", parts[length(parts)])) {
    gene_name <- paste(parts[1:(length(parts)-1)], collapse = ".")
  }
  
  return(gene_name)
}


# Apply the functions to the protein column

results.p.pbc$GeneName <- sapply(results.p.pbc$protein, 
                                           extract_gene_name)

results.p.pbc[1:20,5:8]

results.p.pbc$UniProtID <- sapply(results.p.pbc$protein, 
                                            extract_uniprot_id)

results.p.pbc$GeneName.UniprotID <- results.p.pbc$protein

results.p.pbc[1:20,6:10]

# Make new column name to logFC column (as it is log2FC results)
results.p.pbc[1:3,1:5]
results.p.pbc$log2FC <- results.p.pbc$logFC
results.p.pbc[1:3,c(1,8:11)]



#### Volcano PBC vs healthy

# Identify proteins that meet the filtering criteria

selected_proteins_p_pbc <- results.p.pbc %>%
  filter(adj.P.Val < 0.05)

selected_proteins_p_pbc[60:80,5:9]

# Print the result
# print(selected_accession_names_p_pbc[11:12,1:2])

volcano_plot_p_pbc <- ggplot(results.p.pbc) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), 
                 color = ifelse(logFC <= -0.00001  & adj.P.Val < 0.05, 'blue', 
                                ifelse(logFC >= 0.00001 & adj.P.Val < 0.05, 'red', 'grey')))) +
  geom_text(data = selected_proteins_p_pbc,
            aes(x = logFC, y = -log10(adj.P.Val), 
                label = GeneName),
            size = 3, hjust = 0, vjust = 0) +
  scale_color_identity() +
  labs(x = "log2(FC)", y = "-log10(adj P value)") +
  ggtitle("Differential Plasma Protein Expression, PBC vs. Healthy") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +  # Add horizontal line
  theme(
    axis.title = element_text(size = 20),  # Increase axis title size
    axis.text = element_text(size = 18),   # Increase axis text (tick labels) size
    plot.title = element_text(size = 20)   # Increase plot title size
  )

volcano_plot_p_pbc




#### Check if proteins are up/downregulated in PBC compared to controlgroup ####
# tapply(p$PIGR.P01833,p$group,mean) 
# "checking that up and downregulated are correct (-1,+1)"
# tapply(p$TTR.P02766,p$group, mean)

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/All/volcano_plasma_pbc_healthy_fdr0.05.jpeg",
#         width = 10, height = 10, plot = volcano_plot_p_pbc, device = "jpeg")

# Number of upregulated proteins in PBC:

selected_upregulated_proteins_p_pbc <- results.p.pbc %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_p_pbc) # 45

# Number of downregulated proteins in PBC:
selected_downregulated_proteins_p_pbc <- results.p.pbc %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_p_pbc) # 57

# Sort the upregulated proteins by adjusted p-value
sorted_upregulated_proteins_p_pbc <- 
  arrange(subset(results.p.pbc, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

sorted_upregulated_proteins_p_pbc[1:5,5:10]

# Display the sorted list of upregulated proteins with gene names
print(sorted_upregulated_proteins_p_pbc[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
up_plasma_pbc_vs_healthy <- 
  sorted_upregulated_proteins_p_pbc[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]

# Sort the downregulated proteins by adjusted p-value
sorted_downregulated_proteins_p_pbc <- 
  arrange(subset(results.p.pbc, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of downregulated proteins with gene names
print(sorted_downregulated_proteins_p_pbc[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")])

# Create a data frame with the desired output
down_plasma_pbc_vs_healthy <- 
  sorted_downregulated_proteins_p_pbc[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]



############### PBC vs. healthy - without outliers ################### 

p_pbc_no_outliers <- subset(p_no_outliers, group %in% c("PBC","Healthy"))

levels(p_pbc_no_outliers$group)
table(p_pbc_no_outliers$group)

p_pbc_no_outliers$group <- droplevels(p_pbc_no_outliers$group)

p_pbc_no_outliers[1:3,1:5]

dim(p_pbc_no_outliers)

design_p_pbc_no_outliers <- model.matrix(~ 0 + group, p_pbc_no_outliers)

design_p_pbc_no_outliers[1:20,]

dim(p_pbc_no_outliers)

p_pbc_no_outliers[1:3,1:5]

p_pbc_no_outliers[1:3,584:586]

p_pbc_no_outliers <- p_pbc_no_outliers[,3:584]

dim(p_pbc_no_outliers)

# lmFit expects input array to have structure: protein x sample

fit_p_pbc_no_outliers <- lmFit(t(p_pbc_no_outliers), design_p_pbc_no_outliers)

head(coef(fit_p_pbc_no_outliers))

# Comparisons between groups (log fold-changes) are obtained as contrasts of
# these fitted linear models:

contr_p_pbc_no_outliers <- makeContrasts(
  groupPBC-groupHealthy,
  levels=design_p_pbc_no_outliers
)

contr_p_pbc_no_outliers

# Estimate contrast for each protein

tmp_p_pbc_no_outliers <- contrasts.fit(
  fit_p_pbc_no_outliers,
  contr_p_pbc_no_outliers
)

# Empirical Bayes smoothing of standard errors

tmp_p_pbc_no_outliers <- eBayes(tmp_p_pbc_no_outliers)

tmp_p_pbc_no_outliers


### Diagnostic plots:

# Scatterplot of residual-variances vs average log-expression
# save plot

# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy/Without outliers/plotSA_pbc_vs_healthy_without_outliers_plasma.jpeg")

plotSA(
  tmp_p_pbc_no_outliers,
  main = "Residual variances vs. average log-expression, PBC vs. Healthy without outliers, Plasma"
)

#dev.off()


# Mean Difference plot. Log-intensity ratios (differences, y) versus log-intensity averages (means, x):
# save plot

# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy/Without outliers/MD-plot_pbc_vs_healthy_without_outliers_plasma.jpeg")

plotMD(
  tmp_p_pbc_no_outliers,
  main = "Mean-difference plot, PBC vs. Healthy without outliers, Plasma"
)

# dev.off()

# Extract results

top.table.p.pbc_no_outliers <- topTable(
  tmp_p_pbc_no_outliers,
  coef = 1,
  sort.by = "P",
  n = Inf
)

results.p.pbc_no_outliers <- as.data.frame(top.table.p.pbc_no_outliers)

results.p.pbc_no_outliers$protein <- row.names(results.p.pbc_no_outliers)

names(results.p.pbc_no_outliers)

results.p.pbc_no_outliers[1:20,1:6]


# Apply the functions to the protein column

results.p.pbc_no_outliers$GeneName <- sapply(
  results.p.pbc_no_outliers$protein,
  extract_gene_name
)

results.p.pbc_no_outliers[1:20,5:8]

results.p.pbc_no_outliers$UniProtID <- sapply(
  results.p.pbc_no_outliers$protein,
  extract_uniprot_id
)

results.p.pbc_no_outliers$GeneName.UniprotID <- 
  results.p.pbc_no_outliers$protein

results.p.pbc_no_outliers[1:20,6:10]

# Make new column name to logFC column (as it is log2FC results)

results.p.pbc_no_outliers$log2FC <- 
  results.p.pbc_no_outliers$logFC


#### Volcano PBC vs healthy without outliers ####

selected_proteins_p_pbc_no_outliers <- 
  results.p.pbc_no_outliers %>%
  filter(adj.P.Val < 0.05)

selected_proteins_p_pbc_no_outliers[1:20,5:9]


volcano_plot_p_pbc_no_outliers <- ggplot(results.p.pbc_no_outliers) +
  geom_point(aes(
    x = logFC,
    y = -log10(adj.P.Val), 
    color = ifelse(
      logFC <= -0.00001 & adj.P.Val < 0.05,
      'blue',
      ifelse(
        logFC >= 0.00001 & adj.P.Val < 0.05,
        'red',
        'grey'
      )
    )
  )) +
  geom_text(
    data = selected_proteins_p_pbc_no_outliers,
    aes(
      x = logFC,
      y = -log10(adj.P.Val), 
      label = GeneName
    ),
    size = 3,
    hjust = 0,
    vjust = 0
  ) +
  scale_color_identity() +
  labs(x = "log2(FC)", y = "-log10(adj P value)") +
  ggtitle("Differential Plasma Protein Expression, PBC vs. Healthy without outliers") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "gray"
  ) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 20)
  )

volcano_plot_p_pbc_no_outliers

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/Without outliers/volcano_plasma_pbc_healthy_without_outliers_fdr0.05.jpeg",
#         width = 10, height = 10, plot = volcano_plot_p_pbc_no_outliers, device = "jpeg")

# Number of upregulated proteins in PBC without outliers: #53

selected_upregulated_proteins_p_pbc_no_outliers <- results.p.pbc_no_outliers %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)

nrow(selected_upregulated_proteins_p_pbc_no_outliers)


# Number of downregulated proteins in PBC without outliers: # 56

selected_downregulated_proteins_p_pbc_no_outliers <- results.p.pbc_no_outliers %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)

nrow(selected_downregulated_proteins_p_pbc_no_outliers)


# Sort the upregulated proteins by adjusted p-value

sorted_upregulated_proteins_p_pbc_no_outliers <- 
  arrange(
    subset(results.p.pbc_no_outliers, logFC >= 0.0001 & adj.P.Val < 0.05),
    adj.P.Val
  )

sorted_upregulated_proteins_p_pbc_no_outliers[1:5,5:10]


# Display the sorted list of upregulated proteins with gene names

print(
  sorted_upregulated_proteins_p_pbc_no_outliers[, 
                                                c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]
)


# Create a data frame with the desired output

up_plasma_pbc_vs_healthy_no_outliers <- 
  sorted_upregulated_proteins_p_pbc_no_outliers[, 
                                                c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]


# Sort the downregulated proteins by adjusted p-value

sorted_downregulated_proteins_p_pbc_no_outliers <- 
  arrange(
    subset(results.p.pbc_no_outliers, logFC <= -0.0001 & adj.P.Val < 0.05),
    adj.P.Val
  )


# Display the sorted list of downregulated proteins with gene names

print(
  sorted_downregulated_proteins_p_pbc_no_outliers[, 
                                                  c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]
)


# Create a data frame with the desired output

down_plasma_pbc_vs_healthy_no_outliers <- 
  sorted_downregulated_proteins_p_pbc_no_outliers[, 
                                                  c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]

#################### PBV vs healthy without cirrhosis ################################################

p_pbc_no_c <- subset(p_no_c, group %in% c("PBC","Healthy"))
levels(p_pbc_no_c$group)
table(p_pbc_no_c$group)
p_pbc_no_c$group <- droplevels(p_pbc_no_c$group)
table(p_pbc_no_c$group)

design_p_pbc_no_c <- model.matrix(~ 0 + group, p_pbc_no_c)
design_p_pbc_no_c[1:20,] 

dim(p_pbc_no_c)
p_pbc_no_c[1:3,1:5]
p_pbc_no_c[1:3,584:586]
p_pbc_no_c <- p_pbc_no_c[,3:584] # create only numeric values for p_pbc_no_c, used for model underneath
dim(p_pbc_no_c)

fit_p_pbc_no_c <- lmFit(t(p_pbc_no_c), design_p_pbc_no_c)
head(coef(fit_p_pbc_no_c))

contr_p_pbc_no_c <- makeContrasts(groupPBC-groupHealthy, levels=design_p_pbc_no_c)
contr_p_pbc_no_c

tmp_p_pbc_no_c <- contrasts.fit(fit_p_pbc_no_c, contr_p_pbc_no_c)
tmp_p_pbc_no_c <- eBayes(tmp_p_pbc_no_c)



### Diagnostic plots:
# Scatterplot of residual-variances vs average log-expression
# save plot
# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy/No cirrhosis/plotSA_pbc_vs_healthy_plasma_no.cirrhosis.jpeg")

#plotSA(tmp_p_pbc_no_c, main = "Residual variances vs. average log-expression, non-cirrhotic PBC vs. Healthy, Plasma")
#dev.off()

# Mean Difference plot. Log-intensity ratios (differences, y) versus log-intensity averages (means, x):
# save plot
# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy/No cirrhosis/MD-plot_pbc_vs_healthy_plasma_no.cirrhosis.jpeg")

#plotMD(tmp_p_pbc_no_c, main = "Mean-difference plot, non-cirrhotic PBC vs. Healthy, Plasma") 
#dev.off()

# Extract results
top.table.p.pbc.no.c <- topTable(tmp_p_pbc_no_c, coef = 1, sort.by = "P", n = Inf) #default is "BH" which is alias of "fdr"
results.p.pbc.no.c <-  as.data.frame(top.table.p.pbc.no.c)
results.p.pbc.no.c$protein <- row.names(results.p.pbc.no.c)
names(results.p.pbc.no.c)
results.p.pbc.no.c[60:80,1:7]
results.p.pbc.no.c$GeneName <- sapply(results.p.pbc.no.c$protein, 
                                 extract_gene_name)
results.p.pbc.no.c[60:80,5:8]
results.p.pbc.no.c$UniProtID <- sapply(results.p.pbc.no.c$protein, 
                                  extract_uniprot_id)
results.p.pbc.no.c$GeneName.UniprotID <- results.p.pbc.no.c$protein
results.p.pbc.no.c[60:80,6:10]
results.p.pbc.no.c$log2FC <- results.p.pbc.no.c$logFC
results.p.pbc.no.c[1:3,c(1,8:11)]





############## Differential expression PBC vs MASLD ###################

p_pbc_masld <- subset(p, group %in% c("PBC","NAFLD"))
levels(p_pbc_masld$group)
table(p_pbc_masld$group)
p_pbc_masld$group <- droplevels(p_pbc_masld$group)
p_pbc_masld[1:3,1:5]
dim(p_pbc_masld)

design_p_pbc_masld <- model.matrix(~ 0 + group, p_pbc_masld)
design_p_pbc_masld[1:20,]

dim(p_pbc_masld)
p_pbc_masld[1:3,1:5]
p_pbc_masld[1:3,583:587]
p_pbc_masld <- p_pbc_masld[,3:584] # create only numeric values
dim(p_pbc_masld)

fit_p_pbc_masld <- lmFit(t(p_pbc_masld), design_p_pbc_masld)
head(coef(fit_p_pbc_masld))

contr_p_pbc_masld <- makeContrasts(groupPBC-groupNAFLD, levels=design_p_pbc_masld)
contr_p_pbc_masld

tmp_p_pbc_masld <- contrasts.fit(fit_p_pbc_masld, contr_p_pbc_masld)
tmp_p_pbc_masld <- eBayes(tmp_p_pbc_masld)
tmp_p_pbc_masld

top.table.p.pbc.masld <- topTable(tmp_p_pbc_masld, coef = 1, sort.by = "P", n = Inf)
results.p.pbc.masld <- as.data.frame(top.table.p.pbc.masld)

results.p.pbc.masld$protein <- row.names(results.p.pbc.masld)
names(results.p.pbc.masld)
results.p.pbc.masld[1:20,1:6]

results.p.pbc.masld$GeneName <- sapply(results.p.pbc.masld$protein, 
                                       extract_gene_name)
results.p.pbc.masld[1:20,5:8]

results.p.pbc.masld$UniProtID <- sapply(results.p.pbc.masld$protein, 
                                        extract_uniprot_id)

results.p.pbc.masld$GeneName.UniprotID <- results.p.pbc.masld$protein
results.p.pbc.masld[1:20,6:10]

results.p.pbc.masld[1:3,1:5]
results.p.pbc.masld$log2FC <- results.p.pbc.masld$logFC
results.p.pbc.masld[1:3,c(1,8:11)]

#### Volcano PBC

selected_proteins_p_pbc_masld <- results.p.pbc.masld %>%
  filter(adj.P.Val < 0.05)

selected_proteins_p_pbc_masld[60:80,5:9]

volcano_plot_p_pbc_masld <- ggplot(results.p.pbc.masld) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), 
                 color = ifelse(logFC <= -0.00001  & adj.P.Val < 0.05, 'blue', 
                                ifelse(logFC >= 0.00001 & adj.P.Val < 0.05, 'red', 'grey')))) +
  geom_text(data = selected_proteins_p_pbc_masld,
            aes(x = logFC, y = -log10(adj.P.Val), 
                label = GeneName),
            size = 3, hjust = 0, vjust = 0) +
  scale_color_identity() +
  labs(x = "log2(FC)", y = "-log10(adj P value)") +
  ggtitle("Differential Plasma Protein Expression, PBC vs. MASLD") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 20)
  )

volcano_plot_p_pbc_masld

# tapply(p$ACY1.Q03154, p$group,mean) # "checking that up and downregulated are correct (-1,+1)"
# tapply(p$RBP4.P02753, p$group, mean)

# Example check:
# results.p.pbc.masld[results.p.pbc.masld$GeneName == "ICAM1",c(1,5,7:10)]
# 2^(results.p.pbc.masld[results.p.pbc.masld$GeneName == "ICAM1","logFC"])

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_MASLD/All/volcano_plasma_pbc_masld_fdr0.05.jpeg",
#         width = 10, height = 10, plot = volcano_plot_p_pbc_masld, device = "jpeg")

# Number of upregulated proteins in PBC:

selected_upregulated_proteins_p_pbc_masld <- results.p.pbc.masld %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)

nrow(selected_upregulated_proteins_p_pbc_masld)

# Number of downregulated proteins in PBC:

selected_downregulated_proteins_p_pbc_masld <- results.p.pbc.masld %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)

nrow(selected_downregulated_proteins_p_pbc_masld)

# Sort the upregulated proteins by adjusted p-value

sorted_upregulated_proteins_p_pbc_masld <- 
  arrange(subset(results.p.pbc.masld, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

sorted_upregulated_proteins_p_pbc_masld[1:5,5:10]

# Display the sorted list of upregulated proteins with gene names

print(sorted_upregulated_proteins_p_pbc_masld[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output

up_plasma_pbc_vs_masld <- 
  sorted_upregulated_proteins_p_pbc_masld[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]

# Sort the downregulated proteins by adjusted p-value

sorted_downregulated_proteins_p_pbc_masld <- 
  arrange(subset(results.p.pbc.masld, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of downregulated proteins with gene names

print(sorted_downregulated_proteins_p_pbc_masld[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")])

# Create a data frame with the desired output

down_plasma_pbc_vs_masld <- 
  sorted_downregulated_proteins_p_pbc_masld[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]



#################### PBC vs MASLD without cirrhosis ################################################

p_pbc_masld_no_c <- subset(p_no_c, group %in% c("PBC","NAFLD"))
levels(p_pbc_masld_no_c$group)
table(p_pbc_masld_no_c$group)
p_pbc_masld_no_c$group <- droplevels(p_pbc_masld_no_c$group)
table(p_pbc_masld_no_c$group)

design_p_pbc_masld_no_c <- model.matrix(~ 0 + group, p_pbc_masld_no_c)
design_p_pbc_masld_no_c[1:20,] 

dim(p_pbc_masld_no_c)
p_pbc_masld_no_c[1:3,1:5]
p_pbc_masld_no_c[1:3,583:587]
p_pbc_masld_no_c <- p_pbc_masld_no_c[,3:584] # create only numeric values
dim(p_pbc_masld_no_c)

fit_p_pbc_masld_no_c <- lmFit(t(p_pbc_masld_no_c), design_p_pbc_masld_no_c)
head(coef(fit_p_pbc_masld_no_c))

contr_p_pbc_masld_no_c <- makeContrasts(groupPBC-groupNAFLD, levels=design_p_pbc_masld_no_c)
contr_p_pbc_masld_no_c

tmp_p_pbc_masld_no_c <- contrasts.fit(fit_p_pbc_masld_no_c, contr_p_pbc_masld_no_c)
tmp_p_pbc_masld_no_c <- eBayes(tmp_p_pbc_masld_no_c)

top.table.p.pbc.masld.no.c <- topTable(tmp_p_pbc_masld_no_c, coef = 1, sort.by = "P", n = Inf)
results.p.pbc.masld.no.c <- as.data.frame(top.table.p.pbc.masld.no.c)

results.p.pbc.masld.no.c$protein <- row.names(results.p.pbc.masld.no.c)
names(results.p.pbc.masld.no.c)
results.p.pbc.masld.no.c[60:80,1:7]

results.p.pbc.masld.no.c$GeneName <- sapply(results.p.pbc.masld.no.c$protein, 
                                            extract_gene_name)
results.p.pbc.masld.no.c[60:80,5:8]

results.p.pbc.masld.no.c$UniProtID <- sapply(results.p.pbc.masld.no.c$protein, 
                                             extract_uniprot_id)

results.p.pbc.masld.no.c$GeneName.UniprotID <- results.p.pbc.masld.no.c$protein
results.p.pbc.masld.no.c[60:80,6:10]

results.p.pbc.masld.no.c$log2FC <- results.p.pbc.masld.no.c$logFC
results.p.pbc.masld.no.c[1:3,c(1,8:11)]

#### Volcano PBC

selected_proteins_p_pbc_masld_no_c <- results.p.pbc.masld.no.c %>%
  filter(adj.P.Val < 0.05)

selected_proteins_p_pbc_masld_no_c[1:5,1:5]

volcano_plot_p_pbc_masld_no_c <- ggplot(results.p.pbc.masld.no.c) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), 
                 color = ifelse(logFC <= -0.00001  & adj.P.Val < 0.05, 'blue', 
                                ifelse(logFC >= 0.00001 & adj.P.Val < 0.05, 'red', 'grey')))) +
  geom_text(data = selected_proteins_p_pbc_masld_no_c,
            aes(x = logFC, y = -log10(adj.P.Val), 
                label = GeneName),
            size = 3, hjust = 0, vjust = 0) +
  scale_color_identity() +
  labs(x = "log2(FC)", y = "-log10(adj P value)") +
  ggtitle("Differential Plasma Protein Expression, non-cirrhotic PBC vs. MASLD")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray")+
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 20)
  )

volcano_plot_p_pbc_masld_no_c

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_MASLD/No cirrhosis/volcano_plasma_pbc_masld_fdr0.05_no.cirrhosis.jpeg",
#         width = 10, height = 10, plot = volcano_plot_p_pbc_masld_no_c, device = "jpeg")

# Number of upregulated proteins in PBC:

selected_upregulated_proteins_p_pbc_masld_no_c <- results.p.pbc.masld.no.c %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)

nrow(selected_upregulated_proteins_p_pbc_masld_no_c)

# Number of downregulated proteins in PBC:

selected_downregulated_proteins_p_pbc_masld_no_c <- results.p.pbc.masld.no.c %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)

nrow(selected_downregulated_proteins_p_pbc_masld_no_c)

# Sort the upregulated proteins by adjusted p-value

sorted_upregulated_proteins_p_pbc_masld_no_c <- 
  arrange(subset(results.p.pbc.masld.no.c, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of upregulated proteins with gene names

print(sorted_upregulated_proteins_p_pbc_masld_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")])

# Create a data frame with the desired output

up_plasma_pbc_vs_masld_no_c <- 
  sorted_upregulated_proteins_p_pbc_masld_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]

# Sort the downregulated proteins by adjusted p-value

sorted_downregulated_proteins_p_pbc_masld_no_c <- 
  arrange(subset(results.p.pbc.masld.no.c, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of downregulated proteins with gene names

print(sorted_downregulated_proteins_p_pbc_masld_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")])

# Create a data frame with the desired output

down_plasma_pbc_vs_masld_no_c <- 
  sorted_downregulated_proteins_p_pbc_masld_no_c[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]


############## DE PBC without MASLD ################################

p_pbc_no_masld <- subset(p_no_masld, group %in% c("PBC","Healthy"))
levels(p_pbc_no_masld$group)
table(p_pbc_no_masld$group)
p_pbc_no_masld$group <- droplevels(p_pbc_no_masld$group)
table(p_pbc_no_masld$group)

design_p_pbc_no_masld <- model.matrix(~ 0 + group, p_pbc_no_masld)
design_p_pbc_no_masld[1:20,] 

dim(p_pbc_no_masld)
p_pbc_no_masld[1:3,1:5]
p_pbc_no_masld[1:3,584:586]
p_pbc_no_masld <- p_pbc_no_masld[,3:584] # create only numeric values for p_pbc_no_masld, used for model underneath
dim(p_pbc_no_masld)

fit_p_pbc_no_masld <- lmFit(t(p_pbc_no_masld), design_p_pbc_no_masld)
head(coef(fit_p_pbc_no_masld))

contr_p_pbc_no_masld <- makeContrasts(groupPBC-groupHealthy, levels=design_p_pbc_no_masld)
contr_p_pbc_no_masld

tmp_p_pbc_no_masld <- contrasts.fit(fit_p_pbc_no_masld, contr_p_pbc_no_masld)
tmp_p_pbc_no_masld <- eBayes(tmp_p_pbc_no_masld)



### Diagnostic plots:
# Scatterplot of residual-variances vs average log-expression
# save plot
# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy/No MASLD/plotSA_pbc_vs_healthy_plasma_no.masld.jpeg")
# plotSA(tmp_p_pbc_no_masld, main = "Residual variances vs. average log-expression, PBC without MASLD vs. Healthy, Plasma")
# dev.off()

# Mean Difference plot. Log-intensity ratios (differences, y) versus log-intensity averages (means, x):
# save plot
# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy/No MASLD/MD-plot_pbc_vs_healthy_plasma_no.masld.jpeg")
# plotMD(tmp_p_pbc_no_masld, main = "Mean-difference plot, PBC without MASLD vs. Healthy, Plasma") 
# dev.off()

# Extract results
top.table.p.pbc.no.masld <- topTable(tmp_p_pbc_no_masld, coef = 1, sort.by = "P", n = Inf) #default is "BH" which is alias of "fdr"
results.p.pbc.no.masld <-  as.data.frame(top.table.p.pbc.no.masld)
results.p.pbc.no.masld$protein <- row.names(results.p.pbc.no.masld)
names(results.p.pbc.no.masld)
results.p.pbc.no.masld[60:80,1:7]
results.p.pbc.no.masld$GeneName <- sapply(results.p.pbc.no.masld$protein, 
                                      extract_gene_name)
results.p.pbc.no.masld[60:80,5:8]
results.p.pbc.no.masld$UniProtID <- sapply(results.p.pbc.no.masld$protein, 
                                       extract_uniprot_id)
results.p.pbc.no.masld$GeneName.UniprotID <- results.p.pbc.no.masld$protein
results.p.pbc.no.masld[60:80,6:10]
results.p.pbc.no.masld$log2FC <- results.p.pbc.no.masld$logFC
results.p.pbc.no.masld[1:3,c(1,8:11)]


#### Volcano PBC without cirrhosis  

# Identify proteins that meet the filtering criteria

selected_proteins_p_pbc_no_c <- results.p.pbc.no.c %>%
  filter(adj.P.Val < 0.05)

selected_proteins_p_pbc_no_c[1:5,1:5]

volcano_plot_p_pbc_no_c <- ggplot(results.p.pbc.no.c) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), 
                 color = ifelse(logFC <= -0.00001  & adj.P.Val < 0.05, 'blue', 
                                ifelse(logFC >= 0.00001 & adj.P.Val < 0.05, 'red', 'grey')))) +
  geom_text(data = selected_proteins_p_pbc_no_c,
            aes(x = logFC, y = -log10(adj.P.Val), 
                label = GeneName),
            size = 3, hjust = 0, vjust = 0) +
  scale_color_identity() +
  labs(x = "log2(FC)", y = "-log10(adj P value)") +
  ggtitle("Differential Plasma Protein Expression, non-cirrhotic PBC vs. Healthy")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray")+  # Add horizontal line
  theme(
    axis.title = element_text(size = 20),  # Increase axis title size
    axis.text = element_text(size = 18),   # Increase axis text (tick labels) size
    plot.title = element_text(size = 20)   # Increase plot title size
  )

volcano_plot_p_pbc_no_c

#save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/No cirrhosis/volcano_plasma_pbc_healthy_fdr0.05_no.cirrhosis.jpeg",
#         width = 10, height = 10, plot = volcano_plot_p_pbc_no_c, device = "jpeg")

# Number of upregulated proteins in PBC without cirrhosis vs Healthy:

selected_upregulated_proteins_p_pbc_no_c <- results.p.pbc.no.c %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_p_pbc_no_c) # 38

# Number of downregulated proteins in PBC:
selected_downregulated_proteins_p_pbc_no_c <- results.p.pbc.no.c %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_p_pbc_no_c) # 51

# Sort the upregulated proteins by adjusted p-value
sorted_upregulated_proteins_p_pbc_no_c <- 
  arrange(subset(results.p.pbc.no.c, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of upregulated proteins with gene names
print(sorted_upregulated_proteins_p_pbc_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")])

# Create a data frame with the desired output
up_plasma_pbc_vs_healthy_no_c <- 
  sorted_upregulated_proteins_p_pbc_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]

# Sort the downregulated proteins by adjusted p-value
sorted_downregulated_proteins_p_pbc_no_c <- 
  arrange(subset(results.p.pbc.no.c, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of downregulated proteins with gene names
print(sorted_downregulated_proteins_p_pbc_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")])

# Create a data frame with the desired output
down_plasma_pbc_vs_healthy_no_c <- 
  sorted_downregulated_proteins_p_pbc_no_c[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]



#### Volcano PBC without MASLD 

# Identify proteins that meet the filtering criteria

selected_proteins_p_pbc_no_masld <- results.p.pbc.no.masld %>%
  filter(adj.P.Val < 0.05)

selected_proteins_p_pbc_no_masld[1:5,1:5]

volcano_plot_p_pbc_no_masld <- ggplot(results.p.pbc.no.masld) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), 
                 color = ifelse(logFC <= -0.00001  & adj.P.Val < 0.05, 'blue', 
                                ifelse(logFC >= 0.00001 & adj.P.Val < 0.05, 'red', 'grey')))) +
  geom_text(data = selected_proteins_p_pbc_no_masld,
            aes(x = logFC, y = -log10(adj.P.Val), 
                label = GeneName),
            size = 3, hjust = 0, vjust = 0) +
  scale_color_identity() +
  labs(x = "log2(FC)", y = "-log10(adj P value)") +
  ggtitle("Differential Plasma Protein Expression, PBC without MASLD vs. Healthy")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray")+  # Add horizontal line
  theme(
    axis.title = element_text(size = 20),  # Increase axis title size
    axis.text = element_text(size = 18),   # Increase axis text (tick labels) size
    plot.title = element_text(size = 20)   # Increase plot title size
  )

volcano_plot_p_pbc_no_masld

#save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/No MASLD/volcano_plasma_pbc_healthy_fdr0.05_no.masld.jpeg",
#         width = 10, height = 10, plot = volcano_plot_p_pbc_no_masld, device = "jpeg")


# Number of upregulated proteins in PBC without MASLD vs Healthy:

selected_upregulated_proteins_p_pbc_no_masld <- results.p.pbc.no.masld %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_p_pbc_no_masld) # 41

# Number of downregulated proteins in PBC:
selected_downregulated_proteins_p_pbc_no_masld <- results.p.pbc.no.masld %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_p_pbc_no_masld) # 40

# Sort the upregulated proteins by adjusted p-value
sorted_upregulated_proteins_p_pbc_no_masld <- 
  arrange(subset(results.p.pbc.no.masld, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of upregulated proteins with gene names
print(sorted_upregulated_proteins_p_pbc_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")])

# Create a data frame with the desired output
up_plasma_pbc_vs_healthy_no_masld <- 
  sorted_upregulated_proteins_p_pbc_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]

# Sort the downregulated proteins by adjusted p-value
sorted_downregulated_proteins_p_pbc_no_masld <- 
  arrange(subset(results.p.pbc.no.masld, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of downregulated proteins with gene names
print(sorted_downregulated_proteins_p_pbc_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")])

# Create a data frame with the desired output
down_plasma_pbc_vs_healthy_no_masld <- 
  sorted_downregulated_proteins_p_pbc_no_masld[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]



################## PBC vs Healthy and MASLD - all individuals

# Check the new group counts
table(p$group_1)

### DE analysis of PBC compared with Healthy/MASLD

p_pbc_h.m <- subset(p, group_1 %in% c("PBC","Control"))

table(p_pbc_h.m$group_1)

design_p_pbc_h.m <- model.matrix(~ 0 + group_1, p_pbc_h.m)
design_p_pbc_h.m[1:20,]

dim(p_pbc_h.m)
p_pbc_h.m[1:3,1:5]
p_pbc_h.m[1:3,584:586]

p_pbc_h.m <- p_pbc_h.m[,3:584] # create only numeric values for p_pbc_h.m, used for model underneath
dim(p_pbc_h.m)

fit_p_pbc_h.m <- lmFit(t(p_pbc_h.m), design_p_pbc_h.m)
head(coef(fit_p_pbc_h.m))

contr_p_pbc_h.m <- makeContrasts(group_1PBC-group_1Control, levels=design_p_pbc_h.m)
contr_p_pbc_h.m

tmp_p_pbc_h.m <- contrasts.fit(fit_p_pbc_h.m, contr_p_pbc_h.m)
tmp_p_pbc_h.m <- eBayes(tmp_p_pbc_h.m)
tmp_p_pbc_h.m

### Diagnostic plots:
# Scatterplot of residual-variances vs average log-expression
# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/All/plotSA_pbc_vs_healthy_and_masld_plasma.jpeg")
# plotSA(tmp_p_pbc_h.m, main = "Residual variances vs. average log-expression, PBC vs. Healthy and MASLD, Plasma")
# dev.off()

# Mean Difference plot. Log-intensity ratios (differences, y) versus log-intensity averages (means, x):
# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/All/MD-plot_pbc_vs_healthy_and_masld_plasma.jpeg")# plotMD(tmp_p_pbc_h.m, main = "Mean-difference plot, PBC vs. Healthy and MASLD, Plasma") 
# dev.off()

# Extract results
top.table.p.pbc.h.m <- topTable(tmp_p_pbc_h.m, coef = 1, sort.by = "P", n = Inf) 
results.p.pbc.h.m <-  as.data.frame(top.table.p.pbc.h.m)
results.p.pbc.h.m$protein <- row.names(results.p.pbc.h.m)
names(results.p.pbc.h.m)
results.p.pbc.h.m[1:20,1:7]
results.p.pbc.h.m$GeneName <- sapply(results.p.pbc.h.m$protein, 
                                 extract_gene_name)
results.p.pbc.h.m[1:20,5:8]
results.p.pbc.h.m$UniProtID <- sapply(results.p.pbc.h.m$protein, 
                                  extract_uniprot_id)
results.p.pbc.h.m$GeneName.UniprotID <- results.p.pbc.h.m$protein
results.p.pbc.h.m[1:20,6:10]
results.p.pbc.h.m$log2FC <- results.p.pbc.h.m$logFC
results.p.pbc.h.m[1:3,c(1,8:11)]

#### Create volcano plot

# Identify proteins that meet the filtering criteria

selected_proteins_p_pbc_h.m <- results.p.pbc.h.m %>%
  filter(adj.P.Val < 0.05)

selected_proteins_p_pbc_h.m[1:20,5:10]

volcano_plot_p_pbc_h.m <- ggplot(results.p.pbc.h.m) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), 
                 color = ifelse(logFC <= -0.00001  & adj.P.Val < 0.05, 'blue', 
                                ifelse(logFC >= 0.00001 & adj.P.Val < 0.05, 'red', 'grey')))) +
  geom_text(data = selected_proteins_p_pbc_h.m,
            aes(x = logFC, y = -log10(adj.P.Val), 
                label = GeneName),
            size = 3, hjust = 0, vjust = 0) +
  scale_color_identity() +
  labs(x = "log2 fold change", y = "-log10 adj p-value") +
  ggtitle("Differential Plasma Protein Expression, PBC vs. Healthy and MASLD")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray")  # Add horizontal line

volcano_plot_p_pbc_h.m

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/All/volcano_plasma_pbc_healthy_masld_fdr0.05.jpeg",
#         width = 10, height = 10, plot = volcano_plot_p_pbc_h.m, device = "jpeg")

# Number of upregulated proteins in PBC vs Healthy and MASLD:

selected_upregulated_proteins_p_pbc_h.m <- results.p.pbc.h.m %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_p_pbc_h.m) # 69

# Number of downregulated proteins in PBC:
selected_downregulated_proteins_p_pbc_h.m <- results.p.pbc.h.m %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_p_pbc_h.m) # 52

# Sort the upregulated proteins by adjusted p-value
sorted_upregulated_proteins_p_pbc_h.m <- 
  arrange(subset(results.p.pbc.h.m, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of upregulated proteins with gene names
print(sorted_upregulated_proteins_p_pbc_h.m[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")])

# Create a data frame with the desired output
up_plasma_pbc_vs_h.m <- 
  sorted_upregulated_proteins_p_pbc_h.m[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]

# Sort the downregulated proteins by adjusted p-value
sorted_downregulated_proteins_p_pbc_h.m <- 
  arrange(subset(results.p.pbc.h.m, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of downregulated proteins with gene names
print(sorted_downregulated_proteins_p_pbc_h.m[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
down_plasma_pbc_vs_h.m <- 
  sorted_downregulated_proteins_p_pbc_h.m[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]



################## PBC vs Healthy and MASLD - without outliers ##################

# Check the new group counts
table(p_no_outliers$group_1)

### DE analysis of PBC compared with Healthy/MASLD

p_pbc_h.m_no_outliers <- subset(p_no_outliers, group_1 %in% c("PBC","Control"))

table(p_pbc_h.m_no_outliers$group_1)

design_p_pbc_h.m_no_outliers <- model.matrix(~ 0 + group_1, p_pbc_h.m_no_outliers)

design_p_pbc_h.m_no_outliers[1:20,]

dim(p_pbc_h.m_no_outliers)

p_pbc_h.m_no_outliers[1:3,1:5]

p_pbc_h.m_no_outliers[1:3,584:586]

p_pbc_h.m_no_outliers <- p_pbc_h.m_no_outliers[,3:584]

dim(p_pbc_h.m_no_outliers)

fit_p_pbc_h.m_no_outliers <- lmFit(
  t(p_pbc_h.m_no_outliers),
  design_p_pbc_h.m_no_outliers
)

head(coef(fit_p_pbc_h.m_no_outliers))

contr_p_pbc_h.m_no_outliers <- makeContrasts(
  group_1PBC-group_1Control,
  levels=design_p_pbc_h.m_no_outliers
)

contr_p_pbc_h.m_no_outliers

tmp_p_pbc_h.m_no_outliers <- contrasts.fit(
  fit_p_pbc_h.m_no_outliers,
  contr_p_pbc_h.m_no_outliers
)

tmp_p_pbc_h.m_no_outliers <- eBayes(tmp_p_pbc_h.m_no_outliers)

tmp_p_pbc_h.m_no_outliers


### Diagnostic plots

# Scatterplot of residual-variances vs average log-expression
# save plot

# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/Without outliers/plotSA_pbc_vs_healthy_and_masld_without_outliers_plasma.jpeg")

plotSA(
  tmp_p_pbc_h.m_no_outliers,
  main = "Residual variances vs. average log-expression, PBC vs. Healthy and MASLD without outliers, Plasma"
)

#dev.off()


# Mean Difference plot. Log-intensity ratios versus log-intensity averages
# save plot

# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/Without outliers/MD-plot_pbc_vs_healthy_and_masld_without_outliers_plasma.jpeg")

plotMD(
  tmp_p_pbc_h.m_no_outliers,
  main = "Mean-difference plot, PBC vs. Healthy and MASLD without outliers, Plasma"
)

# dev.off()
# Extract results

top.table.p.pbc.h.m_no_outliers <- topTable(
  tmp_p_pbc_h.m_no_outliers,
  coef = 1,
  sort.by = "P",
  n = Inf
)

results.p.pbc.h.m_no_outliers <- as.data.frame(top.table.p.pbc.h.m_no_outliers)

results.p.pbc.h.m_no_outliers$protein <- row.names(results.p.pbc.h.m_no_outliers)

names(results.p.pbc.h.m_no_outliers)

results.p.pbc.h.m_no_outliers[1:20,1:7]

results.p.pbc.h.m_no_outliers$GeneName <- sapply(
  results.p.pbc.h.m_no_outliers$protein,
  extract_gene_name
)

results.p.pbc.h.m_no_outliers[1:20,5:8]

results.p.pbc.h.m_no_outliers$UniProtID <- sapply(
  results.p.pbc.h.m_no_outliers$protein,
  extract_uniprot_id
)

results.p.pbc.h.m_no_outliers$GeneName.UniprotID <- 
  results.p.pbc.h.m_no_outliers$protein

results.p.pbc.h.m_no_outliers[1:20,6:10]

results.p.pbc.h.m_no_outliers$log2FC <- 
  results.p.pbc.h.m_no_outliers$logFC

results.p.pbc.h.m_no_outliers[1:3,c(1,8:11)]


#### Create volcano plot ####

selected_proteins_p_pbc_h.m_no_outliers <- 
  results.p.pbc.h.m_no_outliers %>%
  filter(adj.P.Val < 0.05)

selected_proteins_p_pbc_h.m_no_outliers[1:20,5:10]

volcano_plot_p_pbc_h.m_no_outliers <- ggplot(results.p.pbc.h.m_no_outliers) +
  geom_point(aes(
    x = logFC,
    y = -log10(adj.P.Val), 
    color = ifelse(
      logFC <= -0.00001 & adj.P.Val < 0.05,
      'blue',
      ifelse(
        logFC >= 0.00001 & adj.P.Val < 0.05,
        'red',
        'grey'
      )
    )
  )) +
  geom_text(
    data = selected_proteins_p_pbc_h.m_no_outliers,
    aes(
      x = logFC,
      y = -log10(adj.P.Val), 
      label = GeneName
    ),
    size = 3,
    hjust = 0,
    vjust = 0
  ) +
  scale_color_identity() +
  labs(x = "log2 fold change", y = "-log10 adj p-value") +
  ggtitle("Differential Plasma Protein Expression, PBC vs. Healthy and MASLD without outliers") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "gray"
  )

volcano_plot_p_pbc_h.m_no_outliers


# Save plot
# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/Without outliers/volcano_plasma_pbc_healthy_masld_without_outliers_fdr0.05.jpeg",
#         width = 10, height = 10, plot = volcano_plot_p_pbc_h.m_no_outliers, device = "jpeg")


# Number of upregulated proteins in PBC vs Healthy and MASLD without outliers: #38

selected_upregulated_proteins_p_pbc_h.m_no_outliers <- 
  results.p.pbc.h.m_no_outliers %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)

nrow(selected_upregulated_proteins_p_pbc_h.m_no_outliers)


# Number of downregulated proteins in PBC vs healthy and MASLD without outliers: 31

selected_downregulated_proteins_p_pbc_h.m_no_outliers <- 
  results.p.pbc.h.m_no_outliers %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)

nrow(selected_downregulated_proteins_p_pbc_h.m_no_outliers)


# Sort the upregulated proteins by adjusted p-value

sorted_upregulated_proteins_p_pbc_h.m_no_outliers <- 
  arrange(
    subset(results.p.pbc.h.m_no_outliers, logFC >= 0.0001 & adj.P.Val < 0.05),
    adj.P.Val
  )


# Display the sorted list of upregulated proteins with gene names

print(
  sorted_upregulated_proteins_p_pbc_h.m_no_outliers[, 
                                                    c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]
)


# Create a data frame with the desired output

up_plasma_pbc_vs_h.m_no_outliers <- 
  sorted_upregulated_proteins_p_pbc_h.m_no_outliers[, 
                                                    c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]


# Sort the downregulated proteins by adjusted p-value

sorted_downregulated_proteins_p_pbc_h.m_no_outliers <- 
  arrange(
    subset(results.p.pbc.h.m_no_outliers, logFC <= -0.0001 & adj.P.Val < 0.05),
    adj.P.Val
  )


# Display the sorted list of downregulated proteins with gene names

print(
  sorted_downregulated_proteins_p_pbc_h.m_no_outliers[, 
                                                      c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]
)


# Create a data frame with the desired output

down_plasma_pbc_vs_h.m_no_outliers <- 
  sorted_downregulated_proteins_p_pbc_h.m_no_outliers[, 
                                                      c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]





#### PBC vs. Healthy and MASLD, no cirrhosis

table(p_no_c$group_1)

p_pbc_h.m_no_c <- subset(p_no_c, group_1 %in% c("PBC","Control"))

table(p_pbc_h.m_no_c$group_1)

design_p_pbc_h.m_no_c <- model.matrix(~ 0 + group_1, p_pbc_h.m_no_c)
design_p_pbc_h.m_no_c[1:20,]

dim(p_pbc_h.m_no_c)
p_pbc_h.m_no_c[1:3,1:5]
p_pbc_h.m_no_c[1:3,584:586]
p_pbc_h.m_no_c <- p_pbc_h.m_no_c[,3:584] # create only numeric values 
dim(p_pbc_h.m_no_c)

fit_p_pbc_h.m_no_c <- lmFit(t(p_pbc_h.m_no_c), design_p_pbc_h.m_no_c)
head(coef(fit_p_pbc_h.m_no_c))

contr_p_pbc_h.m_no_c <- makeContrasts(group_1PBC-group_1Control, levels=design_p_pbc_h.m_no_c)
contr_p_pbc_h.m_no_c
tmp_p_pbc_h.m_no_c <- contrasts.fit(fit_p_pbc_h.m_no_c, contr_p_pbc_h.m_no_c)
tmp_p_pbc_h.m_no_c <- eBayes(tmp_p_pbc_h.m_no_c)
tmp_p_pbc_h.m_no_c
top.table.p.pbc.h.m.no.c <- topTable(tmp_p_pbc_h.m_no_c,coef = 1, sort.by = "P", n = Inf)
results.p.pbc.h.m.no.c <-  as.data.frame(top.table.p.pbc.h.m.no.c)
results.p.pbc.h.m.no.c$protein <- row.names(results.p.pbc.h.m.no.c)
names(results.p.pbc.h.m.no.c)
results.p.pbc.h.m.no.c[1:20,1:7]
results.p.pbc.h.m.no.c$GeneName <- sapply(results.p.pbc.h.m.no.c$protein, 
                                      extract_gene_name)
results.p.pbc.h.m.no.c[1:20,5:8]
results.p.pbc.h.m.no.c$UniProtID <- sapply(results.p.pbc.h.m.no.c$protein, 
                                       extract_uniprot_id)
results.p.pbc.h.m.no.c$GeneName.UniprotID <- results.p.pbc.h.m.no.c$protein
results.p.pbc.h.m.no.c[1:20,6:10]
results.p.pbc.h.m.no.c$log2FC <- results.p.pbc.h.m.no.c$logFC
results.p.pbc.h.m.no.c[1:3,c(1,8:11)]


# volcano PBC vs Healthy and MASLD no cirrhosis

selected_proteins_p_pbc_h.m_no_c <- results.p.pbc.h.m.no.c %>%
  filter(adj.P.Val < 0.05)

selected_proteins_p_pbc_h.m_no_c[1:20,1:5]

volcano_plot_p_pbc_h.m_no_c <- ggplot(results.p.pbc.h.m.no.c) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), 
                 color = ifelse(logFC <= -0.00001  & adj.P.Val < 0.05, 'blue', 
                                ifelse(logFC >= 0.00001 & adj.P.Val < 0.05, 'red', 'grey')))) +
  geom_text(data = selected_proteins_p_pbc_h.m_no_c,
            aes(x = logFC, y = -log10(adj.P.Val), 
                label = GeneName),
            size = 3, hjust = 0, vjust = 0) +
  scale_color_identity() +
  labs(x = "log2 fold change", y = "-log10 adj p-value") +
  ggtitle("Differential Plasma Protein Expression, non-cirrhotic PBC vs. Healthy and MASLD")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray")  # Add horizontal line


volcano_plot_p_pbc_h.m_no_c

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/No cirrhosis/volcano_plasma_pbc_healthy_masld_fdr0.05_no_c.jpeg",
#         width = 10, height = 10, plot = volcano_plot_p_pbc_h.m_no_c, device = "jpeg")

# Number of upregulated proteins in PBC no cirrhosis:

selected_upregulated_proteins_p_pbc_h.m_no_c <- results.p.pbc.h.m.no.c %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_p_pbc_h.m_no_c) # 56

# Number of downregulated proteins in PBC:
selected_downregulated_proteins_p_pbc_h.m_no_c <- results.p.pbc.h.m.no.c %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_p_pbc_h.m_no_c) # 51

# Sort the upregulated proteins by adjusted p-value
sorted_upregulated_proteins_p_pbc_h.m_no_c <- 
  arrange(subset(results.p.pbc.h.m.no.c, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of upregulated proteins with gene names
print(sorted_upregulated_proteins_p_pbc_h.m_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")])

# Create a data frame with the desired output
up_plasma_pbc_vs_h.m_no_c <- 
  sorted_upregulated_proteins_p_pbc_h.m_no_c[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]

# Sort the downregulated proteins by adjusted p-value
sorted_downregulated_proteins_p_pbc_h.m_no_c <- 
  arrange(subset(results.p.pbc.h.m.no.c, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of downregulated proteins with gene names
print(sorted_downregulated_proteins_p_pbc_h.m_no_c[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
down_plasma_pbc_vs_h.m_no_c <- 
  sorted_downregulated_proteins_p_pbc_h.m_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]






####### PBC without MASLD vs Healthy and MASLD ################

table(p_no_masld$group_1)

p_pbc_h.m_no_masld <- subset(p_no_masld, group_1 %in% c("PBC","Control"))

table(p_pbc_h.m_no_masld$group_1)

design_p_pbc_h.m_no_masld <- model.matrix(~ 0 + group_1, p_pbc_h.m_no_masld)
design_p_pbc_h.m_no_masld[1:20,]

dim(p_pbc_h.m_no_masld)
p_pbc_h.m_no_masld[1:3,1:5]
p_pbc_h.m_no_masld[1:3,584:586]
p_pbc_h.m_no_masld <- p_pbc_h.m_no_masld[,3:584] # create only numeric values 
dim(p_pbc_h.m_no_masld)

fit_p_pbc_h.m_no_masld <- lmFit(t(p_pbc_h.m_no_masld), design_p_pbc_h.m_no_masld)
head(coef(fit_p_pbc_h.m_no_masld))

contr_p_pbc_h.m_no_masld <- makeContrasts(group_1PBC-group_1Control, levels=design_p_pbc_h.m_no_masld)
contr_p_pbc_h.m_no_masld
tmp_p_pbc_h.m_no_masld <- contrasts.fit(fit_p_pbc_h.m_no_masld, contr_p_pbc_h.m_no_masld)
tmp_p_pbc_h.m_no_masld <- eBayes(tmp_p_pbc_h.m_no_masld)
tmp_p_pbc_h.m_no_masld
top.table.p.pbc.h.m.no.masld <- topTable(tmp_p_pbc_h.m_no_masld,coef = 1, sort.by = "P", n = Inf)
results.p.pbc.h.m.no.masld <-  as.data.frame(top.table.p.pbc.h.m.no.masld)
results.p.pbc.h.m.no.masld$protein <- row.names(results.p.pbc.h.m.no.masld)
names(results.p.pbc.h.m.no.masld)
results.p.pbc.h.m.no.masld[1:20,1:7]
results.p.pbc.h.m.no.masld$GeneName <- sapply(results.p.pbc.h.m.no.masld$protein, 
                                          extract_gene_name)
results.p.pbc.h.m.no.masld[1:20,5:8]
results.p.pbc.h.m.no.masld$UniProtID <- sapply(results.p.pbc.h.m.no.masld$protein, 
                                           extract_uniprot_id)
results.p.pbc.h.m.no.masld$GeneName.UniprotID <- results.p.pbc.h.m.no.masld$protein
results.p.pbc.h.m.no.masld[1:20,6:10]
results.p.pbc.h.m.no.masld$log2FC <- results.p.pbc.h.m.no.masld$logFC
results.p.pbc.h.m.no.masld[1:3,c(1,8:11)]

# volcano PBC without MASLD vs Healthy and MASLD

selected_proteins_p_pbc_h.m_no_masld <- results.p.pbc.h.m.no.masld %>%
  filter(adj.P.Val < 0.05)

selected_proteins_p_pbc_h.m_no_masld[1:20,1:5]

volcano_plot_p_pbc_h.m_no_masld <- ggplot(results.p.pbc.h.m.no.masld) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), 
                 color = ifelse(logFC <= -0.00001  & adj.P.Val < 0.05, 'blue', 
                                ifelse(logFC >= 0.00001 & adj.P.Val < 0.05, 'red', 'grey')))) +
  geom_text(data = selected_proteins_p_pbc_h.m_no_masld,
            aes(x = logFC, y = -log10(adj.P.Val), 
                label = GeneName),
            size = 3, hjust = 0, vjust = 0) +
  scale_color_identity() +
  labs(x = "log2 fold change", y = "-log10 adj p-value") +
  ggtitle("Differential Plasma Protein Expression,  PBC without MASLD vs. Healthy and MASLD")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray")  # Add horizontal line


volcano_plot_p_pbc_h.m_no_masld

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/No MASLD/volcano_plasma_pbc_h.m_no_masld_fdr0.05.jpeg",
#         width = 10, height = 10, plot = volcano_plot_p_pbc_h.m_no_masld, device = "jpeg")

# Number of upregulated proteins in PBC no MASLD:

selected_upregulated_proteins_p_pbc_h.m_no_masld <- results.p.pbc.h.m.no.masld %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_p_pbc_h.m_no_masld) # 56

# Number of downregulated proteins in PBC no MASLD:
selected_downregulated_proteins_p_pbc_h.m_no_masld <- results.p.pbc.h.m.no.masld %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_p_pbc_h.m_no_masld) # 42

# Sort the upregulated proteins by adjusted p-value
sorted_upregulated_proteins_p_pbc_h.m_no_masld <- 
  arrange(subset(results.p.pbc.h.m.no.masld, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of upregulated proteins with gene names
print(sorted_upregulated_proteins_p_pbc_h.m_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")])

# Create a data frame with the desired output
up_plasma_pbc_vs_h.m_no_masld <- 
  sorted_upregulated_proteins_p_pbc_h.m_no_masld[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]

# Sort the downregulated proteins by adjusted p-value
sorted_downregulated_proteins_p_pbc_h.m_no_masld <- 
  arrange(subset(results.p.pbc.h.m.no.masld, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of downregulated proteins with gene names
print(sorted_downregulated_proteins_p_pbc_h.m_no_masld[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
down_plasma_pbc_vs_h.m_no_masld <- 
  sorted_downregulated_proteins_p_pbc_h.m_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]



############################ Liver #########################################################

###### PCA ######

#### from http://www.sthda.com/english/articles/22-principal-component-methods-videos/65-pca-in-r-using-factominer-quick-scripts-and-videos/ ###

# import baseline characteristics from "merged" document created with "pbc_baseline_characteristics.R"

pca_t <- read.csv("../data/liver_characteristics.csv")

pca_t$group <- factor(pca_t$group)
table(pca_t$group)
pca_t$group <- factor(pca_t$group, labels=c("Healthy","MASLD","PBC"))
table(pca_t$group)



# Split PBC into PBC with MASLD and PBC without MASLD
pca_t$group_1 <- pca_t$group 

pca_t$group_1 <- ifelse(pca_t$group == "PBC" & pca_t$con_masld == "Yes", "PBC with MASLD", 
                      ifelse(pca_t$group == "PBC" & pca_t$con_masld == "No", "PBC without MASLD", 
                             ifelse(pca_t$group == "Healthy","Healthy",
                                    ifelse(pca_t$group=="MASLD","MASLD",
                                           pca_t$group))))

pca_t$group_1 <- factor(pca_t$group_1)
table(pca_t$group_1)

# Split PBC by Ursochol (urso = yes/no)
pca_t$group_urso <- pca_t$group

pca_t$group_urso <- ifelse(
  pca_t$group == "PBC" & pca_t$urso == "yes", "PBC_urso_yes",
  ifelse(
    pca_t$group == "PBC" & pca_t$urso == "no",  "PBC_urso_no",
    as.character(pca_t$group)
  )
)

pca_t$group_urso <- factor(pca_t$group_urso)
table(pca_t$group_urso, useNA = "ifany")


# Subset dataframe
dim(pca_t)
pca_t[1:3,45:55] # proteins from column 52:7683
pca_t[1:3,7680:7685]
pca_t <- pca_t %>% 
  select(all_of(c("subject", "age", "sex", "bmi", "cirrose", "con_masld", "group", "group_1", "pbc_antibody", "group_urso")), 52:7683)

# Correct rest of structure in PCA 
str(pca_t)
pca_t$sex <- factor(pca_t$sex)
pca_t$cirrose <- factor(pca_t$cirrose)
pca_t$con_masld <- factor(pca_t$con_masld)
pca_t$group_1   <- factor(pca_t$group_1)
pca_t$group_urso<- factor(pca_t$group_urso)



# Now to PCA
res.pca.t <- PCA(
  pca_t,
  graph = FALSE,
  quali.sup = c("subject", "sex", "cirrose", "con_masld", 
                "group", "group_1", "pbc_antibody", "group_urso"),
  quanti.sup = c("bmi", "age")
)

  # Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component
eig.val.t <- res.pca.t$eig

# n.pc <- min(10, nrow(eig.val.t))

#barplot(
#  eig.val.t[1:n.pc, 2],
#  names.arg = 1:n.pc,
#  main = "Variances Explained by PCs (%)",
#  xlab = "Principal Components",
#  ylab = "Percentage of variance",
#  col = "steelblue"
#)

# Visualize the graph of individuals. Individuals with a similar profile are grouped together.

# By diagnosis
# jpeg(file = "../output/PCA/pca_liver_diagnosis_pbc_masld_healthy.jpeg")

fviz_pca_ind(res.pca.t, 
             geom.ind = "point", # default is point + text and text is 1:n samples
             col.ind = pca_t$group, # color by group
             palette = c("#58aadb","grey","#AA0233"),
             pointsize = 3, # increased size
             pointshape = 19, # makes circles
             invisible = "quali", # removes centroids (group mean dots/points)
             legend.title = "Diagnosis") +
  ggtitle("Individuals - PCA liver") +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),  # Increase legend text size
    axis.text = element_text(size = 14),    # Increase x- and y-axis text size
    axis.title = element_text(size = 16),   # Increase x- and y-axis title size
    plot.title = element_text(size = 18),   # Increase plot title size
    legend.position = "bottom"              # Move legend to the bottom
  )

# dev.off()

# save PCA plot with sample IDs
# jpeg(file = "../output/PCA/pca_liver_diagnosis_with_IDs.jpeg",
#      width = 3000, height = 2500, res = 300)

fviz_pca_ind(
  res.pca.t,
  geom.ind = "point",
  col.ind = pca_t$group,
  palette = c("#58aadb", "grey", "#AA0233"),
  pointsize = 3,
  pointshape = 19,
  invisible = "quali",
  legend.title = "Diagnosis"
) +
  geom_text_repel(
    aes(label = pca_t$subject),
    size = 3,
    max.overlaps = Inf
  ) +
  ggtitle("Individuals - PCA liver with sample IDs") +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.position = "bottom"
  )

dev.off()


###### 3D PCA ######

# Extract PCA coordinates
pca_3d_t <- as.data.frame(res.pca.t$ind$coord)

# Add metadata
pca_3d_t$subject <- pca_t$subject
pca_3d_t$group <- pca_t$group

# Create 3D PCA plot
pca_3d_plot_t <- plot_ly(
  data = pca_3d_t,
  x = ~Dim.1,
  y = ~Dim.2,
  z = ~Dim.3,
  color = ~group,
  colors = c("#58aadb", "grey", "#AA0233"),
  text = ~subject,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 5)
) %>%
  layout(
    title = "3D PCA liver",
    scene = list(
      xaxis = list(title = paste0("PC1 (", round(eig.val.t[1, 2], 1), "%)")),
      yaxis = list(title = paste0("PC2 (", round(eig.val.t[2, 2], 1), "%)")),
      zaxis = list(title = paste0("PC3 (", round(eig.val.t[3, 2], 1), "%)"))
    ),
    legend = list(title = list(text = "Diagnosis"))
  )

# View plot
pca_3d_plot_t


# save plot
# htmlwidgets::saveWidget(pca_3d_plot_t, "../output/PCA/pca_3d_liver_diagnosis.html",
#                         selfcontained = TRUE)


###### 3D PCA without subject IDs ######

pca_3d_plot_t_no_ids <- plot_ly(
  data = pca_3d_t,
  x = ~Dim.1,
  y = ~Dim.2,
  z = ~Dim.3,
  color = ~group,
  colors = c("#58aadb", "grey", "#AA0233"),
  text = ~paste("Group:", group),
  hoverinfo = "text",
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 5)
) %>%
  layout(
    title = "3D PCA liver",
    scene = list(
      xaxis = list(title = paste0("PC1 (", round(eig.val.t[1, 2], 1), "%)")),
      yaxis = list(title = paste0("PC2 (", round(eig.val.t[2, 2], 1), "%)")),
      zaxis = list(title = paste0("PC3 (", round(eig.val.t[3, 2], 1), "%)"))
    ),
    legend = list(title = list(text = "Diagnosis"))
  )

pca_3d_plot_t_no_ids

# save plot
# htmlwidgets::saveWidget(pca_3d_plot_t_no_ids, "../output/PCA/pca_3d_liver_diagnosis_no_IDs.html",
#                         selfcontained = TRUE)


# Extract PCA coordinates (first 5 PCs)
pca_coord_t <- as.data.frame(res.pca.t$ind$coord[, 1:5])

# Calculate Mahalanobis distance
md_t <- mahalanobis(
  pca_coord_t,
  center = colMeans(pca_coord_t),
  cov = cov(pca_coord_t)
)

# Add Mahalanobis distance and subject IDs
pca_coord_t$Mahalanobis <- md_t
pca_coord_t$subject <- pca_t$subject

# Define outlier threshold (99% confidence interval)
threshold_t <- qchisq(0.99, df = 5)

# Identify potential outliers
outliers_t <- pca_coord_t[pca_coord_t$Mahalanobis > threshold_t, ]

# View potential outliers # FALL_75, PLS_49, PLS_51
outliers_t

#Dim.1     Dim.2      Dim.3     Dim.4      Dim.5 Mahalanobis subject
#27   4.871669 104.76611 -33.916941 -3.748867  32.953226    17.72373 FALL_75
#37 225.579974  24.97052   6.193744 40.638660 -24.713131    31.85935  PLS_49
#38  -5.873517  32.39725 106.379378 -2.483169   8.207143    22.57845  PLS_51

# Define liver outliers
liver_outlier_ids <- c("FALL_75", "PLS_49", "PLS_51")


###### PCA without outliers ######

# Remove outliers
pca_t_no_outliers <- pca_t %>%
  filter(!subject %in% liver_outlier_ids)

# Check dimensions
dim(pca_t)
dim(pca_t_no_outliers)

# Check that outliers are removed
subset(pca_t_no_outliers, subject %in% liver_outlier_ids)


# Now to PCA
res.pca.t.no_outliers <- PCA(
  pca_t_no_outliers,
  graph = FALSE,
  quali.sup = c("subject", "sex", "cirrose", "con_masld",
                "group", "group_1", "pbc_antibody", "group_urso"),
  quanti.sup = c("bmi", "age")
)

eig.val.t.no_outliers <- res.pca.t.no_outliers$eig

# n.pc.t.no_outliers <- min(10, nrow(eig.val.t.no_outliers))

#barplot(
#  eig.val.t.no_outliers[1:n.pc.t.no_outliers, 2],
#  names.arg = 1:n.pc.t.no_outliers,
#  main = "Variances Explained by PCs (%)",
#  xlab = "Principal Components",
#  ylab = "Percentage of variance",
#  col = "steelblue"
#)

# Visualize the graph of individuals. Individuals with a similar profile are grouped together.

# By diagnosis

fviz_pca_ind(
  res.pca.t.no_outliers,
  geom.ind = "point",
  col.ind = pca_t_no_outliers$group,
  palette = c("#58aadb","grey","#AA0233"),
  pointsize = 3,
  pointshape = 19,
  invisible = "quali",
  legend.title = "Diagnosis"
) +
  ggtitle("Individuals - PCA liver without outliers") +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.position = "bottom"
  )

# Save plot

# ggsave("../output/PCA/pca_liver_diagnosis_without_outliers.jpeg",
#         width = 10, height = 8, device = "jpeg")

# By diagnosis and concurrent MASLD
# jpeg(file = "../output/PCA/pca_liver_diagnosis_pbc_con.masld_masld_healthy.jpeg")

fviz_pca_ind(res.pca.t, 
             geom.ind = "point", # default is point + text and text is 1:n samples
             col.ind = pca_t$group_1, # color by group
             palette = c("#dba458","#58aadb","grey","#AA0233"),             
             pointsize = 3, # increased size
             pointshape = 19, # makes circles
             invisible = "quali", # removes centroids (group mean dots/points)
             legend.title = "Diagnosis") + 
  ggtitle("Individuals - PCA liver")

# dev.off()

# By sex
# jpeg(file = "../output/PCA/pca_liver_sex_pbc_masld_healthy.jpeg")

fviz_pca_ind(res.pca.t, 
             geom.ind = "point", # default er c("point","text")
             col.ind = pca_t$sex,
             pointsize = 3,
             pointshape = 19,
             invisible = "quali",
             legend.title = "Sex") + 
  ggtitle("Individuals - PCA liver")

# dev.off()

# By cirrhosis
# jpeg(file = "../output/PCA/pca_liver_cirrhosis_pbc_masld_healthy.jpeg")

fviz_pca_ind(res.pca.t, 
             geom.ind = "point", # default er c("point","text")
             col.ind = pca_t$cirrose,
             pointsize = 3,
             pointshape = 19,
             invisible = "quali",
             legend.title = "Cirrhosis") +
  ggtitle("Individuals - PCA liver") +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),  # Increase legend text size
    axis.text = element_text(size = 14),    # Increase x- and y-axis text size
    axis.title = element_text(size = 16),   # Increase x- and y-axis title size
    plot.title = element_text(size = 18),   # Increase plot title size
    legend.position = "bottom"              # Move legend to the bottom
  )

# dev.off()

# Cirrhosis cluster together

# subset only pbc and healthy
pca_t_pbc <- subset(pca_t, group %in% c("PBC","Healthy"))
table(pca_t_pbc$group)
pca_t_pbc$group <- droplevels(pca_t_pbc$group)
table(pca_t_pbc$group)

# Now to PCA

res.pca.t.pbc <- PCA(pca_t_pbc, 
               graph = FALSE, 
               quali.sup = c("subject", "age", "sex", "bmi", "cirrose", "con_masld", "group", "group_1", "pbc_antibody", "group_urso"),
               quanti.sup = c("bmi","age"))

# Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component
eig.val.t.pbc <- res.pca.t.pbc$eig
#
#barplot(eig.val.t.pbc[1:10, 2], 
#        names.arg = 1:10, 
#        main = "Variances Explained by PCs (%)",
#        xlab = "Principal Components",
#        ylab = "Percentage of variances",
#        col ="steelblue")

# Visualize the graph of individuals. Individuals with a similar profile are grouped together.

# By diagnosis
# jpeg(file = "../output/PCA/pca_liver_diagnosis_pbc_healthy.jpeg")

fviz_pca_ind(res.pca.t.pbc, 
             geom.ind = "point", # default is point + text and text is 1:n samples
             col.ind = pca_t_pbc$group, # color by group
             palette = c("#58aadb","grey"),
             pointsize = 3, # increased size
             pointshape = 19, # makes circles
             invisible = "quali", # removes centroids (group mean dots/points)
             legend.title = "Diagnosis") +
  ggtitle("Individuals - PCA liver")

# dev.off()

# By sex
# jpeg(file = "../output/PCA/pca_liver_sex_pbc_healthy.jpeg")

fviz_pca_ind(res.pca.t.pbc, 
             geom.ind = "point", # default er c("point","text") 
             col.ind = pca_t_pbc$sex,
             pointsize = 3,
             pointshape = 19,
             invisible = "quali",
             legend.title = "Sex") +
  ggtitle("Individuals - PCA liver")

# dev.off()

# By cirrhosis
# jpeg(file = "../output/PCA/pca_liver_cirrhosis_pbc_healthy.jpeg")

fviz_pca_ind(res.pca.t.pbc, 
             geom.ind = "point", # default er c("point","text") 
             col.ind = pca_t_pbc$cirrose,
             pointsize = 3,
             pointshape = 19,
             invisible = "quali",
             legend.title = "Cirrhosis") +
  ggtitle("Individuals - PCA liver")

# dev.off()

# Cirrhosis cluster together


################## Differential expression analysis liver ########################################################

t <- read.csv("../data/liver_minprob_data_for_proteomics.tsv", sep = "\t")

t[1:10, 1:10] #instead of head because the df is too big
dim(t)
str(t[1:5])
t$group <- as.factor(t$group)
table(t$group)

# Check if there are there any missing values
sum(is.na(t))

# Check if any neg values
any(t[,3:7634] < 0)

# Preparing for other group comparisons

t$group_1 <- t$group

t <- t %>%
  mutate(group_1 = case_when(
    group_1 %in% c("Healthy", "NAFLD") ~ "Control",    # Healthy and MASLD as a single group
    TRUE ~ group_1
  ))

# Check the new group counts
table(t$group_1)

t$group_2 <- t$group

t <- t %>%
  mutate(group_2 = case_when(
    group_2 %in% c("Healthy", "NAFLD","AIH","PSC") ~ "Others",   # All groups except PBC as a single group
    TRUE ~ group_2
  ))

# Check the new group counts
table(t$group_2)

# Dataset without outliers
t_no_outliers <- subset(t, !subject %in% liver_outlier_ids)

# Check dimensions
dim(t)
dim(t_no_outliers)

# Check that outliers are removed
subset(t_no_outliers, subject %in% liver_outlier_ids)


# Dataset without cirrhosis

pca_t$subject[pca_t$cirrose=="Yes"]

t_no_c <- subset(t, ! subject %in% c("FALL_6","FALL_75"))
dim(t)
dim(t_no_c)

# Dataset PBC without MASLD 

t_no_masld <- t %>%
  filter(!(subject %in% pca_t$subject[pca_t$group_1 == "PBC with MASLD"]))

dim(t)
dim(t_no_masld)


# Differential expression PBC vs Healthy 

t_pbc <- subset(t, group %in% c("PBC","Healthy"))
levels(t_pbc$group)
table(t_pbc$group)
t_pbc$group <- droplevels(t_pbc$group)
table(t_pbc$group)
t_pbc[1:20,1:5]
dim(t_pbc)

design_t_pbc <- model.matrix(~ 0 + group, t_pbc)
design_t_pbc[1:20,]

dim(t_pbc)
t_pbc[1:3,1:5]
t_pbc[1:3,7632:7636]

t_pbc <- t_pbc[,3:7634] # create only numeric values
dim(t_pbc)

# lmFit expects input array to have structure: protein x sample
# lmFit fits a linear model using weighted least squares for each protein:
fit_t_pbc <- lmFit(t(t_pbc), design_t_pbc)
head(coef(fit_t_pbc))

# Comparisons between groups (log fold-changes) are obtained as contrasts of
# these fitted linear models:
# Samples are grouped based on experimental condition
# The variability of protein expression is compared between these groups
contr_t_pbc <- makeContrasts(groupPBC-groupHealthy, levels=design_t_pbc)

contr_t_pbc

# Estimate contrast for each protein
tmp_t_pbc <- contrasts.fit(fit_t_pbc, contr_t_pbc)

# Empirical Bayes smoothing of standard errors (shrinks standard errors
# that are much larger or smaller than those from other proteins towards the average standard error)
tmp_t_pbc <- eBayes(tmp_t_pbc)
# could also set robust = T and trend = T in eBayes (gives somewhat same result). 
# Check: https://support.bioconductor.org/p/56560/

### Diagnostic plots:

# Scatterplot of residual-variances vs average log-expression
# Save plot
# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy/All/plotSA_pbc_vs_healthy_liver.jpeg")
# plotSA(tmp_t_pbc, main = "Residual variances vs average log-expression, PBC vs. Healthy, Liver") 
# dev.off()

# Mean Difference plot. Log-intensity ratios (differences, y) versus log-intensity averages (means, x):
# Save plot
# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy/All/MD-plot_pbc_vs_healthy_liver.jpeg")
# plotMD(tmp_t_pbc, main = "Mean-difference plot, PBC vs. Healthy, Liver") 
# dev.off()

# Extract results
top.table.t.pbc <- topTable(tmp_t_pbc, coef = 1, sort.by = "P", n = Inf) #default is "BH" which is alias of "fdr"
results.t.pbc <-  as.data.frame(top.table.t.pbc)
results.t.pbc$protein <- row.names(results.t.pbc)
names(results.t.pbc)
results.t.pbc[170:190,1:7]
results.t.pbc$GeneName <- sapply(results.t.pbc$protein, 
                                 extract_gene_name)
results.t.pbc[170:190,5:8]
results.t.pbc$UniProtID <- sapply(results.t.pbc$protein, 
                                  extract_uniprot_id)
results.t.pbc$GeneName.UniprotID <- results.t.pbc$protein
results.t.pbc[170:190,6:10]
results.t.pbc$log2FC <- results.t.pbc$logFC
results.t.pbc[1:3,c(1,8:11)]

#### Create plot PBC tissue

# Identify proteins that meet the filtering criteria

selected_proteins_t_pbc <- results.t.pbc %>%
  filter(adj.P.Val < 0.05)

# As volcano plot cannot print all the names and becomes impossible to view, I filter only those proteins
# that are regulated more than 1 or less than -1
selected_proteins_t_pbc_volcano <- results.t.pbc %>%
  filter(logFC >= 1 & adj.P.Val < 0.05 | logFC <= -1 & adj.P.Val < 0.05)

# Volcano plot

volcano_plot_t_pbc <- ggplot(results.t.pbc) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val),
                 color = ifelse(logFC <= -0.0001 & adj.P.Val < 0.05, 'blue', 
                ifelse(logFC >= 0.0001 & adj.P.Val < 0.05, 'red', 'grey')))) +
  geom_text(data = selected_proteins_t_pbc_volcano,
            aes(x = logFC, 
                y = -log10(adj.P.Val), 
                label = GeneName),
            size = 3, hjust = 0, vjust = 0) +
  scale_color_identity() +
  labs(x = "log2(FC)", y = "-log10(adj P value)") +
  ggtitle("Differential Liver Protein Expression, PBC vs. Healthy")+
  annotate("segment", x = -Inf, xend = Inf, y = -log10(0.05), yend = -log10(0.05),
           linetype = "dashed", color = "gray") +  # Add horizontal line
  theme(
    axis.title = element_text(size = 20),  # Increase axis title size
    axis.text = element_text(size = 18),   # Increase axis text (tick labels) size
    plot.title = element_text(size = 20)   # Increase plot title size
  )

volcano_plot_t_pbc

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/All/volcano_liver_pbc_healthy_fdr0.05.jpeg",
#         width = 10, height = 10, plot = volcano_plot_t_pbc, device = "jpeg")

# Number of upregulated proteins in PBC:

selected_upregulated_proteins_t_pbc <- results.t.pbc %>%
  filter((logFC >= 0.0001 & adj.P.Val < 0.05))
nrow(selected_upregulated_proteins_t_pbc) # 1341

# Number of downregulated proteins in PBC:
selected_downregulated_proteins_t_pbc <- results.t.pbc %>%
  filter(logFC <= -0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_t_pbc) # 1180

# Sort the upregulated proteins by adjusted p-value
sorted_upregulated_proteins_t_pbc <- 
  arrange(subset(results.t.pbc, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of upregulated proteins with gene names
print(sorted_upregulated_proteins_t_pbc[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
up_liver_pbc_vs_healthy <- 
  sorted_upregulated_proteins_t_pbc[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]

# Sort the downregulated proteins by adjusted p-value
sorted_downregulated_proteins_t_pbc <- 
  arrange(subset(results.t.pbc, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of downregulated proteins with gene names
print(sorted_downregulated_proteins_t_pbc[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
down_liver_pbc_vs_healthy <- 
  sorted_downregulated_proteins_t_pbc[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]




################## Differential expression PBC vs Healthy - without outliers ##################

t_pbc_no_outliers <- subset(t_no_outliers, group %in% c("PBC","Healthy"))

levels(t_pbc_no_outliers$group)

table(t_pbc_no_outliers$group)

t_pbc_no_outliers$group <- droplevels(t_pbc_no_outliers$group)

table(t_pbc_no_outliers$group)

t_pbc_no_outliers[1:20,1:5]

dim(t_pbc_no_outliers)


design_t_pbc_no_outliers <- model.matrix(~ 0 + group, t_pbc_no_outliers)

design_t_pbc_no_outliers[1:20,]


dim(t_pbc_no_outliers)

t_pbc_no_outliers[1:3,1:5]

t_pbc_no_outliers[1:3,7632:7636]


t_pbc_no_outliers <- t_pbc_no_outliers[,3:7634] # create only numeric values

dim(t_pbc_no_outliers)


# lmFit expects input array to have structure: protein x sample
# lmFit fits a linear model using weighted least squares for each protein:

fit_t_pbc_no_outliers <- lmFit(
  t(t_pbc_no_outliers),
  design_t_pbc_no_outliers
)

head(coef(fit_t_pbc_no_outliers))


# Comparisons between groups (log fold-changes) are obtained as contrasts of
# these fitted linear models:
# Samples are grouped based on experimental condition
# The variability of protein expression is compared between these groups

contr_t_pbc_no_outliers <- makeContrasts(
  groupPBC-groupHealthy,
  levels=design_t_pbc_no_outliers
)

contr_t_pbc_no_outliers


# Estimate contrast for each protein

tmp_t_pbc_no_outliers <- contrasts.fit(
  fit_t_pbc_no_outliers,
  contr_t_pbc_no_outliers
)


# Empirical Bayes smoothing of standard errors

tmp_t_pbc_no_outliers <- eBayes(tmp_t_pbc_no_outliers)

tmp_t_pbc_no_outliers


### Diagnostic plots:

# Scatterplot of residual-variances vs average log-expression
# Save plot

# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy/Without outliers/plotSA_pbc_vs_healthy_without_outliers_liver.jpeg")

#plotSA(
#  tmp_t_pbc_no_outliers,
#  main = "Residual variances vs average log-expression, PBC vs. Healthy without outliers, Liver"
#)

#dev.off()


# Mean Difference plot. Log-intensity ratios (differences, y) versus log-intensity averages (means, x):
# Save plot

# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy/Without outliers/MD-plot_pbc_vs_healthy_without_outliers_liver.jpeg")

#plotMD(
#  tmp_t_pbc_no_outliers,
#  main = "Mean-difference plot, PBC vs. Healthy without outliers, Liver"
#)

#dev.off()


# Extract results

top.table.t.pbc.no_outliers <- topTable(
  tmp_t_pbc_no_outliers,
  coef = 1,
  sort.by = "P",
  n = Inf
)

results.t.pbc.no_outliers <- as.data.frame(top.table.t.pbc.no_outliers)

results.t.pbc.no_outliers$protein <- row.names(results.t.pbc.no_outliers)

names(results.t.pbc.no_outliers)

results.t.pbc.no_outliers[170:190,1:7]

results.t.pbc.no_outliers$GeneName <- sapply(
  results.t.pbc.no_outliers$protein,
  extract_gene_name
)

results.t.pbc.no_outliers[170:190,5:8]

results.t.pbc.no_outliers$UniProtID <- sapply(
  results.t.pbc.no_outliers$protein,
  extract_uniprot_id
)

results.t.pbc.no_outliers$GeneName.UniprotID <- 
  results.t.pbc.no_outliers$protein

results.t.pbc.no_outliers[170:190,6:10]

results.t.pbc.no_outliers$log2FC <- 
  results.t.pbc.no_outliers$logFC

results.t.pbc.no_outliers[1:3,c(1,8:11)]


#### Create plot PBC tissue ####

# Identify proteins that meet the filtering criteria

selected_proteins_t_pbc_no_outliers <- 
  results.t.pbc.no_outliers %>%
  filter(adj.P.Val < 0.05)


# As volcano plot cannot print all the names and becomes impossible to view,
# I filter only those proteins that are regulated more than 1 or less than -1

selected_proteins_t_pbc_volcano_no_outliers <- 
  results.t.pbc.no_outliers %>%
  filter(logFC >= 1 & adj.P.Val < 0.05 |
           logFC <= -1 & adj.P.Val < 0.05)


# Volcano plot

volcano_plot_t_pbc_no_outliers <- ggplot(results.t.pbc.no_outliers) +
  geom_point(aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = ifelse(
      logFC <= -0.0001 & adj.P.Val < 0.05,
      'blue',
      ifelse(
        logFC >= 0.0001 & adj.P.Val < 0.05,
        'red',
        'grey'
      )
    )
  )) +
  geom_text(
    data = selected_proteins_t_pbc_volcano_no_outliers,
    aes(
      x = logFC,
      y = -log10(adj.P.Val),
      label = GeneName
    ),
    size = 3,
    hjust = 0,
    vjust = 0
  ) +
  scale_color_identity() +
  labs(x = "log2(FC)", y = "-log10(adj P value)") +
  ggtitle("Differential Liver Protein Expression, PBC vs. Healthy without outliers") +
  annotate(
    "segment",
    x = -Inf,
    xend = Inf,
    y = -log10(0.05),
    yend = -log10(0.05),
    linetype = "dashed",
    color = "gray"
  ) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 20)
  )

volcano_plot_t_pbc_no_outliers


# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/Without outliers/volcano_liver_pbc_healthy_without_outliers_fdr0.05.jpeg",
#         width = 10, height = 10, plot = volcano_plot_t_pbc_no_outliers, device = "jpeg")


# Number of upregulated proteins in PBC without outliers: #1216

selected_upregulated_proteins_t_pbc_no_outliers <- 
  results.t.pbc.no_outliers %>%
  filter((logFC >= 0.0001 & adj.P.Val < 0.05))

nrow(selected_upregulated_proteins_t_pbc_no_outliers)


# Number of downregulated proteins in PBC without outliers: #980

selected_downregulated_proteins_t_pbc_no_outliers <- 
  results.t.pbc.no_outliers %>%
  filter(logFC <= -0.0001 & adj.P.Val < 0.05)

nrow(selected_downregulated_proteins_t_pbc_no_outliers)


# Sort the upregulated proteins by adjusted p-value

sorted_upregulated_proteins_t_pbc_no_outliers <- 
  arrange(
    subset(results.t.pbc.no_outliers,
           logFC >= 0.0001 & adj.P.Val < 0.05),
    adj.P.Val
  )


# Display the sorted list of upregulated proteins with gene names

print(
  sorted_upregulated_proteins_t_pbc_no_outliers[, 
                                                c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]
)


# Create a data frame with the desired output

up_liver_pbc_vs_healthy_no_outliers <- 
  sorted_upregulated_proteins_t_pbc_no_outliers[, 
                                                c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]


# Sort the downregulated proteins by adjusted p-value

sorted_downregulated_proteins_t_pbc_no_outliers <- 
  arrange(
    subset(results.t.pbc.no_outliers,
           logFC <= -0.0001 & adj.P.Val < 0.05),
    adj.P.Val
  )


# Display the sorted list of downregulated proteins with gene names

print(
  sorted_downregulated_proteins_t_pbc_no_outliers[, 
                                                  c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]
)


# Create a data frame with the desired output

down_liver_pbc_vs_healthy_no_outliers <- 
  sorted_downregulated_proteins_t_pbc_no_outliers[, 
                                                  c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]




###### PBC vs. Healthy, no cirrhosis

table(t_no_c$group)

t_pbc_no_c <- subset(t_no_c, group %in% c("PBC","Healthy"))
levels(t_pbc_no_c$group)
table(t_pbc_no_c$group)
t_pbc_no_c$group <- droplevels(t_pbc_no_c$group)
table(t_pbc_no_c$group)

design_t_pbc_no_c <- model.matrix(~ 0 + group, t_pbc_no_c)
design_t_pbc_no_c[1:20,] 

dim(t_pbc_no_c)
t_pbc_no_c[1:3,1:5]
t_pbc_no_c[1:3,7632:7636]
t_pbc_no_c <- t_pbc_no_c[,3:7634] # create only numeric values
dim(t_pbc_no_c)

# lmFit expects input array to have structure: protein x sample
# lmFit fits a linear model using weighted least squares for each protein:
fit_t_pbc_no_c <- lmFit(t(t_pbc_no_c), design_t_pbc_no_c)
head(coef(fit_t_pbc_no_c))

# Comparisons between groups (log fold-changes) are obtained as contrasts of
# these fitted linear models:
# Samples are grouped based on experimental condition
# The variability of protein expression is compared between these groups
contr_t_pbc_no_c <- makeContrasts(groupPBC-groupHealthy, levels=design_t_pbc_no_c)

contr_t_pbc_no_c

# Estimate contrast for each protein
tmp_t_pbc_no_c <- contrasts.fit(fit_t_pbc_no_c, contr_t_pbc_no_c)

# Empirical Bayes smoothing of standard errors (shrinks standard errors
# that are much larger or smaller than those from other proteins towards the average standard error)
tmp_t_pbc_no_c <- eBayes(tmp_t_pbc_no_c)
# could also set robust = T and trend = T in eBayes (gives somewhat same result). 
# Check: https://support.bioconductor.org/p/56560/

### Diagnostic plots:
# Scatterplot of residual-variances vs average log-expression
# save plot
# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy/No cirrhosis/plotSA_pbc_vs_healthy_liver_no.cirrhosis.jpeg")

# plotSA(tmp_t_pbc_no_c, main = "Residual variances vs average log-expression, non-cirrhotic PBC vs. Healthy, Liver") 
# dev.off()

# Mean Difference plot. Log-intensity ratios (differences, y) versus log-intensity averages (means, x):
# save plot
# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy/No cirrhosis/MD-plot_pbc_vs_healthy_liver_no.cirrhosis.jpeg")
# plotMD(tmp_t_pbc_no_c, main = "Mean-difference plot, non-cirrhotic PBC vs. Healthy, Liver") 
# dev.off()

# Extract results
top.table.t.pbc.no.c <- topTable(tmp_t_pbc_no_c, coef = 1, sort.by = "P", n = Inf) #default is "BH" which is alias of "fdr"
results.t.pbc.no.c <-  as.data.frame(top.table.t.pbc.no.c)
results.t.pbc.no.c$protein <- row.names(results.t.pbc.no.c)
names(results.t.pbc.no.c)
results.t.pbc.no.c[1:20,1:7]
results.t.pbc.no.c$GeneName <- sapply(results.t.pbc.no.c$protein, 
                                      extract_gene_name)
results.t.pbc.no.c[1:20,5:8]
results.t.pbc.no.c$UniProtID <- sapply(results.t.pbc.no.c$protein, 
                                       extract_uniprot_id)
results.t.pbc.no.c$GeneName.UniprotID <- results.t.pbc.no.c$protein
results.t.pbc.no.c[1:20,6:10]
results.t.pbc.no.c$log2FC <- results.t.pbc.no.c$logFC
results.t.pbc.no.c[1:3,c(1,8:11)]





########### Differential expression PBC vs MASLD ###################


t_pbc_masld <- subset(t, group %in% c("PBC","NAFLD"))
levels(t_pbc_masld$group)
table(t_pbc_masld$group)
t_pbc_masld$group <- droplevels(t_pbc_masld$group)
table(t_pbc_masld$group)
t_pbc_masld[1:20,1:5]
dim(t_pbc_masld)

design_t_pbc_masld <- model.matrix(~ 0 + group, t_pbc_masld)
design_t_pbc_masld[1:20,]

dim(t_pbc_masld)
t_pbc_masld[1:3,1:5]
t_pbc_masld[1:3,7632:7636]

t_pbc_masld <- t_pbc_masld[,3:7634] # create only numeric values
dim(t_pbc_masld)

fit_t_pbc_masld <- lmFit(t(t_pbc_masld), design_t_pbc_masld)
head(coef(fit_t_pbc_masld))

contr_t_pbc_masld <- makeContrasts(groupPBC-groupNAFLD, levels=design_t_pbc_masld)
contr_t_pbc_masld

tmp_t_pbc_masld <- contrasts.fit(fit_t_pbc_masld, contr_t_pbc_masld)
tmp_t_pbc_masld <- eBayes(tmp_t_pbc_masld)

top.table.t.pbc.masld <- topTable(tmp_t_pbc_masld, coef = 1, sort.by = "P", n = Inf)
results.t.pbc.masld <- as.data.frame(top.table.t.pbc.masld)

results.t.pbc.masld$protein <- row.names(results.t.pbc.masld)
names(results.t.pbc.masld)
results.t.pbc.masld[170:190,1:7]

results.t.pbc.masld$GeneName <- sapply(results.t.pbc.masld$protein, 
                                       extract_gene_name)
results.t.pbc.masld[170:190,5:8]

results.t.pbc.masld$UniProtID <- sapply(results.t.pbc.masld$protein, 
                                        extract_uniprot_id)

results.t.pbc.masld$GeneName.UniprotID <- results.t.pbc.masld$protein
results.t.pbc.masld[170:190,6:10]

results.t.pbc.masld$log2FC <- results.t.pbc.masld$logFC
results.t.pbc.masld[1:3,c(1,8:11)]

#### Create plot PBC tissue

selected_proteins_t_pbc_masld <- results.t.pbc.masld %>%
  filter(adj.P.Val < 0.05)

selected_proteins_t_pbc_masld_volcano <- results.t.pbc.masld %>%
  filter(logFC >= 1 & adj.P.Val < 0.05 | logFC <= -1 & adj.P.Val < 0.05)

volcano_plot_t_pbc_masld <- ggplot(results.t.pbc.masld) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val),
                 color = ifelse(logFC <= -0.0001 & adj.P.Val < 0.05, 'blue', 
                                ifelse(logFC >= 0.0001 & adj.P.Val < 0.05, 'red', 'grey')))) +
  geom_text(data = selected_proteins_t_pbc_masld_volcano,
            aes(x = logFC, 
                y = -log10(adj.P.Val), 
                label = GeneName),
            size = 3, hjust = 0, vjust = 0) +
  scale_color_identity() +
  labs(x = "log2(FC)", y = "-log10(adj P value)") +
  ggtitle("Differential Liver Protein Expression, PBC vs. MASLD")+
  annotate("segment", x = -Inf, xend = Inf, y = -log10(0.05), yend = -log10(0.05),
           linetype = "dashed", color = "gray") +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 20)
  )

volcano_plot_t_pbc_masld

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_MASLD/All/volcano_liver_pbc_masld_fdr0.05.jpeg",
#         width = 10, height = 10, plot = volcano_plot_t_pbc_masld, device = "jpeg")

# Number of upregulated proteins in PBC:

selected_upregulated_proteins_t_pbc_masld <- results.t.pbc.masld %>%
  filter((logFC >= 0.0001 & adj.P.Val < 0.05))
nrow(selected_upregulated_proteins_t_pbc_masld)

# Number of downregulated proteins in PBC:

selected_downregulated_proteins_t_pbc_masld <- results.t.pbc.masld %>%
  filter(logFC <= -0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_t_pbc_masld)

# Sort the upregulated proteins by adjusted p-value

sorted_upregulated_proteins_t_pbc_masld <- 
  arrange(subset(results.t.pbc.masld, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

print(sorted_upregulated_proteins_t_pbc_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")])

up_liver_pbc_vs_masld <- 
  sorted_upregulated_proteins_t_pbc_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]

# Sort the downregulated proteins by adjusted p-value

sorted_downregulated_proteins_t_pbc_masld <- 
  arrange(subset(results.t.pbc.masld, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

print(sorted_downregulated_proteins_t_pbc_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")])

down_liver_pbc_vs_masld <- 
  sorted_downregulated_proteins_t_pbc_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]


# write.xlsx(list("downreg_plasma_pbc_vs_masld" = down_plasma_pbc_vs_masld,
#                 "upreg_plasma_pbc_vs_masld" = up_plasma_pbc_vs_masld,
#                 "downreg_liver_pbc_vs_masld" = down_liver_pbc_vs_masld,
#                 "upreg_liver_pbc_vs_masld" = up_liver_pbc_vs_masld),
#            "../output/FDR0.05_all/PBC_vs_MASLD/All/DE_proteins_fdr0.05_pbc_vs_masld.xlsx")

###### PBC vs. MASLD, no cirrhosis

table(t_no_c$group)

t_pbc_masld_no_c <- subset(t_no_c, group %in% c("PBC","NAFLD"))
levels(t_pbc_masld_no_c$group)
table(t_pbc_masld_no_c$group)
t_pbc_masld_no_c$group <- droplevels(t_pbc_masld_no_c$group)
table(t_pbc_masld_no_c$group)

design_t_pbc_masld_no_c <- model.matrix(~ 0 + group, t_pbc_masld_no_c)
design_t_pbc_masld_no_c[1:20,]

dim(t_pbc_masld_no_c)
t_pbc_masld_no_c[1:3,1:5]
t_pbc_masld_no_c[1:3,7632:7636]
t_pbc_masld_no_c <- t_pbc_masld_no_c[,3:7634] # create only numeric values
dim(t_pbc_masld_no_c)

fit_t_pbc_masld_no_c <- lmFit(t(t_pbc_masld_no_c), design_t_pbc_masld_no_c)
head(coef(fit_t_pbc_masld_no_c))

contr_t_pbc_masld_no_c <- makeContrasts(groupPBC-groupNAFLD, levels=design_t_pbc_masld_no_c)
contr_t_pbc_masld_no_c
tmp_t_pbc_masld_no_c <- contrasts.fit(fit_t_pbc_masld_no_c, contr_t_pbc_masld_no_c)
tmp_t_pbc_masld_no_c <- eBayes(tmp_t_pbc_masld_no_c)
top.table.t.pbc.masld.no.c <- topTable(tmp_t_pbc_masld_no_c, coef = 1, sort.by = "P", n = Inf)
results.t.pbc.masld.no.c <-  as.data.frame(top.table.t.pbc.masld.no.c)
results.t.pbc.masld.no.c$protein <- row.names(results.t.pbc.masld.no.c)
names(results.t.pbc.masld.no.c)
results.t.pbc.masld.no.c[1:20,1:7]
results.t.pbc.masld.no.c$GeneName <- sapply(results.t.pbc.masld.no.c$protein, 
                                      extract_gene_name)
results.t.pbc.masld.no.c[1:20,5:8]
results.t.pbc.masld.no.c$UniProtID <- sapply(results.t.pbc.masld.no.c$protein, 
                                       extract_uniprot_id)
results.t.pbc.masld.no.c$GeneName.UniprotID <- results.t.pbc.masld.no.c$protein
results.t.pbc.masld.no.c[1:20,6:10]
results.t.pbc.masld.no.c$log2FC <- results.t.pbc.masld.no.c$logFC
results.t.pbc.masld.no.c[1:3,c(1,8:11)]

#### volcano PBC vs. MASLD tissue

selected_proteins_t_pbc_masld_no_c <- results.t.pbc.masld.no.c %>%
  filter(adj.P.Val < 0.05)

# As volcano plot cannot print all the names and becomes impossible to view, I filter only those proteins
# that are regulated more than 1 or less than -1
selected_proteins_t_pbc_masld_volcano_no_c <- results.t.pbc.masld.no.c %>%
  filter(logFC >= 1 & adj.P.Val < 0.05 | logFC <= -1 & adj.P.Val < 0.05)

volcano_plot_t_pbc_masld_no_c <- ggplot(results.t.pbc.masld.no.c) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val),
                 color = ifelse(logFC <= -0.0001 & adj.P.Val < 0.05, 'blue', 
                                ifelse(logFC >= 0.0001 & adj.P.Val < 0.05, 'red', 'grey')))) +
  geom_text(data = selected_proteins_t_pbc_masld_volcano_no_c,
            aes(x = logFC, 
                y = -log10(adj.P.Val), 
                label = GeneName),
            size = 3, hjust = 0, vjust = 0) +
  scale_color_identity() +
  labs(x = "log2(FC)", y = "-log10(adj P value)") +
  ggtitle("Differential Liver Protein Expression, non-cirrhotic PBC vs. MASLD")+
  annotate("segment", x = -Inf, xend = Inf, y = -log10(0.05), yend = -log10(0.05),
           linetype = "dashed", color = "gray") +  # Add horizontal line
  theme(
    axis.title = element_text(size = 20),  # Increase axis title size
    axis.text = element_text(size = 18),   # Increase axis text (tick labels) size
    plot.title = element_text(size = 20)   # Increase plot title size
  )

volcano_plot_t_pbc_masld_no_c

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_MASLD/No cirrhosis/volcano_liver_pbc_masld_fdr0.05_no.cirrhosis.jpeg",
#         width = 10, height = 10, plot = volcano_plot_t_pbc_masld_no_c, device = "jpeg")

# Number of upregulated proteins in PBC:

selected_upregulated_proteins_t_pbc_masld_no_c <- results.t.pbc.masld.no.c %>%
  filter((logFC >= 0.0001 & adj.P.Val < 0.05))
nrow(selected_upregulated_proteins_t_pbc_masld_no_c)

# Number of downregulated proteins in PBC:

selected_downregulated_proteins_t_pbc_masld_no_c <- results.t.pbc.masld.no.c %>%
  filter(logFC <= -0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_t_pbc_masld_no_c)

# Sort the upregulated proteins by adjusted p-value
sorted_upregulated_proteins_t_pbc_masld_no_c <- 
  arrange(subset(results.t.pbc.masld.no.c, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of upregulated proteins with gene names
print(sorted_upregulated_proteins_t_pbc_masld_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
up_liver_pbc_vs_masld_no_c <- 
  sorted_upregulated_proteins_t_pbc_masld_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]

# Sort the downregulated proteins by adjusted p-value
sorted_downregulated_proteins_t_pbc_masld_no_c <- 
  arrange(subset(results.t.pbc.masld.no.c, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of downregulated proteins with gene names
print(sorted_downregulated_proteins_t_pbc_masld_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
down_liver_pbc_vs_masld_no_c <- 
  sorted_downregulated_proteins_t_pbc_masld_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]

# write.xlsx(list("down_plasma_pbc_masld_no_c" = down_plasma_pbc_vs_masld_no_c,
#                 "up_plasma_pbc_masld_no_c" = up_plasma_pbc_vs_masld_no_c,
#                 "down_liver_pbc_masld_no_c" = down_liver_pbc_vs_masld_no_c,
#                 "up_liver_pbc_masld_no_c" = up_liver_pbc_vs_masld_no_c),
#            "../output/FDR0.05_all/PBC_vs_MASLD/No cirrhosis/DE_proteins_fdr0.05_pbc_vs_masld_no_c.xlsx")

###### PBC without MASLD vs. Healthy

table(t_no_masld$group)

t_pbc_no_masld <- subset(t_no_masld, group %in% c("PBC","Healthy"))
levels(t_pbc_no_masld$group)
table(t_pbc_no_masld$group)
t_pbc_no_masld$group <- droplevels(t_pbc_no_masld$group)
table(t_pbc_no_masld$group)

design_t_pbc_no_masld <- model.matrix(~ 0 + group, t_pbc_no_masld)
design_t_pbc_no_masld[1:17,] 

dim(t_pbc_no_masld)
t_pbc_no_masld[1:3,1:5]
t_pbc_no_masld[1:3,7632:7636]
t_pbc_no_masld <- t_pbc_no_masld[,3:7634] # create only numeric values
dim(t_pbc_no_masld)

# lmFit expects input array to have structure: protein x sample
# lmFit fits a linear model using weighted least squares for each protein:
fit_t_pbc_no_masld <- lmFit(t(t_pbc_no_masld), design_t_pbc_no_masld)
head(coef(fit_t_pbc_no_masld))

# Comparisons between groups (log fold-changes) are obtained as contrasts of
# these fitted linear models:
# Samples are grouped based on experimental condition
# The variability of protein expression is compared between these groups
contr_t_pbc_no_masld <- makeContrasts(groupPBC-groupHealthy, levels=design_t_pbc_no_masld)

contr_t_pbc_no_masld

# Estimate contrast for each protein
tmp_t_pbc_no_masld <- contrasts.fit(fit_t_pbc_no_masld, contr_t_pbc_no_masld)

# Empirical Bayes smoothing of standard errors (shrinks standard errors
# that are much larger or smaller than those from other proteins towards the average standard error)
tmp_t_pbc_no_masld <- eBayes(tmp_t_pbc_no_masld)
# could also set robust = T and trend = T in eBayes (gives somewhat same result). 
# Check: https://support.bioconductor.org/p/56560/

### Diagnostic plots:
# Scatterplot of residual-variances vs average log-expression
# save plot
# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy/No MASLD/plotSA_pbc_vs_healthy_liver_no.masld.jpeg")
#plotSA(tmp_t_pbc_no_masld, main = "Residual variances vs average log-expression, PBC without vs. Healthy, Liver") 
#dev.off()

# Mean Difference plot. Log-intensity ratios (differences, y) versus log-intensity averages (means, x):
# save plot
# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy/No MASLD/MD-plot_pbc_vs_healthy_liver_no.masld.jpeg")
#plotMD(tmp_t_pbc_no_masld, main = "Mean-difference plot, PBC without MASLD vs. Healthy, Liver") 
#dev.off()

# Extract results
top.table.t.pbc.no.masld <- topTable(tmp_t_pbc_no_masld, coef = 1, sort.by = "P", n = Inf) #default is "BH" which is alias of "fdr"
results.t.pbc.no.masld <-  as.data.frame(top.table.t.pbc.no.masld)
results.t.pbc.no.masld$protein <- row.names(results.t.pbc.no.masld)
names(results.t.pbc.no.masld)
results.t.pbc.no.masld[1:20,1:7]
results.t.pbc.no.masld$GeneName <- sapply(results.t.pbc.no.masld$protein, 
                                      extract_gene_name)
results.t.pbc.no.masld[1:20,5:8]
results.t.pbc.no.masld$UniProtID <- sapply(results.t.pbc.no.masld$protein, 
                                       extract_uniprot_id)
results.t.pbc.no.masld$GeneName.UniprotID <- results.t.pbc.no.masld$protein
results.t.pbc.no.masld[1:20,6:10]
results.t.pbc.no.masld$log2FC <- results.t.pbc.no.masld$logFC
results.t.pbc.no.masld[1:3,c(1,8:11)]




#### volcano PBC vs. healthy tissue - no cirrhosis  

# Identify proteins that meet the filtering criteria

selected_proteins_t_pbc_no_c <- results.t.pbc.no.c %>%
  filter(adj.P.Val < 0.05)

# As volcano plot cannot print all the names and becomes impossible to view, I filter only those proteins
# that are regulated more than 1 or less than -1
selected_proteins_t_pbc_volcano_no_c <- results.t.pbc.no.c %>%
  filter(logFC >= 1 & adj.P.Val < 0.05 | logFC <= -1 & adj.P.Val < 0.05)

volcano_plot_t_pbc_no_c <- ggplot(results.t.pbc.no.c) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val),
                 color = ifelse(logFC <= -0.0001 & adj.P.Val < 0.05, 'blue', 
                                ifelse(logFC >= 0.0001 & adj.P.Val < 0.05, 'red', 'grey')))) +
  geom_text(data = selected_proteins_t_pbc_volcano_no_c,
            aes(x = logFC, 
                y = -log10(adj.P.Val), 
                label = GeneName),
            size = 3, hjust = 0, vjust = 0) +
  scale_color_identity() +
  labs(x = "log2(FC)", y = "-log10(adj P value)") +
  ggtitle("Differential Liver Protein Expression, non-cirrhotic PBC vs. Healthy")+
  annotate("segment", x = -Inf, xend = Inf, y = -log10(0.05), yend = -log10(0.05),
           linetype = "dashed", color = "gray") +  # Add horizontal line
  theme(
    axis.title = element_text(size = 20),  # Increase axis title size
    axis.text = element_text(size = 18),   # Increase axis text (tick labels) size
    plot.title = element_text(size = 20)   # Increase plot title size
  )

volcano_plot_t_pbc_no_c

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/No cirrhosis/volcano_liver_pbc_healthy_fdr0.05_no.cirrhosis.jpeg",
#         width = 10, height = 10, plot = volcano_plot_t_pbc_no_c, device = "jpeg")

# Number of upregulated proteins in PBC:

selected_upregulated_proteins_t_pbc_no_c <- results.t.pbc.no.c %>%
  filter((logFC >= 0.0001 & adj.P.Val < 0.05))
nrow(selected_upregulated_proteins_t_pbc_no_c) # 1039

# Number of downregulated proteins in PBC:
selected_downregulated_proteins_t_pbc_no_c <- results.t.pbc.no.c %>%
  filter(logFC <= -0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_t_pbc_no_c) # 1272

# Sort the upregulated proteins by adjusted p-value
sorted_upregulated_proteins_t_pbc_no_c <- 
  arrange(subset(results.t.pbc.no.c, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of upregulated proteins with gene names
print(sorted_upregulated_proteins_t_pbc_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
up_liver_pbc_vs_healthy_no_c <- 
  sorted_upregulated_proteins_t_pbc_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]

# Sort the downregulated proteins by adjusted p-value
sorted_downregulated_proteins_t_pbc_no_c <- 
  arrange(subset(results.t.pbc.no.c, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of downregulated proteins with gene names
print(sorted_downregulated_proteins_t_pbc_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
down_liver_pbc_vs_healthy_no_c <- 
  sorted_downregulated_proteins_t_pbc_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]





#### volcano PBC without MASLD vs. healthy liver

# Identify proteins that meet the filtering criteria

selected_proteins_t_pbc_no_masld <- results.t.pbc.no.masld %>%
  filter(adj.P.Val < 0.05)

# As volcano plot cannot print all the names and becomes impossible to view, I filter only those proteins
# that are regulated more than 1 or less than -1
selected_proteins_t_pbc_volcano_no_masld <- results.t.pbc.no.masld %>%
  filter(logFC >= 1 & adj.P.Val < 0.05 | logFC <= -1 & adj.P.Val < 0.05)

volcano_plot_t_pbc_no_masld <- ggplot(results.t.pbc.no.masld) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val),
                 color = ifelse(logFC <= -0.0001 & adj.P.Val < 0.05, 'blue', 
                                ifelse(logFC >= 0.0001 & adj.P.Val < 0.05, 'red', 'grey')))) +
  geom_text(data = selected_proteins_t_pbc_volcano_no_masld,
            aes(x = logFC, 
                y = -log10(adj.P.Val), 
                label = GeneName),
            size = 3, hjust = 0, vjust = 0) +
  scale_color_identity() +
  labs(x = "log2(FC)", y = "-log10(adj P value)") +
  ggtitle("Differential Liver Protein Expression, PBC without MASLD vs. Healthy")+
  annotate("segment", x = -Inf, xend = Inf, y = -log10(0.05), yend = -log10(0.05),
           linetype = "dashed", color = "gray") +  # Add horizontal line
  theme(
    axis.title = element_text(size = 20),  # Increase axis title size
    axis.text = element_text(size = 18),   # Increase axis text (tick labels) size
    plot.title = element_text(size = 20)   # Increase plot title size
  )

volcano_plot_t_pbc_no_masld

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/No MASLD/volcano_liver_pbc_healthy_fdr0.05_no.masld.jpeg",
#         width = 10, height = 10, plot = volcano_plot_t_pbc_no_masld, device = "jpeg")

# Number of upregulated proteins in PBC:

selected_upregulated_proteins_t_pbc_no_masld <- results.t.pbc.no.masld %>%
  filter((logFC >= 0.0001 & adj.P.Val < 0.05))
nrow(selected_upregulated_proteins_t_pbc_no_masld) # 1039

# Number of downregulated proteins in PBC:
selected_downregulated_proteins_t_pbc_no_masld <- results.t.pbc.no.masld %>%
  filter(logFC <= -0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_t_pbc_no_masld) # 1272

# Sort the upregulated proteins by adjusted p-value
sorted_upregulated_proteins_t_pbc_no_masld <- 
  arrange(subset(results.t.pbc.no.masld, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of upregulated proteins with gene names
print(sorted_upregulated_proteins_t_pbc_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
up_liver_pbc_vs_healthy_no_masld <- 
  sorted_upregulated_proteins_t_pbc_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]

# Sort the downregulated proteins by adjusted p-value
sorted_downregulated_proteins_t_pbc_no_masld <- 
  arrange(subset(results.t.pbc.no.masld, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of downregulated proteins with gene names
print(sorted_downregulated_proteins_t_pbc_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
down_liver_pbc_vs_healthy_no_masld <- 
  sorted_downregulated_proteins_t_pbc_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]


########## PBC vs Healthy and MASLD

# Check the new group counts
table(t$group_1)

### DE analysis of PBC compared with Healthy/MASLD

t_pbc_h.m <- subset(t, group_1 %in% c("PBC","Control"))

table(t_pbc_h.m$group_1)

design_t_pbc_h.m <- model.matrix(~ 0 + group_1, t_pbc_h.m)
design_t_pbc_h.m[1:20,]

dim(t_pbc_h.m)
t_pbc_h.m[1:3,1:5]
t_pbc_h.m[1:3,7632:7636]

t_pbc_h.m <- t_pbc_h.m[,3:7634] # create only numeric values
dim(t_pbc_h.m)

fit_t_pbc_h.m <- lmFit(t(t_pbc_h.m), design_t_pbc_h.m)
head(coef(fit_t_pbc_h.m))

contr_t_pbc_h.m <- makeContrasts(group_1PBC-group_1Control, levels=design_t_pbc_h.m)
contr_t_pbc_h.m

tmp_t_pbc_h.m <- contrasts.fit(fit_t_pbc_h.m, contr_t_pbc_h.m)
tmp_t_pbc_h.m <- eBayes(tmp_t_pbc_h.m)
tmp_t_pbc_h.m

### Diagnostic plots:
# Scatterplot of residual-variances vs average log-expression
# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/All/plotSA_pbc_vs_healthy_and_masld_liver.jpeg")
#plotSA(tmp_t_pbc_h.m, main = "Residual variances vs. average log-expression, PBC vs. Healthy and MASLD, Liver")
# dev.off()

# Mean Difference plot. Log-intensity ratios (differences, y) versus log-intensity averages (means, x):
# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/All/MD-plot_pbc_vs_healthy_and_masld_liver.jpeg")
#plotMD(tmp_t_pbc_h.m, main = "Mean-difference plot, PBC vs. Healthy and MASLD, Liver") 
# dev.off()

# Extract results
top.table.t.pbc.h.m <- topTable(tmp_t_pbc_h.m, coef = 1, sort.by = "P", n = Inf) 
results.t.pbc.h.m <-  as.data.frame(top.table.t.pbc.h.m)
results.t.pbc.h.m$protein <- row.names(results.t.pbc.h.m)
names(results.t.pbc.h.m)
results.t.pbc.h.m[1:20,1:7]
results.t.pbc.h.m$GeneName <- sapply(results.t.pbc.h.m$protein, 
                                 extract_gene_name)
results.t.pbc.h.m[1:20,5:8]
results.t.pbc.h.m$UniProtID <- sapply(results.t.pbc.h.m$protein, 
                                  extract_uniprot_id)
results.t.pbc.h.m$GeneName.UniprotID <- results.t.pbc.h.m$protein
results.t.pbc.h.m[1:20,6:10]
results.t.pbc.h.m$log2FC <- results.t.pbc.h.m$logFC
results.t.pbc.h.m[1:3,c(1,8:11)]

#### Create volcano plot

# Identify proteins that meet the filtering criteria

selected_proteins_t_pbc_h.m <- results.t.pbc.h.m %>%
  filter(adj.P.Val < 0.05)

selected_proteins_t_pbc_h.m[1:20,1:5]

# As volcano plot cannot print all the names and becomes impossible to view, I filter only those proteins
# that are regulated more than 1 or less than -1
selected_proteins_t_pbc_h.m_volcano <- results.t.pbc.h.m %>%
  filter(logFC >= 1 & adj.P.Val < 0.05 | logFC <= -1 & adj.P.Val < 0.05)

volcano_plot_t_pbc_h.m <- ggplot(results.t.pbc.h.m) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), 
                 color = ifelse(logFC <= -0.00001  & adj.P.Val < 0.05, 'blue', 
                                ifelse(logFC >= 0.00001 & adj.P.Val < 0.05, 'red', 'grey')))) +
  geom_text(data = selected_proteins_t_pbc_h.m_volcano,
            aes(x = logFC, y = -log10(adj.P.Val), 
                label = GeneName),
            size = 3, hjust = 0, vjust = 0) +
  scale_color_identity() +
  labs(x = "log2 fold change", y = "-log10 adj p-value") +
  ggtitle("Differential Liver Protein Expression, PBC vs. Healthy and MASLD")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray")  # Add horizontal line

volcano_plot_t_pbc_h.m

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/All/volcano_liver_pbc_healthy_and_masld_fdr0.05.jpeg",
#         width = 10, height = 10, plot = volcano_plot_t_pbc_h.m, device = "jpeg")

# Number of upregulated proteins in PBC:

selected_upregulated_proteins_t_pbc_h.m <- results.t.pbc.h.m %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_t_pbc_h.m) # 899

# Number of downregulated proteins in PBC:
selected_downregulated_proteins_t_pbc_h.m <- results.t.pbc.h.m %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_t_pbc_h.m) # 411

# Sort the upregulated proteins by adjusted p-value
sorted_upregulated_proteins_t_pbc_h.m <- 
  arrange(subset(results.t.pbc.h.m, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of upregulated proteins with gene names
print(sorted_upregulated_proteins_t_pbc_h.m[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
up_liver_pbc_vs_h.m <- 
  sorted_upregulated_proteins_t_pbc_h.m[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]

# Sort the downregulated proteins by adjusted p-value
sorted_downregulated_proteins_t_pbc_h.m <- 
  arrange(subset(results.t.pbc.h.m, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of downregulated proteins with gene names
print(sorted_downregulated_proteins_t_pbc_h.m[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
down_liver_pbc_vs_h.m <- 
  sorted_downregulated_proteins_t_pbc_h.m[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]



################### Plasma PBC vs Healthy - single covariate sensitivity analyses ###################

p_match_single <- pca_p

p_match_single <- p_match_single %>%
  select(-any_of(c("cirrose", "con_masld", "group_1", "pbc_antibody", "group_urso")))

p_match_single <- subset(p_match_single, group %in% c("PBC","Healthy"))
levels(p_match_single$group)
table(p_match_single$group)
p_match_single$group <- droplevels(p_match_single$group)

dim(p_match_single)
p_match_single[1:3,1:6]
p_match_single[1:3,583:587]

# Here metadata columns are:
# 1 subject
# 2 age
# 3 sex
# 4 bmi
# 5 group
# proteins start at column 6

p_match_single_numeric <- p_match_single[,6:587]
dim(p_match_single_numeric)

################### plasma age only ###################

design_p_match_age <- model.matrix(~ 0 + group + age, p_match_single)
design_p_match_age[1:20,]

fit_p_match_age <- lmFit(t(p_match_single_numeric), design_p_match_age)
contr_p_match_age <- makeContrasts(groupPBC-groupHealthy, levels=design_p_match_age)

tmp_p_match_age <- contrasts.fit(fit_p_match_age, contr_p_match_age)
tmp_p_match_age <- eBayes(tmp_p_match_age)

results.p.match.age <- as.data.frame(topTable(tmp_p_match_age, coef = 1, sort.by = "P", n = Inf))
results.p.match.age$protein <- row.names(results.p.match.age)

results.p.match.age$GeneName <- sapply(results.p.match.age$protein, extract_gene_name)
results.p.match.age$UniProtID <- sapply(results.p.match.age$protein, extract_uniprot_id)
results.p.match.age$GeneName.UniprotID <- results.p.match.age$protein
results.p.match.age$log2FC <- results.p.match.age$logFC

selected_upregulated_proteins_p_match_age <- results.p.match.age %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_p_match_age)

selected_downregulated_proteins_p_match_age <- results.p.match.age %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_p_match_age)

sorted_upregulated_proteins_p_match_age <- 
  arrange(subset(results.p.match.age, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

up_plasma_pbc_vs_healthy_age_match <- 
  sorted_upregulated_proteins_p_match_age[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]

sorted_downregulated_proteins_p_match_age <- 
  arrange(subset(results.p.match.age, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

down_plasma_pbc_vs_healthy_age_match <- 
  sorted_downregulated_proteins_p_match_age[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]


################### plama sex only ###################

design_p_match_sex <- model.matrix(~ 0 + group + sex, p_match_single)
design_p_match_sex[1:20,]

fit_p_match_sex <- lmFit(t(p_match_single_numeric), design_p_match_sex)
contr_p_match_sex <- makeContrasts(groupPBC-groupHealthy, levels=design_p_match_sex)

tmp_p_match_sex <- contrasts.fit(fit_p_match_sex, contr_p_match_sex)
tmp_p_match_sex <- eBayes(tmp_p_match_sex)

results.p.match.sex <- as.data.frame(topTable(tmp_p_match_sex, coef = 1, sort.by = "P", n = Inf))
results.p.match.sex$protein <- row.names(results.p.match.sex)

results.p.match.sex$GeneName <- sapply(results.p.match.sex$protein, extract_gene_name)
results.p.match.sex$UniProtID <- sapply(results.p.match.sex$protein, extract_uniprot_id)
results.p.match.sex$GeneName.UniprotID <- results.p.match.sex$protein
results.p.match.sex$log2FC <- results.p.match.sex$logFC

selected_upregulated_proteins_p_match_sex <- results.p.match.sex %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_p_match_sex)

selected_downregulated_proteins_p_match_sex <- results.p.match.sex %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_p_match_sex)

sorted_upregulated_proteins_p_match_sex <- 
  arrange(subset(results.p.match.sex, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

up_plasma_pbc_vs_healthy_sex_match <- 
  sorted_upregulated_proteins_p_match_sex[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]

sorted_downregulated_proteins_p_match_sex <- 
  arrange(subset(results.p.match.sex, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

down_plasma_pbc_vs_healthy_sex_match <- 
  sorted_downregulated_proteins_p_match_sex[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]


################### plama BMI only ###################

design_p_match_bmi <- model.matrix(~ 0 + group + bmi, p_match_single)
design_p_match_bmi[1:20,]

fit_p_match_bmi <- lmFit(t(p_match_single_numeric), design_p_match_bmi)
contr_p_match_bmi <- makeContrasts(groupPBC-groupHealthy, levels=design_p_match_bmi)

tmp_p_match_bmi <- contrasts.fit(fit_p_match_bmi, contr_p_match_bmi)
tmp_p_match_bmi <- eBayes(tmp_p_match_bmi)

results.p.match.bmi <- as.data.frame(topTable(tmp_p_match_bmi, coef = 1, sort.by = "P", n = Inf))
results.p.match.bmi$protein <- row.names(results.p.match.bmi)

results.p.match.bmi$GeneName <- sapply(results.p.match.bmi$protein, extract_gene_name)
results.p.match.bmi$UniProtID <- sapply(results.p.match.bmi$protein, extract_uniprot_id)
results.p.match.bmi$GeneName.UniprotID <- results.p.match.bmi$protein
results.p.match.bmi$log2FC <- results.p.match.bmi$logFC

selected_upregulated_proteins_p_match_bmi <- results.p.match.bmi %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_p_match_bmi)

selected_downregulated_proteins_p_match_bmi <- results.p.match.bmi %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_p_match_bmi)

sorted_upregulated_proteins_p_match_bmi <- 
  arrange(subset(results.p.match.bmi, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

up_plasma_pbc_vs_healthy_bmi_match <- 
  sorted_upregulated_proteins_p_match_bmi[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]

sorted_downregulated_proteins_p_match_bmi <- 
  arrange(subset(results.p.match.bmi, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

down_plasma_pbc_vs_healthy_bmi_match <- 
  sorted_downregulated_proteins_p_match_bmi[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]


################### plasma cirrhosis only ###################

p_match_cirrhosis <- pca_p

p_match_cirrhosis <- p_match_cirrhosis %>%
  select(-any_of(c("con_masld", "group_1", "pbc_antibody", "group_urso")))

p_match_cirrhosis <- subset(p_match_cirrhosis, group %in% c("PBC","Healthy"))
p_match_cirrhosis$group <- droplevels(p_match_cirrhosis$group)

# Metadata columns:
# 1 subject
# 2 age
# 3 sex
# 4 bmi
# 5 cirrose
# 6 group
# proteins start at column 7

p_match_cirrhosis_numeric <- p_match_cirrhosis[,7:588]

design_p_match_cirrhosis <- model.matrix(~ 0 + group + cirrose, p_match_cirrhosis)
design_p_match_cirrhosis[1:20,]

fit_p_match_cirrhosis <- lmFit(t(p_match_cirrhosis_numeric), design_p_match_cirrhosis)

contr_p_match_cirrhosis <- makeContrasts(groupPBC-groupHealthy,
                                         levels = design_p_match_cirrhosis)

tmp_p_match_cirrhosis <- contrasts.fit(fit_p_match_cirrhosis,
                                       contr_p_match_cirrhosis)
tmp_p_match_cirrhosis <- eBayes(tmp_p_match_cirrhosis)

results.p.match.cirrhosis <- as.data.frame(
  topTable(tmp_p_match_cirrhosis, coef = 1, sort.by = "P", n = Inf)
)

results.p.match.cirrhosis$protein <- row.names(results.p.match.cirrhosis)

results.p.match.cirrhosis$GeneName <- sapply(results.p.match.cirrhosis$protein, extract_gene_name)
results.p.match.cirrhosis$UniProtID <- sapply(results.p.match.cirrhosis$protein, extract_uniprot_id)
results.p.match.cirrhosis$GeneName.UniprotID <- results.p.match.cirrhosis$protein
results.p.match.cirrhosis$log2FC <- results.p.match.cirrhosis$logFC

selected_upregulated_proteins_p_match_cirrhosis <- results.p.match.cirrhosis %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_p_match_cirrhosis)

selected_downregulated_proteins_p_match_cirrhosis <- results.p.match.cirrhosis %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_p_match_cirrhosis)


################### Age + sex matching PBC vs Healthy #########################################################


### Plasma

dim(pca_p)
pca_p[1:5, 1:10] 
pca_p[1:5, 589:591]

p_match <- pca_p
str(p_match[1:9])
table(p_match$group)

# Exclude metadata not used in model

p_match <- p_match %>%
  select(-any_of(c("bmi", "cirrose", "con_masld", "group_1", "pbc_antibody", "group_urso")))

# Dataset without cirrhosis

pca_p$subject[pca_p$cirrose=="Yes"]

p_match_no_c <- subset(p_match, ! subject %in% c("FALL_6", "FALL_75"))

dim(p_match)
dim(p_match_no_c)

p_match <- subset(p_match, group %in% c("PBC","Healthy"))
levels(p_match$group)
table(p_match$group)
p_match$group <- droplevels(p_match$group)
p_match[1:3,1:5]
dim(p_match)

design_p_match <- model.matrix(~ 0 + group + age + sex, p_match)
design_p_match[1:20,]

dim(p_match)
p_match[1:3,1:5]
p_match[1:3,582:586]
p_match <- p_match[,5:586] # create only numeric protein values
dim(p_match)

fit_p_match <- lmFit(t(p_match), design_p_match)
head(coef(fit_p_match))

contr_p_match <- makeContrasts(groupPBC-groupHealthy, levels=design_p_match)
contr_p_match

tmp_p_match <- contrasts.fit(fit_p_match, contr_p_match)
tmp_p_match <- eBayes(tmp_p_match)

top.table.p.match <- topTable(tmp_p_match, coef = 1, sort.by = "P", n = Inf) 
results.p.match <- as.data.frame(top.table.p.match)
results.p.match$protein <- row.names(results.p.match)

results.p.match$GeneName <- sapply(results.p.match$protein, extract_gene_name)
results.p.match$UniProtID <- sapply(results.p.match$protein, extract_uniprot_id)
results.p.match$GeneName.UniprotID <- results.p.match$protein
results.p.match$log2FC <- results.p.match$logFC

selected_upregulated_proteins_p_match <- results.p.match %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_p_match)

selected_downregulated_proteins_p_match <- results.p.match %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_p_match)

sorted_upregulated_proteins_p_match <- 
  arrange(subset(results.p.match, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

up_plasma_pbc_vs_healthy_match <- 
  sorted_upregulated_proteins_p_match[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]

sorted_downregulated_proteins_p_match <- 
  arrange(subset(results.p.match, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

down_plasma_pbc_vs_healthy_match <- 
  sorted_downregulated_proteins_p_match[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]


################### plasma cirrhosis only ###################

p_match_cirrhosis <- pca_p

p_match_cirrhosis <- p_match_cirrhosis %>%
  select(-any_of(c("con_masld", "group_1", "pbc_antibody", "group_urso")))

p_match_cirrhosis <- subset(p_match_cirrhosis, group %in% c("PBC","Healthy"))
levels(p_match_cirrhosis$group)
table(p_match_cirrhosis$group)
p_match_cirrhosis$group <- droplevels(p_match_cirrhosis$group)

dim(p_match_cirrhosis)
p_match_cirrhosis[1:3,1:7]
p_match_cirrhosis[1:3,584:588]

# Here metadata columns are:
# 1 subject
# 2 age
# 3 sex
# 4 bmi
# 5 cirrose
# 6 group
# proteins start at column 7

p_match_cirrhosis_numeric <- p_match_cirrhosis[,7:588]
dim(p_match_cirrhosis_numeric)

design_p_match_cirrhosis <- model.matrix(~ 0 + group + cirrose, p_match_cirrhosis)
design_p_match_cirrhosis[1:20,]

fit_p_match_cirrhosis <- lmFit(t(p_match_cirrhosis_numeric), design_p_match_cirrhosis)
contr_p_match_cirrhosis <- makeContrasts(groupPBC-groupHealthy, levels=design_p_match_cirrhosis)

tmp_p_match_cirrhosis <- contrasts.fit(fit_p_match_cirrhosis, contr_p_match_cirrhosis)
tmp_p_match_cirrhosis <- eBayes(tmp_p_match_cirrhosis)

results.p.match.cirrhosis <- as.data.frame(topTable(tmp_p_match_cirrhosis, coef = 1, sort.by = "P", n = Inf))
results.p.match.cirrhosis$protein <- row.names(results.p.match.cirrhosis)

results.p.match.cirrhosis$GeneName <- sapply(results.p.match.cirrhosis$protein, extract_gene_name)
results.p.match.cirrhosis$UniProtID <- sapply(results.p.match.cirrhosis$protein, extract_uniprot_id)
results.p.match.cirrhosis$GeneName.UniprotID <- results.p.match.cirrhosis$protein
results.p.match.cirrhosis$log2FC <- results.p.match.cirrhosis$logFC

selected_upregulated_proteins_p_match_cirrhosis <- results.p.match.cirrhosis %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_p_match_cirrhosis)

selected_downregulated_proteins_p_match_cirrhosis <- results.p.match.cirrhosis %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_p_match_cirrhosis)

sorted_upregulated_proteins_p_match_cirrhosis <- 
  arrange(subset(results.p.match.cirrhosis, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

up_plasma_pbc_vs_healthy_cirrhosis_match <- 
  sorted_upregulated_proteins_p_match_cirrhosis[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]

sorted_downregulated_proteins_p_match_cirrhosis <- 
  arrange(subset(results.p.match.cirrhosis, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

down_plasma_pbc_vs_healthy_cirrhosis_match <- 
  sorted_downregulated_proteins_p_match_cirrhosis[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]



### Plasma without cirrhosis

p_match_no_c <- subset(p_match_no_c, group %in% c("PBC","Healthy"))
levels(p_match_no_c$group)
table(p_match_no_c$group)
p_match_no_c$group <- droplevels(p_match_no_c$group)
table(p_match_no_c$group)

design_p_match_no_c <- model.matrix(~ 0 + group + age + sex, p_match_no_c)
design_p_match_no_c[1:20,] 

dim(p_match_no_c)
p_match_no_c[1:3,1:5]
p_match_no_c[1:3,582:586]

p_match_no_c <- p_match_no_c[,5:586] # create only numeric protein values
dim(p_match_no_c)

fit_p_match_no_c <- lmFit(t(p_match_no_c), design_p_match_no_c)
head(coef(fit_p_match_no_c))

contr_p_match_no_c <- makeContrasts(groupPBC-groupHealthy, levels=design_p_match_no_c)
contr_p_match_no_c

tmp_p_match_no_c <- contrasts.fit(fit_p_match_no_c, contr_p_match_no_c)
tmp_p_match_no_c <- eBayes(tmp_p_match_no_c)

top.table.p.match.no.c <- topTable(tmp_p_match_no_c, coef = 1, sort.by = "P", n = Inf)
results.p.match.no.c <- as.data.frame(top.table.p.match.no.c)
results.p.match.no.c$protein <- row.names(results.p.match.no.c)

results.p.match.no.c$GeneName <- sapply(results.p.match.no.c$protein, extract_gene_name)
results.p.match.no.c$UniProtID <- sapply(results.p.match.no.c$protein, extract_uniprot_id)
results.p.match.no.c$GeneName.UniprotID <- results.p.match.no.c$protein
results.p.match.no.c$log2FC <- results.p.match.no.c$logFC

selected_upregulated_proteins_p_match_no_c <- results.p.match.no.c %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_p_match_no_c)

selected_downregulated_proteins_p_match_no_c <- results.p.match.no.c %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_p_match_no_c)

sorted_upregulated_proteins_p_match_no_c <- 
  arrange(subset(results.p.match.no.c, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

up_plasma_pbc_vs_healthy_no_c_match <- 
  sorted_upregulated_proteins_p_match_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]

sorted_downregulated_proteins_p_match_no_c <- 
  arrange(subset(results.p.match.no.c, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

down_plasma_pbc_vs_healthy_no_c_match <- 
  sorted_downregulated_proteins_p_match_no_c[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]



################### Summary plasma ###################

data.frame(
  Model = c(
    "Unadjusted",
    "Age only",
    "Sex only",
    "BMI only",
    "Cirrhosis only",
    "Age + sex",
    "Age + sex (no cirrhosis)"
  ),
  Upregulated = c(
    nrow(selected_upregulated_proteins_p_pbc),
    nrow(selected_upregulated_proteins_p_match_age),
    nrow(selected_upregulated_proteins_p_match_sex),
    nrow(selected_upregulated_proteins_p_match_bmi),
    nrow(selected_upregulated_proteins_p_match_cirrhosis),
    nrow(selected_upregulated_proteins_p_match),
    nrow(selected_upregulated_proteins_p_match_no_c)
  ),
  Downregulated = c(
    nrow(selected_downregulated_proteins_p_pbc),
    nrow(selected_downregulated_proteins_p_match_age),
    nrow(selected_downregulated_proteins_p_match_sex),
    nrow(selected_downregulated_proteins_p_match_bmi),
    nrow(selected_downregulated_proteins_p_match_cirrhosis),
    nrow(selected_downregulated_proteins_p_match),
    nrow(selected_downregulated_proteins_p_match_no_c)
  )
)




################### Liver PBC vs Healthy - single covariate sensitivity analyses ###################

t_match_single <- pca_t

t_match_single <- t_match_single %>%
  select(-any_of(c("cirrose", "con_masld", "group_1", "pbc_antibody", "group_urso")))

t_match_single <- subset(t_match_single, group %in% c("PBC","Healthy"))
levels(t_match_single$group)
table(t_match_single$group)
t_match_single$group <- droplevels(t_match_single$group)

dim(t_match_single)
t_match_single[1:3,1:6]
t_match_single[1:3,7633:7637]

# Here metadata columns are:
# 1 subject
# 2 age
# 3 sex
# 4 bmi
# 5 group
# proteins start at column 6

t_match_single_numeric <- t_match_single[,6:7637]
dim(t_match_single_numeric)

################### liver age only ###################

design_t_match_age <- model.matrix(~ 0 + group + age, t_match_single)
design_t_match_age[1:20,]

fit_t_match_age <- lmFit(t(t_match_single_numeric), design_t_match_age)
contr_t_match_age <- makeContrasts(groupPBC-groupHealthy, levels=design_t_match_age)

tmp_t_match_age <- contrasts.fit(fit_t_match_age, contr_t_match_age)
tmp_t_match_age <- eBayes(tmp_t_match_age)

results.t.match.age <- as.data.frame(topTable(tmp_t_match_age, coef = 1, sort.by = "P", n = Inf))
results.t.match.age$protein <- row.names(results.t.match.age)

results.t.match.age$GeneName <- sapply(results.t.match.age$protein, extract_gene_name)
results.t.match.age$UniProtID <- sapply(results.t.match.age$protein, extract_uniprot_id)
results.t.match.age$GeneName.UniprotID <- results.t.match.age$protein
results.t.match.age$log2FC <- results.t.match.age$logFC

selected_upregulated_proteins_t_match_age <- results.t.match.age %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_t_match_age)

selected_downregulated_proteins_t_match_age <- results.t.match.age %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_t_match_age)

sorted_upregulated_proteins_t_match_age <- 
  arrange(subset(results.t.match.age, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

up_liver_pbc_vs_healthy_age_match <- 
  sorted_upregulated_proteins_t_match_age[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]

sorted_downregulated_proteins_t_match_age <- 
  arrange(subset(results.t.match.age, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

down_liver_pbc_vs_healthy_age_match <- 
  sorted_downregulated_proteins_t_match_age[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]



################### liver sex only ###################

design_t_match_sex <- model.matrix(~ 0 + group + sex, t_match_single)
design_t_match_sex[1:20,]

fit_t_match_sex <- lmFit(t(t_match_single_numeric), design_t_match_sex)
contr_t_match_sex <- makeContrasts(groupPBC-groupHealthy, levels=design_t_match_sex)

tmp_t_match_sex <- contrasts.fit(fit_t_match_sex, contr_t_match_sex)
tmp_t_match_sex <- eBayes(tmp_t_match_sex)

results.t.match.sex <- as.data.frame(topTable(tmp_t_match_sex, coef = 1, sort.by = "P", n = Inf))
results.t.match.sex$protein <- row.names(results.t.match.sex)

results.t.match.sex$GeneName <- sapply(results.t.match.sex$protein, extract_gene_name)
results.t.match.sex$UniProtID <- sapply(results.t.match.sex$protein, extract_uniprot_id)
results.t.match.sex$GeneName.UniprotID <- results.t.match.sex$protein
results.t.match.sex$log2FC <- results.t.match.sex$logFC


selected_upregulated_proteins_t_match_sex <- results.t.match.sex %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_t_match_sex)

selected_downregulated_proteins_t_match_sex <- results.t.match.sex %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_t_match_sex)

sorted_upregulated_proteins_t_match_sex <- 
  arrange(subset(results.t.match.sex, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

up_liver_pbc_vs_healthy_sex_match <- 
  sorted_upregulated_proteins_t_match_sex[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]

sorted_downregulated_proteins_t_match_sex <- 
  arrange(subset(results.t.match.sex, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

down_liver_pbc_vs_healthy_sex_match <- 
  sorted_downregulated_proteins_t_match_sex[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]


################### liver BMI only ###################

design_t_match_bmi <- model.matrix(~ 0 + group + bmi, t_match_single)
design_t_match_bmi[1:20,]

fit_t_match_bmi <- lmFit(t(t_match_single_numeric), design_t_match_bmi)
contr_t_match_bmi <- makeContrasts(groupPBC-groupHealthy, levels=design_t_match_bmi)

tmp_t_match_bmi <- contrasts.fit(fit_t_match_bmi, contr_t_match_bmi)
tmp_t_match_bmi <- eBayes(tmp_t_match_bmi)

results.t.match.bmi <- as.data.frame(topTable(tmp_t_match_bmi, coef = 1, sort.by = "P", n = Inf))
results.t.match.bmi$protein <- row.names(results.t.match.bmi)

results.t.match.bmi$GeneName <- sapply(results.t.match.bmi$protein, extract_gene_name)
results.t.match.bmi$UniProtID <- sapply(results.t.match.bmi$protein, extract_uniprot_id)
results.t.match.bmi$GeneName.UniprotID <- results.t.match.bmi$protein
results.t.match.bmi$log2FC <- results.t.match.bmi$logFC

selected_upregulated_proteins_t_match_bmi <- results.t.match.bmi %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_t_match_bmi)

selected_downregulated_proteins_t_match_bmi <- results.t.match.bmi %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_t_match_bmi)

sorted_upregulated_proteins_t_match_bmi <- 
  arrange(subset(results.t.match.bmi, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

up_liver_pbc_vs_healthy_bmi_match <- 
  sorted_upregulated_proteins_t_match_bmi[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]

sorted_downregulated_proteins_t_match_bmi <- 
  arrange(subset(results.t.match.bmi, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

down_liver_pbc_vs_healthy_bmi_match <- 
  sorted_downregulated_proteins_t_match_bmi[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]


################### liver cirrhosis only ###################

t_match_cirrhosis <- pca_t

t_match_cirrhosis <- t_match_cirrhosis %>%
  select(-any_of(c("con_masld", "group_1", "pbc_antibody", "group_urso")))

t_match_cirrhosis <- subset(t_match_cirrhosis, group %in% c("PBC","Healthy"))
levels(t_match_cirrhosis$group)
table(t_match_cirrhosis$group)
t_match_cirrhosis$group <- droplevels(t_match_cirrhosis$group)

dim(t_match_cirrhosis)
t_match_cirrhosis[1:3,1:7]
t_match_cirrhosis[1:3,7634:7638]

# Here metadata columns are:
# 1 subject
# 2 age
# 3 sex
# 4 bmi
# 5 cirrose
# 6 group
# proteins start at column 7

t_match_cirrhosis_numeric <- t_match_cirrhosis[,7:7638]
dim(t_match_cirrhosis_numeric)

design_t_match_cirrhosis <- model.matrix(~ 0 + group + cirrose, t_match_cirrhosis)
design_t_match_cirrhosis[1:20,]

fit_t_match_cirrhosis <- lmFit(t(t_match_cirrhosis_numeric), design_t_match_cirrhosis)
contr_t_match_cirrhosis <- makeContrasts(groupPBC-groupHealthy, levels=design_t_match_cirrhosis)

tmp_t_match_cirrhosis <- contrasts.fit(fit_t_match_cirrhosis, contr_t_match_cirrhosis)
tmp_t_match_cirrhosis <- eBayes(tmp_t_match_cirrhosis)

results.t.match.cirrhosis <- as.data.frame(topTable(tmp_t_match_cirrhosis, coef = 1, sort.by = "P", n = Inf))
results.t.match.cirrhosis$protein <- row.names(results.t.match.cirrhosis)

results.t.match.cirrhosis$GeneName <- sapply(results.t.match.cirrhosis$protein, extract_gene_name)
results.t.match.cirrhosis$UniProtID <- sapply(results.t.match.cirrhosis$protein, extract_uniprot_id)
results.t.match.cirrhosis$GeneName.UniprotID <- results.t.match.cirrhosis$protein
results.t.match.cirrhosis$log2FC <- results.t.match.cirrhosis$logFC

selected_upregulated_proteins_t_match_cirrhosis <- results.t.match.cirrhosis %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_t_match_cirrhosis)

selected_downregulated_proteins_t_match_cirrhosis <- results.t.match.cirrhosis %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_t_match_cirrhosis)

sorted_upregulated_proteins_t_match_cirrhosis <- 
  arrange(subset(results.t.match.cirrhosis, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

up_liver_pbc_vs_healthy_cirrhosis_match <- 
  sorted_upregulated_proteins_t_match_cirrhosis[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]

sorted_downregulated_proteins_t_match_cirrhosis <- 
  arrange(subset(results.t.match.cirrhosis, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

down_liver_pbc_vs_healthy_cirrhosis_match <- 
  sorted_downregulated_proteins_t_match_cirrhosis[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]



################## Liver age + sex matching PBC vs healthy #####

dim(pca_t)
pca_t[1:5, 1:10] 
pca_t[1:5, 7638:7642]

t_match <- pca_t
str(t_match[1:9])
table(t_match$group)

t_match <- t_match %>%
  select(-any_of(c("bmi", "cirrose", "con_masld", "group_1", "pbc_antibody", "group_urso")))

pca_t$subject[pca_t$cirrose=="Yes"]

t_match_no_c <- subset(t_match, ! subject %in% c("FALL_6", "FALL_75"))

dim(t_match)
dim(t_match_no_c)

t_match <- subset(t_match, group %in% c("PBC","Healthy"))
levels(t_match$group)
table(t_match$group)
t_match$group <- droplevels(t_match$group)
t_match[1:3,1:5]
dim(t_match)

design_t_match <- model.matrix(~ 0 + group + age + sex, t_match)
design_t_match[1:20,]

dim(t_match)
t_match[1:3,1:5]
t_match[1:3, (ncol(t_match)-2):ncol(t_match)]
t_match <- t_match[,5:ncol(t_match)]
dim(t_match)

fit_t_match <- lmFit(t(t_match), design_t_match)
head(coef(fit_t_match))

contr_t_match <- makeContrasts(groupPBC-groupHealthy, levels=design_t_match)
contr_t_match

tmp_t_match <- contrasts.fit(fit_t_match, contr_t_match)
tmp_t_match <- eBayes(tmp_t_match)

top.table.t.match <- topTable(tmp_t_match, coef = 1, sort.by = "P", n = Inf) 
results.t.match <- as.data.frame(top.table.t.match)
results.t.match$protein <- row.names(results.t.match)

results.t.match$GeneName <- sapply(results.t.match$protein, extract_gene_name)
results.t.match$UniProtID <- sapply(results.t.match$protein, extract_uniprot_id)
results.t.match$GeneName.UniprotID <- results.t.match$protein
results.t.match$log2FC <- results.t.match$logFC

selected_upregulated_proteins_t_match <- results.t.match %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_t_match)

selected_downregulated_proteins_t_match <- results.t.match %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_t_match)

sorted_upregulated_proteins_t_match <- 
  arrange(subset(results.t.match, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

up_liver_pbc_vs_healthy_match <- 
  sorted_upregulated_proteins_t_match[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]

sorted_downregulated_proteins_t_match <- 
  arrange(subset(results.t.match, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

down_liver_pbc_vs_healthy_match <- 
  sorted_downregulated_proteins_t_match[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]



### Liver without cirrhosis

t_match_no_c <- subset(t_match_no_c, group %in% c("PBC","Healthy"))
levels(t_match_no_c$group)
table(t_match_no_c$group)
t_match_no_c$group <- droplevels(t_match_no_c$group)
table(t_match_no_c$group)

design_t_match_no_c <- model.matrix(~ 0 + group + age + sex, t_match_no_c)
design_t_match_no_c[1:20,] 

dim(t_match_no_c)
t_match_no_c[1:3,1:5]
t_match_no_c[1:3, (ncol(t_match_no_c)-2):ncol(t_match_no_c)]
t_match_no_c <- t_match_no_c[,5:ncol(t_match_no_c)]
dim(t_match_no_c)

fit_t_match_no_c <- lmFit(t(t_match_no_c), design_t_match_no_c)
head(coef(fit_t_match_no_c))

contr_t_match_no_c <- makeContrasts(groupPBC-groupHealthy, levels=design_t_match_no_c)
contr_t_match_no_c

tmp_t_match_no_c <- contrasts.fit(fit_t_match_no_c, contr_t_match_no_c)
tmp_t_match_no_c <- eBayes(tmp_t_match_no_c)

top.table.t.match.no.c <- topTable(tmp_t_match_no_c, coef = 1, sort.by = "P", n = Inf)
results.t.match.no.c <- as.data.frame(top.table.t.match.no.c)
results.t.match.no.c$protein <- row.names(results.t.match.no.c)

results.t.match.no.c$GeneName <- sapply(results.t.match.no.c$protein, extract_gene_name)
results.t.match.no.c$UniProtID <- sapply(results.t.match.no.c$protein, extract_uniprot_id)
results.t.match.no.c$GeneName.UniprotID <- results.t.match.no.c$protein
results.t.match.no.c$log2FC <- results.t.match.no.c$logFC

selected_upregulated_proteins_t_match_no_c <- results.t.match.no.c %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_t_match_no_c)

selected_downregulated_proteins_t_match_no_c <- results.t.match.no.c %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_t_match_no_c)

sorted_upregulated_proteins_t_match_no_c <- 
  arrange(subset(results.t.match.no.c, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

up_liver_pbc_vs_healthy_no_c_match <- 
  sorted_upregulated_proteins_t_match_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC", "adj.P.Val")]

sorted_downregulated_proteins_t_match_no_c <- 
  arrange(subset(results.t.match.no.c, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

down_liver_pbc_vs_healthy_no_c_match <- 
  sorted_downregulated_proteins_t_match_no_c[, c("GeneName.UniprotID","GeneName","UniProtID", "log2FC", "adj.P.Val")]



################### Summary liver ###################

data.frame(
  Model = c(
    "Unadjusted",
    "Age only",
    "Sex only",
    "BMI only",
    "Cirrhosis only",
    "Age + sex",
    "Age + sex (no cirrhosis)"
  ),
  Upregulated = c(
    nrow(selected_upregulated_proteins_t_pbc),
    nrow(selected_upregulated_proteins_t_match_age),
    nrow(selected_upregulated_proteins_t_match_sex),
    nrow(selected_upregulated_proteins_t_match_bmi),
    nrow(selected_upregulated_proteins_t_match_cirrhosis),
    nrow(selected_upregulated_proteins_t_match),
    nrow(selected_upregulated_proteins_t_match_no_c)
  ),
  Downregulated = c(
    nrow(selected_downregulated_proteins_t_pbc),
    nrow(selected_downregulated_proteins_t_match_age),
    nrow(selected_downregulated_proteins_t_match_sex),
    nrow(selected_downregulated_proteins_t_match_bmi),
    nrow(selected_downregulated_proteins_t_match_cirrhosis),
    nrow(selected_downregulated_proteins_t_match),
    nrow(selected_downregulated_proteins_t_match_no_c)
  )
)



# ######################## Save in an excel file ##########
# write.xlsx(list(
#   "plasma_age_up" = up_plasma_pbc_vs_healthy_age_match, "plasma_age_down" = down_plasma_pbc_vs_healthy_age_match,
#   "plasma_sex_up" = up_plasma_pbc_vs_healthy_sex_match, "plasma_sex_down" = down_plasma_pbc_vs_healthy_sex_match,
#   "plasma_bmi_up" = up_plasma_pbc_vs_healthy_bmi_match, "plasma_bmi_down" = down_plasma_pbc_vs_healthy_bmi_match,
#   "plasma_cirrhosis_up" = up_plasma_pbc_vs_healthy_cirrhosis_match, "plasma_cirrhosis_down" = down_plasma_pbc_vs_healthy_cirrhosis_match,
#   "plasma_age_sex_up" = up_plasma_pbc_vs_healthy_match, "plasma_age_sex_down" = down_plasma_pbc_vs_healthy_match,
#   "plasma_age_sex_nc_up" = up_plasma_pbc_vs_healthy_no_c_match, "plasma_age_sex_nc_down" = down_plasma_pbc_vs_healthy_no_c_match,
#   "liver_age_up" = up_liver_pbc_vs_healthy_age_match, "liver_age_down" = down_liver_pbc_vs_healthy_age_match,
#   "liver_sex_up" = up_liver_pbc_vs_healthy_sex_match, "liver_sex_down" = down_liver_pbc_vs_healthy_sex_match,
#   "liver_bmi_up" = up_liver_pbc_vs_healthy_bmi_match, "liver_bmi_down" = down_liver_pbc_vs_healthy_bmi_match,
#   "liver_cirrhosis_up" = up_liver_pbc_vs_healthy_cirrhosis_match, "liver_cirrhosis_down" = down_liver_pbc_vs_healthy_cirrhosis_match,
#   "liver_age_sex_up" = up_liver_pbc_vs_healthy_match, "liver_age_sex_down" = down_liver_pbc_vs_healthy_match,
#   "liver_age_sex_nc_up" = up_liver_pbc_vs_healthy_no_c_match, "liver_age_sex_nc_down" = down_liver_pbc_vs_healthy_no_c_match),
#   "../output/FDR0.05_all/PBC_vs_Healthy/All/sensitivity_adjusted_DE_proteins_pbc_vs_healthy.xlsx")


######################## Save in an excel file age + sex matching #####################################

# write.xlsx(list("down_plasma_pbc_healthy" = down_plasma_pbc_vs_healthy_match,
#                 "up_plasma_pbc_healthy" = up_plasma_pbc_vs_healthy_match,
#                 "down_liver_pbc_healthy" = down_liver_pbc_vs_healthy_match,
#                 "up_liver_pbc_healthy" = up_liver_pbc_vs_healthy_match),
#            "../output/FDR0.05_all/PBC_vs_Healthy/All/age_sex_corr_DE_proteins_fdr0.05_pbc_vs_healthy.xlsx")

# write.xlsx(list("down_plasma_pbc_h_no_c" = down_plasma_pbc_vs_healthy_no_c_match,
#                 "up_plasma_pbc_h_no_c" = up_plasma_pbc_vs_healthy_no_c_match,
#                 "down_liver_pbc_h_no_c" = down_liver_pbc_vs_healthy_no_c_match,
#                 "up_liver_pbc_h_no_c" = up_liver_pbc_vs_healthy_no_c_match),
#            "../output/FDR0.05_all/PBC_vs_Healthy/No cirrhosis/age_sex_corr_DE_proteins_fdr0.05_pbc_vs_healthy_no.cirrhosis.xlsx")


########## PBC vs Healthy and MASLD - without outliers

# Check the new group counts
table(t_no_outliers$group_1)

### DE analysis of PBC compared with Healthy/MASLD

t_pbc_h.m_no_outliers <- subset(
  t_no_outliers,
  group_1 %in% c("PBC","Control")
)

table(t_pbc_h.m_no_outliers$group_1)


design_t_pbc_h.m_no_outliers <- model.matrix(
  ~ 0 + group_1,
  t_pbc_h.m_no_outliers
)

design_t_pbc_h.m_no_outliers[1:20,]


dim(t_pbc_h.m_no_outliers)

t_pbc_h.m_no_outliers[1:3,1:5]

t_pbc_h.m_no_outliers[1:3,7632:7636]


t_pbc_h.m_no_outliers <- t_pbc_h.m_no_outliers[,3:7634] # create only numeric values

dim(t_pbc_h.m_no_outliers)


fit_t_pbc_h.m_no_outliers <- lmFit(
  t(t_pbc_h.m_no_outliers),
  design_t_pbc_h.m_no_outliers
)

head(coef(fit_t_pbc_h.m_no_outliers))


contr_t_pbc_h.m_no_outliers <- makeContrasts(
  group_1PBC-group_1Control,
  levels=design_t_pbc_h.m_no_outliers
)

contr_t_pbc_h.m_no_outliers


tmp_t_pbc_h.m_no_outliers <- contrasts.fit(
  fit_t_pbc_h.m_no_outliers,
  contr_t_pbc_h.m_no_outliers
)

tmp_t_pbc_h.m_no_outliers <- eBayes(tmp_t_pbc_h.m_no_outliers)

tmp_t_pbc_h.m_no_outliers


### Diagnostic plots:

# Scatterplot of residual-variances vs average log-expression
# Save plot

# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/Without outliers/plotSA_pbc_vs_healthy_and_masld_without_outliers_liver.jpeg")

#plotSA(
#  tmp_t_pbc_h.m_no_outliers,
#  main = "Residual variances vs. average log-expression, PBC vs. Healthy and MASLD without outliers, Liver"
#)

#dev.off()


# Mean Difference plot. Log-intensity ratios (differences, y) versus log-intensity averages (means, x):
# Save plot

# jpeg(file = "../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/Without outliers/MD-plot_pbc_vs_healthy_and_masld_without_outliers_liver.jpeg")

#plotMD(
#  tmp_t_pbc_h.m_no_outliers,
#  main = "Mean-difference plot, PBC vs. Healthy and MASLD without outliers, Liver"
#)

#dev.off()


# Extract results

top.table.t.pbc.h.m.no_outliers <- topTable(
  tmp_t_pbc_h.m_no_outliers,
  coef = 1,
  sort.by = "P",
  n = Inf
)

results.t.pbc.h.m.no_outliers <- as.data.frame(top.table.t.pbc.h.m.no_outliers)

results.t.pbc.h.m.no_outliers$protein <- row.names(results.t.pbc.h.m.no_outliers)

names(results.t.pbc.h.m.no_outliers)

results.t.pbc.h.m.no_outliers[1:20,1:7]


results.t.pbc.h.m.no_outliers$GeneName <- sapply(
  results.t.pbc.h.m.no_outliers$protein,
  extract_gene_name
)

results.t.pbc.h.m.no_outliers[1:20,5:8]


results.t.pbc.h.m.no_outliers$UniProtID <- sapply(
  results.t.pbc.h.m.no_outliers$protein,
  extract_uniprot_id
)

results.t.pbc.h.m.no_outliers$GeneName.UniprotID <- 
  results.t.pbc.h.m.no_outliers$protein

results.t.pbc.h.m.no_outliers[1:20,6:10]


results.t.pbc.h.m.no_outliers$log2FC <- 
  results.t.pbc.h.m.no_outliers$logFC

results.t.pbc.h.m.no_outliers[1:3,c(1,8:11)]


#### Create volcano plot

# Identify proteins that meet the filtering criteria

selected_proteins_t_pbc_h.m_no_outliers <- 
  results.t.pbc.h.m.no_outliers %>%
  filter(adj.P.Val < 0.05)

selected_proteins_t_pbc_h.m_no_outliers[1:20,1:5]


# As volcano plot cannot print all the names and becomes impossible to view,
# I filter only those proteins that are regulated more than 1 or less than -1

selected_proteins_t_pbc_h.m_volcano_no_outliers <- 
  results.t.pbc.h.m.no_outliers %>%
  filter(logFC >= 1 & adj.P.Val < 0.05 |
           logFC <= -1 & adj.P.Val < 0.05)


volcano_plot_t_pbc_h.m_no_outliers <- 
  ggplot(results.t.pbc.h.m.no_outliers) +
  geom_point(aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = ifelse(
      logFC <= -0.00001 & adj.P.Val < 0.05,
      'blue',
      ifelse(
        logFC >= 0.00001 & adj.P.Val < 0.05,
        'red',
        'grey'
      )
    )
  )) +
  geom_text(
    data = selected_proteins_t_pbc_h.m_volcano_no_outliers,
    aes(
      x = logFC,
      y = -log10(adj.P.Val),
      label = GeneName
    ),
    size = 3,
    hjust = 0,
    vjust = 0
  ) +
  scale_color_identity() +
  labs(x = "log2 fold change", y = "-log10 adj p-value") +
  ggtitle("Differential Liver Protein Expression, PBC vs. Healthy and MASLD without outliers") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "gray"
  )

volcano_plot_t_pbc_h.m_no_outliers


# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/Without outliers/volcano_liver_pbc_healthy_and_masld_without_outliers_fdr0.05.jpeg",
#         width = 10, height = 10, plot = volcano_plot_t_pbc_h.m_no_outliers, device = "jpeg")


# Number of upregulated proteins in PBC: #357 

selected_upregulated_proteins_t_pbc_h.m_no_outliers <- 
  results.t.pbc.h.m.no_outliers %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)

nrow(selected_upregulated_proteins_t_pbc_h.m_no_outliers)


# Number of downregulated proteins in PBC: #115

selected_downregulated_proteins_t_pbc_h.m_no_outliers <- 
  results.t.pbc.h.m.no_outliers %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)

nrow(selected_downregulated_proteins_t_pbc_h.m_no_outliers)


# Sort the upregulated proteins by adjusted p-value

sorted_upregulated_proteins_t_pbc_h.m_no_outliers <- 
  arrange(
    subset(
      results.t.pbc.h.m.no_outliers,
      logFC >= 0.0001 & adj.P.Val < 0.05
    ),
    adj.P.Val
  )


# Display the sorted list of upregulated proteins with gene names

print(
  sorted_upregulated_proteins_t_pbc_h.m_no_outliers[, 
                                                    c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]
)


# Create a data frame with the desired output

up_liver_pbc_vs_h.m_no_outliers <- 
  sorted_upregulated_proteins_t_pbc_h.m_no_outliers[, 
                                                    c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]


# Sort the downregulated proteins by adjusted p-value

sorted_downregulated_proteins_t_pbc_h.m_no_outliers <- 
  arrange(
    subset(
      results.t.pbc.h.m.no_outliers,
      logFC <= -0.0001 & adj.P.Val < 0.05
    ),
    adj.P.Val
  )


# Display the sorted list of downregulated proteins with gene names

print(
  sorted_downregulated_proteins_t_pbc_h.m_no_outliers[, 
                                                      c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]
)


# Create a data frame with the desired output

down_liver_pbc_vs_h.m_no_outliers <- 
  sorted_downregulated_proteins_t_pbc_h.m_no_outliers[, 
                                                      c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]




#### PBC vs. Healthy and MASLD, no cirrhosis

table(t_no_c$group_1)

t_pbc_h.m_no_c <- subset(t_no_c, group_1 %in% c("PBC","Control"))

table(t_pbc_h.m_no_c$group_1)

design_t_pbc_h.m_no_c <- model.matrix(~ 0 + group_1, t_pbc_h.m_no_c)
design_t_pbc_h.m_no_c[1:20,]

dim(t_pbc_h.m_no_c)
t_pbc_h.m_no_c[1:3,1:5]
t_pbc_h.m_no_c[1:3,7632:7636]

t_pbc_h.m_no_c <- t_pbc_h.m_no_c[,3:7634] # create only numeric values 
dim(t_pbc_h.m_no_c)

fit_t_pbc_h.m_no_c <- lmFit(t(t_pbc_h.m_no_c), design_t_pbc_h.m_no_c)
head(coef(fit_t_pbc_h.m_no_c))

contr_t_pbc_h.m_no_c <- makeContrasts(group_1PBC-group_1Control, levels=design_t_pbc_h.m_no_c)
contr_t_pbc_h.m_no_c
tmp_t_pbc_h.m_no_c <- contrasts.fit(fit_t_pbc_h.m_no_c, contr_t_pbc_h.m_no_c)
tmp_t_pbc_h.m_no_c <- eBayes(tmp_t_pbc_h.m_no_c)
tmp_t_pbc_h.m_no_c
top.table.t.pbc.h.m.no.c <- topTable(tmp_t_pbc_h.m_no_c,coef = 1, sort.by = "P", n = Inf)
results.t.pbc.h.m.no.c <-  as.data.frame(top.table.t.pbc.h.m.no.c)
results.t.pbc.h.m.no.c$protein <- row.names(results.t.pbc.h.m.no.c)
names(results.t.pbc.h.m.no.c)
results.t.pbc.h.m.no.c[1:20,1:7]
results.t.pbc.h.m.no.c$GeneName <- sapply(results.t.pbc.h.m.no.c$protein, 
                                     extract_gene_name)
results.t.pbc.h.m.no.c[1:20,5:8]
results.t.pbc.h.m.no.c$UniProtID <- sapply(results.t.pbc.h.m.no.c$protein, 
                                      extract_uniprot_id)
results.t.pbc.h.m.no.c$GeneName.UniprotID <- results.t.pbc.h.m.no.c$protein
results.t.pbc.h.m.no.c[1:20,6:10]
results.t.pbc.h.m.no.c$log2FC <- results.t.pbc.h.m.no.c$logFC
results.t.pbc.h.m.no.c[1:3,c(1,8:11)]

# volcano PBC vs Healthy and MASLD no cirrhosis 

selected_proteins_t_pbc_h.m_no_c <- results.t.pbc.h.m.no.c %>%
  filter(adj.P.Val < 0.05)

selected_proteins_t_pbc_h.m_no_c[1:20,1:5]

# As volcano plot cannot print all the names and becomes impossible to view, I filter only those proteins
# that are regulated more than 1 or less than -1
selected_proteins_t_pbc_h.m_no_c_volcano <- results.t.pbc.h.m.no.c %>%
  filter(logFC >= 1 & adj.P.Val < 0.05 | logFC <= -1 & adj.P.Val < 0.05)

volcano_plot_t_pbc_h.m_no_c <- ggplot(results.t.pbc.h.m.no.c) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), 
                 color = ifelse(logFC <= -0.00001  & adj.P.Val < 0.05, 'blue', 
                                ifelse(logFC >= 0.00001 & adj.P.Val < 0.05, 'red', 'grey')))) +
  geom_text(data = selected_proteins_t_pbc_h.m_no_c_volcano,
            aes(x = logFC, y = -log10(adj.P.Val), 
                label = GeneName),
            size = 3, hjust = 0, vjust = 0) +
  scale_color_identity() +
  labs(x = "log2 fold change", y = "-log10 adj p-value") +
  ggtitle("Differential Liver Protein Expression, non-cirrhotic PBC vs. Healthy and MASLD")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray")  # Add horizontal line


volcano_plot_t_pbc_h.m_no_c

# tapply(t$ANKRD22.Q5VYY1,t$group,mean) # "checking that up and downregulated are correct (-1,+1)"

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/No cirrhosis/volcano_liver_pbc_healthy_and_MASLD_no_cirrhosis_fdr0.05.jpeg",
#         width = 10, height = 10, plot = volcano_plot_t_pbc_h.m_no_c, device = "jpeg")

# Number of upregulated proteins in PBC:

selected_upregulated_proteins_t_pbc_h.m_no_c <- results.t.pbc.h.m.no.c %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_t_pbc_h.m_no_c) # 242

# Number of downregulated proteins in PBC:
selected_downregulated_proteins_t_pbc_h.m_no_c <- results.t.pbc.h.m.no.c %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_t_pbc_h.m_no_c) # 94

# Sort the upregulated proteins by adjusted p-value
sorted_upregulated_proteins_t_pbc_h.m_no_c <- 
  arrange(subset(results.t.pbc.h.m.no.c, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of upregulated proteins with gene names
print(sorted_upregulated_proteins_t_pbc_h.m_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
up_liver_pbc_vs_h.m_no_c <- 
  sorted_upregulated_proteins_t_pbc_h.m_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]

# Sort the downregulated proteins by adjusted p-value
sorted_downregulated_proteins_t_pbc_h.m_no_c <- 
  arrange(subset(results.t.pbc.h.m.no.c, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of downregulated proteins with gene names
print(sorted_downregulated_proteins_t_pbc_h.m_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
down_liver_pbc_vs_h.m_no_c <- 
  sorted_downregulated_proteins_t_pbc_h.m_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]





#### PBC without MASLD vs. Healthy and MASLD

table(t_no_masld$group_1)

t_pbc_h.m_no_masld <- subset(t_no_masld, group_1 %in% c("PBC","Control"))

table(t_pbc_h.m_no_masld$group_1)

design_t_pbc_h.m_no_masld <- model.matrix(~ 0 + group_1, t_pbc_h.m_no_masld)
design_t_pbc_h.m_no_masld[1:20,]

dim(t_pbc_h.m_no_masld)
t_pbc_h.m_no_masld[1:3,1:5]
t_pbc_h.m_no_masld[1:3,7632:7636]

t_pbc_h.m_no_masld <- t_pbc_h.m_no_masld[,3:7634] # create only numeric values 
dim(t_pbc_h.m_no_masld)

fit_t_pbc_h.m_no_masld <- lmFit(t(t_pbc_h.m_no_masld), design_t_pbc_h.m_no_masld)
head(coef(fit_t_pbc_h.m_no_masld))

contr_t_pbc_h.m_no_masld <- makeContrasts(group_1PBC-group_1Control, levels=design_t_pbc_h.m_no_masld)
contr_t_pbc_h.m_no_masld
tmp_t_pbc_h.m_no_masld <- contrasts.fit(fit_t_pbc_h.m_no_masld, contr_t_pbc_h.m_no_masld)
tmp_t_pbc_h.m_no_masld <- eBayes(tmp_t_pbc_h.m_no_masld)
tmp_t_pbc_h.m_no_masld
top.table.t.pbc.h.m.no.masld <- topTable(tmp_t_pbc_h.m_no_masld,coef = 1, sort.by = "P", n = Inf)
results.t.pbc.h.m.no.masld <-  as.data.frame(top.table.t.pbc.h.m.no.masld)
results.t.pbc.h.m.no.masld$protein <- row.names(results.t.pbc.h.m.no.masld)
names(results.t.pbc.h.m.no.masld)
results.t.pbc.h.m.no.masld[1:20,1:7]
results.t.pbc.h.m.no.masld$GeneName <- sapply(results.t.pbc.h.m.no.masld$protein, 
                                          extract_gene_name)
results.t.pbc.h.m.no.masld[1:20,5:8]
results.t.pbc.h.m.no.masld$UniProtID <- sapply(results.t.pbc.h.m.no.masld$protein, 
                                           extract_uniprot_id)
results.t.pbc.h.m.no.masld$GeneName.UniprotID <- results.t.pbc.h.m.no.masld$protein
results.t.pbc.h.m.no.masld[1:20,6:10]
results.t.pbc.h.m.no.masld$log2FC <- results.t.pbc.h.m.no.masld$logFC
results.t.pbc.h.m.no.masld[1:3,c(1,8:11)]

# volcano PBC without MASLD vs Healthy and MASLD 

selected_proteins_t_pbc_h.m_no_masld <- results.t.pbc.h.m.no.masld %>%
  filter(adj.P.Val < 0.05)

selected_proteins_t_pbc_h.m_no_masld[1:20,1:5]

# As volcano plot cannot print all the names and becomes impossible to view, I filter only those proteins
# that are regulated more than 1 or less than -1
selected_proteins_t_pbc_h.m_no_masld_volcano <- results.t.pbc.h.m.no.masld %>%
  filter(logFC >= 1 & adj.P.Val < 0.05 | logFC <= -1 & adj.P.Val < 0.05)

volcano_plot_t_pbc_h.m_no_masld <- ggplot(results.t.pbc.h.m.no.masld) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), 
                 color = ifelse(logFC <= -0.00001  & adj.P.Val < 0.05, 'blue', 
                                ifelse(logFC >= 0.00001 & adj.P.Val < 0.05, 'red', 'grey')))) +
  geom_text(data = selected_proteins_t_pbc_h.m_no_masld_volcano,
            aes(x = logFC, y = -log10(adj.P.Val), 
                label = GeneName),
            size = 3, hjust = 0, vjust = 0) +
  scale_color_identity() +
  labs(x = "log2 fold change", y = "-log10 adj p-value") +
  ggtitle("Differential Liver Protein Expression, PBC without MASLD vs. Healthy and MASLD")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray")  # Add horizontal line


volcano_plot_t_pbc_h.m_no_masld

# tapply(t$ANKRD22.Q5VYY1,t$group,mean) # "checking that up and downregulated are correct (-1,+1)"

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/No MASLD/volcano_liver_pbc_h.m_no_masld_fdr0.05.jpeg",
#         width = 10, height = 10, plot = volcano_plot_t_pbc_h.m_no_masld, device = "jpeg")

# Number of upregulated proteins in PBC:

selected_upregulated_proteins_t_pbc_h.m_no_masld <- results.t.pbc.h.m.no.masld %>%
  filter(logFC >= 0.0001 & adj.P.Val < 0.05)
nrow(selected_upregulated_proteins_t_pbc_h.m_no_masld) # 242

# Number of downregulated proteins in PBC:
selected_downregulated_proteins_t_pbc_h.m_no_masld <- results.t.pbc.h.m.no.masld %>%
  filter(logFC < 0.0001 & adj.P.Val < 0.05)
nrow(selected_downregulated_proteins_t_pbc_h.m_no_masld) # 94

# Sort the upregulated proteins by adjusted p-value
sorted_upregulated_proteins_t_pbc_h.m_no_masld <- 
  arrange(subset(results.t.pbc.h.m.no.masld, logFC >= 0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of upregulated proteins with gene names
print(sorted_upregulated_proteins_t_pbc_h.m_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
up_liver_pbc_vs_h.m_no_masld <- 
  sorted_upregulated_proteins_t_pbc_h.m_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]

# Sort the downregulated proteins by adjusted p-value
sorted_downregulated_proteins_t_pbc_h.m_no_masld <- 
  arrange(subset(results.t.pbc.h.m.no.masld, logFC <= -0.0001 & adj.P.Val < 0.05), adj.P.Val)

# Display the sorted list of downregulated proteins with gene names
print(sorted_downregulated_proteins_t_pbc_h.m_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")])

# Create a data frame with the desired output
down_liver_pbc_vs_h.m_no_masld <- 
  sorted_downregulated_proteins_t_pbc_h.m_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC", "adj.P.Val")]




##################### Paired liver and plasma ###########################################

# Combine Results on GeneTables with Subject IDs, excluding proteins not present in both datasets
combined_results <- merge(results.p.pbc, results.t.pbc,
                          by = c("protein", "GeneName", "UniProtID", "GeneName.UniprotID"),
                          suffixes = c("_plasma", "_liver"),
                          all = FALSE)

nrow(combined_results) # 414 proteins present in both liver and plasma

combined_results[1:3,]
combined_results$FC_liver <- 2^combined_results$log2FC_liver
combined_results$FC_plasma <- 2^combined_results$log2FC_plasma
combined_results[1:3,]

# Correlation Analysis with Paired Samples

correlation_df <- cor.test(combined_results$log2FC_plasma, combined_results$log2FC_liver, method = "pearson", paired = TRUE)
correlation_df

# Values for correlation coefficient, confidence interval, and p-value used for plot underneath
rho_value <- round(correlation_df$estimate,3)
ci_lower <- round(correlation_df$conf.int[1],3)
ci_upper <- round(correlation_df$conf.int[2],3)
p_value <- "< 0.001" # remember to change if another correlation anaylsis with another p-value

correlation_liver_plasma_pbc_healthy <- 
  ggplot(combined_results, aes(x = log2FC_liver, y = log2FC_plasma)) +
  geom_point() +
  labs(x = "Log2 FC Liver", y = "Log2 FC Plasma") +
  ggtitle("Correlation of Protein Regulation in PBC compared with Healthy") +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = "blue") + 
  annotate("text", x = Inf, y = -Inf, 
           label = paste0("italic(r) == ", rho_value, " * ' (95% CI ", ci_lower, "-", ci_upper, ", ' * italic(p) * ' ", p_value, ")'"), 
           parse = TRUE, hjust = 1.1, vjust = -2)

correlation_liver_plasma_pbc_healthy

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/All/scatter_correlation_fc_liver_plasma_pbc_healthy.jpeg",
#         width = 10, height = 10, plot = correlation_liver_plasma_pbc_healthy, device = "jpeg")

# only for significant proteins:

combined_results_sign <- combined_results %>% filter(adj.P.Val_plasma < 0.05 & adj.P.Val_liver < 0.05)

nrow(combined_results_sign) # 26 proteins significantly regulated in both liver and plasma (PBC vs Healthy)

# Correlation Analysis with Paired Samples sigificant proteins
correlation_df_sign <- cor.test(combined_results_sign$log2FC_plasma, combined_results_sign$log2FC_liver, 
                                method = "pearson", paired = TRUE)
correlation_df_sign
combined_results_sign[21:30,]

######## fc scatterplot

correlation_sign_liver_plasma_pbc_healthy <- ggplot(combined_results_sign, aes(x = FC_liver, 
                                                                      y = FC_plasma)) +
  geom_point() +
  geom_text_repel(aes(label = GeneName), 
                  box.padding = 0.5, # Adjust padding around text labels
                  max.overlaps = Inf, # Allow for more attempts to avoid overlaps
                  segment.color = "black", # Line color
                  segment.size = 0.5, # Line width
                  segment.curvature = 0.2, # Line curvature
                  segment.alpha = 0.5, # Line transparency
                  hjust = 0, vjust = 0) + # Text label alignment
    labs(x = "FC Liver", y = "FC Plasma") +
  ggtitle("Differentially expressed proteins in liver and plasma in PBC compared with healthy") 

correlation_sign_liver_plasma_pbc_healthy

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/All/scatter_correlation_fc_sign_liver_plasma_all_proteins_named_pbc_healthy_fdr0.05.jpeg",
#         width = 8, height = 6, plot = correlation_sign_liver_plasma_pbc_healthy, device = "jpeg")

# Correlation analysis
correlation_df_sign <- cor.test(combined_results_sign$log2FC_plasma, 
                                combined_results_sign$log2FC_liver, 
                                method = "pearson", paired = TRUE)



# log2fc scatterplot 

correlation_sign_log2fc_liver_plasma_pbc_healthy <- ggplot(combined_results_sign, 
                                                      aes(x = log2FC_liver, y = logFC_plasma)) +
  geom_point() +
  geom_text_repel(aes(label = GeneName), max.overlaps = Inf) +  # Tilføjet max.overlaps
  labs(x = "log2FC Liver", y = "log2FC Plasma") +
  ggtitle("Differentially expressed proteins in liver and plasma in PBC compared with healthy") +
  theme(
    axis.title = element_text(size = 14),   # Increase axis title size
    axis.text = element_text(size = 14),    # Increase axis text (tick labels) size
    plot.title = element_text(size = 14)    # Increase plot title size
  )


correlation_sign_log2fc_liver_plasma_pbc_healthy

# save plot
# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/All/scatter_correlation_log2fc_sign_liver_plasma_pbc_healthy_fdr0.05.jpeg",
#         width = 8, height = 6, plot = correlation_sign_log2fc_liver_plasma_pbc_healthy, device = "jpeg")


# log2fc scatterplot liver/plasma with colouring according to regulation in liver and plasma
## Create a new column in your dataframe for colouring of only the dots
combined_results_sign$Regulation <- with(combined_results_sign, ifelse(
  log2FC_plasma > 0 & log2FC_liver    > 0, "Upregulated in plasma and liver",
  ifelse(log2FC_plasma > 0 & log2FC_liver < 0, "Upregulated in plasma and downregulated in liver",
         "Downregulated in plasma and liver")
))

# 2) Force factor levels into the exact legend‐order you want
combined_results_sign$Regulation <- factor(
  combined_results_sign$Regulation,
  levels = c(
    "Upregulated in plasma and liver",
    "Upregulated in plasma and downregulated in liver",
    "Downregulated in plasma and liver"
  )
)

correlation_sign_log2fc_liver_plasma_pbc_healthy_dots <- 
  ggplot(combined_results_sign, aes(x = log2FC_liver, y = log2FC_plasma)) +
  geom_point(aes(color = Regulation), size = 2) +
  geom_text_repel(aes(label = GeneName),
                  color = "black", max.overlaps = Inf) +
  scale_color_manual(
    values = c(
      "Upregulated in plasma and liver"                  = "red",
      "Upregulated in plasma and downregulated in liver" = "goldenrod",
      "Downregulated in plasma and liver"                = "blue"
    ),
    breaks = levels(combined_results_sign$Regulation)
  ) +
  labs(
    x     = "log2FC Liver",
    y     = "log2FC Plasma",
    title = "Differentially expressed proteins in liver and plasma in PBC compared with healthy"
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 12),
    plot.title      = element_text(size = 12),
    legend.title    = element_blank(),
    legend.text     = element_text(size = 8),
    legend.position = "bottom",        # <— here
    legend.direction= "horizontal",    # <— here
    legend.key.width= unit(1.5, "lines"),
    legend.spacing.x= unit(0.5, "lines")
  )

correlation_sign_log2fc_liver_plasma_pbc_healthy_dots

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/All/scatter_correlation_log2fc_sign_liver_plasma_pbc_healthy_dots_fdr0.05.jpeg",
#         width = 8, height = 6, plot = correlation_sign_log2fc_liver_plasma_pbc_healthy_dots, device = "jpeg")

# Identify Commonly Regulated Proteins
common_upregulated <- combined_results_sign %>% filter(log2FC_plasma >= 0 & log2FC_liver >= 0)
nrow(common_upregulated) # 9

common_up_pbc_vs_healthy <- 
  common_upregulated[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                         "adj.P.Val_liver")]

common_downregulated <- combined_results_sign %>% filter(log2FC_plasma < 0 & log2FC_liver < 0)
nrow(common_downregulated) # 15
common_down_pbc_vs_healthy <- 
  common_downregulated[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                         "adj.P.Val_liver")]

common_upreg.liver_downreg.plasma <- combined_results_sign %>% filter(log2FC_plasma < 0 & log2FC_liver >= 0)
nrow(common_upreg.liver_downreg.plasma) # 0
common_up.liver_down.plasma_pbc_vs_healthy <- 
  common_upreg.liver_downreg.plasma[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                         "adj.P.Val_liver")]

common_downreg.liver_upreg.plasma <- combined_results_sign %>% filter(log2FC_plasma >= 0 & log2FC_liver < 0)
nrow(common_downreg.liver_upreg.plasma) # 2

common_down.liver_up.plasma_pbc_vs_healthy <- 
  common_downreg.liver_upreg.plasma[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                         "adj.P.Val_liver")]

######################## Save in an excel file

# write.xlsx(list("downreg_plasma_pbc_vs_healthy" = down_plasma_pbc_vs_healthy,
#                 "upreg_plasma_pbc_vs_healthy" = up_plasma_pbc_vs_healthy,
#                 "downreg_liver_pbc_vs_healthy" = down_liver_pbc_vs_healthy,
#                 "upreg_liver_pbc_vs_healthy" = up_liver_pbc_vs_healthy,
#                 "up_liver_and_plasma_pbc" = common_up_pbc_vs_healthy,
#                 "down_liver_and_plasma_pbc" = common_down_pbc_vs_healthy,
#                 "up_liver_down_plasma_pbc" = common_up.liver_down.plasma_pbc_vs_healthy,
#                 "down_liver_up_plasma_pbc" = common_down.liver_up.plasma_pbc_vs_healthy),
#            "../output/FDR0.05_all/PBC_vs_Healthy/All/DE_proteins_fdr0.05_pbc_vs_healthy.xlsx")





##################### Paired liver and plasma - without outliers ###########################################

# Combine Results on GeneTables with Subject IDs, excluding proteins not present in both datasets
combined_results_no_outliers <- merge(
  results.p.pbc_no_outliers,
  results.t.pbc.no_outliers,
  by = c("protein", "GeneName", "UniProtID", "GeneName.UniprotID"),
  suffixes = c("_plasma", "_liver"),
  all = FALSE
)

nrow(combined_results_no_outliers)

combined_results_no_outliers[1:3,]

combined_results_no_outliers$FC_liver <- 2^combined_results_no_outliers$log2FC_liver
combined_results_no_outliers$FC_plasma <- 2^combined_results_no_outliers$log2FC_plasma

combined_results_no_outliers[1:3,]


# Correlation Analysis with Paired Samples

correlation_df_no_outliers <- cor.test(
  combined_results_no_outliers$log2FC_plasma,
  combined_results_no_outliers$log2FC_liver,
  method = "pearson",
  paired = TRUE
)

correlation_df_no_outliers


# Values for correlation coefficient, confidence interval, and p-value used for plot underneath

rho_value_no_outliers <- round(correlation_df_no_outliers$estimate, 3)
ci_lower_no_outliers <- round(correlation_df_no_outliers$conf.int[1], 3)
ci_upper_no_outliers <- round(correlation_df_no_outliers$conf.int[2], 3)
p_value_no_outliers <- "< 0.001"


correlation_liver_plasma_pbc_healthy_no_outliers <- 
  ggplot(combined_results_no_outliers, aes(x = log2FC_liver, y = log2FC_plasma)) +
  geom_point() +
  labs(x = "Log2 FC Liver", y = "Log2 FC Plasma") +
  ggtitle("Correlation of Protein Regulation in PBC compared with Healthy without outliers") +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = "blue") + 
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = paste0(
      "italic(r) == ", rho_value_no_outliers,
      " * ' (95% CI ", ci_lower_no_outliers, "-", ci_upper_no_outliers,
      ", ' * italic(p) * ' ", p_value_no_outliers, ")'"
    ),
    parse = TRUE,
    hjust = 1.1,
    vjust = -2
  )

correlation_liver_plasma_pbc_healthy_no_outliers


# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/Without outliers/scatter_correlation_fc_liver_plasma_pbc_healthy_without_outliers.jpeg",
#         width = 10, height = 10, plot = correlation_liver_plasma_pbc_healthy_no_outliers, device = "jpeg")


# only for significant proteins:

combined_results_sign_no_outliers <- combined_results_no_outliers %>%
  filter(adj.P.Val_plasma < 0.05 & adj.P.Val_liver < 0.05)

nrow(combined_results_sign_no_outliers)


# Correlation Analysis with Paired Samples significant proteins

correlation_df_sign_no_outliers <- cor.test(
  combined_results_sign_no_outliers$log2FC_plasma,
  combined_results_sign_no_outliers$log2FC_liver,
  method = "pearson",
  paired = TRUE
)

correlation_df_sign_no_outliers

combined_results_sign_no_outliers[1:10,]


######## fc scatterplot

correlation_sign_liver_plasma_pbc_healthy_no_outliers <- 
  ggplot(combined_results_sign_no_outliers, aes(x = FC_liver, y = FC_plasma)) +
  geom_point() +
  geom_text_repel(
    aes(label = GeneName),
    box.padding = 0.5,
    max.overlaps = Inf,
    segment.color = "black",
    segment.size = 0.5,
    segment.curvature = 0.2,
    segment.alpha = 0.5,
    hjust = 0,
    vjust = 0
  ) +
  labs(x = "FC Liver", y = "FC Plasma") +
  ggtitle("Differentially expressed proteins in liver and plasma in PBC compared with healthy without outliers")

correlation_sign_liver_plasma_pbc_healthy_no_outliers


# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/Without outliers/scatter_correlation_fc_sign_liver_plasma_all_proteins_named_pbc_healthy_without_outliers_fdr0.05.jpeg",
#         width = 8, height = 6, plot = correlation_sign_liver_plasma_pbc_healthy_no_outliers, device = "jpeg")


# log2fc scatterplot 

correlation_sign_log2fc_liver_plasma_pbc_healthy_no_outliers <- 
  ggplot(
    combined_results_sign_no_outliers,
    aes(x = log2FC_liver, y = logFC_plasma)
  ) +
  geom_point() +
  geom_text_repel(aes(label = GeneName), max.overlaps = Inf) +
  labs(x = "log2FC Liver", y = "log2FC Plasma") +
  ggtitle("Differentially expressed proteins in liver and plasma in PBC compared with healthy without outliers") +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14)
  )

correlation_sign_log2fc_liver_plasma_pbc_healthy_no_outliers


# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/Without outliers/scatter_correlation_log2fc_sign_liver_plasma_pbc_healthy_without_outliers_fdr0.05.jpeg",
#         width = 8, height = 6, plot = correlation_sign_log2fc_liver_plasma_pbc_healthy_no_outliers, device = "jpeg")


# log2fc scatterplot liver/plasma with colouring according to regulation in liver and plasma

combined_results_sign_no_outliers$Regulation <- with(
  combined_results_sign_no_outliers,
  ifelse(
    log2FC_plasma > 0 & log2FC_liver > 0,
    "Upregulated in plasma and liver",
    ifelse(
      log2FC_plasma > 0 & log2FC_liver < 0,
      "Upregulated in plasma and downregulated in liver",
      "Downregulated in plasma and liver"
    )
  )
)

combined_results_sign_no_outliers$Regulation <- factor(
  combined_results_sign_no_outliers$Regulation,
  levels = c(
    "Upregulated in plasma and liver",
    "Upregulated in plasma and downregulated in liver",
    "Downregulated in plasma and liver"
  )
)

correlation_sign_log2fc_liver_plasma_pbc_healthy_dots_no_outliers <- 
  ggplot(combined_results_sign_no_outliers, aes(x = log2FC_liver, y = log2FC_plasma)) +
  geom_point(aes(color = Regulation), size = 2) +
  geom_text_repel(aes(label = GeneName), color = "black", max.overlaps = Inf) +
  scale_color_manual(
    values = c(
      "Upregulated in plasma and liver"                  = "red",
      "Upregulated in plasma and downregulated in liver" = "goldenrod",
      "Downregulated in plasma and liver"                = "blue"
    ),
    breaks = levels(combined_results_sign_no_outliers$Regulation)
  ) +
  labs(
    x     = "log2FC Liver",
    y     = "log2FC Plasma",
    title = "Differentially expressed proteins in liver and plasma in PBC compared with healthy without outliers"
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 12),
    plot.title      = element_text(size = 12),
    legend.title    = element_blank(),
    legend.text     = element_text(size = 8),
    legend.position = "bottom",
    legend.direction= "horizontal",
    legend.key.width= unit(1.5, "lines"),
    legend.spacing.x= unit(0.5, "lines")
  )

correlation_sign_log2fc_liver_plasma_pbc_healthy_dots_no_outliers


# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/Without outliers/scatter_correlation_log2fc_sign_liver_plasma_pbc_healthy_dots_without_outliers_fdr0.05.jpeg",
#         width = 8, height = 6, plot = correlation_sign_log2fc_liver_plasma_pbc_healthy_dots_no_outliers, device = "jpeg")


# Identify Commonly Regulated Proteins

common_upregulated_no_outliers <- combined_results_sign_no_outliers %>%
  filter(log2FC_plasma >= 0 & log2FC_liver >= 0)

nrow(common_upregulated_no_outliers)

common_up_pbc_vs_healthy_no_outliers <- 
  common_upregulated_no_outliers[, c(
    "GeneName.UniprotID",
    "GeneName",
    "UniProtID",
    "log2FC_plasma",
    "adj.P.Val_plasma",
    "log2FC_liver",
    "adj.P.Val_liver"
  )]


common_downregulated_no_outliers <- combined_results_sign_no_outliers %>%
  filter(log2FC_plasma < 0 & log2FC_liver < 0)

nrow(common_downregulated_no_outliers)

common_down_pbc_vs_healthy_no_outliers <- 
  common_downregulated_no_outliers[, c(
    "GeneName.UniprotID",
    "GeneName",
    "UniProtID",
    "log2FC_plasma",
    "adj.P.Val_plasma",
    "log2FC_liver",
    "adj.P.Val_liver"
  )]


common_upreg.liver_downreg.plasma_no_outliers <- combined_results_sign_no_outliers %>%
  filter(log2FC_plasma < 0 & log2FC_liver >= 0)

nrow(common_upreg.liver_downreg.plasma_no_outliers)

common_up.liver_down.plasma_pbc_vs_healthy_no_outliers <- 
  common_upreg.liver_downreg.plasma_no_outliers[, c(
    "GeneName.UniprotID",
    "GeneName",
    "UniProtID",
    "log2FC_plasma",
    "adj.P.Val_plasma",
    "log2FC_liver",
    "adj.P.Val_liver"
  )]


common_downreg.liver_upreg.plasma_no_outliers <- combined_results_sign_no_outliers %>%
  filter(log2FC_plasma >= 0 & log2FC_liver < 0)

nrow(common_downreg.liver_upreg.plasma_no_outliers)

common_down.liver_up.plasma_pbc_vs_healthy_no_outliers <- 
  common_downreg.liver_upreg.plasma_no_outliers[, c(
    "GeneName.UniprotID",
    "GeneName",
    "UniProtID",
    "log2FC_plasma",
    "adj.P.Val_plasma",
    "log2FC_liver",
    "adj.P.Val_liver"
  )]


######################## Save in an excel file

# write.xlsx(list(
#   "downreg_plasma_pbc_vs_healthy" = down_plasma_pbc_vs_healthy_no_outliers,
#   "upreg_plasma_pbc_vs_healthy" = up_plasma_pbc_vs_healthy_no_outliers,
#   "downreg_liver_pbc_vs_healthy" = down_liver_pbc_vs_healthy_no_outliers,
#   "upreg_liver_pbc_vs_healthy" = up_liver_pbc_vs_healthy_no_outliers,
#   "up_liver_and_plasma_pbc" = common_up_pbc_vs_healthy_no_outliers,
#   "down_liver_and_plasma_pbc" = common_down_pbc_vs_healthy_no_outliers,
#   "up_liver_down_plasma_pbc" = common_up.liver_down.plasma_pbc_vs_healthy_no_outliers,
#   "down_liver_up_plasma_pbc" = common_down.liver_up.plasma_pbc_vs_healthy_no_outliers),
#   "../output/FDR0.05_all/PBC_vs_Healthy/Without outliers/DE_proteins_fdr0.05_pbc_vs_healthy_without_outliers.xlsx")




##################### Paired liver and plasma without cirrhosis #######################

# Combine Results Tables with Subject IDs, excluding proteins not present in both datasets
combined_results_no_c <- merge(results.p.pbc.no.c, results.t.pbc.no.c,
                          by = c("protein", "GeneName", "UniProtID", "GeneName.UniprotID"),
                          suffixes = c("_plasma", "_liver"),
                          all = FALSE)

nrow(combined_results_no_c) # 414 proteins present in both liver and plasma
combined_results_no_c[1:3,]
combined_results_no_c$FC_liver <- 2^combined_results_no_c$log2FC_liver
combined_results_no_c$FC_plasma <- 2^combined_results_no_c$log2FC_plasma
combined_results_no_c[1:3,]

# Correlation fold change Analysis with Paired Samples

correlation_df_no_c <- cor.test(combined_results_no_c$log2FC_plasma, 
                                combined_results_no_c$log2FC_liver, method = "pearson", paired = TRUE)
correlation_df_no_c

# Values for correlation coefficient, confidence interval, and p-value used for plot underneath
rho_value_no_c <- round(correlation_df_no_c$estimate,3)
ci_lower_no_c <- round(correlation_df_no_c$conf.int[1],3)
ci_upper_no_c <- round(correlation_df_no_c$conf.int[2],3)
p_value_no_c <- "< 0.001" # remember to change if another correlation anaylsis with another p-value

# Correlation plot

correlation_liver_plasma_pbc_healthy_no_c <- 
  ggplot(combined_results_no_c, aes(x = log2FC_liver, y = log2FC_plasma)) +
  geom_point() +
  labs(x = "Log2 FC Liver", y = "Log2 FC Plasma") +
  ggtitle("Correlation of Protein Regulation in non-cirrhotic PBC compared with Healthy") +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = "blue") + 
  annotate("text", x = Inf, y = -Inf, 
           label = paste0("italic(r) == ", rho_value_no_c, " * ' (95% CI ", ci_lower_no_c, "-", ci_upper_no_c, ", ' * italic(p) * ' ", p_value_no_c, ")'"), 
           parse = TRUE, hjust = 1.1, vjust = -2)

correlation_liver_plasma_pbc_healthy_no_c

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/No cirrhosis/scatter_correlation_fc_liver_plasma_pbc_healthy_no_cirrhosis.jpeg",
#         width = 10, height = 10, plot = correlation_liver_plasma_pbc_healthy_no_c, device = "jpeg")
# only for significant proteins:

combined_results_sign_no_c <- combined_results_no_c %>% filter(adj.P.Val_plasma < 0.05 & adj.P.Val_liver < 0.05)

nrow(combined_results_sign_no_c) # 26 proteins significantly regulated in both liver and plasma

# Correlation Analysis with Paired Samples sigificant proteins
correlation_df_sign_no_c <- cor.test(combined_results_sign_no_c$log2FC_plasma, combined_results_sign_no_c$log2FC_liver, 
                                     method = "pearson", paired = TRUE)
correlation_df_sign_no_c
combined_results_sign_no_c[1:5,]

correlation_sign_liver_plasma_pbc_healthy_no_c <- ggplot(combined_results_sign_no_c, aes(x = FC_liver, 
                                                                                         y = FC_plasma)) +
  geom_point() +
  geom_text_repel(aes(label = GeneName), 
                  box.padding = 0.5, # Adjust padding around text labels
                  max.overlaps = Inf, # Allow for more attempts to avoid overlaps
                  segment.color = "black", # Line color
                  segment.size = 0.5, # Line width
                  segment.curvature = 0.2, # Line curvature
                  segment.alpha = 0.5, # Line transparency
                  hjust = 0, vjust = 0) + # Text label alignment
  labs(x = "FC Liver", y = "FC Plasma") +
  ggtitle("Differentially expressed proteins in liver and plasma in non-cirrhotic PBC compared with healthy") 

correlation_sign_liver_plasma_pbc_healthy_no_c

# Save plot
# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/No cirrhosis/scatter_correlation_fc_sign_liver_plasma_all_proteins_named_pbc_healthy_no_cirrhosis.jpeg",
#         width = 10, height = 10, plot = correlation_sign_liver_plasma_pbc_healthy_no_c, device = "jpeg")

# log2fc scatterplot liver/plasma with colouring according to regulation in liver and plasma
## Create a new column in your dataframe for colouring of only the dots
combined_results_sign_no_c$Regulation <- with(combined_results_sign_no_c, ifelse(
  log2FC_plasma > 0 & log2FC_liver    > 0, "Upregulated in plasma and liver",
  ifelse(log2FC_plasma > 0 & log2FC_liver < 0, "Upregulated in plasma and downregulated in liver",
         "Downregulated in plasma and liver")
))

# 2) Force factor levels into the exact legend‐order you want
combined_results_sign_no_c$Regulation <- factor(
  combined_results_sign_no_c$Regulation,
  levels = c(
    "Upregulated in plasma and liver",
    "Upregulated in plasma and downregulated in liver",
    "Downregulated in plasma and liver"
  )
)

correlation_sign_no_c_log2fc_liver_plasma_pbc_healthy_dots <- 
  ggplot(combined_results_sign_no_c, aes(x = log2FC_liver, y = log2FC_plasma)) +
  geom_point(aes(color = Regulation), size = 2) +
  geom_text_repel(aes(label = GeneName),
                  color = "black", max.overlaps = Inf) +
  scale_color_manual(
    values = c(
      "Upregulated in plasma and liver"                  = "red",
      "Upregulated in plasma and downregulated in liver" = "goldenrod",
      "Downregulated in plasma and liver"                = "blue"
    ),
    breaks = levels(combined_results_sign_no_c$Regulation)
  ) +
  labs(
    x     = "log2FC Liver",
    y     = "log2FC Plasma",
    title = "Differentially expressed proteins in liver and plasma in non-cirrhotic PBC compared with healthy"
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 12),
    plot.title      = element_text(size = 12),
    legend.title    = element_blank(),
    legend.text     = element_text(size = 8),
    legend.position = "bottom",        # <— here
    legend.direction= "horizontal",    # <— here
    legend.key.width= unit(1.5, "lines"),
    legend.spacing.x= unit(0.5, "lines")
  )

correlation_sign_no_c_log2fc_liver_plasma_pbc_healthy_dots

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/No cirrhosis/scatter_correlation_log2fc_sign_no_cirrhosis_liver_plasma_pbc_healthy_dots_fdr0.05.jpeg",
#         width = 8, height = 6, plot = correlation_sign_no_c_log2fc_liver_plasma_pbc_healthy_dots, device = "jpeg")
 
# Identify Commonly Regulated Proteins
common_upregulated_no_c <- combined_results_sign_no_c %>% filter(log2FC_plasma >= 0 & log2FC_liver >= 0)
nrow(common_upregulated_no_c) # 19

common_up_pbc_vs_healthy_no_c <- 
  common_upregulated_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                         "adj.P.Val_liver")]

common_downregulated_no_c <- combined_results_sign_no_c %>% filter(log2FC_plasma < 0 & log2FC_liver < 0)
nrow(common_downregulated_no_c) # 9
common_down_pbc_vs_healthy_no_c <- 
  common_downregulated_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                           "adj.P.Val_liver")]

common_upreg.liver_downreg.plasma.no.c <- combined_results_sign_no_c %>% filter(log2FC_plasma < 0 & log2FC_liver >= 0)
nrow(common_upreg.liver_downreg.plasma.no.c) # 1
common_up.liver_down.plasma_pbc_vs_healthy_no_c <- 
  common_upreg.liver_downreg.plasma.no.c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                        "adj.P.Val_liver")]

common_downreg.liver_upreg.plasma.no.c <- combined_results_sign_no_c %>% filter(log2FC_plasma >= 0 & log2FC_liver < 0)
nrow(common_downreg.liver_upreg.plasma.no.c) # 5

common_down.liver_up.plasma_pbc_vs_healthy_no_c <- 
  common_downreg.liver_upreg.plasma.no.c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                        "adj.P.Val_liver")]

######################## Save in an excel file

# write.xlsx(list("downreg_plasma_pbc_vs_healthy" = down_plasma_pbc_vs_healthy_no_c,
#                 "upreg_plasma_pbc_vs_healthy" = up_plasma_pbc_vs_healthy_no_c,
#                 "downreg_liver_pbc_vs_healthy" = down_liver_pbc_vs_healthy_no_c,
#                 "upreg_liver_pbc_vs_healthy" = up_liver_pbc_vs_healthy_no_c,
#                 "up_liver_and_plasma_pbc_no" = common_up_pbc_vs_healthy_no_c,
#                 "down_liver_and_plasma_pbc" = common_down_pbc_vs_healthy_no_c,
#                 "up_liver_down_plasma_pbc" = common_up.liver_down.plasma_pbc_vs_healthy_no_c,
#                 "down_liver_up_plasma_pbc" = common_down.liver_up.plasma_pbc_vs_healthy_no_c),
#            "../output/FDR0.05_all/PBC_vs_Healthy/No cirrhosis/DE_proteins_fdr0.05_pbc_vs_healthy_no.cirrhosis.xlsx")

##################### Paired liver and plasma without MASLD in PBC group #######################

# Combine Results Tables with Subject IDs, excluding proteins not present in both datasets
combined_results_no_masld <- merge(results.p.pbc.no.masld, results.t.pbc.no.masld,
                               by = c("protein", "GeneName", "UniProtID", "GeneName.UniprotID"),
                               suffixes = c("_plasma", "_liver"),
                               all = FALSE)

nrow(combined_results_no_masld) # 414 proteins present in both liver and plasma
combined_results_no_masld[1:3,]
combined_results_no_masld$FC_liver <- 2^combined_results_no_masld$log2FC_liver
combined_results_no_masld$FC_plasma <- 2^combined_results_no_masld$log2FC_plasma
combined_results_no_masld[1:3,]

# Correlation fold change Analysis with Paired Samples

correlation_df_no_masld <- cor.test(combined_results_no_masld$log2FC_plasma, 
                                combined_results_no_masld$log2FC_liver, method = "pearson", paired = TRUE)
correlation_df_no_masld

# Values for correlation coefficient, confidence interval, and p-value used for plot underneath
rho_value_no_masld <- round(correlation_df_no_masld$estimate,3)
ci_lower_no_masld <- round(correlation_df_no_masld$conf.int[1],3)
ci_upper_no_masld <- round(correlation_df_no_masld$conf.int[2],3)
p_value_no_masld <- "< 0.001" # remember to change if another correlation anaylsis with another p-value

# Correlation plot

correlation_liver_plasma_pbc_healthy_no_masld <- 
  ggplot(combined_results_no_masld, aes(x = log2FC_liver, y = log2FC_plasma)) +
  geom_point() +
  labs(x = "Log2 FC Liver", y = "Log2 FC Plasma") +
  ggtitle("Correlation of Protein Regulation in PBC without MASLD compared with Healthy") +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = "blue") + 
  annotate("text", x = Inf, y = -Inf, 
           label = paste0("italic(r) == ", rho_value_no_masld, " * ' (95% CI ", ci_lower_no_masld, "-", ci_upper_no_masld, ", ' * italic(p) * ' ", p_value_no_masld, ")'"), 
           parse = TRUE, hjust = 1.1, vjust = -2)

correlation_liver_plasma_pbc_healthy_no_masld

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/No MASLD/scatter_correlation_fc_liver_plasma_pbc_healthy_no_masld.jpeg",
#         width = 10, height = 10, plot = correlation_liver_plasma_pbc_healthy_no_masld, device = "jpeg")

# only for significant proteins:

combined_results_sign_no_masld <- combined_results_no_masld %>% filter(adj.P.Val_plasma < 0.05 & adj.P.Val_liver < 0.05)

nrow(combined_results_sign_no_masld) # ONLY 11 proteins significantly regulated in both liver and plasma

# Correlation Analysis with Paired Samples sigificant proteins
correlation_df_sign_no_masld <- cor.test(combined_results_sign_no_masld$log2FC_plasma, combined_results_sign_no_masld$log2FC_liver, 
                                     method = "pearson", paired = TRUE)
correlation_df_sign_no_masld
combined_results_sign_no_masld[1:5,]

correlation_sign_liver_plasma_pbc_healthy_no_masld <- ggplot(combined_results_sign_no_masld, aes(x = FC_liver, 
                                                                                         y = FC_plasma)) +
  geom_point() +
  geom_text_repel(aes(label = GeneName), 
                  box.padding = 0.5, # Adjust padding around text labels
                  max.overlaps = Inf, # Allow for more attempts to avoid overlaps
                  segment.color = "black", # Line color
                  segment.size = 0.5, # Line width
                  segment.curvature = 0.2, # Line curvature
                  segment.alpha = 0.5, # Line transparency
                  hjust = 0, vjust = 0) + # Text label alignment
  labs(x = "FC Liver", y = "FC Plasma") +
  ggtitle("Differentially expressed proteins in liver and plasma in PBC without MASLD compared with healthy") 

correlation_sign_liver_plasma_pbc_healthy_no_masld


# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/No MASLD/scatter_correlation_fc_sign_liver_plasma_all_proteins_named_pbc_healthy_no_masld.jpeg",
#         width = 10, height = 10, plot = correlation_sign_liver_plasma_pbc_healthy_no_masld, device = "jpeg")


correlation_sign_liver_plasma_pbc_healthy_no_masld_1 <- ggplot(combined_results_sign_no_masld, aes(x = FC_liver, 
                                                                                                   y = FC_plasma)) +
  geom_point() +
  geom_text_repel(aes(label = GeneName)) + # Text label alignment
  labs(x = "FC Liver", y = "FC Plasma") +
  ggtitle("Differentially expressed proteins in liver and plasma in PBC without MASLD compared with healthy") +
  theme(
    axis.title = element_text(size = 16),   # Increase axis title size
    axis.text = element_text(size = 14),    # Increase axis text (tick labels) size
    plot.title = element_text(size = 18)    # Increase plot title size
  ) 

correlation_sign_liver_plasma_pbc_healthy_no_masld_1

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/No MASLD/scatter_correlation_fc_sign_liver_plasma_pbc_healthy_no_masld.fdr0.05.jpeg",
#         width = 8, height = 6, plot = correlation_sign_liver_plasma_pbc_healthy_no_masld_1, device = "jpeg")


# log2fc scatterplot liver/plasma with colouring according to regulation in liver and plasma
## Create a new column in your dataframe for colouring of only the dots
combined_results_sign_no_masld$Regulation <- with(combined_results_sign_no_masld, ifelse(
  log2FC_plasma > 0 & log2FC_liver    > 0, "Upregulated in plasma and liver",
  ifelse(log2FC_plasma > 0 & log2FC_liver < 0, "Upregulated in plasma and downregulated in liver",
         "Downregulated in plasma and liver")
))

# 2) Force factor levels into the exact legend‐order you want
combined_results_sign_no_masld$Regulation <- factor(
  combined_results_sign_no_masld$Regulation,
  levels = c(
    "Upregulated in plasma and liver",
    "Upregulated in plasma and downregulated in liver",
    "Downregulated in plasma and liver"
  )
)

correlation_sign_no_masld_log2fc_liver_plasma_pbc_healthy_dots <- 
  ggplot(combined_results_sign_no_masld, aes(x = log2FC_liver, y = log2FC_plasma)) +
  geom_point(aes(color = Regulation), size = 2) +
  geom_text_repel(aes(label = GeneName),
                  color = "black", max.overlaps = Inf) +
  scale_color_manual(
    values = c(
      "Upregulated in plasma and liver"                  = "red",
      "Upregulated in plasma and downregulated in liver" = "goldenrod",
      "Downregulated in plasma and liver"                = "blue"
    ),
    breaks = levels(combined_results_sign_no_masld$Regulation)
  ) +
  labs(
    x     = "log2FC Liver",
    y     = "log2FC Plasma",
    title = "Differentially expressed proteins in liver and plasma in PBC without MASLD compared with healthy"
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 12),
    plot.title      = element_text(size = 12),
    legend.title    = element_blank(),
    legend.text     = element_text(size = 8),
    legend.position = "bottom",        # <— here
    legend.direction= "horizontal",    # <— here
    legend.key.width= unit(1.5, "lines"),
    legend.spacing.x= unit(0.5, "lines")
  )

correlation_sign_no_masld_log2fc_liver_plasma_pbc_healthy_dots

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/No MASLD/scatter_correlation_log2fc_sign_no_masld_liver_plasma_pbc_healthy_dots_fdr0.05.jpeg",
#         width = 8, height = 6, plot = correlation_sign_no_masld_log2fc_liver_plasma_pbc_healthy_dots, device = "jpeg")

# Identify Commonly Regulated Proteins
common_upregulated_no_masld <- combined_results_sign_no_masld %>% filter(log2FC_plasma >= 0 & log2FC_liver >= 0)
nrow(common_upregulated_no_masld) # 19

common_up_pbc_vs_healthy_no_masld <- 
  common_upregulated_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                              "adj.P.Val_liver")]

common_downregulated_no_masld <- combined_results_sign_no_masld %>% filter(log2FC_plasma < 0 & log2FC_liver < 0)
nrow(common_downregulated_no_masld) # 9
common_down_pbc_vs_healthy_no_masld <- 
  common_downregulated_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                "adj.P.Val_liver")]

common_upreg.liver_downreg.plasma.no.c <- combined_results_sign_no_masld %>% filter(log2FC_plasma < 0 & log2FC_liver >= 0)
nrow(common_upreg.liver_downreg.plasma.no.c) # 1
common_up.liver_down.plasma_pbc_vs_healthy_no_masld <- 
  common_upreg.liver_downreg.plasma.no.c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                             "adj.P.Val_liver")]

common_downreg.liver_upreg.plasma.no.c <- combined_results_sign_no_masld %>% filter(log2FC_plasma >= 0 & log2FC_liver < 0)
nrow(common_downreg.liver_upreg.plasma.no.c) # 5

common_down.liver_up.plasma_pbc_vs_healthy_no_masld <- 
  common_downreg.liver_upreg.plasma.no.c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                             "adj.P.Val_liver")]

######################## Save in an excel file

# write.xlsx(list("downreg_plasma_pbc_vs_healthy" = down_plasma_pbc_vs_healthy_no_masld,
#                 "upreg_plasma_pbc_vs_healthy" = up_plasma_pbc_vs_healthy_no_masld,
#                 "downreg_liver_pbc_vs_healthy" = down_liver_pbc_vs_healthy_no_masld,
#                 "upreg_liver_pbc_vs_healthy" = up_liver_pbc_vs_healthy_no_masld,
#                 "up_liver_and_plasma_pbc_no" = common_up_pbc_vs_healthy_no_masld,
#                 "down_liver_and_plasma_pbc" = common_down_pbc_vs_healthy_no_masld,
#                 "up_liver_down_plasma_pbc" = common_up.liver_down.plasma_pbc_vs_healthy_no_masld,
#                 "down_liver_up_plasma_pbc" = common_down.liver_up.plasma_pbc_vs_healthy_no_masld),
#            "../output/FDR0.05_all/PBC_vs_Healthy/No MASLD/DE_proteins_fdr0.05_pbc_vs_healthy_no.masld.xlsx")




##################### Paired liver and plasma - PBC vs Healthy and MASLD

# Combine Results Tables with Subject IDs, excluding proteins not present in both datasets
combined_results_h.m <- merge(results.p.pbc.h.m, results.t.pbc.h.m,
                               by = c("protein", "GeneName", "UniProtID", "GeneName.UniprotID"),
                               suffixes = c("_plasma", "_liver"),
                               all = FALSE)

nrow(combined_results_h.m) # 414 proteins present in both liver and plasma
combined_results_h.m[1:3,]
combined_results_h.m$FC_liver <- 2^combined_results_h.m$log2FC_liver
combined_results_h.m$FC_plasma <- 2^combined_results_h.m$log2FC_plasma
combined_results_h.m[1:3,]

# Correlation fold change Analysis with Paired Samples

correlation_df_h.m <- cor.test(combined_results_h.m$log2FC_plasma, combined_results_h.m$log2FC_liver, 
                               method = "pearson", paired = TRUE)
correlation_df_h.m

# Values for correlation coefficient, confidence interval, and p-value used for plot underneath
rho_value_h.m <- round(correlation_df_h.m$estimate,3)
ci_lower_h.m <- round(correlation_df_h.m$conf.int[1],3)
ci_upper_h.m <- round(correlation_df_h.m$conf.int[2],3)
p_value_h.m <- "< 0.001" # remember to change if another correlation anaylsis with another p-value

correlation_liver_plasma_pbc_h.m <- 
  ggplot(combined_results_h.m, aes(x = log2FC_liver, y = log2FC_plasma)) +
  geom_point() +
  labs(x = "Log2 FC Liver", y = "Log2 FC Plasma") +
  ggtitle("Correlation of Protein Regulation in PBC compared with Healthy and MASLD") +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = "blue") + 
  annotate("text", x = Inf, y = -Inf, 
           label = paste0("italic(r) == ", rho_value_h.m, " * ' (95% CI ", ci_lower_h.m, "-", ci_upper_h.m, ", ' * italic(p) * ' ", p_value_h.m, ")'"), 
           parse = TRUE, hjust = 1.1, vjust = -2)

correlation_liver_plasma_pbc_h.m

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/All/scatter_correlation_fc_liver_plasma_pbc_healthy_masld.jpeg",
#         width = 10, height = 10, plot = correlation_liver_plasma_pbc_h.m, device = "jpeg")

# only for significant proteins:

combined_results_sign_h.m <- combined_results_h.m %>% filter(adj.P.Val_plasma < 0.05 & adj.P.Val_liver < 0.05)

nrow(combined_results_sign_h.m) # 58 proteins significantly regulated in both liver and plasma

# Correlation Analysis with Paired Samples sigificant proteins
correlation_df_sign_h.m <- cor.test(combined_results_sign_h.m$log2FC_plasma, 
                                    combined_results_sign_h.m$log2FC_liver, 
                                    method = "pearson", paired = TRUE)

correlation_df_sign_h.m
combined_results_sign_h.m[1:5,]

correlation_sign_liver_plasma_pbc_h.m <- ggplot(combined_results_sign_h.m, aes(x = FC_liver, 
                                                                               y = FC_plasma)) +
  geom_point() +
  geom_text_repel(aes(label = GeneName), 
                  box.padding = 0.5, # Adjust padding around text labels
                  max.overlaps = Inf, # Allow for more attempts to avoid overlaps
                  segment.color = "black", # Line color
                  segment.size = 0.5, # Line width
                  segment.curvature = 0.2, # Line curvature
                  segment.alpha = 0.5, # Line transparency
                  hjust = 0, vjust = 0) + # Text label alignment
  labs(x = "FC Liver", y = "FC Plasma") +
  ggtitle("Significantly Differentially Expressed Proteins in Liver and Plasma in PBC compared with Healthy and MASLD") 

correlation_sign_liver_plasma_pbc_h.m

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/All/scatter_correlation_fc_sign_liver_plasma_all_proteins_named_pbc_healthy_masld.jpeg",
#         width = 10, height = 10, plot = correlation_sign_liver_plasma_pbc_h.m, device = "jpeg")


correlation_sign_liver_plasma_pbc_h.m_1 <- ggplot(combined_results_sign_h.m, aes(x = FC_liver, 
                                                                                y = FC_plasma)) +
   geom_point() +
   geom_text_repel(aes(label = GeneName)) + 
   labs(x = "FC Liver", y = "FC Plasma") +
   ggtitle("Significantly Differentially Expressed Proteins in Liver and Plasma in PBC compared with Healthy and MASLD") 
 
correlation_sign_liver_plasma_pbc_h.m_1
 
# Save plot
 
# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/All/scatter_correlation_fc_sign_liver_plasma_pbc_healthy_masld.jpeg",
#         width = 10, height = 10, plot = correlation_sign_liver_plasma_pbc_h.m_1, device = "jpeg")

# log2fc scatterplot liver/plasma with colouring according to regulation in liver and plasma
## Create a new column in your dataframe for colouring of only the dots
combined_results_sign_h.m$Regulation <- with(combined_results_sign_h.m, ifelse(
  log2FC_plasma > 0 & log2FC_liver    > 0, "Upregulated in plasma and liver",
  ifelse(log2FC_plasma > 0 & log2FC_liver < 0, "Upregulated in plasma and downregulated in liver",
         "Downregulated in plasma and liver")
))

# 2) Force factor levels into the exact legend‐order you want
combined_results_sign_h.m$Regulation <- factor(
  combined_results_sign_h.m$Regulation,
  levels = c(
    "Upregulated in plasma and liver",
    "Upregulated in plasma and downregulated in liver",
    "Downregulated in plasma and liver"
  )
)

correlation_sign_h.m_log2fc_liver_plasma_pbc_healthy_dots <- 
  ggplot(combined_results_sign_h.m, aes(x = log2FC_liver, y = log2FC_plasma)) +
  geom_point(aes(color = Regulation), size = 2) +
  geom_text_repel(aes(label = GeneName),
                  color = "black", max.overlaps = Inf) +
  scale_color_manual(
    values = c(
      "Upregulated in plasma and liver"                  = "red",
      "Upregulated in plasma and downregulated in liver" = "goldenrod",
      "Downregulated in plasma and liver"                = "blue"
    ),
    breaks = levels(combined_results_sign_h.m$Regulation)
  ) +
  labs(
    x     = "log2FC Liver",
    y     = "log2FC Plasma",
    title = "Differentially expressed proteins in liver and plasma in PBC compared with healthy and MASLD"
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 12),
    plot.title      = element_text(size = 12),
    legend.title    = element_blank(),
    legend.text     = element_text(size = 8),
    legend.position = "bottom",        # <— here
    legend.direction= "horizontal",    # <— here
    legend.key.width= unit(1.5, "lines"),
    legend.spacing.x= unit(0.5, "lines")
  )

correlation_sign_h.m_log2fc_liver_plasma_pbc_healthy_dots

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/All/scatter_correlation_log2fc_sign_h.m_liver_plasma_pbc_healthy_dots_fdr0.05.jpeg",
#         width = 8, height = 6, plot = correlation_sign_h.m_log2fc_liver_plasma_pbc_healthy_dots, device = "jpeg")

# Identify Commonly Regulated Proteins
common_upregulated_h.m <- combined_results_sign_h.m %>% filter(log2FC_plasma >= 0 & log2FC_liver >= 0)
nrow(common_upregulated_h.m) # 9
common_up_pbc_vs_h.m <- 
  common_upregulated_h.m[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                "adj.P.Val_liver")]

common_downregulated_h.m <- combined_results_sign_h.m %>% filter(log2FC_plasma < 0 & log2FC_liver < 0)
nrow(common_downregulated_h.m) # 3
common_down_pbc_vs_h.m <- 
  common_downregulated_h.m[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                  "adj.P.Val_liver")]

common_upreg.liver_downreg.plasma_h.m <- combined_results_sign_h.m %>% filter(log2FC_plasma < 0 & log2FC_liver >= 0)
nrow(common_upreg.liver_downreg.plasma_h.m) # 0
common_up.liver_down.plasma_pbc_vs_h.m <- 
  common_upreg.liver_downreg.plasma_h.m[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                               "adj.P.Val_liver")]

common_downreg.liver_upreg.plasma_h.m <- combined_results_sign_h.m %>% filter(log2FC_plasma >= 0 & log2FC_liver < 0)
nrow(common_downreg.liver_upreg.plasma_h.m) # 0

common_down.liver_up.plasma_pbc_vs_h.m <- 
  common_downreg.liver_upreg.plasma_h.m[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                               "adj.P.Val_liver")]

######################## Save in an excel file

# write.xlsx(list("downreg_plasma_pbc_vs_h_and_m" = down_plasma_pbc_vs_h.m,
#                 "upreg_plasma_pbc_vs_h_and_m" = up_plasma_pbc_vs_h.m,
#                 "downreg_liver_pbc_vs_h_and_m" = down_liver_pbc_vs_h.m,
#                 "upreg_liver_pbc_vs_h_and_m" = up_liver_pbc_vs_h.m,
#                 "up_liver_and_plasma_pbc" = common_up_pbc_vs_h.m,
#                 "down_liver_and_plasma_pbc" = common_down_pbc_vs_h.m,
#                 "up_liver_down_plasma_pbc" = common_up.liver_down.plasma_pbc_vs_h.m,
#                 "down_liver_up_plasma_pbc" = common_down.liver_up.plasma_pbc_vs_h.m),
#            "../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/All/DE_proteins_fdr0.05_pbc_vs_healthy_and_masld.xlsx")


##################### Paired liver and plasma - PBC vs Healthy and MASLD - without outliers ####################


# Combine Results Tables with Subject IDs, excluding proteins not present in both datasets
combined_results_h.m_no_outliers <- merge(results.p.pbc.h.m_no_outliers, results.t.pbc.h.m.no_outliers,
                                          by = c("protein", "GeneName", "UniProtID", "GeneName.UniprotID"),
                                          suffixes = c("_plasma", "_liver"),
                                          all = FALSE)

nrow(combined_results_h.m_no_outliers)
combined_results_h.m_no_outliers[1:3,]
combined_results_h.m_no_outliers$FC_liver <- 2^combined_results_h.m_no_outliers$log2FC_liver
combined_results_h.m_no_outliers$FC_plasma <- 2^combined_results_h.m_no_outliers$log2FC_plasma
combined_results_h.m_no_outliers[1:3,]

# Correlation fold change Analysis with Paired Samples

correlation_df_h.m_no_outliers <- cor.test(combined_results_h.m_no_outliers$log2FC_plasma,
                                           combined_results_h.m_no_outliers$log2FC_liver, 
                                           method = "pearson", paired = TRUE)
correlation_df_h.m_no_outliers

# Values for correlation coefficient, confidence interval, and p-value used for plot underneath
rho_value_h.m_no_outliers <- round(correlation_df_h.m_no_outliers$estimate,3)
ci_lower_h.m_no_outliers <- round(correlation_df_h.m_no_outliers$conf.int[1],3)
ci_upper_h.m_no_outliers <- round(correlation_df_h.m_no_outliers$conf.int[2],3)
p_value_h.m_no_outliers <- "< 0.001"

correlation_liver_plasma_pbc_h.m_no_outliers <- 
  ggplot(combined_results_h.m_no_outliers, aes(x = log2FC_liver, y = log2FC_plasma)) +
  geom_point() +
  labs(x = "Log2 FC Liver", y = "Log2 FC Plasma") +
  ggtitle("Correlation of Protein Regulation in PBC compared with Healthy and MASLD without outliers") +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = "blue") + 
  annotate("text", x = Inf, y = -Inf, 
           label = paste0("italic(r) == ", rho_value_h.m_no_outliers, " * ' (95% CI ", ci_lower_h.m_no_outliers, "-", ci_upper_h.m_no_outliers, ", ' * italic(p) * ' ", p_value_h.m_no_outliers, ")'"), 
           parse = TRUE, hjust = 1.1, vjust = -2)

correlation_liver_plasma_pbc_h.m_no_outliers

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/Without outliers/scatter_correlation_fc_liver_plasma_pbc_healthy_masld_without_outliers.jpeg",
#         width = 10, height = 10, plot = correlation_liver_plasma_pbc_h.m_no_outliers, device = "jpeg")

# only for significant proteins:

combined_results_sign_h.m_no_outliers <- combined_results_h.m_no_outliers %>% filter(adj.P.Val_plasma < 0.05 & adj.P.Val_liver < 0.05)

nrow(combined_results_sign_h.m_no_outliers)

# Correlation Analysis with Paired Samples sigificant proteins
correlation_df_sign_h.m_no_outliers <- cor.test(combined_results_sign_h.m_no_outliers$log2FC_plasma, 
                                                combined_results_sign_h.m_no_outliers$log2FC_liver, 
                                                method = "pearson", paired = TRUE)

correlation_df_sign_h.m_no_outliers
combined_results_sign_h.m_no_outliers[1:5,]

correlation_sign_liver_plasma_pbc_h.m_no_outliers <- ggplot(combined_results_sign_h.m_no_outliers, aes(x = FC_liver, 
                                                                                                       y = FC_plasma)) +
  geom_point() +
  geom_text_repel(aes(label = GeneName), 
                  box.padding = 0.5,
                  max.overlaps = Inf,
                  segment.color = "black",
                  segment.size = 0.5,
                  segment.curvature = 0.2,
                  segment.alpha = 0.5,
                  hjust = 0, vjust = 0) +
  labs(x = "FC Liver", y = "FC Plasma") +
  ggtitle("Significantly Differentially Expressed Proteins in Liver and Plasma in PBC compared with Healthy and MASLD without outliers") 

correlation_sign_liver_plasma_pbc_h.m_no_outliers

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/Without outliers/scatter_correlation_fc_sign_liver_plasma_all_proteins_named_pbc_healthy_masld_without_outliers.jpeg",
#         width = 10, height = 10, plot = correlation_sign_liver_plasma_pbc_h.m_no_outliers, device = "jpeg")

correlation_sign_liver_plasma_pbc_h.m_1_no_outliers <- ggplot(combined_results_sign_h.m_no_outliers, aes(x = FC_liver, 
                                                                                                         y = FC_plasma)) +
  geom_point() +
  geom_text_repel(aes(label = GeneName)) + 
  labs(x = "FC Liver", y = "FC Plasma") +
  ggtitle("Significantly Differentially Expressed Proteins in Liver and Plasma in PBC compared with Healthy and MASLD without outliers") 

correlation_sign_liver_plasma_pbc_h.m_1_no_outliers

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/Without outliers/scatter_correlation_fc_sign_liver_plasma_pbc_healthy_masld_without_outliers.jpeg",
#         width = 10, height = 10, plot = correlation_sign_liver_plasma_pbc_h.m_1_no_outliers, device = "jpeg")

# log2fc scatterplot liver/plasma with colouring according to regulation in liver and plasma

combined_results_sign_h.m_no_outliers$Regulation <- with(combined_results_sign_h.m_no_outliers, ifelse(
  log2FC_plasma > 0 & log2FC_liver > 0, "Upregulated in plasma and liver",
  ifelse(log2FC_plasma > 0 & log2FC_liver < 0, "Upregulated in plasma and downregulated in liver",
         "Downregulated in plasma and liver")
))

combined_results_sign_h.m_no_outliers$Regulation <- factor(
  combined_results_sign_h.m_no_outliers$Regulation,
  levels = c(
    "Upregulated in plasma and liver",
    "Upregulated in plasma and downregulated in liver",
    "Downregulated in plasma and liver"
  )
)

correlation_sign_h.m_log2fc_liver_plasma_pbc_healthy_dots_no_outliers <- 
  ggplot(combined_results_sign_h.m_no_outliers, aes(x = log2FC_liver, y = log2FC_plasma)) +
  geom_point(aes(color = Regulation), size = 2) +
  geom_text_repel(aes(label = GeneName),
                  color = "black", max.overlaps = Inf) +
  scale_color_manual(
    values = c(
      "Upregulated in plasma and liver"                  = "red",
      "Upregulated in plasma and downregulated in liver" = "goldenrod",
      "Downregulated in plasma and liver"                = "blue"
    ),
    breaks = levels(combined_results_sign_h.m_no_outliers$Regulation)
  ) +
  labs(
    x     = "log2FC Liver",
    y     = "log2FC Plasma",
    title = "Differentially expressed proteins in liver and plasma in PBC compared with healthy and MASLD without outliers"
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 12),
    plot.title      = element_text(size = 12),
    legend.title    = element_blank(),
    legend.text     = element_text(size = 8),
    legend.position = "bottom",
    legend.direction= "horizontal",
    legend.key.width= unit(1.5, "lines"),
    legend.spacing.x= unit(0.5, "lines")
  )

correlation_sign_h.m_log2fc_liver_plasma_pbc_healthy_dots_no_outliers

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/Without outliers/scatter_correlation_log2fc_sign_h.m_liver_plasma_pbc_healthy_dots_without_outliers_fdr0.05.jpeg",
#         width = 8, height = 6, plot = correlation_sign_h.m_log2fc_liver_plasma_pbc_healthy_dots_no_outliers, device = "jpeg")

# Identify Commonly Regulated Proteins

common_upregulated_h.m_no_outliers <- combined_results_sign_h.m_no_outliers %>% filter(log2FC_plasma >= 0 & log2FC_liver >= 0)
nrow(common_upregulated_h.m_no_outliers)

common_up_pbc_vs_h.m_no_outliers <- 
  common_upregulated_h.m_no_outliers[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                         "adj.P.Val_liver")]

common_downregulated_h.m_no_outliers <- combined_results_sign_h.m_no_outliers %>% filter(log2FC_plasma < 0 & log2FC_liver < 0)
nrow(common_downregulated_h.m_no_outliers)

common_down_pbc_vs_h.m_no_outliers <- 
  common_downregulated_h.m_no_outliers[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                           "adj.P.Val_liver")]

common_upreg.liver_downreg.plasma_h.m_no_outliers <- combined_results_sign_h.m_no_outliers %>% filter(log2FC_plasma < 0 & log2FC_liver >= 0)
nrow(common_upreg.liver_downreg.plasma_h.m_no_outliers)

common_up.liver_down.plasma_pbc_vs_h.m_no_outliers <- 
  common_upreg.liver_downreg.plasma_h.m_no_outliers[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                                        "adj.P.Val_liver")]

common_downreg.liver_upreg.plasma_h.m_no_outliers <- combined_results_sign_h.m_no_outliers %>% filter(log2FC_plasma >= 0 & log2FC_liver < 0)
nrow(common_downreg.liver_upreg.plasma_h.m_no_outliers)

common_down.liver_up.plasma_pbc_vs_h.m_no_outliers <- 
  common_downreg.liver_upreg.plasma_h.m_no_outliers[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                                        "adj.P.Val_liver")]


######################## Save in an excel file

# write.xlsx(list("downreg_plasma_pbc_vs_h_and_m" = down_plasma_pbc_vs_h.m_no_outliers,
#                 "upreg_plasma_pbc_vs_h_and_m" = up_plasma_pbc_vs_h.m_no_outliers,
#                 "downreg_liver_pbc_vs_h_and_m" = down_liver_pbc_vs_h.m_no_outliers,
#                 "upreg_liver_pbc_vs_h_and_m" = up_liver_pbc_vs_h.m_no_outliers,
#                 "up_liver_and_plasma_pbc" = common_up_pbc_vs_h.m_no_outliers,
#                 "down_liver_and_plasma_pbc" = common_down_pbc_vs_h.m_no_outliers,
#                 "up_liver_down_plasma_pbc" = common_up.liver_down.plasma_pbc_vs_h.m_no_outliers,
#                 "down_liver_up_plasma_pbc" = common_down.liver_up.plasma_pbc_vs_h.m_no_outliers),
#            "../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/Without outliers/DE_proteins_fdr0.05_pbc_vs_healthy_and_masld_without_outliers.xlsx")



#################### Paired liver and plasma without cirrhosis #######################

# Combine Results Tables with Subject IDs, excluding proteins not present in both datasets
combined_results_h.m_no_c <- merge(results.p.pbc.h.m.no.c, results.t.pbc.h.m.no.c,
                              by = c("protein", "GeneName", "UniProtID", "GeneName.UniprotID"),
                              suffixes = c("_plasma", "_liver"),
                              all = FALSE)

nrow(combined_results_h.m_no_c) # 414 proteins present in both liver and plasma
combined_results_h.m_no_c[1:3,]
combined_results_h.m_no_c$FC_liver <- 2^combined_results_h.m_no_c$log2FC_liver
combined_results_h.m_no_c$FC_plasma <- 2^combined_results_h.m_no_c$log2FC_plasma
combined_results_h.m_no_c[1:3,]

# Correlation fold change Analysis with Paired Samples

correlation_df_h.m_no_c <- cor.test(combined_results_h.m_no_c$log2FC_plasma, 
                                    combined_results_h.m_no_c$log2FC_liver, method = "pearson", paired = TRUE)
correlation_df_h.m_no_c

# Values for correlation coefficient, confidence interval, and p-value used for plot underneath
rho_value_h.m_no_c <- round(correlation_df_h.m_no_c$estimate,3)
ci_lower_h.m_no_c <- round(correlation_df_h.m_no_c$conf.int[1],3)
ci_upper_h.m_no_c <- round(correlation_df_h.m_no_c$conf.int[2],3)
p_value_h.m_no_c <- "< 0.001" # remember to change if another correlation anaylsis with another p-value

# Correlation plot

correlation_liver_plasma_pbc_h.m_no_c <- 
  ggplot(combined_results_h.m_no_c, aes(x = log2FC_liver, y = log2FC_plasma)) +
  geom_point() +
  labs(x = "Log2 FC Liver", y = "Log2 FC Plasma") +
  ggtitle("Correlation of Protein Regulation in non-cirrhotic PBC compared with Healthy and MASLD") +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = "blue") + 
  annotate("text", x = Inf, y = -Inf, 
           label = paste0("italic(r) == ", rho_value_h.m_no_c, " * ' (95% CI ", ci_lower_h.m_no_c, "-", ci_upper_h.m_no_c, ", ' * italic(p) * ' ", p_value_h.m_no_c, ")'"), 
           parse = TRUE, hjust = 1.1, vjust = -2)

correlation_liver_plasma_pbc_h.m_no_c

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/No cirrhosis/scatter_correlation_fc_liver_plasma_all_proteins_named_pbc_healthy_masld_no_cirrhosis.jpeg",
#         width = 10, height = 10, plot = correlation_liver_plasma_pbc_h.m_no_c, device = "jpeg")

# Only for significant proteins:

combined_results_sign_h.m_no_c <- combined_results_h.m_no_c %>% filter(adj.P.Val_plasma < 0.05 & adj.P.Val_liver < 0.05)

nrow(combined_results_sign_h.m_no_c) # 12 proteins significantly regulated in both liver and plasma

# Correlation Analysis with Paired Samples sigificant proteins
correlation_df_sign_h.m_no_c <- cor.test(combined_results_sign_h.m_no_c$log2FC_plasma, 
                                         combined_results_sign_h.m_no_c$log2FC_liver, 
                                         method = "pearson", paired = TRUE)
correlation_df_sign_h.m_no_c
combined_results_sign_h.m_no_c[1:5,]

correlation_sign_liver_plasma_pbc_h.m_no_c <- ggplot(combined_results_sign_h.m_no_c, aes(x = FC_liver, 
                                                                                         y = FC_plasma)) +
  geom_point() +
  geom_text_repel(aes(label = GeneName), 
                  box.padding = 0.5, # Adjust padding around text labels
                  max.overlaps = Inf, # Allow for more attempts to avoid overlaps
                  segment.color = "black", # Line color
                  segment.size = 0.5, # Line width
                  segment.curvature = 0.2, # Line curvature
                  segment.alpha = 0.5, # Line transparency
                  hjust = 0, vjust = 0) + # Text label alignment
  labs(x = "FC Liver", y = "FC Plasma") +
  ggtitle("Significantly Differentially Expressed Proteins in Liver and Plasma in non-cirrhotic PBC compared with Healthy and MASLD") 

correlation_sign_liver_plasma_pbc_h.m_no_c

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/No cirrhosis/scatter_correlation_fc_sign_liver_plasma_all_proteins_named_pbc_healthy_masld_no_cirrhosis.jpeg",
#         width = 10, height = 10, plot = correlation_sign_liver_plasma_pbc_h.m_no_c, device = "jpeg")

correlation_sign_liver_plasma_pbc_h.m_no_c_1 <- ggplot(combined_results_sign_h.m_no_c, aes(x = FC_liver, 
                                                                                            y = FC_plasma)) +
   geom_point() +
   geom_text_repel(aes(label = GeneName)) + 
   labs(x = "FC Liver", y = "FC Plasma") +
   ggtitle("Significantly Differentially Expressed Proteins in Liver and Plasma in non-cirrhotic PBC compared with Healthy and MASLD") 
 
correlation_sign_liver_plasma_pbc_h.m_no_c_1
 
# Save plot
 
# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/No cirrhosis/scatter_correlation_fc_sign_liver_plasma_pbc_healthy_masld_no_cirrhosis.jpeg",
#         width = 10, height = 10, plot = correlation_sign_liver_plasma_pbc_h.m_no_c_1, device = "jpeg")

# log2fc scatterplot liver/plasma with colouring according to regulation in liver and plasma
## Create a new column in your dataframe for colouring of only the dots
combined_results_sign_h.m_no_c$Regulation <- with(combined_results_sign_h.m_no_c, ifelse(
  log2FC_plasma > 0 & log2FC_liver    > 0, "Upregulated in plasma and liver",
  ifelse(log2FC_plasma > 0 & log2FC_liver < 0, "Upregulated in plasma and downregulated in liver",
         "Downregulated in plasma and liver")
))

# 2) Force factor levels into the exact legend‐order you want
combined_results_sign_h.m_no_c$Regulation <- factor(
  combined_results_sign_h.m_no_c$Regulation,
  levels = c(
    "Upregulated in plasma and liver",
    "Upregulated in plasma and downregulated in liver",
    "Downregulated in plasma and liver"
  )
)

correlation_sign_h.m_no_c_log2fc_liver_plasma_pbc_healthy_dots <- 
  ggplot(combined_results_sign_h.m_no_c, aes(x = log2FC_liver, y = log2FC_plasma)) +
  geom_point(aes(color = Regulation), size = 2) +
  geom_text_repel(aes(label = GeneName),
                  color = "black", max.overlaps = Inf) +
  scale_color_manual(
    values = c(
      "Upregulated in plasma and liver"                  = "red",
      "Upregulated in plasma and downregulated in liver" = "goldenrod",
      "Downregulated in plasma and liver"                = "blue"
    ),
    breaks = levels(combined_results_sign_h.m_no_c$Regulation)
  ) +
  labs(
    x     = "log2FC Liver",
    y     = "log2FC Plasma",
    title = "Differentially expressed proteins in liver and plasma in non-cirrhotic PBC compared with healthy and MASLD"
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 12),
    plot.title      = element_text(size = 12),
    legend.title    = element_blank(),
    legend.text     = element_text(size = 8),
    legend.position = "bottom",        # <— here
    legend.direction= "horizontal",    # <— here
    legend.key.width= unit(1.5, "lines"),
    legend.spacing.x= unit(0.5, "lines")
  )

correlation_sign_h.m_no_c_log2fc_liver_plasma_pbc_healthy_dots

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/No cirrhosis/scatter_correlation_log2fc_sign_h.m_no_c_liver_plasma_pbc_healthy_dots_fdr0.05.jpeg",
#         width = 8, height = 6, plot = correlation_sign_h.m_no_c_log2fc_liver_plasma_pbc_healthy_dots, device = "jpeg")

# Identify Commonly Regulated Proteins
common_upregulated_h.m_no_c <- combined_results_sign_h.m_no_c %>% filter(log2FC_plasma >= 0 & log2FC_liver >= 0)
nrow(common_upregulated_h.m_no_c) # 5
common_up_pbc_vs_h.m_no_c <- 
  common_upregulated_h.m_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                             "adj.P.Val_liver")]

common_downregulated_h.m_no_c <- combined_results_sign_h.m_no_c %>% filter(log2FC_plasma < 0 & log2FC_liver < 0)
nrow(common_downregulated_h.m_no_c) # 6
common_down_pbc_vs_h.m_no_c <- 
  common_downregulated_h.m_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                               "adj.P.Val_liver")]

common_upreg.liver_downreg.plasma_h.m_no_c <- combined_results_sign_h.m_no_c %>% filter(log2FC_plasma < 0 & log2FC_liver >= 0)
nrow(common_upreg.liver_downreg.plasma_h.m_no_c) # 0
common_up.liver_down.plasma_pbc_vs_h.m_no_c <- 
  common_upreg.liver_downreg.plasma_h.m_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID","log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                            "adj.P.Val_liver")]

common_downreg.liver_upreg.plasma_h.m_no_c <- combined_results_sign_h.m_no_c %>% filter(log2FC_plasma >= 0 & log2FC_liver < 0)
nrow(common_downreg.liver_upreg.plasma_h.m_no_c) # 0

common_down.liver_up.plasma_pbc_vs_h.m_no_c <- 
  common_downreg.liver_upreg.plasma_h.m_no_c[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                            "adj.P.Val_liver")]


######################## Save in an excel file

# write.xlsx(list("downreg_plasma_pbc_vs_h_and_m" = down_plasma_pbc_vs_h.m_no_c,
#                 "upreg_plasma_pbc_vs_h_and_m" = up_plasma_pbc_vs_h.m_no_c,
#                 "downreg_liver_pbc_vs_h_and_m" = down_liver_pbc_vs_h.m_no_c,
#                 "upreg_liver_pbc_vs_h_and_m" = up_liver_pbc_vs_h.m_no_c,
#                 "up_liver_and_plasma_pbc" = common_up_pbc_vs_h.m_no_c,
#                 "down_liver_and_plasma_pbc" = common_down_pbc_vs_h.m_no_c,
#                 "up_liver_down_plasma_pbc" = common_up.liver_down.plasma_pbc_vs_h.m_no_c,
#                 "down_liver_up_plasma_pbc" = common_down.liver_up.plasma_pbc_vs_h.m_no_c),
#            "../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/No cirrhosis/DE_proteins_fdr0.05_pbc_vs_healthy_and_masld_no_cir.xlsx")


##################### paired liver and plasma PBC without MASLD vs. Healthy and MASLD ################ 



# Combine Results Tables with Subject IDs, excluding proteins not present in both datasets
combined_results_h.m_no_masld <- merge(results.p.pbc.h.m.no.masld, results.t.pbc.h.m.no.masld,
                                   by = c("protein", "GeneName", "UniProtID", "GeneName.UniprotID"),
                                   suffixes = c("_plasma", "_liver"),
                                   all = FALSE)

nrow(combined_results_h.m_no_masld) # 414 proteins present in both liver and plasma
combined_results_h.m_no_masld[1:3,]
combined_results_h.m_no_masld$FC_liver <- 2^combined_results_h.m_no_masld$log2FC_liver
combined_results_h.m_no_masld$FC_plasma <- 2^combined_results_h.m_no_masld$log2FC_plasma
combined_results_h.m_no_masld[1:3,]

# Correlation fold change Analysis with Paired Samples

correlation_df_h.m_no_masld <- cor.test(combined_results_h.m_no_masld$log2FC_plasma, 
                                    combined_results_h.m_no_masld$log2FC_liver, method = "pearson", paired = TRUE)
correlation_df_h.m_no_masld

# Values for correlation coefficient, confidence interval, and p-value used for plot underneath
rho_value_h.m_no_masld <- round(correlation_df_h.m_no_masld$estimate,3)
ci_lower_h.m_no_masld <- round(correlation_df_h.m_no_masld$conf.int[1],3)
ci_upper_h.m_no_masld <- round(correlation_df_h.m_no_masld$conf.int[2],3)
p_value_h.m_no_masld <- "< 0.001" # remember to change if another correlation anaylsis with another p-value

# Correlation plot

correlation_liver_plasma_pbc_h.m_no_masld <- 
  ggplot(combined_results_h.m_no_masld, aes(x = log2FC_liver, y = log2FC_plasma)) +
  geom_point() +
  labs(x = "Log2 FC Liver", y = "Log2 FC Plasma") +
  ggtitle("Correlation of Protein Regulation in PBC without MASLD compared with Healthy and MASLD") +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = "blue") + 
  annotate("text", x = Inf, y = -Inf, 
           label = paste0("italic(r) == ", rho_value_h.m_no_masld, " * ' (95% CI ", ci_lower_h.m_no_masld, "-", ci_upper_h.m_no_masld, ", ' * italic(p) * ' ", p_value_h.m_no_masld, ")'"), 
           parse = TRUE, hjust = 1.1, vjust = -2)

correlation_liver_plasma_pbc_h.m_no_masld

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/No MASLD/scatter_correlation_fc_liver_plasma_pbc_h.m_no_masld.jpeg",
#         width = 10, height = 10, plot = correlation_liver_plasma_pbc_healthy_h.m_no_masld, device = "jpeg")

# only for significant proteins:

combined_results_sign_h.m_no_masld <- combined_results_h.m_no_masld %>% filter(adj.P.Val_plasma < 0.05 & adj.P.Val_liver < 0.05)

nrow(combined_results_sign_h.m_no_masld) # ONLY 1 proteins significantly regulated in both liver and plasma

# Correlation Analysis with Paired Samples sigificant proteins
# correlation_df_sign_h.m_no_masld <- cor.test(combined_results_sign_h.m_no_masld$log2FC_plasma, combined_results_sign_h.m_no_masld$log2FC_liver, 
#                                         method = "pearson", paired = TRUE)
# correlation_df_sign_h.m_no_masld
combined_results_sign_h.m_no_masld[1:5,]

correlation_sign_liver_plasma_pbc_h.m_no_masld <- ggplot(combined_results_sign_h.m_no_masld, aes(x = FC_liver, 
                                                                                                 y = FC_plasma)) +
  geom_point() +
  geom_text_repel(aes(label = GeneName), 
                  box.padding = 0.5, # Adjust padding around text labels
                  max.overlaps = Inf, # Allow for more attempts to avoid overlaps
                  segment.color = "black", # Line color
                  segment.size = 0.5, # Line width
                  segment.curvature = 0.2, # Line curvature
                  segment.alpha = 0.5, # Line transparency
                  hjust = 0, vjust = 0) + # Text label alignment
  labs(x = "FC Liver", y = "FC Plasma") +
  ggtitle("Differentially expressed proteins in liver and plasma in PBC without MASLD compared with Healthy and MASLD") 

correlation_sign_liver_plasma_pbc_h.m_no_masld


# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/No MASLD/scatter_correlation_fc_sign_liver_plasma_all_proteins_named_pbc_h.m_no_masld.jpeg",
#         width = 10, height = 10, plot = correlation_sign_liver_plasma_pbc_h.m_no_masld, device = "jpeg")


correlation_sign_liver_plasma_pbc_h.m_no_masld_1 <- ggplot(combined_results_sign_h.m_no_masld, aes(x = FC_liver, 
                                                                                                   y = FC_plasma)) +
  geom_point() +
  geom_text_repel(aes(label = GeneName)) + # Text label alignment
  labs(x = "FC Liver", y = "FC Plasma") +
  ggtitle("Differentially expressed proteins in liver and plasma in PBC without MASLD compared with Healthy and MASLD") +
  theme(
    axis.title = element_text(size = 16),   # Increase axis title size
    axis.text = element_text(size = 14),    # Increase axis text (tick labels) size
    plot.title = element_text(size = 18)    # Increase plot title size
  ) 

correlation_sign_liver_plasma_pbc_h.m_no_masld_1

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/No MASLD/scatter_correlation_fc_sign_liver_plasma_pbc_h.m_no_masld.fdr0.05.jpeg",
#         width = 8, height = 6, plot = correlation_sign_liver_plasma_pbc_h.m_no_masld_1, device = "jpeg")

# log2fc scatterplot liver/plasma with colouring according to regulation in liver and plasma
## Create a new column in your dataframe for colouring of only the dots
combined_results_sign_h.m_no_masld$Regulation <- with(combined_results_sign_h.m_no_masld, ifelse(
  log2FC_plasma > 0 & log2FC_liver    > 0, "Upregulated in plasma and liver",
  ifelse(log2FC_plasma > 0 & log2FC_liver < 0, "Upregulated in plasma and downregulated in liver",
         "Downregulated in plasma and liver")
))

# 2) Force factor levels into the exact legend‐order you want
combined_results_sign_h.m_no_masld$Regulation <- factor(
  combined_results_sign_h.m_no_masld$Regulation,
  levels = c(
    "Upregulated in plasma and liver",
    "Upregulated in plasma and downregulated in liver",
    "Downregulated in plasma and liver"
  )
)

correlation_sign_log2fc_liver_plasma_pbc_h.m_no_masld_dots <- 
  ggplot(combined_results_sign_h.m_no_masld, aes(x = log2FC_liver, y = log2FC_plasma)) +
  geom_point(aes(color = Regulation), size = 2) +
  geom_text_repel(aes(label = GeneName),
                  color = "black", max.overlaps = Inf) +
  scale_color_manual(
    values = c(
      "Upregulated in plasma and liver"                  = "red",
      "Upregulated in plasma and downregulated in liver" = "goldenrod",
      "Downregulated in plasma and liver"                = "blue"
    ),
    breaks = levels(combined_results_sign_h.m_no_masld$Regulation)
  ) +
  labs(
    x     = "log2FC Liver",
    y     = "log2FC Plasma",
    title = "Differentially expressed proteins in liver and plasma in PBC without MASLD compared with healthy"
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 12),
    plot.title      = element_text(size = 12),
    legend.title    = element_blank(),
    legend.text     = element_text(size = 8),
    legend.position = "bottom",        # <— here
    legend.direction= "horizontal",    # <— here
    legend.key.width= unit(1.5, "lines"),
    legend.spacing.x= unit(0.5, "lines")
  )

correlation_sign_log2fc_liver_plasma_pbc_h.m_no_masld_dots

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/No MASLD/correlation_sign_log2fc_liver_plasma_pbc_h.m_no_masld_dots_fdr0.05.jpeg",
#         width = 8, height = 6, plot = correlation_sign_log2fc_liver_plasma_pbc_h.m_no_masld_dots, device = "jpeg")


# Identify Commonly Regulated Proteins
common_upregulated_h.m_no_masld <- combined_results_sign_h.m_no_masld %>% filter(log2FC_plasma >= 0 & log2FC_liver >= 0)
nrow(common_upregulated_h.m_no_masld) # 1

common_up_pbc_vs_h.m_no_masld <- 
  common_upregulated_h.m_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                  "adj.P.Val_liver")]

common_downregulated_h.m_no_masld <- combined_results_sign_h.m_no_masld %>% filter(log2FC_plasma < 0 & log2FC_liver < 0)
nrow(common_downregulated_h.m_no_masld) # 9

common_down_pbc_vs_h.m_no_masld <- 
  common_downregulated_h.m_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                    "adj.P.Val_liver")]

common_upreg.liver_downreg.plasma_h.m_no_masld <- combined_results_sign_h.m_no_masld %>% filter(log2FC_plasma < 0 & log2FC_liver >= 0)
nrow(common_upreg.liver_downreg.plasma_h.m_no_masld ) # 1

common_up.liver_down.plasma_pbc_vs_h.m_no_masld <- 
  common_upreg.liver_downreg.plasma_h.m_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                             "adj.P.Val_liver")]

common_downreg.liver_upreg.plasma_h.m_no_masld <- combined_results_sign_h.m_no_masld %>% filter(log2FC_plasma >= 0 & log2FC_liver < 0)
nrow(common_downreg.liver_upreg.plasma_h.m_no_masld) # 5

common_down.liver_up.plasma_pbc_vs_h.m_no_masld <- 
  common_downreg.liver_upreg.plasma_h.m_no_masld[, c("GeneName.UniprotID","GeneName", "UniProtID", "log2FC_plasma", "adj.P.Val_plasma","log2FC_liver",
                                             "adj.P.Val_liver")]

# write.xlsx(list("down_plasma_pbc_no_masld" = down_plasma_pbc_vs_h.m_no_masld,
#                 "up_plasma_pbc_no_masld" = up_plasma_pbc_vs_h.m_no_masld,
#                 "down_liver_pbc_no_masld" = down_liver_pbc_vs_h.m_no_masld,
#                 "up_liver_pbc_no_masld" = up_liver_pbc_vs_h.m_no_masld,
#                 "up_both_pbc_no_masld" = common_up_pbc_vs_h.m_no_masld,
#                 "down_both_pbc_no_masld" = common_down_pbc_vs_h.m_no_masld,
#                 "up_l_down_p_pbc_no_masld" = common_up.liver_down.plasma_pbc_vs_h.m_no_masld,
#                 "down_l_up_p_pbc_no_masld" = common_down.liver_up.plasma_pbc_vs_h.m_no_masld),
#            "../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/No MASLD/DE_proteins_fdr0.05_pbc_vs_h.m_no.masld.xlsx")

### PBC vs. Healthy

# Downregulated plasma
# Get a list of unique proteins in down_plasma_pbc_vs_healthy not present in down_plasma_pbc_vs_healthy_no_c
# Get the unique proteins from down_plasma_pbc_vs_healthy that are not present in down_plasma_pbc_vs_healthy_no_c
unique_down_p_pbc <- setdiff(down_plasma_pbc_vs_healthy$GeneName.UniprotID, 
                             down_plasma_pbc_vs_healthy_no_c$GeneName.UniprotID)

# Retrieve rows from down_plasma_pbc_vs_healthy (with unique proteins) keeping log2FC and adj.P.Val
filtered_down_p_pbc <- down_plasma_pbc_vs_healthy %>%
  filter(GeneName.UniprotID %in% unique_down_p_pbc) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

# Retrieve rows from results.p.pbc.no.c based on unique protein list and keep log2FC and adj.P.Val
filtered_down_p_pbc_no_c <- results.p.pbc.no.c %>%
  filter(GeneName.UniprotID %in% unique_down_p_pbc) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

# Merge the two datasets (from step 2 and step 3) based on GeneName.UniprotID
# Rename columns to indicate the source of log2FC and adj.P.Val (down_plasma and results.p.pbc.no.c)
merged_results_down_p_pbc <- filtered_down_p_pbc %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(filtered_down_p_pbc_no_c %>%
               rename(log2FC_no_cir = log2FC, adj.P.Val_no_cir = adj.P.Val),
             by = c("GeneName.UniprotID", "GeneName", "UniProtID"))

# Upregulated plasma
unique_up_p_pbc <- setdiff(up_plasma_pbc_vs_healthy$GeneName.UniprotID, 
                             up_plasma_pbc_vs_healthy_no_c$GeneName.UniprotID)
filtered_up_p_pbc <- up_plasma_pbc_vs_healthy %>%
  filter(GeneName.UniprotID %in% unique_up_p_pbc) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)
filtered_up_p_pbc_no_c <- results.p.pbc.no.c %>%
  filter(GeneName.UniprotID %in% unique_up_p_pbc) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)
merged_results_up_p_pbc <- filtered_up_p_pbc %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(filtered_up_p_pbc_no_c %>%
               rename(log2FC_no_cir = log2FC, adj.P.Val_no_cir = adj.P.Val),
             by = c("GeneName.UniprotID", "GeneName", "UniProtID"))

# Downregulated liver
unique_down_t_pbc <- setdiff(down_liver_pbc_vs_healthy$GeneName.UniprotID, 
                           down_liver_pbc_vs_healthy_no_c$GeneName.UniprotID)
filtered_down_t_pbc <- down_liver_pbc_vs_healthy %>%
  filter(GeneName.UniprotID %in% unique_down_t_pbc) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)
filtered_down_t_pbc_no_c <- results.t.pbc.no.c %>%
  filter(GeneName.UniprotID %in% unique_down_t_pbc) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)
merged_results_down_t_pbc <- filtered_down_t_pbc %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(filtered_down_t_pbc_no_c %>%
               rename(log2FC_no_cir = log2FC, adj.P.Val_no_cir = adj.P.Val),
             by = c("GeneName.UniprotID", "GeneName", "UniProtID"))

# Upregulated liver
unique_up_t_pbc <- setdiff(up_liver_pbc_vs_healthy$GeneName.UniprotID, 
                             up_liver_pbc_vs_healthy_no_c$GeneName.UniprotID)
filtered_up_t_pbc <- up_liver_pbc_vs_healthy %>%
  filter(GeneName.UniprotID %in% unique_up_t_pbc) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)
filtered_up_t_pbc_no_c <- results.t.pbc.no.c %>%
  filter(GeneName.UniprotID %in% unique_up_t_pbc) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)
merged_results_up_t_pbc <- filtered_up_t_pbc %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(filtered_up_t_pbc_no_c %>%
               rename(log2FC_no_cir = log2FC, adj.P.Val_no_cir = adj.P.Val),
             by = c("GeneName.UniprotID", "GeneName", "UniProtID"))

# Save in an excel file

# write.xlsx(list("down_plasma" = merged_results_down_p_pbc,
#                 "up_plasma" = merged_results_up_p_pbc,
#                 "down_liver" = merged_results_down_t_pbc,
#                 "up_liver" = merged_results_up_t_pbc),
#            "../output/FDR0.05_all/PBC_vs_Healthy/proteins_no_longer_diff_exp_when_excluding_cirrhosis_pbc_vs_healthy.xlsx")


################ Create Excel file with proteins affected by excluding outliers

### PBC vs. Healthy

#### Proteins no longer significant after excluding outliers ####

# Downregulated plasma
lost_down_p_pbc_outliers <- setdiff(
  down_plasma_pbc_vs_healthy$GeneName.UniprotID,
  down_plasma_pbc_vs_healthy_no_outliers$GeneName.UniprotID
)

filtered_lost_down_p_pbc_total <- down_plasma_pbc_vs_healthy %>%
  filter(GeneName.UniprotID %in% lost_down_p_pbc_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

filtered_lost_down_p_pbc_no_outliers <- results.p.pbc_no_outliers %>%
  filter(GeneName.UniprotID %in% lost_down_p_pbc_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

merged_lost_down_p_pbc_outliers <- filtered_lost_down_p_pbc_total %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(
    filtered_lost_down_p_pbc_no_outliers %>%
      rename(log2FC_without_outliers = log2FC,
             adj.P.Val_without_outliers = adj.P.Val),
    by = c("GeneName.UniprotID", "GeneName", "UniProtID")
  )


# Upregulated plasma
lost_up_p_pbc_outliers <- setdiff(
  up_plasma_pbc_vs_healthy$GeneName.UniprotID,
  up_plasma_pbc_vs_healthy_no_outliers$GeneName.UniprotID
)

filtered_lost_up_p_pbc_total <- up_plasma_pbc_vs_healthy %>%
  filter(GeneName.UniprotID %in% lost_up_p_pbc_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

filtered_lost_up_p_pbc_no_outliers <- results.p.pbc_no_outliers %>%
  filter(GeneName.UniprotID %in% lost_up_p_pbc_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

merged_lost_up_p_pbc_outliers <- filtered_lost_up_p_pbc_total %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(
    filtered_lost_up_p_pbc_no_outliers %>%
      rename(log2FC_without_outliers = log2FC,
             adj.P.Val_without_outliers = adj.P.Val),
    by = c("GeneName.UniprotID", "GeneName", "UniProtID")
  )


# Downregulated liver
lost_down_t_pbc_outliers <- setdiff(
  down_liver_pbc_vs_healthy$GeneName.UniprotID,
  down_liver_pbc_vs_healthy_no_outliers$GeneName.UniprotID
)

filtered_lost_down_t_pbc_total <- down_liver_pbc_vs_healthy %>%
  filter(GeneName.UniprotID %in% lost_down_t_pbc_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

filtered_lost_down_t_pbc_no_outliers <- results.t.pbc.no_outliers %>%
  filter(GeneName.UniprotID %in% lost_down_t_pbc_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

merged_lost_down_t_pbc_outliers <- filtered_lost_down_t_pbc_total %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(
    filtered_lost_down_t_pbc_no_outliers %>%
      rename(log2FC_without_outliers = log2FC,
             adj.P.Val_without_outliers = adj.P.Val),
    by = c("GeneName.UniprotID", "GeneName", "UniProtID")
  )


# Upregulated liver
lost_up_t_pbc_outliers <- setdiff(
  up_liver_pbc_vs_healthy$GeneName.UniprotID,
  up_liver_pbc_vs_healthy_no_outliers$GeneName.UniprotID
)

filtered_lost_up_t_pbc_total <- up_liver_pbc_vs_healthy %>%
  filter(GeneName.UniprotID %in% lost_up_t_pbc_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

filtered_lost_up_t_pbc_no_outliers <- results.t.pbc.no_outliers %>%
  filter(GeneName.UniprotID %in% lost_up_t_pbc_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

merged_lost_up_t_pbc_outliers <- filtered_lost_up_t_pbc_total %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(
    filtered_lost_up_t_pbc_no_outliers %>%
      rename(log2FC_without_outliers = log2FC,
             adj.P.Val_without_outliers = adj.P.Val),
    by = c("GeneName.UniprotID", "GeneName", "UniProtID")
  )


#### Proteins significant after excluding outliers ####

# Downregulated plasma
new_down_p_pbc_outliers <- setdiff(
  down_plasma_pbc_vs_healthy_no_outliers$GeneName.UniprotID,
  down_plasma_pbc_vs_healthy$GeneName.UniprotID
)

filtered_new_down_p_pbc_no_outliers <- down_plasma_pbc_vs_healthy_no_outliers %>%
  filter(GeneName.UniprotID %in% new_down_p_pbc_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

filtered_new_down_p_pbc_total <- results.p.pbc %>%
  filter(GeneName.UniprotID %in% new_down_p_pbc_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

merged_new_down_p_pbc_outliers <- filtered_new_down_p_pbc_total %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(
    filtered_new_down_p_pbc_no_outliers %>%
      rename(log2FC_without_outliers = log2FC,
             adj.P.Val_without_outliers = adj.P.Val),
    by = c("GeneName.UniprotID", "GeneName", "UniProtID")
  )


# Upregulated plasma
new_up_p_pbc_outliers <- setdiff(
  up_plasma_pbc_vs_healthy_no_outliers$GeneName.UniprotID,
  up_plasma_pbc_vs_healthy$GeneName.UniprotID
)

filtered_new_up_p_pbc_no_outliers <- up_plasma_pbc_vs_healthy_no_outliers %>%
  filter(GeneName.UniprotID %in% new_up_p_pbc_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

filtered_new_up_p_pbc_total <- results.p.pbc %>%
  filter(GeneName.UniprotID %in% new_up_p_pbc_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

merged_new_up_p_pbc_outliers <- filtered_new_up_p_pbc_total %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(
    filtered_new_up_p_pbc_no_outliers %>%
      rename(log2FC_without_outliers = log2FC,
             adj.P.Val_without_outliers = adj.P.Val),
    by = c("GeneName.UniprotID", "GeneName", "UniProtID")
  )


# Downregulated liver
new_down_t_pbc_outliers <- setdiff(
  down_liver_pbc_vs_healthy_no_outliers$GeneName.UniprotID,
  down_liver_pbc_vs_healthy$GeneName.UniprotID
)

filtered_new_down_t_pbc_no_outliers <- down_liver_pbc_vs_healthy_no_outliers %>%
  filter(GeneName.UniprotID %in% new_down_t_pbc_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

filtered_new_down_t_pbc_total <- results.t.pbc %>%
  filter(GeneName.UniprotID %in% new_down_t_pbc_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

merged_new_down_t_pbc_outliers <- filtered_new_down_t_pbc_total %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(
    filtered_new_down_t_pbc_no_outliers %>%
      rename(log2FC_without_outliers = log2FC,
             adj.P.Val_without_outliers = adj.P.Val),
    by = c("GeneName.UniprotID", "GeneName", "UniProtID")
  )


# Upregulated liver
new_up_t_pbc_outliers <- setdiff(
  up_liver_pbc_vs_healthy_no_outliers$GeneName.UniprotID,
  up_liver_pbc_vs_healthy$GeneName.UniprotID
)

filtered_new_up_t_pbc_no_outliers <- up_liver_pbc_vs_healthy_no_outliers %>%
  filter(GeneName.UniprotID %in% new_up_t_pbc_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

filtered_new_up_t_pbc_total <- results.t.pbc %>%
  filter(GeneName.UniprotID %in% new_up_t_pbc_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

merged_new_up_t_pbc_outliers <- filtered_new_up_t_pbc_total %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(
    filtered_new_up_t_pbc_no_outliers %>%
      rename(log2FC_without_outliers = log2FC,
             adj.P.Val_without_outliers = adj.P.Val),
    by = c("GeneName.UniprotID", "GeneName", "UniProtID")
  )


# write.xlsx(list("lost_down_plasma" = merged_lost_down_p_pbc_outliers,
#                 "lost_up_plasma" = merged_lost_up_p_pbc_outliers,
#                 "lost_down_liver" = merged_lost_down_t_pbc_outliers,
#                 "lost_up_liver" = merged_lost_up_t_pbc_outliers,
#                 "new_down_plasma" = merged_new_down_p_pbc_outliers,
#                 "new_up_plasma" = merged_new_up_p_pbc_outliers,
#                 "new_down_liver" = merged_new_down_t_pbc_outliers,
#                 "new_up_liver" = merged_new_up_t_pbc_outliers),
#            "../output/FDR0.05_all/PBC_vs_Healthy/proteins_changed_significance_when_excluding_outliers_pbc_vs_healthy.xlsx")



# PBC vs. Healthy/MASLD

# Downregulated plasma
unique_down_p_pbc_h.m <- setdiff(down_plasma_pbc_vs_h.m$GeneName.UniprotID, 
                             down_plasma_pbc_vs_h.m_no_c$GeneName.UniprotID)
filtered_down_p_pbc_h.m <- down_plasma_pbc_vs_h.m %>%
  filter(GeneName.UniprotID %in% unique_down_p_pbc_h.m) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)
filtered_down_p_pbc_h.m_no_c <- results.p.pbc.h.m.no.c %>%
  filter(GeneName.UniprotID %in% unique_down_p_pbc_h.m) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)
merged_results_down_p_pbc_h.m <- filtered_down_p_pbc_h.m %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(filtered_down_p_pbc_h.m_no_c %>%
               rename(log2FC_no_cir = log2FC, adj.P.Val_no_cir = adj.P.Val),
             by = c("GeneName.UniprotID", "GeneName", "UniProtID"))

# Upregulated plasma
unique_up_p_pbc_h.m <- setdiff(up_plasma_pbc_vs_h.m$GeneName.UniprotID, 
                           up_plasma_pbc_vs_h.m_no_c$GeneName.UniprotID)
filtered_up_p_pbc_h.m <- up_plasma_pbc_vs_h.m %>%
  filter(GeneName.UniprotID %in% unique_up_p_pbc_h.m) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)
filtered_up_p_pbc_h.m_no_c <- results.p.pbc.h.m.no.c %>%
  filter(GeneName.UniprotID %in% unique_up_p_pbc_h.m) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)
merged_results_up_p_pbc_h.m <- filtered_up_p_pbc_h.m %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(filtered_up_p_pbc_h.m_no_c %>%
               rename(log2FC_no_cir = log2FC, adj.P.Val_no_cir = adj.P.Val),
             by = c("GeneName.UniprotID", "GeneName", "UniProtID"))

# Downregulated liver
unique_down_t_pbc_h.m <- setdiff(down_liver_pbc_vs_h.m$GeneName.UniprotID, 
                             down_liver_pbc_vs_h.m_no_c$GeneName.UniprotID)
filtered_down_t_pbc_h.m <- down_liver_pbc_vs_h.m %>%
  filter(GeneName.UniprotID %in% unique_down_t_pbc_h.m) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)
filtered_down_t_pbc_h.m_no_c <- results.t.pbc.h.m.no.c %>%
  filter(GeneName.UniprotID %in% unique_down_t_pbc_h.m) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)
merged_results_down_t_pbc_h.m <- filtered_down_t_pbc_h.m %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(filtered_down_t_pbc_h.m_no_c %>%
               rename(log2FC_no_cir = log2FC, adj.P.Val_no_cir = adj.P.Val),
             by = c("GeneName.UniprotID", "GeneName", "UniProtID"))

# Upregulated liver
unique_up_t_pbc_h.m <- setdiff(up_liver_pbc_vs_h.m$GeneName.UniprotID, 
                           up_liver_pbc_vs_h.m_no_c$GeneName.UniprotID)
filtered_up_t_pbc_h.m <- up_liver_pbc_vs_h.m %>%
  filter(GeneName.UniprotID %in% unique_up_t_pbc_h.m) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)
filtered_up_t_pbc_h.m_no_c <- results.t.pbc.h.m.no.c %>%
  filter(GeneName.UniprotID %in% unique_up_t_pbc_h.m) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)
merged_results_up_t_pbc_h.m <- filtered_up_t_pbc_h.m %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(filtered_up_t_pbc_h.m_no_c %>%
               rename(log2FC_no_cir = log2FC, adj.P.Val_no_cir = adj.P.Val),
             by = c("GeneName.UniprotID", "GeneName", "UniProtID"))

# Save in an excel file

# write.xlsx(list("down_plasma" = merged_results_down_p_pbc_h.m,
#                 "up_plasma" = merged_results_up_p_pbc_h.m,
#                 "down_liver" = merged_results_down_t_pbc_h.m,
#                 "up_liver" = merged_results_up_t_pbc_h.m),
#            "../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/proteins_no_longer_diff_exp_when_excluding_cirrhosis_pbc_vs_healthy_and_masld.xlsx")



################ Create Excel file with proteins affected by excluding outliers

### PBC vs. Healthy/MASLD


#### Proteins no longer significant after excluding outliers ####

# Downregulated plasma
lost_down_p_pbc_h.m_outliers <- setdiff(
  down_plasma_pbc_vs_h.m$GeneName.UniprotID,
  down_plasma_pbc_vs_h.m_no_outliers$GeneName.UniprotID
)

filtered_lost_down_p_pbc_h.m_total <- down_plasma_pbc_vs_h.m %>%
  filter(GeneName.UniprotID %in% lost_down_p_pbc_h.m_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

filtered_lost_down_p_pbc_h.m_no_outliers <- results.p.pbc.h.m_no_outliers %>%
  filter(GeneName.UniprotID %in% lost_down_p_pbc_h.m_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

merged_lost_down_p_pbc_h.m_outliers <- filtered_lost_down_p_pbc_h.m_total %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(
    filtered_lost_down_p_pbc_h.m_no_outliers %>%
      rename(log2FC_without_outliers = log2FC,
             adj.P.Val_without_outliers = adj.P.Val),
    by = c("GeneName.UniprotID", "GeneName", "UniProtID")
  )


# Upregulated plasma
lost_up_p_pbc_h.m_outliers <- setdiff(
  up_plasma_pbc_vs_h.m$GeneName.UniprotID,
  up_plasma_pbc_vs_h.m_no_outliers$GeneName.UniprotID
)

filtered_lost_up_p_pbc_h.m_total <- up_plasma_pbc_vs_h.m %>%
  filter(GeneName.UniprotID %in% lost_up_p_pbc_h.m_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

filtered_lost_up_p_pbc_h.m_no_outliers <- results.p.pbc.h.m_no_outliers %>%
  filter(GeneName.UniprotID %in% lost_up_p_pbc_h.m_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

merged_lost_up_p_pbc_h.m_outliers <- filtered_lost_up_p_pbc_h.m_total %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(
    filtered_lost_up_p_pbc_h.m_no_outliers %>%
      rename(log2FC_without_outliers = log2FC,
             adj.P.Val_without_outliers = adj.P.Val),
    by = c("GeneName.UniprotID", "GeneName", "UniProtID")
  )


# Downregulated liver
lost_down_t_pbc_h.m_outliers <- setdiff(
  down_liver_pbc_vs_h.m$GeneName.UniprotID,
  down_liver_pbc_vs_h.m_no_outliers$GeneName.UniprotID
)

filtered_lost_down_t_pbc_h.m_total <- down_liver_pbc_vs_h.m %>%
  filter(GeneName.UniprotID %in% lost_down_t_pbc_h.m_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

filtered_lost_down_t_pbc_h.m_no_outliers <- results.t.pbc.h.m.no_outliers %>%
  filter(GeneName.UniprotID %in% lost_down_t_pbc_h.m_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

merged_lost_down_t_pbc_h.m_outliers <- filtered_lost_down_t_pbc_h.m_total %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(
    filtered_lost_down_t_pbc_h.m_no_outliers %>%
      rename(log2FC_without_outliers = log2FC,
             adj.P.Val_without_outliers = adj.P.Val),
    by = c("GeneName.UniprotID", "GeneName", "UniProtID")
  )


# Upregulated liver
lost_up_t_pbc_h.m_outliers <- setdiff(
  up_liver_pbc_vs_h.m$GeneName.UniprotID,
  up_liver_pbc_vs_h.m_no_outliers$GeneName.UniprotID
)

filtered_lost_up_t_pbc_h.m_total <- up_liver_pbc_vs_h.m %>%
  filter(GeneName.UniprotID %in% lost_up_t_pbc_h.m_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

filtered_lost_up_t_pbc_h.m_no_outliers <- results.t.pbc.h.m.no_outliers %>%
  filter(GeneName.UniprotID %in% lost_up_t_pbc_h.m_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

merged_lost_up_t_pbc_h.m_outliers <- filtered_lost_up_t_pbc_h.m_total %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(
    filtered_lost_up_t_pbc_h.m_no_outliers %>%
      rename(log2FC_without_outliers = log2FC,
             adj.P.Val_without_outliers = adj.P.Val),
    by = c("GeneName.UniprotID", "GeneName", "UniProtID")
  )


#### Proteins significant after excluding outliers ####

# Downregulated plasma
new_down_p_pbc_h.m_outliers <- setdiff(
  down_plasma_pbc_vs_h.m_no_outliers$GeneName.UniprotID,
  down_plasma_pbc_vs_h.m$GeneName.UniprotID
)

filtered_new_down_p_pbc_h.m_no_outliers <- down_plasma_pbc_vs_h.m_no_outliers %>%
  filter(GeneName.UniprotID %in% new_down_p_pbc_h.m_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

filtered_new_down_p_pbc_h.m_total <- results.p.pbc.h.m %>%
  filter(GeneName.UniprotID %in% new_down_p_pbc_h.m_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

merged_new_down_p_pbc_h.m_outliers <- filtered_new_down_p_pbc_h.m_total %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(
    filtered_new_down_p_pbc_h.m_no_outliers %>%
      rename(log2FC_without_outliers = log2FC,
             adj.P.Val_without_outliers = adj.P.Val),
    by = c("GeneName.UniprotID", "GeneName", "UniProtID")
  )


# Upregulated plasma
new_up_p_pbc_h.m_outliers <- setdiff(
  up_plasma_pbc_vs_h.m_no_outliers$GeneName.UniprotID,
  up_plasma_pbc_vs_h.m$GeneName.UniprotID
)

filtered_new_up_p_pbc_h.m_no_outliers <- up_plasma_pbc_vs_h.m_no_outliers %>%
  filter(GeneName.UniprotID %in% new_up_p_pbc_h.m_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

filtered_new_up_p_pbc_h.m_total <- results.p.pbc.h.m %>%
  filter(GeneName.UniprotID %in% new_up_p_pbc_h.m_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

merged_new_up_p_pbc_h.m_outliers <- filtered_new_up_p_pbc_h.m_total %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(
    filtered_new_up_p_pbc_h.m_no_outliers %>%
      rename(log2FC_without_outliers = log2FC,
             adj.P.Val_without_outliers = adj.P.Val),
    by = c("GeneName.UniprotID", "GeneName", "UniProtID")
  )


# Downregulated liver
new_down_t_pbc_h.m_outliers <- setdiff(
  down_liver_pbc_vs_h.m_no_outliers$GeneName.UniprotID,
  down_liver_pbc_vs_h.m$GeneName.UniprotID
)

filtered_new_down_t_pbc_h.m_no_outliers <- down_liver_pbc_vs_h.m_no_outliers %>%
  filter(GeneName.UniprotID %in% new_down_t_pbc_h.m_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

filtered_new_down_t_pbc_h.m_total <- results.t.pbc.h.m %>%
  filter(GeneName.UniprotID %in% new_down_t_pbc_h.m_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

merged_new_down_t_pbc_h.m_outliers <- filtered_new_down_t_pbc_h.m_total %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(
    filtered_new_down_t_pbc_h.m_no_outliers %>%
      rename(log2FC_without_outliers = log2FC,
             adj.P.Val_without_outliers = adj.P.Val),
    by = c("GeneName.UniprotID", "GeneName", "UniProtID")
  )


# Upregulated liver
new_up_t_pbc_h.m_outliers <- setdiff(
  up_liver_pbc_vs_h.m_no_outliers$GeneName.UniprotID,
  up_liver_pbc_vs_h.m$GeneName.UniprotID
)

filtered_new_up_t_pbc_h.m_no_outliers <- up_liver_pbc_vs_h.m_no_outliers %>%
  filter(GeneName.UniprotID %in% new_up_t_pbc_h.m_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

filtered_new_up_t_pbc_h.m_total <- results.t.pbc.h.m %>%
  filter(GeneName.UniprotID %in% new_up_t_pbc_h.m_outliers) %>%
  select(GeneName.UniprotID, GeneName, UniProtID, log2FC, adj.P.Val)

merged_new_up_t_pbc_h.m_outliers <- filtered_new_up_t_pbc_h.m_total %>%
  rename(log2FC_total_cohort = log2FC, adj.P.Val_total_cohort = adj.P.Val) %>%
  inner_join(
    filtered_new_up_t_pbc_h.m_no_outliers %>%
      rename(log2FC_without_outliers = log2FC,
             adj.P.Val_without_outliers = adj.P.Val),
    by = c("GeneName.UniprotID", "GeneName", "UniProtID")
  )


# write.xlsx(list("lost_down_plasma" = merged_lost_down_p_pbc_h.m_outliers,
#                 "lost_up_plasma" = merged_lost_up_p_pbc_h.m_outliers,
#                 "lost_down_liver" = merged_lost_down_t_pbc_h.m_outliers,
#                 "lost_up_liver" = merged_lost_up_t_pbc_h.m_outliers,
#                 "new_down_plasma" = merged_new_down_p_pbc_h.m_outliers,
#                 "new_up_plasma" = merged_new_up_p_pbc_h.m_outliers,
#                 "new_down_liver" = merged_new_down_t_pbc_h.m_outliers,
#                 "new_up_liver" = merged_new_up_t_pbc_h.m_outliers),
#            "../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/proteins_changed_significance_when_excluding_outliers_pbc_vs_healthy_and_masld.xlsx")

#################### Venn diagram

### Total protein expression:

input_total <- list("Plasma proteome" = results.p.pbc[,7], 
                     "Liver proteome" = results.t.pbc[,7])

venn_plot_total <- ggvenn(input_total, fill_color = c("royal blue","brown"), 
                           stroke_size = 0.5, set_name_size = 10, show_percentage = F, text_size = 10) +  
   ggtitle("Venn Diagram of Total Protein Detection") +
   theme_minimal() + # Different theme possible
   theme(
     axis.title = element_blank(),  # Hide axis titles
     axis.text = element_blank(),   # Hide axis text
     axis.ticks = element_blank(),   # Hide axis ticks
     plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),  # Center and adjust main title size
     legend.position = "none",  # Remove legend
     panel.background = element_blank(),  # Remove panel background
     plot.background = element_blank(),  # Remove plot background
     panel.grid = element_blank()  # Remove panel grid lines
   )
 
venn_plot_total
 
# ggsave("../output/venn_total.jpeg",
#         width = 12, height = 12, plot = venn_plot_total, device = "jpeg")

### Venn PBC vs Healthy
 
input_pbc <- list(Plasma = selected_proteins_p_pbc[,7], 
                  Liver = selected_proteins_t_pbc[,7])

venn_plot_pbc <- ggvenn(input_pbc, fill_color = c("royal blue","brown"), 
                        stroke_size = 0.5, set_name_size = 10, show_percentage = F, text_size = 10) +
  ggtitle("Venn Diagram of\nSignificant Differential Protein Expression in\nPBC vs. Healthy") +
  theme_minimal() + 
  theme(
    axis.title = element_blank(),  
    axis.text = element_blank(),   
    axis.ticks = element_blank(),   
    plot.title = element_text(hjust = 0.5, size = 30,face="bold"),  
    legend.position = "none", 
    panel.background = element_blank(), 
    plot.background = element_blank(),  
    panel.grid = element_blank()  
  )

venn_plot_pbc

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/All/venn_pbc_healthy_fdr0.05.jpeg",
#         width = 12, height = 12, plot = venn_plot_pbc, device = "jpeg")

########## Venn diagram without cirrhosis

input_pbc_no_c <- list(Plasma = selected_proteins_p_pbc_no_c[,7], 
                  Liver = selected_proteins_t_pbc_no_c[,7])

venn_plot_pbc_no_c <- ggvenn(input_pbc_no_c, fill_color = c("royal blue","brown"), 
                        stroke_size = 0.5, set_name_size = 10, show_percentage = F, text_size = 10) +
  ggtitle("Venn Diagram of\nSignificant Differential Protein Expression in\nnon-cirrhotic PBC vs. Healthy") +
  theme_minimal() + 
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(),   
    axis.ticks = element_blank(),  
    plot.title = element_text(hjust = 0.5, size = 30,face="bold"),  
    legend.position = "none", 
    panel.background = element_blank(), 
    plot.background = element_blank(),  
    panel.grid = element_blank()  
  )

venn_plot_pbc_no_c

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/No cirrhosis/venn_pbc_healthy_fdr0.05_no_cirrhosis.jpeg",
#         width = 12, height = 12, plot = venn_plot_pbc_no_c, device = "jpeg")


###### Venn Diagram PBC without MASLD vs. healthy

input_pbc_no_masld <- list(Plasma = selected_proteins_p_pbc_no_masld[,7], 
                       Liver = selected_proteins_t_pbc_no_masld[,7])

venn_plot_pbc_no_masld <- ggvenn(input_pbc_no_masld, fill_color = c("royal blue","brown"), 
                             stroke_size = 0.5, set_name_size = 10, show_percentage = F, text_size = 10) +
  ggtitle("Venn Diagram of\nSignificant Differential Protein Expression in\nPBC without MASLD vs. Healthy") +
  theme_minimal() + 
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(),   
    axis.ticks = element_blank(),  
    plot.title = element_text(hjust = 0.5, size = 30,face="bold"),  
    legend.position = "none", 
    panel.background = element_blank(), 
    plot.background = element_blank(),  
    panel.grid = element_blank()  
  )

venn_plot_pbc_no_masld 

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/No MASLD/venn_pbc_healthy_fdr0.05_no_masld.jpeg",
#         width = 12, height = 12, plot = venn_plot_pbc_no_masld, device = "jpeg")

### Venn PBC vs MASLD

input_pbc_masld <- list(Plasma = selected_proteins_p_pbc_masld[,7], 
                        Liver = selected_proteins_t_pbc_masld[,7])

venn_plot_pbc_masld <- ggvenn(input_pbc_masld, fill_color = c("royal blue","brown"), 
                              stroke_size = 0.5, set_name_size = 10, 
                              show_percentage = F, text_size = 10) +
  ggtitle("Venn Diagram of\nSignificant Differential Protein Expression in\nPBC vs. MASLD") +
  theme_minimal() + 
  theme(
    axis.title = element_blank(),  
    axis.text = element_blank(),   
    axis.ticks = element_blank(),   
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),  
    legend.position = "none", 
    panel.background = element_blank(), 
    plot.background = element_blank(),  
    panel.grid = element_blank()  
  )

venn_plot_pbc_masld

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_MASLD/All/venn_pbc_masld_fdr0.05.jpeg",
#         width = 10, height = 10, plot = venn_plot_pbc_masld, device = "jpeg")

### Venn PBC vs MASLD, no cirrhosis

input_pbc_masld_no_c <- list(Plasma = selected_proteins_p_pbc_masld_no_c[,7], 
                             Liver = selected_proteins_t_pbc_masld_no_c[,7])

venn_plot_pbc_masld_no_c <- ggvenn(input_pbc_masld_no_c, fill_color = c("royal blue","brown"), 
                                   stroke_size = 0.5, set_name_size = 10, 
                                   show_percentage = F, text_size = 10) +
  ggtitle("Venn Diagram of\nSignificant Differential Protein Expression in\nnon-cirrhotic PBC vs. MASLD") +
  theme_minimal() + 
  theme(
    axis.title = element_blank(),  
    axis.text = element_blank(),   
    axis.ticks = element_blank(),   
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),  
    legend.position = "none", 
    panel.background = element_blank(), 
    plot.background = element_blank(),  
    panel.grid = element_blank()  
  )

venn_plot_pbc_masld_no_c

# Save plot
# ggsave("../output/FDR0.05_all/PBC_vs_MASLD/No cirrhosis/venn_pbc_masld_fdr0.05_no.cirrhosis.jpeg",
#         width = 10, height = 10, plot = venn_plot_pbc_masld_no_c, device = "jpeg")


### PBC vs Healthy and MASLD 

input_pbc_h.m <- list(Plasma = selected_proteins_p_pbc_h.m[,7], 
                         Liver = selected_proteins_t_pbc_h.m[,7])

venn_plot_pbc_h.m <- ggvenn(input_pbc_h.m, fill_color = c("royal blue","brown"), 
                               stroke_size = 0.5, set_name_size = 10, show_percentage = F, text_size = 10) +
  ggtitle("Venn Diagram of\nSignificant Differential Protein Expression in\nPBC vs. Healthy and MASLD") +
  theme_minimal() +  
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(),   
    axis.ticks = element_blank(),   
    plot.title = element_text(hjust = 0.5, size = 30,face="bold"),  
    legend.position = "none", 
    panel.background = element_blank(),  
    plot.background = element_blank(),  
    panel.grid = element_blank()  
  )

venn_plot_pbc_h.m

# Save plot
# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/All/venn_pbc_h.m_fdr0.05.jpeg",
#         width = 15, height = 12, plot = venn_plot_pbc_h.m, device = "jpeg")

########## Venn diagram without cirrhosis

input_pbc_h.m_no_c <- list(Plasma = selected_proteins_p_pbc_h.m_no_c[,7], 
                              Liver = selected_proteins_t_pbc_h.m_no_c[,7])

venn_plot_pbc_h.m_no_c <- ggvenn(input_pbc_h.m_no_c, fill_color = c("royal blue","brown"), 
                                    stroke_size = 0.5, set_name_size = 10, show_percentage = F, text_size = 10) +
  ggtitle("Venn Diagram of\nSignificant Differential Protein Expression in\nnon-cirrhotic PBC vs. Healthy and MASLD") +
  theme_minimal() +  
  theme(
    axis.title = element_blank(),  
    axis.text = element_blank(),   
    axis.ticks = element_blank(),  
    plot.title = element_text(hjust = 0.5, size = 30,face="bold"),  
    legend.position = "none", 
    panel.background = element_blank(),  
    plot.background = element_blank(),  
    panel.grid = element_blank() 
  )

venn_plot_pbc_h.m_no_c

# Save plot

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy_and_MASLD/No cirrhosis/venn_pbc_h.m_fdr0.05_no_cirrhosis.jpeg",
#         width = 16, height = 10, plot = venn_plot_pbc_h.m_no_c, device = "jpeg")


#### GO enrichment ####

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

enrich_ont <- function(genes, ont) {
  enrichGO(gene          = genes,
           OrgDb         = org.Hs.eg.db,
           keyType       = "ENTREZID",
           ont           = ont,
           pAdjustMethod = "BH",
           pvalueCutoff  = 1,
           qvalueCutoff  = 0.20,
           readable      = TRUE) |>
    as.data.frame() |>
    mutate(ONTOLOGY = ont)
}

pal <- c("GO Biological Process" = "#3B82F6",
         "GO Cellular Component" = "#F59E0B",
         "GO Molecular Function" = "#10B981")



#### GO Plasma ####

ids_uniprot <- unique(na.omit(selected_proteins_p_pbc$UniProtID))

gene_ids_p_pbc <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = ids_uniprot,
                                        columns = "ENTREZID",
                                        keytype = "UNIPROT") |>
  distinct(ENTREZID, .keep_all = TRUE) |>
  filter(!is.na(ENTREZID))

genes_entrez <- unique(gene_ids_p_pbc$ENTREZID)

res_bp <- enrich_ont(genes_entrez, "BP")
res_cc <- enrich_ont(genes_entrez, "CC")
res_mf <- enrich_ont(genes_entrez, "MF")

go_all <- bind_rows(res_bp, res_cc, res_mf) |>
  mutate(Count = str_count(geneID, "/") + 1) |>
  filter(p.adjust <= 0.05) |>
  mutate(ONTOLOGY = recode(ONTOLOGY,
                           BP = "GO Biological Process",
                           CC = "GO Cellular Component",
                           MF = "GO Molecular Function"))

top_n_per_ont <- 12

plot_df <- go_all |>
  group_by(ONTOLOGY) |>
  arrange(desc(Count), p.adjust, .by_group = TRUE) |>
  slice_head(n = top_n_per_ont) |>
  ungroup() |>
  mutate(
    Term = fct_reorder(Description, Count),
    padj_lab = label_scientific(digits = 2)(p.adjust)
  )

max_count <- max(plot_df$Count, na.rm = TRUE)
pad <- max(2, ceiling(max_count * 0.20))

p_plasma <- ggplot(plot_df, aes(x = Count, y = Term, fill = ONTOLOGY)) +
  geom_col(width = 0.8) +
  geom_text(aes(x = max_count + pad * 0.55, label = padj_lab),
            hjust = 0, size = 3.1) +
  facet_grid(rows = vars(ONTOLOGY), scales = "free_y", space = "free_y") +
  scale_fill_manual(values = pal, guide = guide_legend(title = NULL)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)),
                     limits = c(0, max_count + pad)) +
  labs(title = "Plasma — Gene Ontology Enrichment Analysis",
       x = "Protein count", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.text.y = element_text(face = "bold"),
    legend.position = "top",
    panel.spacing.y = unit(10, "pt"),
    plot.margin = margin(10, 22, 10, 10),
    axis.text.y = element_text(size = 9)
  ) +
  annotate("text",
           x = max_count + pad * 0.55, y = Inf, label = "adj. P-value",
           vjust = -0.5, hjust = 0, size = 3.1)

print(p_plasma)

# ggsave("../output/GO_enrichment_plasma_pbc_vs_healthy.jpeg",
#         p_plasma, width = 8, height = 10, dpi = 300)


#### GO Liver ####

ids_uniprot_t <- unique(na.omit(selected_proteins_t_pbc$UniProtID))

gene_ids_t_pbc <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = ids_uniprot_t,
                                        columns = "ENTREZID",
                                        keytype = "UNIPROT") |>
  distinct(ENTREZID, .keep_all = TRUE) |>
  filter(!is.na(ENTREZID))

genes_entrez_t <- unique(gene_ids_t_pbc$ENTREZID)

res_bp_t <- enrich_ont(genes_entrez_t, "BP")
res_cc_t <- enrich_ont(genes_entrez_t, "CC")
res_mf_t <- enrich_ont(genes_entrez_t, "MF")

go_all_t <- bind_rows(res_bp_t, res_cc_t, res_mf_t) |>
  mutate(Count = str_count(geneID, "/") + 1) |>
  filter(p.adjust <= 0.05) |>
  mutate(ONTOLOGY = recode(ONTOLOGY,
                           BP = "GO Biological Process",
                           CC = "GO Cellular Component",
                           MF = "GO Molecular Function"))

plot_df_t <- go_all_t |>
  group_by(ONTOLOGY) |>
  arrange(desc(Count), p.adjust, .by_group = TRUE) |>
  slice_head(n = top_n_per_ont) |>
  ungroup() |>
  mutate(
    Term = fct_reorder(Description, Count),
    padj_lab = label_scientific(digits = 2)(p.adjust)
  )

max_count_t <- max(plot_df_t$Count, na.rm = TRUE)
pad_t <- max(2, ceiling(max_count_t * 0.20))

p_liver <- ggplot(plot_df_t, aes(x = Count, y = Term, fill = ONTOLOGY)) +
  geom_col(width = 0.8) +
  geom_text(aes(x = max_count_t + pad_t * 0.55, label = padj_lab),
            hjust = 0, size = 3.1) +
  facet_grid(rows = vars(ONTOLOGY), scales = "free_y", space = "free_y") +
  scale_fill_manual(values = pal, guide = guide_legend(title = NULL)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)),
                     limits = c(0, max_count_t + pad_t)) +
  labs(title = "Liver — Gene Ontology Enrichment Analysis",
       x = "Protein count", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.text.y = element_text(face = "bold"),
    legend.position = "top",
    panel.spacing.y = unit(10, "pt"),
    plot.margin = margin(10, 22, 10, 10),
    axis.text.y = element_text(size = 9)
  ) +
  annotate("text",
           x = max_count_t + pad_t * 0.55, y = Inf, label = "adj. P-value",
           vjust = -0.5, hjust = 0, size = 3.1)

print(p_liver)

# ggsave("GO_enrichment_liver_pbc_vs_healthy.jpeg",
#        p_liver, width = 8, height = 10, dpi = 300)

################################# Enrichment analysis GOBP ############################################

library("BiocManager")
library("clusterProfiler")
library("org.Hs.eg.db") # For human gene annotations, adjust for other species if needed
library("AnnotationDbi")
library(dplyr)

# 1. Provide a list of gene symbols or Entrez IDs which we're interested in analyzing for enrichment.
# 2. The org.Hs.eg.db package contains annotation mappings for genes in the specified organism (in this case, human). 
#   It includes mappings between gene symbols, Entrez IDs, and various other identifiers, as well as associated 
#   annotation data like GO terms.
# 3. enrichGO uses the provided gene list and the organism database to determine whether certain GO terms are 
#   overrepresented in your gene list compared to what would be expected by chance. It calculates enrichment scores 
#   and associated statistics, such as adjusted p-values (q-values) using methods like the Benjamini-Hochberg (BH) 
#   correction for multiple testing.
# 4. The result of enrichGO is a data frame containing enriched GO terms along with associated statistics like p-values, 
#   adjusted p-values (q-values), and gene counts. This output helps identify which biological processes, 
#   molecular functions, or cellular components are overrepresented in gene list. So, even though we don't 
#   explicitly provide p-values or significance values, enrichGO calculates these statistics based on the gene 
#   list we provided and the annotations available in the organism database.

#################### PBC vs. Healthy plasma  #######################################################################

dim(selected_proteins_p_pbc) # those proteins that have an adj p value < 0.05
selected_proteins_p_pbc[1:5,]

# Map gene symbols to Entrez IDs (required for clusterProfiler)http://127.0.0.1:17963/graphics/plot_zoom_png?width=1872&height=900
gene_ids_p_pbc <- AnnotationDbi::select(org.Hs.eg.db, 
                                        keys = selected_proteins_p_pbc$UniProtID, 
                                        columns = c("ENTREZID"), 
                                        keytype = "UNIPROT")

# Remove any NA values that might be present
dim(gene_ids_p_pbc)
gene_ids_p_pbc[1:10,]
gene_ids_p_pbc <- na.omit(gene_ids_p_pbc)
dim(gene_ids_p_pbc)

# Perform GO enrichment analysis - first biological procesess
go_enrich_p_pbc_bp <- enrichGO(gene = gene_ids_p_pbc$ENTREZID, 
                               OrgDb = org.Hs.eg.db, 
                               keyType = "ENTREZID", 
                               ont = "BP", # Can be "BP", "MF", or "CC" for Biological Process, Molecular Function, or Cellular Component
                               pAdjustMethod = "BH", 
                               pvalueCutoff = 0.05, 
                               qvalueCutoff = 0.05)

# View the results
head(go_enrich_p_pbc_bp)

# Visualize the results
barplot(go_enrich_p_pbc_bp, showCategory = 20)

go_enrich_p_pbc_bp_dotplot <- dotplot(go_enrich_p_pbc_bp, showCategory = 20) +
  theme(axis.text.y = element_text(size = 10)) # Adjust text size
go_enrich_p_pbc_bp_dotplot

# Save the plot using ggsave

# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/All/enrich_bp_plasma_pbc_healthy_fdr0.05.jpeg",
#         width = 12, height = 12, plot = go_enrich_p_pbc_bp_dotplot, dpi = 300, device = "jpeg")

# Now Molecular function
go_enrich_p_pbc_mf <- enrichGO(gene = gene_ids_p_pbc$ENTREZID, 
                               OrgDb = org.Hs.eg.db, 
                               keyType = "ENTREZID", 
                               ont = "MF", # Can be "BP", "MF", or "CC" for Biological Process, Molecular Function, or Cellular Component
                               pAdjustMethod = "BH", 
                               pvalueCutoff = 0.05, 
                               qvalueCutoff = 0.05)

# View the results
head(go_enrich_p_pbc_mf)

# Visualize the results
barplot(go_enrich_p_pbc_mf, showCategory = 20)

go_enrich_p_pbc_mf_dotplot <- dotplot(go_enrich_p_pbc_mf, showCategory = 20) +
  theme(axis.text.y = element_text(size = 10)) # Adjust text size
go_enrich_p_pbc_mf_dotplot

# Save the plot using ggsave
# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/All/enrich_mf_plasma_pbc_healthy_fdr0.05.jpeg",
#         width = 12, height = 12, plot = go_enrich_p_pbc_mf_dotplot, dpi = 300, device = "jpeg")

# Now cellular component
go_enrich_p_pbc_cc <- enrichGO(gene = gene_ids_p_pbc$ENTREZID, 
                               OrgDb = org.Hs.eg.db, 
                               keyType = "ENTREZID", 
                               ont = "CC", # Can be "BP", "MF", or "CC" for Biological Process, Molecular Function, or Cellular Component
                               pAdjustMethod = "BH", 
                               pvalueCutoff = 0.05, 
                               qvalueCutoff = 0.05)

# View the results
head(go_enrich_p_pbc_cc)

# Visualize the results
barplot(go_enrich_p_pbc_cc, showCategory = 20)

go_enrich_p_pbc_cc_dotplot <- dotplot(go_enrich_p_pbc_cc, showCategory = 20) +
  theme(axis.text.y = element_text(size = 10)) # Adjust text size
go_enrich_p_pbc_cc_dotplot

# Save the plot using ggsave
# ggsave("../output/FDR0.05_all/PBC_vs_Healthy/All/enrich_cc_plasma_pbc_healthy_fdr0.05.jpeg",
#         width = 12, height = 12, plot = go_enrich_p_pbc_cc_dotplot, dpi = 300, device = "jpeg")



