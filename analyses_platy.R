# loading packages
library(dartR)
library(stringr)
library(ape)
library(data.table)
library(hierfstat)
library(reshape2)
library(broom)
library(related)
library(plotly)
################################################################################
########################## DATA INPUT SNP DATA #################################
################################################################################
# Read data
gl.data_original <-
  gl.read.dart("Report_DPla19-4171_SNP_mapping_2.csv",
               ind.metafile = "ID_Pop_Platypus.csv")
# Assigning populations
gl.data_original <- gl.recode.pop(gl.data_original,
                                  pop.recode = "new_pop_assignments.csv")
table(pop(gl.data_original))
# Formatting chromosome names
gl.data_original$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1 <-
  str_extract(gl.data_original$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1,
              "[^_]*_[^_]*")
# Assigning chromosome name
gl.data_original$chromosome <-
  as.factor(gl.data_original$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1)
# Assigning SNP position
gl.data_original$position <-
  gl.data_original$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
# Amending data
# T42 was a recapture of T28
# T3 (reported as male) and T5 (reported as female) are the same sample.
# T3 has more missing data.
gl.report.callrate(gl.data_original, method = "ind")
# E32 and V34 have missing data in more than 50% of loci
# Dropping individuals
gl.data <-
  gl.drop.ind(gl.data_original, ind.list = c("T42", "T3", "E32", "V34"))
# V32 and V30 were mislabeled, they belong to the opposite population.
V32 <- which(gl.data$ind.names == "V32")
V30 <- which(gl.data$ind.names == "V30")
gps_temp_1 <- gl.data$other$latlong[V32, ]
gps_temp_2 <- gl.data$other$latlong[V30, ]
gl.data$other$latlong[V32, ] <- gps_temp_2
gl.data$other$latlong[V30, ] <- gps_temp_1
gl.data$pop[V32] <- "OVENS"
gl.data$pop[V30]  <- "MITTA_ABOVE"

################################################################################
#################### DATA INPUT MICROSATELLITE DATA ############################
################################################################################
# read genepop
platy_micro <- read.genepop(file = "platy_microsats.gen", ncode = 2)
ind_info <- read.table("ind_microsat.csv", sep = ",", header = TRUE)
platy_micro@pop <-  as.factor(ind_info$Population)
strata(platy_micro) <- data.frame(pop = ind_info$Population)
platy_micro_hierfstat <- genind2hierfstat(platy_micro)

below_above <- read.genepop(file = "below_above.gen", ncode = 2)
above_unreg <- read.genepop(file = "above_unreg.gen", ncode = 2)
below_unreg <- read.genepop(file = "below_unreg.gen", ncode = 2)
#just the top and bottom locations
unreg_below_above <-
  read.genepop(file = "unreg_below_above.gen", ncode = 2)

################################################################################
########################### FILTERING SNP DATA #################################
################################################################################

# Reproducibility (RepAvg)
gl.report.reproducibility(gl.data)
platy <-  gl.filter.reproducibility(gl.data, threshold = 1)
# Retain only one SNP per read
gl.report.secondaries(platy)
platy <- gl.filter.secondaries(platy, method = "best")
# Departure from Hardy-Weinberg proportions
gl.report.hwe(
  platy,
  multi_comp = TRUE,
  multi_comp_method = "bonferroni",
  alpha_val = 0.05
)
platy <- gl.filter.hwe(
  platy,
  multi_comp = TRUE,
  multi_comp_method = "bonferroni",
  alpha_val = 0.05
)
# Mapped to chromosome
gl.report.locmetric(platy, metric = "AlnCnt_Platypus_Chrom_NCBIv1")
platy <- gl.filter.locmetric(
  platy,
  metric = "AlnCnt_Platypus_Chrom_NCBIv1",
  lower = 1,
  upper = max(platy$other$loc.metrics$AlnCnt_Platypus_Chrom_NCBIv1)
)
# BLAST alignment E-value
gl.report.locmetric(platy, metric = "AlnEvalue_Platypus_Chrom_NCBIv1")
platy <- gl.filter.locmetric(
  platy,
  metric = "AlnEvalue_Platypus_Chrom_NCBIv1",
  lower = min(platy$other$loc.metrics$AlnEvalue_Platypus_Chrom_NCBIv1),
  upper = 1e-20
)

# Dataset used to analyse variation between groups
# Missing data by site
gl.report.callrate(platy)
platy_between <- gl.filter.callrate(platy, threshold = 0.9)
# Minor allele count (MAC)
gl.report.maf(platy_between)
platy_between <- gl.filter.maf(platy_between, threshold = 3)
# Linkage disequilibrium (r2)
res <- gl.report.ld.map(platy_between, ld_max_pairwise = 1e7)
gl.ld.distance(res)
platy_between <-
  gl.filter.ld(platy_between, ld_report = res , threshold = 0.2)
# Remove sites located within coding regions
gff <- read.gff("ref_mOrnAna1.p.v1_top_level.gff3")
genes <- gff[gff$type == "gene", ]

genes_by_chr <- split(genes, genes$seqid)
loci_by_chr <- data.frame(
  chrom = platy_between$chromosome,
  loc_name = locNames(platy_between),
  loc_pos = platy_between$position
)
loci_by_chr <- split(loci_by_chr, loci_by_chr$chrom)

vec_loci_no_genes <- NULL
for (i in 1:length(loci_by_chr)) {
  loci <- loci_by_chr[[i]]
  gene_chr <- which(names(genes_by_chr) == names(loci_by_chr[i]))
  gene_chr_2 <- genes_by_chr[[gene_chr]]
  y  <- list(gene_chr_2$start, gene_chr_2$end)
  loci_no_genes <- loci[!loci$loc_pos %inrange% y, ]
  vec_loci_no_genes <- c(vec_loci_no_genes, loci_no_genes$loc_name)
}

platy_between <-
  gl.keep.loc(platy_between, loc.list = vec_loci_no_genes)

# Dataset used to analyse variation within groups
# Missing data by site
platy_within <- gl.filter.callrate(platy, threshold = 1)
# Remove sites located within sex chromosomes
chrom_platy <- read.csv("chrom_platypus.csv")
sex_chrom <- chrom_platy[22:31, "RefSeq_mOrnAna_v1"]
sex_chrom_platy <-
  locNames(platy_within)[which(platy_within$chromosome %in%
                                 sex_chrom)]
platy_within <- gl.drop.loc(platy_within, loc.list = sex_chrom_platy)

################################################################################
################## Genetic variation within groups #############################
################################################################################
################## Microsatellite data #########################################
# standard error is the standard deviation divided by the square root of the
# sample size.
std <- function(x) {
  x <- x[!is.na(x)]
  se <- sd(x) / sqrt(length(x))
  return(se)
}

stats_micro <- basic.stats(platy_micro_hierfstat)
ho_micro <- round(colMeans(stats_micro$Ho, na.rm = TRUE), 3)
ho_micro_se <- round(apply(stats_micro$Ho, 2, std), 3)
he_micro <- round(colMeans(stats_micro$Hs, na.rm = TRUE), 3)
he_micro_se <- round(apply(stats_micro$Hs, 2, std), 3)
fis_micro <- round(colMeans(stats_micro$Fis, na.rm = TRUE), 3)
fis_micro_se <- round(apply(stats_micro$Fis, 2, std), 3)
Ar_micro <-
  round(colMeans(allelic.richness(platy_micro_hierfstat)$Ar,
                 na.rm = TRUE), 3)

# testing normality
het_groups_micro <- as.data.frame(stats_micro$Hs)
plot(density(het_groups_micro$above))
shapiro.test(het_groups_micro$above)
shapiro.test(het_groups_micro$below)
shapiro.test(het_groups_micro$unreg)

# testing whether He was significantly different between groups
##### NOT significant
below_above_test <-
  wilcox.test(het_groups_micro$above, het_groups_micro$below)
##### NOT significant
below_unreg_test <-
  wilcox.test(het_groups_micro$below, het_groups_micro$unreg)
##### NOT significant
above_unreg_test <-
  wilcox.test(het_groups_micro$above, het_groups_micro$unreg)

################## SNP data ####################################################
stats_snp <- gl.basic.stats(platy_within)
ho_snp <- round(colMeans(stats_snp$Ho, na.rm = TRUE), 3)
ho_snp_se <- round(apply(stats_snp$Ho, 2, std), 3)
he_snp <- round(colMeans(stats_snp$Hs, na.rm = TRUE), 3)
he_snp_se <- round(apply(stats_snp$Hs, 2, std), 3)
fis_snp <- round(colMeans(stats_snp$Fis, na.rm = TRUE), 3)
fis_snp_se <- round(apply(stats_snp$Fis, 2, std), 3)

platy_genind <- gl2gi(platy_within)
platy_snp_hierfstat <- genind2hierfstat(platy_genind)
Ar_snp <-
  round(colMeans(allelic.richness(platy_snp_hierfstat)$Ar, na.rm = TRUE), 3)

# Mean SNP genetic variation across all rivers (expected heterozygosity)
mean_het_snp <- stats_snp$overall["Hs"]

# testing normality
# data did not conform to a normal distribution
shapiro_test <- apply(stats_snp$Hs, 2, shapiro.test)

# testing whether He was significantly different between groups
het_groups_snp <- as.data.frame(stats_snp$Hs)
##### significant
snowy_thredbo <-
  wilcox.test(het_groups_snp$SNOWY, het_groups_snp$THREDBO)
##### significant
snowy_eucumbene_above <- wilcox.test(het_groups_snp$SNOWY,
                                     het_groups_snp$EUCUMBENE_ABOVE)
##### significant
snowy_eucumbene_below <- wilcox.test(het_groups_snp$SNOWY,
                                     het_groups_snp$EUCUMBENE_BELOW)
##### significant
thredbo_eucumbene_below <- wilcox.test(het_groups_snp$THREDBO,
                                       het_groups_snp$EUCUMBENE_BELOW)
##### significant
thredbo_eucumbene_above <- wilcox.test(het_groups_snp$THREDBO,
                                       het_groups_snp$EUCUMBENE_ABOVE)
##### significant
eucumbene_below_eucumbene_above <-
  wilcox.test(het_groups_snp$EUCUMBENE_BELOW,
              het_groups_snp$EUCUMBENE_ABOVE)
##### significant
ovens_mitta_below <- wilcox.test(het_groups_snp$OVENS,
                                 het_groups_snp$MITTA_BELOW)
##### significant
ovens_mitta_above <- wilcox.test(het_groups_snp$OVENS,
                                 het_groups_snp$MITTA_ABOVE)
##### significant
mitta_above_mitta_below <- wilcox.test(het_groups_snp$MITTA_ABOVE,
                                       het_groups_snp$MITTA_BELOW)
##### significant
tenter_severn_above <- wilcox.test(het_groups_snp$TENTERFIELD,
                                   het_groups_snp$SEVERN_ABOVE)
##### significant
tenter_severn_below <- wilcox.test(het_groups_snp$TENTERFIELD,
                                   het_groups_snp$SEVERN_BELOW)
##### NOT significant
severn_below_severn_above <-
  wilcox.test(het_groups_snp$SEVERN_BELOW,
              het_groups_snp$SEVERN_ABOVE)

# testing whether He was significantly different between regions
diversity_regions <- gl.recode.pop(platy_within,
                             pop.recode = "new_pop_assignments_regions.csv")
table(pop(diversity_regions))
regions_stats <- gl.basic.stats(diversity_regions)
regions_expected_he <- colMeans(regions_stats$Hs, na.rm = TRUE)

het_groups_regions <- as.data.frame(regions_stats$Hs)

# significant
border_u_murray <- wilcox.test(het_groups_regions$BORDER,
                               het_groups_regions$U_MURRAY)
# significant
border_snowy <- wilcox.test(het_groups_regions$BORDER,
                            het_groups_regions$SNOWY)
# NOT significant
snowy_u_murray <- wilcox.test(het_groups_regions$SNOWY,
                              het_groups_regions$U_MURRAY)

################################################################################
################## Genetic variation between groups ############################
################################################################################
# testing whether groups separated by major dams are more genetically different
# than otherwise
dam_snp <- platy_between
pop(dam_snp) <- dam_snp$other$ind.metrics$pop
dam_snp <-
  gl.recode.pop(dam_snp, pop.recode = "new_pop_assignments_dam.csv")
table(pop(dam_snp))
# Dividing sampling sites of each pair of rivers into comparable upstream and
# downstream groups.
# For unregulated rivers (Wingecarribee, Tenterfield and Ovens), the division
# point was chosen at a comparable position to the dam in the paired regulated
# river.
# For regulated rivers (Nepean, Severn and Mitta-Mitta), the dam, ignoring the
# reservoir, was used as reference point for the division.

severn_above_tenter_dam <- gl.keep.pop(platy_between,
                                  pop.list = c("SEVERN_ABOVE", "TENTERFIELD"))
severn_below_tenter_dam <- gl.keep.pop(platy_between,
                                 pop.list = c("SEVERN_BELOW", "TENTERFIELD"))
severn_dam <- gl.keep.pop(platy_between,
                          pop.list = c("SEVERN_BELOW", "SEVERN_ABOVE"))
tenter_dam <-
  gl.keep.pop(dam_snp, pop.list = c("tenter_below", "tenter_above"))
mitta_above_ovens_dam <- gl.keep.pop(platy_between,
                                     pop.list = c("MITTA_ABOVE", "OVENS"))
mitta_below_ovens_dam <- gl.keep.pop(platy_between,
                                     pop.list = c("MITTA_BELOW", "OVENS"))
mitta_dam <-
  gl.keep.pop(platy_between, pop.list = c("MITTA_BELOW", "MITTA_ABOVE"))
ovens_dam <-
  gl.keep.pop(dam_snp, pop.list = c("ovens_below", "ovens_above"))
snowy_thredbo_dam <- gl.keep.pop(platy_between,
                                 pop.list = c("THREDBO", "SNOWY"))
snowy_eucumbene_above_dam <- gl.keep.pop(platy_between,
                                    pop.list = c("SNOWY", "EUCUMBENE_ABOVE"))
snowy_eucumbene_below_dam <- gl.keep.pop(platy_between,
                                   pop.list = c("SNOWY", "EUCUMBENE_BELOW"))
thredbo_eucumbene_above_dam <- gl.keep.pop(platy_between,
                                  pop.list = c("THREDBO", "EUCUMBENE_ABOVE"))
thredbo_eucumbene_below_dam <- gl.keep.pop(platy_between,
                                  pop.list = c("THREDBO", "EUCUMBENE_BELOW"))
eucumbene_dam <- gl.keep.pop(platy_between,
                             pop.list = c("EUCUMBENE_BELOW", "EUCUMBENE_ABOVE"))

pop_list <- list(
  severn_above_tenter = severn_above_tenter_dam,
  severn_below_tenter = severn_below_tenter_dam,
  severn = severn_dam,
  tenter = tenter_dam,
  mitta_above_ovens = mitta_above_ovens_dam,
  mitta_below_ovens = mitta_below_ovens_dam,
  mitta = mitta_dam,
  ovens = ovens_dam,
  snowy_thredbo = snowy_thredbo_dam,
  snowy_eucumbene_above = snowy_eucumbene_above_dam,
  snowy_eucumbene_below = snowy_eucumbene_below_dam,
  thredbo_eucumbene_above = thredbo_eucumbene_above_dam,
  thredbo_eucumbene_below = thredbo_eucumbene_below_dam,
  eucumbene = eucumbene_dam,
  winge_nepean_above = above_unreg,
  winge_nepean_below = below_unreg,
  nepean = below_above,
  winge = unreg_below_above
)

# x is a list of genind and/or genlight objects
gen_diff <- function(x) {
  mutual_information <- function(mat) {
    EstMLEFun <- function(mat) {
      # MLE
      entropyFun <- function(p) {
        p <- p[p > 0]
        out <- -sum(p * log(p))
        return(out)
      }
      n <- sum(mat)
      prob.hat <- mat / n
      px.hat <- apply(prob.hat, 1, sum)
      py.hat <- apply(prob.hat, 2, sum)
      I.hat <-
        entropyFun(px.hat) + entropyFun(py.hat) - entropyFun(prob.hat)
      # MLE of Mutual Information!
      return(I.hat)
    }
    mydata <- as.matrix(mat)
    est <- EstMLEFun(mydata)
    return(est)
  }
  
  res <- as.data.frame(matrix(nrow = length(x), ncol = 7))
  colnames(res) <-
    c("Group", "FST", "FST_SE", "MI", "MI_SE", "Dest", "Dest_SE")
  
  for (i in 1:length(x)) {
    tmp <- x[[i]]
    
    if (is.genind(tmp)) {
      tmp_hierf <- hierfstat::genind2hierfstat(tmp)
    } else{
      tmp_genind <- dartR::gl2gi(tmp)
      tmp_hierf <-  hierfstat::genind2hierfstat(tmp_genind)
    }
    
    tmp_res <- hierfstat::basic.stats(tmp_hierf)
    res[i, "Group"] <- names(x[i])
    res[i, "FST"] <- round(tmp_res$overall["Fstp"], 3)
    res[i, "FST_SE"] <-  round(std(tmp_res$perloc$Fstp), 3)
    res[i, "Dest"] <- round(tmp_res$overall["Dest"], 3)
    res[i, "Dest_SE"] <- round(std(tmp_res$perloc$Dest), 3)
    
    tmp_allele_matrix <- allele.count(tmp_hierf)
    tmp_MI <- lapply(tmp_allele_matrix, mutual_information)
    res[i, "MI"] <- round(mean(unlist(tmp_MI)), 3)
    res[i, "MI_SE"] <- round(std(unlist(tmp_MI)), 3)
    
  }
  
  return(res)
  
}

fin_res <- gen_diff(pop_list)

# Testing the significance of the difference of FST values between dammed and
# unregulated rivers
severn_dam_res <- gl.basic.stats(severn_dam)
tenter_res <- gl.basic.stats(tenter_dam)
mitta_res <- gl.basic.stats(mitta_dam)
ovens_res <- gl.basic.stats(ovens_dam)
below_above_hierf <- genind2hierfstat(below_above)
below_above_res <- basic.stats(below_above_hierf)
unreg_below_above_hierf <- genind2hierfstat(unreg_below_above)
unreg_below_above_res <- basic.stats(unreg_below_above_hierf)

# NOT significant
wilcoxon_mitta_ovens <- wilcox.test(mitta_res$perloc$Fstp,
                                    ovens_res$perloc$Fstp)
# significant
wilcoxon_tenter_severn <- wilcox.test(severn_dam_res$perloc$Fstp,
                                      tenter_res$perloc$Fstp)
# significant
wilcoxon_winge_nepean <- wilcox.test(below_above_res$perloc$Fstp,
                                     unreg_below_above_res$perloc$Fstp)

# testing whether the number of platypus generations since the building of the
# dams can predict the genetic differentiation

dams_diff <- read.csv("fst_dams.csv")
fst_gen <- glance(lm(dams_diff$gen_barrier ~ dams_diff$FST))
MI_gen <- glance(lm(dams_diff$gen_barrier ~ dams_diff$MI))
Dest_gen <- glance(lm(dams_diff$gen_barrier ~ dams_diff$Dest))

reg_res <- rbind(fst_gen, MI_gen, Dest_gen)
reg_res <-
  data.frame(gen_diff = c("FST", "Mutual Information", "Jost's D"), reg_res)

ggplot(dams_diff, aes(x = gen_barrier, y = FST)) +
  geom_point(color = "deeppink", size = 2) +
  geom_smooth(method = "lm") +
  ylab("FST") +
  xlab("Generations after dam building") +
  theme_bw()

ggsave(
  "diff_dam_fst.pdf",
  width = 3,
  height = 3,
  units = "in",
  dpi = "retina",
  bg = "transparent"
)

dams_diff_2 <-
  melt(dams_diff[, c("River", "FST", "MI", "Dest", "gen_barrier"), ],
       id = c("River", "gen_barrier"))

ggplot(dams_diff_2, aes(x = gen_barrier, y = value, color = variable)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm") +
  facet_wrap( ~ variable, scales = "free_y") +
  ylab("Genetic differentiation") +
  xlab("Generations after dam building") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(
  "diff_dam_all.pdf",
  width = 8,
  height = 4,
  units = "in",
  dpi = "retina",
  bg = "transparent"
)

################## PCA  ########################################################

plot_pca <- function(x,
                     pt.colors,
                     xaxis = 1,
                     yaxis = 2,
                     zaxis = NULL,
                     pt.size = 6,
                     axis.label.size = 14) {
  # Create a dataframe to hold the required scores
  # the tab	slot is a matrix containing the retained principal components of PCA
  
  if (is.null(zaxis)) {
    m <- cbind(x$tab[, xaxis],x$tab[, yaxis])
  } else {
    m <- cbind(x$tab[, xaxis], x$tab[, yaxis],x$tab[, zaxis])
  }
  df <- data.frame(m)
  
  # Convert the eigenvalues to percentages
  s <- sum(x$pca.eig[x$pca.eig >= 0])
  e <- round(x$pca.eig * 100 / s, 1)

  xlab <- paste("PCA Axis", xaxis, "(", e[xaxis], "%)")
  ylab <- paste("PCA Axis", yaxis, "(", e[yaxis], "%)")
  if (!is.null(zaxis)) {
    zlab <- paste("PCA Axis", zaxis, "(", e[zaxis], "%)")
  }
  ind <- rownames(x$ind.coord)
  pop <- factor(x$grp)
  df <- cbind(df, ind, pop)
  
  # Labels for the axes and points
  if (is.null(zaxis)) {
    colnames(df) <- c("PCoAx", "PCoAy", "ind", "pop")
  } else {
    colnames(df) <- c("PCoAx", "PCoAy", "PCoAz", "ind", "pop")
  }

  # Plot
  plott <- ggplot(df, aes(x = PCoAx,
                          y = PCoAy,
                          color = pop)) +
    geom_point(size = pt.size, aes(color = pop), alpha = 0.5) +
    theme(
      axis.title = element_text(
        face = "bold.italic",
        size = axis.label.size,
        color = "black"
      ),
      axis.text.x = element_text(
        face = "bold",
        angle = 0,
        vjust = 0.5,
        size = axis.label.size
      ),
      axis.text.y = element_text(
        face = "bold",
        angle = 0,
        vjust = 0.5,
        size = axis.label.size
      )
    ) +
    labs(x = xlab, y = ylab) +
    scale_color_manual(values = pt.colors) +
    geom_hline(yintercept = 0) +
     geom_vline(xintercept = 0)  +
     theme(legend.position = "none")
  
  if (!is.null(zaxis)) {
    
    plott <-
      plotly::plot_ly(
        df,
        x = ~ PCoAx,
        y = ~ PCoAy,
        z = ~ PCoAz,
        marker = list(size = pt.size),
        colors = pt.colors,
        text = ind
      ) %>%
      plotly::add_markers(color = ~ pop) %>%
      plotly::layout(
        legend = list(title = list(text = "Populations")),
        scene = list(
          xaxis = list(
            title = xlab,
            titlefont = list(size = axis.label.size)
          ),
          yaxis = list(
            title = ylab,
            titlefont = list(size = axis.label.size)
          ),
          zaxis = list(
            title = zlab,
            titlefont = list(size = axis.label.size)
          )
        )
      )
    show(plott)
  }
  
  print(plott)
  
}

# separating populations
snowy <- gl.keep.pop(platy_between,
                     pop.list = c("THREDBO" , "SNOWY",
                                  "EUCUMBENE_BELOW",
                                  "EUCUMBENE_ABOVE"))
murray <-
  gl.keep.pop(platy_between,
              pop.list = c("MITTA_ABOVE" , "MITTA_BELOW",
                           "OVENS"))
border <-
  gl.keep.pop(platy_between,
              pop.list = c("TENTERFIELD", "SEVERN_BELOW",
                           "SEVERN_ABOVE"))

myCol <- c("deeppink", "darkslategray4", "darkgoldenrod1", "blue")
################### Central NSW Rivers #########################################
pcoa_central <- dapc(platy_micro, n.pca = 10, n.da = 3)
# 2D plot
plot_pca(pcoa_central, pt.colors = myCol)

ggsave(
  "pcoa_central.pdf",
  width = 4,
  height = 4,
  units = "in",
  dpi = "retina",
  bg = "transparent"
)
# 3D plot
plot_pca(pcoa_central, pt.colors = myCol,zaxis = 3)
################### Snowy Rivers ###############################################
pcoa_snowy <- dapc(snowy, n.pca = 10, n.da = 3)
# 2D plot
plot_pca(pcoa_snowy, pt.colors = myCol)

ggsave(
  "pcoa_snowy.pdf",
  width = 4,
  height = 4,
  units = "in",
  dpi = "retina",
  bg = "transparent"
)
# 3D plot
plot_pca(pcoa_snowy, pt.colors = myCol,zaxis = 3)
################### Upper Murray Rivers #######################################
pcoa_murray <- dapc(murray, n.pca = 10, n.da = 3)
# 2D plot
plot_pca(pcoa_murray, pt.colors = myCol)

ggsave(
  "pcoa_murray.pdf",
  width = 4,
  height = 4,
  units = "in",
  dpi = "retina",
  bg = "transparent"
)

# 3D plot
plot_pca(pcoa_murray, pt.colors = myCol,zaxis = 3)
################### Border Rivers ##############################################
pcoa_border <- dapc(border, n.pca = 10, n.da = 3)
# 2D plot
plot_pca(pcoa_border, pt.colors = myCol)

ggsave(
  "pcoa_border.pdf",
  width = 4,
  height = 4,
  units = "in",
  dpi = "retina",
  bg = "transparent"
)

# 3D plot
plot_pca(pcoa_border, pt.colors = myCol,zaxis = 3)

################################################################################
######## Testing contrasting patterns in samples V30 and V32 ###################
################################################################################

# Two samples, each collected in a different river (V30 in Ovens and V32 in
# Mitta Mitta), showed contrasting genetic patterns relative to samples
#collected in the same river.

murray_test <-
  gl.keep.pop(gl.data_original,
              pop.list = c("MITTA_ABOVE",
                           "MITTA_BELOW",
                           "OVENS"))
# filtering on missing data
murray_test <- gl.filter.callrate(murray_test, threshold = 1)
pcoa_murray <- gl.pcoa(murray_test)
gl.pcoa.plot(pcoa_murray, murray_test)

ggsave(
  "V30_V32_PCoA.pdf",
  width = 4,
  height = 4,
  units = "in",
  dpi = "retina",
  bg = "transparent"
)

################## Relatedness analyses ########################################
platy_related <-  gl.filter.callrate(gl.data_original, threshold = 1)
df_rel_platypus <- gl2related(platy_related, save = FALSE)
res_coancestry_platypus <- coancestry(df_rel_platypus, wang = 1)
res_coancestry_2_platypus <-
  res_coancestry_platypus$relatedness[, c(2, 3, 6)]
colnames(res_coancestry_2_platypus) <-
  c("ind1", "ind2", "relatedness")

ids_tmp <- as.data.frame(cbind(indNames(platy_related),
                               as.character(pop(platy_related))))
ids_1 <- ids_tmp
colnames(ids_1) <- c("ind1", "pop1")
ids_2 <- ids_tmp
colnames(ids_2) <- c("ind2", "pop2")

relatedness <- merge(res_coancestry_2_platypus, ids_1, by = "ind1")
relatedness <- merge(relatedness, ids_2, by = "ind2")
relatedness <-
  relatedness[order(relatedness$relatedness, decreasing = TRUE), ]
relatedness <- data.frame(
  ind1 = relatedness$ind1,
  pop1 = relatedness$pop1,
  ind2 = relatedness$ind2,
  pop2 = relatedness$pop2,
  relatedness = relatedness$relatedness
)
# Relatedness analyses identified two pairs of samples in which each pair
# was collected from the same individual (i.e., recaptures; samples T3-T5 and
# T28-T42
head(relatedness)
# Relatedness analyses revealed samples V30 and V32 had closer relatives in the
# opposite river
relatedness_2 <- relatedness[relatedness$relatedness > 0, ]
relatedness_2 <-
  relatedness_2[which(relatedness_2$pop1 != relatedness_2$pop2), ]
V30_sample <- relatedness_2[which(relatedness_2$ind2 == "V30" |
                                    relatedness_2$ind1 == "V30"), ]
V32_sample <- relatedness_2[which(relatedness_2$ind2 == "V32" |
                                    relatedness_2$ind1 == "V32"), ]
head(V30_sample)
head(V32_sample)

