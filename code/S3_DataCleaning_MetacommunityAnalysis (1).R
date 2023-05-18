library("metabaR")
library("tidyverse")
library("reshape2")
library("ggpubr")
library("vegan")
library("adespatial")
library("ggrepel")


## R code for manuscript "Metabarcoding metacommunities: time, space, and 
## land use interact to structure aquatic macroinvertebrate communities in streams"

### PART ONE - DATA CLEANING ###### -----

## Note: Part One Code was adapted from a metabaR tutorial on github: 
## https://github.com/metabaRfactory/metabaR/blob/master/vignettes/metabaRF-vignette.Rmd


## 1.0 Data import -----

inverts <- tabfiles_to_metabarlist(file_reads = "SI_1A_ReadTable.csv",
                                   file_motus = "SI_1B_OTU_taxonomy.csv",
                                   file_pcrs = "SI_1C_PCRs.csv",
                                   file_samples = "SI_1D_SampleMetadata.csv",
                                   files_sep = ",")

## will get a warning about not having index primer sequences but not necessary for downstream analysis


inverts_samples <- subset_metabarlist(inverts,
                                      table = "pcrs",
                                      indices =  is.na(inverts$pcrs$control_type))

inverts_neg <- subset_metabarlist(inverts,
                                  table = "pcrs",
                                  indices =  !is.na(inverts$pcrs$control_type))


summary_metabarlist(inverts)
summary_metabarlist(inverts_samples)
summary_metabarlist(inverts_neg)

## 240 samples, 318 pcrs
## 51334969 reads, 2276
## AVG reads 161430.7 +/- 69285.96 
## 71.53459 +/-  33.0555

## 2.0 Basic data set stats & visualization ---- 

## number of sequence reads per pcr
inverts$pcrs$nb_reads <- rowSums(inverts$reads)

##number of motus (richness) per pcr
inverts$pcrs$nb_motus <- rowSums(inverts$reads>0)

tibble(inverts$pcrs)

summary_metabarlist(inverts)

## 2.1 Checking reads/OTUs in negative controls ----

check1 <- melt(inverts$pcrs[,c("control_type", "nb_reads", "nb_motus")])
check1$Type <- check1$control_type %>% replace_na("samples")

check1$Type = factor(check1$Type, levels = c("extraction", "pcr", "sequencing", "samples"))
labels <- c(nb_reads = "Sequence Reads", nb_motus = "OTUs")

ggplot(data <- check1, aes(x=Type, y=value, color=Type)) + 
  geom_boxplot() + theme_classic(14) + 
  geom_jitter(alpha=0.2) + 
  scale_color_manual(values = c("steelblue", "orange", "darkred", "black")) +
  facet_wrap(~variable, labeller = labeller(variable = labels), scales = "free_y") + 
  theme(axis.text.x = element_text(angle=45, h=1)) +
  ylab("Total reads/OTUs") 

## the negative controls have few sequences & OTUs compared to samples (good! low contamination!)

## 2.2 Assessing sequencing depth ----

## number of reads vs. number of OTUs for controls and samples, lack of correlation 
## between reads and OTUs in samples suggests sequencing depth is sufficient for coverage

inverts$pcrs$Type <- inverts$pcrs$control_type
inverts$pcrs$Type <- inverts$pcrs$Type %>% replace_na("samples")
inverts$pcrs$Type <- factor(inverts$pcrs$Type, levels = c("extraction", "pcr", "sequencing", "samples"))

ggplot(inverts$pcrs, aes(x=nb_reads, y=nb_motus, color = Type)) + 
  geom_point(alpha = 0.5) + theme_classic(14) + 
  scale_y_log10() + scale_x_log10() + 
  scale_color_manual(values = c("steelblue", "orange", "darkred", "black")) +
  xlab("Sequence reads") +
  ylab("OTUs")

##visualization in PCR plate

## 2.3 Visualization of reads in samples and controls in PCR plates ----

inverts$pcrs$plate_no <- factor(inverts$pcrs$plate_no, levels = c("M2019", "J2019", "S2019", "redo"))
labels = c(M2019 = "Plate 1 - May", J2019 = "Plate 2 - July", S2019 = "Plate 3 - Sept", redo = "Plate 4 - Redo")

ggpcrplate(inverts, FUN = function(m){rowSums(m$reads)}, legend_title = "# of reads per PCR") +
  theme_bw(14) +
  scale_fill_manual(values = c("steelblue", "orange", "darkred"), na.translate = F) +
  facet_wrap(~plate_no, labeller = labeller(plate_no = labels), scales = "free_y") +
  guides(fill = guide_legend(title="Control Type"))

## can see low reads in most controls (good!) 

## 2.4 Rarefaction curves for sequencing depth ----

inverts.raref = hill_rarefaction(inverts, nboot = 20, nsteps = 10)
head(inverts.raref$hill_table)

labels = c(D0 = "Species richness", D1 = "Shannon index", D2 = "Inverted Simpson's index", coverage = "Good's coverage index")

gghill_rarefaction(inverts.raref) +
  geom_ribbon(aes(ymin = .data$value - .data$value.sd,
                  ymax = .data$value + .data$value.sd), alpha = 0.1) +
  xlab("Sequence reads") + 
  ylab("Diversity/Coverage index estimate") +
  facet_wrap(~.data$variable, labeller = labeller(variable = labels), ncol = 4, scales = "free_y")

## everything flattens out indicating sequencing depth was sufficient

## 3.0 FLAGGING SPURIOUS SIGNAL ----

## 3.0 FLAGGING SPURIOUS SIGNAL ----

## 3.1 Flagging contaminants ----

## First must subset dataset into plates to evaluate respective neg controls

## 3.1.1 May subset ----

inverts_may <- subset_metabarlist(inverts, 
                                  table = "pcrs",
                                  indices = inverts$pcrs$plate_no == "M2019")


summary_metabarlist(inverts_may) ## data check

inverts_may <- contaslayer(inverts_may,
                           control_types = c("pcr", "sequencing", "extraction"),
                           output_col = "not_a_conta")


##OTU flagged as a contaminant (FALSE in not_a_conta) if relative abundance of OTU across whole dataset is highest in neg controls

table(inverts_may$motus$not_a_conta) ## only 1 OTU flagged as a contaminant!

dt <- inverts_may$motus[!inverts_may$motus$not_a_conta,
                        c("seq_count", "occurrence_count", "Similarity", "Genus_species", "sequence")]

colnames(dt) <- c("total # reads", "total occurrence", "% similarity", "taxonomy", "sequence")
rownames(df) <- NULL

## this should order top ten based on frequency 
## here only one contaminant and is a maple tree (Acer saccharinum)

kable(dt[order(dt[,1], decreasing = TRUE)[1:10],], row.names=F) %>%
  kable_styling(bootstrap_options= c("striped", "hover", "condensed"), 
                font_size = 8, full_width = F)

# Identify the most common contaminant (here is only one, Silver maple!) out of 1301 OTUs
# get contaminant ids
id <- !inverts_may$motus$not_a_conta
max.conta <- rownames(inverts_may$motus[id,])[which.max(inverts_may$motus[id, "seq_count"])]

#and its distribution and relative abundance in each pcr
ggpcrplate(inverts_may, legend_title = "#reads of most \nabundant contaminant",
           FUN = function(m) {m$reads[, max.conta]/rowSums(m$reads)})

# Compute relative abundance of all pcr contaminants together - cannot do if only 1
## a <- data.frame(conta.relab = rowSums(inverts_may$reads[,!inverts_may$motus$not_a_conta]) / 
##                  rowSums(inverts_may$reads))


test <- inverts_may$reads[,!inverts_may$motus$not_a_conta, drop = FALSE] / 
  rowSums(inverts_may$reads)

conta.relab <- rowSums(test)

a <- data.frame(conta.relab)

# flag pcrs with total contaminant relative abundance > 10% of reads)
inverts_may$pcrs$low_contamination_level <- 
  ifelse(a$conta.relab[match(rownames(inverts_may$pcrs), rownames(a))]>1e-1,  F, T)

# Proportion of potentially functional (TRUE) vs. failed (FALSE) pcrs
# (controls included) based on this criterion
table(inverts_may$pcrs$low_contamination_level) / nrow(inverts_may$pcrs)

colnames(inverts_may$pcrs) # data check

## 3.1.2 July subset ----

inverts_july <- subset_metabarlist(inverts, 
                                   table = "pcrs",
                                   indices = inverts$pcrs$plate_no == "J2019")

summary_metabarlist(inverts_july) ## data check

inverts_july <- contaslayer(inverts_july,
                            control_types = c("pcr", "extraction", "sequencing"),
                            output_col = "not_a_conta")

table(inverts_july$motus$not_a_conta) ## NO OTUs flagged out of 1215

dt <- inverts_july$motus[!inverts_july$motus$not_a_conta,
                         c("seq_count", "occurrence_count", "Similarity", "Genus_species", "sequence")]

colnames(dt) <- c("total # reads", "total occurrence", "% similarity", "taxonomy", "sequence")
rownames(df) <- NULL

## this should order top ten based on frequency in a bigger list
## here there are none

kable(dt[order(dt[,1], decreasing = TRUE)[1:10],], row.names=F) %>%
  kable_styling(bootstrap_options= c("striped", "hover", "condensed"), 
                font_size = 8, full_width = F)

# Compute relative abundance of all pcr contaminants together - cannot do if 0
#a <- data.frame(conta.relab = rowSums(inverts_july$reads[,!inverts_july$motus$not_a_conta]) / 
#  rowSums(inverts_july$reads))

test2 <- inverts_july$reads[,!inverts_july$motus$not_a_conta, drop = FALSE] / 
  rowSums(inverts_july$reads)

conta.relab2 <- rowSums(test2)

b <- data.frame(conta.relab2)

# flag pcrs with total contaminant relative abundance > 10% of reads)
inverts_july$pcrs$low_contamination_level <- 
  ifelse(b$conta.relab2[match(rownames(inverts_july$pcrs), rownames(b))]>1e-1,  F, T)

# Proportion of potentially functional (TRUE) vs. failed (FALSE) pcrs
# (controls included) based on this criterion
table(inverts_july$pcrs$low_contamination_level) / nrow(inverts_july$pcrs)

colnames(inverts_july$pcrs)

## 3.1.3 Sept subset ----

inverts_sept <- subset_metabarlist(inverts, 
                                   table = "pcrs",
                                   indices = inverts$pcrs$plate_no == "S2019")


summary_metabarlist(inverts_sept) # data check

##contaminants in pcr controls

inverts_sept <- contaslayer(inverts_sept,
                            control_types = c("pcr", "extraction"),
                            output_col = "not_a_conta")

table(inverts_sept$motus$not_a_conta) ## 3 out of 1375 OTUs flagged

dt <- inverts_sept$motus[!inverts_sept$motus$not_a_conta,
                         c("seq_count", "occurrence_count", "Similarity", "Phylum", "Genus_species", "sequence")]

colnames(dt) <- c("total # reads", "total occurrence", "% similarity", "Phylum", "taxonomy", "sequence")
rownames(df) <- NULL

## OTUs flagged as contaminants match to human, dog and diatom
kable(dt[order(dt[,1], decreasing = TRUE)[1:10],], row.names=F) %>%
  kable_styling(bootstrap_options= c("striped", "hover", "condensed"), 
                font_size = 8, full_width = F)

# Compute relative abundance of all pcr contaminants together
c <- data.frame(conta.relab3 = rowSums(inverts_sept$reads[,!inverts_sept$motus$not_a_conta]) / rowSums(inverts_sept$reads))

# flag pcrs with total contaminant relative abundance > 10% of reads)
inverts_sept$pcrs$low_contamination_level <- 
  ifelse(c$conta.relab3[match(rownames(inverts_sept$pcrs), rownames(c))]>1e-1,  F, T)

# Proportion of potentially functional (TRUE) vs. failed (FALSE) pcrs
# (controls included) based on this criterion
table(inverts_sept$pcrs$low_contamination_level) / nrow(inverts_sept$pcrs)

## 3.1.4 Redo subset ----

inverts_redo <- subset_metabarlist(inverts, 
                                   table = "pcrs",
                                   indices = inverts$pcrs$plate_no == "redo")

summary_metabarlist(inverts_redo) ##data check

##contaminants in controls

inverts_redo <- contaslayer(inverts_redo,
                            control_types = c("pcr", "extraction", "sequencing"),
                            output_col = "not_a_conta")

table(inverts_redo$motus$not_a_conta) ## 2 out of 734 OTUs flagged

dt <- inverts_redo$motus[!inverts_redo$motus$not_a_conta,
                         c("seq_count", "occurrence_count", "Similarity", "Phylum", "Genus_species", "sequence")]

colnames(dt) <- c("total # reads", "total occurrence", "% similarity", "Phylum", "taxonomy", "sequence")
rownames(df) <- NULL

## this should order top ten based on frequency in a bigger list
## here only 2 contaminants (Ascomycota)

kable(dt[order(dt[,1], decreasing = TRUE)[1:10],], row.names=F) %>%
  kable_styling(bootstrap_options= c("striped", "hover", "condensed"), 
                font_size = 8, full_width = F)

# Compute relative abundance of all pcr contaminants together - cannot do if only 1
d <- data.frame(conta.relab4 = rowSums(inverts_redo$reads[,!inverts_redo$motus$not_a_conta]) / rowSums(inverts_redo$reads))

# flag pcrs with total contaminant relative abundance > 10% of reads)
inverts_redo$pcrs$low_contamination_level <- 
  ifelse(d$conta.relab4[match(rownames(inverts_redo$pcrs), rownames(d))]>1e-1,  F, T)

# Proportion of potentially functional (TRUE) vs. failed (FALSE) pcrs
# (controls included) based on this criterion
table(inverts_redo$pcrs$low_contamination_level) / nrow(inverts_redo$pcrs)

## 3.1.5 Combine back together ----

## no sample pcrs were flagged as contaminanted so can combine this way

inverts$pcrs$low_contamination_level <- NA
inverts$pcrs[rownames(inverts_may$pcrs),"low_contamination_level"] <- inverts_may$pcrs$low_contamination_level

inverts$pcrs$low_contamination_level ## it worked, repeat for july/sept/redo

inverts$pcrs[rownames(inverts_july$pcrs),"low_contamination_level"] <- inverts_july$pcrs$low_contamination_level
inverts$pcrs[rownames(inverts_sept$pcrs),"low_contamination_level"] <- inverts_sept$pcrs$low_contamination_level
inverts$pcrs[rownames(inverts_redo$pcrs),"low_contamination_level"] <- inverts_redo$pcrs$low_contamination_level

## repeat for OTU dataframe

inverts$motus$not_a_conta <- NA
inverts$motus[rownames(inverts_may$motus),"not_a_conta"] <- inverts_may$motus$not_a_conta
inverts$motus[rownames(inverts_july$motus),"not_a_conta"] <- inverts_july$motus$not_a_conta
inverts$motus[rownames(inverts_sept$motus),"not_a_conta"] <- inverts_sept$motus$not_a_conta
inverts$motus[rownames(inverts_redo$motus),"not_a_conta"] <- inverts_redo$motus$not_a_conta

inverts$motus$not_a_conta

## 3.2 Flagging non-target OTUs & low quality matches ----

colnames(inverts$pcrs)

inverts$motus$target_taxon <- grepl("Arthropoda|Mollusca|Annelida", inverts$motus$Phylum)

# Proportion of each of these over total number of MOTUs
table(inverts$motus$target_taxon) / nrow(inverts$motus) ### 87.67 OTUs are benthic inverts

inverts$motus$Similarity <- as.numeric(inverts$motus$Similarity)
inverts$motus$Similarity[is.na(inverts$motus$Similarity)] <- 0 ## replace NAs (OTUs with no match) as 0% similarity

# Plot the unweighted distribution OTU sequence similarity scores 
a <- ggplot(inverts$motus, aes(x=Similarity)) + 
  geom_histogram(color="black", fill="grey", bins=20) + 
  geom_vline(xintercept = 90, col="red", lty=2) + ## using 90% but can change
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(x="Percent similarity against best match", y="Number of OTUs") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1250))

# Same for the weighted distribution based on total sequence reads of OTU
b <- ggplot(inverts$motus, 
            aes(x=Similarity, weight = seq_count)) + 
  geom_histogram(color="black", fill="grey", bins=20) + 
  geom_vline(xintercept = 90, col="red", lty=2) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(x="Percent similarity against best match", y="Number of sequence reads") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 50000000))

ggarrange(a,b)


inverts$motus$not_degraded <-
  ifelse(inverts$motus$Similarity < 90, F, T) ### flagging OTUs with under 90% match

# Proportion of each of these over total number of MOTUs
table(inverts$motus$not_degraded) / nrow(inverts$motus) ### 78.7% of OTUs had a match of at least 90

## Intersection with target taxa and sequence matching
table(inverts$motus$target_taxon, 
      inverts$motus$not_degraded)

## 3.3 Flagging pcr outliers based on sequencing depth ----

##flagged based on sequencing depth (includes controls)

ggplot(inverts$pcrs, aes(nb_reads)) +
  geom_histogram(bins=40, color="black", fill="grey") + 
  geom_vline(xintercept = 87344, lty=2, color="red") + ## threshold is average sequencing depth minus 1 SD of first 3 plates
  scale_x_log10() + 
  labs(x="Number of sequence reads", 
       y="Number of samples") +
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 120))


inverts$pcrs$seqdepth_ok <- ifelse(inverts$pcrs$nb_reads < 87344, F, T) ## using one less than SD, 87344, 87.1% of pcrs pass

# Proportion of each of these over total number of pcrs, control excluded
table(inverts$pcrs$seqdepth_ok[inverts$pcrs$type=="sample"]) /
  nrow(inverts$pcrs[inverts$pcrs$type=="sample",]) 

# Flag pcrs with an acceptable sequencing depth (TRUE) or inacceptable one (FALSE) 
# dataset in package set at 1000 (seems way too low for my data, picking 50,000 here?)

inverts$pcrs$seqdepth_ok <- ifelse(inverts$pcrs$nb_reads <   87344, F, T) ## using one less than SD, 87344, 87.1% of pcrs pass

# Proportion of each of these over total number of pcrs, control excluded
table(inverts$pcrs$seqdepth_ok[inverts$pcrs$type=="sample"]) /
  nrow(inverts$pcrs[inverts$pcrs$type=="sample",]) 


## 3.4 Assessing reproducibility via technical replicates ----

## first need to remove samples without tech reps, pcrs with no reads and neg controls
##subset data


inverts$pcrs$tech_rep

inverts_sub <- subset_metabarlist(inverts, 
                                  table = "pcrs",
                                  indices = !is.na(inverts$pcrs$tech_rep))

inverts_sub$pcrs$tech_rep

inverts_sub$pcrs$sample_id

inverts_sub$pcrs$nb_motus

pcr_within_between(
  inverts_sub,
  replicates = inverts_sub$pcrs$sample_id,
  FUN = FUN_pcrdist_bray_freq,
  method = "centroid"
)

pcrslayer(
  inverts_sub,
  replicates = NULL,
  method = "centroid",
  FUN = FUN_pcrdist_bray_freq,
  thresh.method = "intersect",
  output_col = "functional_pcr",
  plot = T
)

#visualization

comp1 <- pcr_within_between(inverts_sub) 
## detecting PCR replicate outliers based on PCR similarity/reproducibility
## uses bray distance to compare dissimilarities within sample and between samples

check_pcr_thresh(comp1) +
  xlab("Bray-Curtis dissimilarity") +
  ylab("Density") +
  scale_color_manual(values = c("steelblue", "orange")) +
  geom_line(size = 1) +
  theme_classic(14)

### very low dissimilarity/distance within samples indicates that tech reps were very consistent


## flags PCR as failed if PCR replicates are outliers by comparing dissimilarities 
## in taxonomic composition within PCR reps vs between samples
## threshold for elimination having distance within samples greater than threshold of intersection

inverts_sub <- pcrslayer(inverts_sub, output_col = "replicating_pcr", method = "pairwise", plot = F) ## flags PCR as failed 

table(inverts_sub$pcrs$replicating_pcr) /
  nrow(inverts_sub$pcrs)

## report flagging in initial metabarlist 
inverts$pcrs$replicating_pcr <- TRUE
inverts$pcrs[rownames(inverts_sub$pcrs),"replicating_pcr"] <- inverts_sub$pcrs$replicating_pcr ## would replace with false if any

inverts$pcrs$replicating_pcr #check


## 3.4 Lowering tag jumps ----

thresholds <- c(0,1e-4,1e-3, 1e-2, 3e-2, 5e-2) 

tests <- lapply(thresholds, function(x) tagjumpslayer(inverts,x))
names(tests) <- paste("t_", thresholds, sep="")

tmp <- melt(as.matrix(do.call("rbind", lapply(tests, function(x) rowSums(x$reads)))))
colnames(tmp) <- c("threshold", "sample", "Sequences")

tmp$OTUs <-
  melt(as.matrix(do.call("rbind", lapply(tests, function(x) {
    rowSums(x$reads > 0)
  }))))$value

# Add control type information on pcrs and make data curation threshold numeric
tmp$controls <- inverts$pcrs$control_type[match(tmp$sample, rownames(inverts$pcrs))]
tmp$threshold <- as.numeric(gsub("t_", "", tmp$threshold))

tmp2 <- melt(tmp, id.vars=colnames(tmp)[-grep("Sequences|OTUs", colnames(tmp))])

unique(tmp2$variable)

tmp2$controls <- tmp2$controls %>% replace_na("samples")
tmp2$controls <- factor(tmp2$controls, levels = c("extraction", "pcr", "sequencing", "samples"))


ggplot(tmp2, aes(x=as.factor(threshold), y=value)) + 
  geom_boxplot(color="black") + 
  geom_vline(xintercept = which(levels(as.factor(tmp2$threshold)) == "0.001"), col="red", lty=2) + 
  geom_jitter(aes(color=controls), width = 0.2, alpha=0.5) + 
  scale_color_manual(values = c("steelblue", "orange", "darkred", "darkgrey")) +
  facet_wrap(~variable + controls, scale="free", ncol=4) + 
  theme_bw() + 
  labs(x="Filtering threshold", y="Reads or OTUs") + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle=40, h=1), 
        legend.position = "none")

## 4.0 Summarizing noise in dataset ----

colnames(inverts$motus)

motus.qual <- !inverts$motus[,c("not_a_conta", "target_taxon", "not_degraded")]
colnames(motus.qual) <- c("Contaminant", "Untargeted", "LowSimilarity")

prop.table(table(apply(motus.qual, 1, sum) > 0)) ###25% OTUs not high quality and/or target taxa

# Proportion of MOTUs and reads artifcat
apply(motus.qual, 2, sum) / nrow(motus.qual)
apply(motus.qual, 2, function(x) sum(inverts$motus$count[x])/sum(inverts$motus$count))

tmp.motus <- 
  apply(sapply(1:ncol(motus.qual), function(x) {
    ifelse(motus.qual[,x]==T, colnames(motus.qual)[x], NA)}), 1, function(x) {
      paste(sort(unique(x)), collapse = "|")
    })
tmp.motus <- as.data.frame(gsub("^$", "Retained", tmp.motus))
colnames(tmp.motus) <-  "artefact_type"

unique(tmp.motus$artefact_type)
count(tmp.motus$artefact_type$Contaminant&LowSimilarity&Untargeted)

ggplot(tmp.motus, aes(x = 1, fill=artefact_type)) +
  geom_bar() + 
  labs(fill="Artifact type") + 
  scale_fill_manual(values = c("black", "black", "black", "darkred", "orange", "steelblue", "darkgrey")) +
  theme_classic(14) + 
  xlim(0, 2) +
  ylab("OTU count") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 2300)) +
  theme(legend.position = "right") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Create a table of pcrs quality criteria 
pcrs.qual <- !inverts$pcrs[,c("low_contamination_level", "seqdepth_ok", "replicating_pcr")]
colnames(pcrs.qual) <- c("high_contamination_level", "low_seqdepth", "techrep_outliers")

colnames(inverts$pcrs)
# Proportion of pcrs potentially artifactual (TRUE) based on the criteria used
# excluding controls
prop.table(table(apply(pcrs.qual[inverts$pcrs$type=="sample",], 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(pcrs.qual[inverts$pcrs$type=="sample",], 2, sum) / nrow(pcrs.qual[inverts$pcrs$type=="sample",])

tmp.pcrs <- 
  apply(sapply(1:ncol(pcrs.qual), function(x) {
    ifelse(pcrs.qual[inverts$pcrs$type=="sample",x]==T, 
           colnames(pcrs.qual)[x], NA)}), 1, function(x) {
             paste(sort(unique(x)), collapse = "|")
           })
tmp.pcrs <- as.data.frame(gsub("^$", "Retained", tmp.pcrs))

colnames(tmp.pcrs) <- "artefact_type"

ggplot(tmp.pcrs, aes(x = 1, fill=artefact_type)) +
  geom_bar() + 
  labs(fill="Artifact type") + 
  scale_fill_manual(values = c("darkred", "steelblue")) +
  theme_classic(14) + 
  xlim(0, 2) +
  ylab("PCRs in dataset") +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "right") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

## 5.0 Data cleaning and aggregation  ----
## 5.1 Removing spurious signals (flags) ----

##subsetting data to remove negs

inverts_s <- subset_metabarlist(inverts, 
                                table = "pcrs",
                                indices =  inverts$pcrs$nb_reads > 0 & (
                                  is.na(inverts$pcrs$control_type)
                                ))

# Subset on MOTUs: we keep motus that are defined as TRUE following the 
# three criteria below (sum of three TRUE is equal to 3 with the rowSums function)

tmp <- tests[["t_0.001"]] ## tagjump threshold we selected

inverts_s <- subset_metabarlist(tmp, 
                                table = "pcrs",
                                indices =  tmp$pcrs$nb_reads > 0 & (
                                  is.na(tmp$pcrs$control_type)
                                ))

tmp$motus$target_taxon

tmp2 <- subset_metabarlist(inverts_s, "motus", 
                           indices = rowSums(inverts_s$motus[,c("not_a_conta", "target_taxon",
                                                                "not_degraded")]) == 3)
unique(tmp2$motus$Phylum)

summary_metabarlist(tmp) ## sanity check

inverts_clean <- subset_metabarlist(tmp2, "pcrs", 
                                    indices = rowSums(tmp2$pcrs[,c("low_contamination_level", 
                                                                   "seqdepth_ok", "replicating_pcr")]) == 3)


summary_metabarlist(inverts_clean) ## removed 3 samples due to sequencing depth

unique(inverts_clean$motus$Phylum)

#check is this led to empty pcrs or otus
if(sum(colSums(inverts_clean$reads)==0)>0){print("empty motus present")}
if(sum(colSums(inverts_clean$reads)==0)>0){print("empty pcrs present")}

#update parameters with removed motus or reduced read counts
inverts_clean$pcrs$nb_reads_postmetabaR = rowSums(inverts_clean$reads)
inverts_clean$pcrs$nb_motus_postmetabaR = rowSums(ifelse(inverts_clean$reads>0, T, F))

#compare basic characteristics before and after data curation
check <- melt(inverts_clean$pcrs[,c("nb_reads", "nb_reads_postmetabaR", 
                                    "nb_motus", "nb_motus_postmetabaR")])
check$type <- ifelse(grepl("motus", check$variable), "OTUs", "Sequences")

tibble(check)

ggplot(data = check, aes(x = variable, y = value)) +
  geom_boxplot(colour = "black") +
  geom_jitter(alpha=0.5, aes(colour = variable)) +
  scale_colour_manual(values = c("orange", "steelblue", "orange", "steelblue")) +
  theme_bw(14) +
  ylab("Number of OTUs or Sequences") +
  facet_wrap(~type, scales = "free") +
  theme(axis.text.x = element_text(angle=45, h=1)) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("nb_motus" = "Pre", "nb_motus_postmetabaR" = "Post", "nb_reads" = "Pre", "nb_reads_postmetabaR" = "Post")) +
  theme(axis.title.x=element_blank())
## not a big difference here

## clean sequences per plate to check number of sequences per run

may_clean <- subset_metabarlist(inverts_clean, 
                                table = "pcrs",
                                indices = inverts_clean$pcrs$plate_no == "M2019")

summary_metabarlist(may_clean)

july_clean <- subset_metabarlist(inverts_clean, 
                                 table = "pcrs",
                                 indices = inverts_clean$pcrs$plate_no == "J2019")

summary_metabarlist(july_clean)

sept_clean <- subset_metabarlist(inverts_clean, 
                                 table = "pcrs",
                                 indices = inverts_clean$pcrs$plate_no == "S2019")

summary_metabarlist(sept_clean)

redo_clean <- subset_metabarlist(inverts_clean, 
                                 table = "pcrs",
                                 indices = inverts_clean$pcrs$plate_no == "redo")

summary_metabarlist(redo_clean)

## 5.2 Combining tech reps (sum otus across reps) and save final dataset ----

inverts_clean$pcrs$nb_motus_postmetabaR

inverts_all <- aggregate_pcrs(inverts_clean, FUN = FUN_agg_pcrs_mean)

inverts_all$pcrs$nb_motus_postmetabaR

tibble(inverts_all$pcrs)

## number of sequence reads per pcr
inverts_all$pcrs$total_reads <- rowSums(inverts_all$reads)

##number of motus (richness) per pcr
inverts_all$pcrs$total_motus <- rowSums(inverts_all$reads>0)

inverts_all$pcrs$nb_reads_postmetabaR
inverts_all$pcrs$total_reads

inverts_all$pcrs$nb_motus_postmetabaR
inverts_all$pcrs$total_motus

summary_metabarlist(inverts_all) 



### 1681 total OTUs (avg 65.05063 +/ 29.35445)
### 41978040 total reads (avg 177122.5 +/ 46880.63)

saveRDS(inverts_all, file = "Inverts_clean_dataset.rds")

### PART TWO - METACOMMUNITY ANALYSIS ###### -----

### 1.0 LOADING DATA ----

##inverts_all <- readRDS("Inverts_clean_dataset.rds") ## read in end files from Part One or continue
env_all <- read.csv("SI_2_Env_metadata.csv") ## environmetal metadata

summary_metabarlist(inverts_all) ## summary of sequences/OTUs

unique(inverts_all$pcrs$control_type) ## double checking no controls in dataset for analysis


## 1.1 Invert OTUs aggregation of bio reps ----

## field reps need to be combined

inverts_all$samples$Site_month <- paste(inverts_all$samples$Site_code, inverts_all$samples$Month)

?aggregate

inverts_all_ag <- aggregate(inverts_all$reads, list(inverts_all$samples$Site_month), FUN = mean)
tibble(inverts_all_ag)

inverts_all_ag <- column_to_rownames(inverts_all_ag, var = "Group.1") #cleanup
inverts_all_agPA <- decostand(inverts_all_ag, "pa") ##convert to presence/absence

tibble(inverts_all_ag)

#1681 OTUs total

### Subset by month

inverts_may <- subset_metabarlist(inverts_all, 
                                  table = "samples",
                                  indices = inverts_all$samples$Month == "May")

tibble(inverts_may$reads)
## 20 sites, 1077 OTUs

inverts_may_ag <- aggregate(inverts_may$reads, list(inverts_may$samples$Site_code), FUN = mean)
tibble(inverts_may_ag)

inverts_july <- subset_metabarlist(inverts_all, 
                                   table = "samples",
                                   indices = inverts_all$samples$Month == "July")


tibble(inverts_july$reads)


inverts_july_ag <- aggregate(inverts_july$reads, list(inverts_july$samples$Site_code), FUN = mean)
tibble(inverts_july_ag)

## 20 sites, 968 OTUs

inverts_sept <- subset_metabarlist(inverts_all, 
                                   table = "samples",
                                   indices = inverts_all$samples$Month == "Sept")


tibble(inverts_sept$reads)

inverts_sept_ag <- aggregate(inverts_sept$reads, list(inverts_sept$samples$Site_code), FUN = mean)
tibble(inverts_sept_ag)

## 20 sites, 1040 OTUs

## 1.2 Data transformations & richness -----

inverts_may_ag <- column_to_rownames(inverts_may_ag, var = "Group.1")
inverts_may_agPA <- decostand(inverts_may_ag, "pa")

tibble(inverts_may_ag)

r1 <- rowSums(inverts_may_agPA)
mean(r1) ## 148.65 avg OTU per site
sd(r1) ## 58.26

inverts_july_ag <- column_to_rownames(inverts_july_ag,  var = "Group.1")
inverts_july_agPA <- decostand(inverts_july_ag, "pa")

tibble(inverts_july_agPA)
r2 <- rowSums(inverts_july_agPA)
mean(r2) ## 144.75 avg OTU per site
sd(r2) ## 45.7

inverts_sept_ag <- column_to_rownames(inverts_sept_ag, var = "Group.1")
inverts_sept_agPA <- decostand(inverts_sept_ag, "pa")

tibble(inverts_sept_agPA)
r3 <- rowSums(inverts_sept_agPA)
mean(r3) ## 162 avg OTU per site
sd(r3) ##42.24

## Env vectors by month

env_may <- env_all %>% filter(Month == "May")
env_july <- env_all %>% filter(Month == "July")
env_sept <- env_all %>% filter(Month == "Sept")


### 2.0 RICHNESS OVER TIME ----

r <- cbind(r1, r2, r3)
r <- as.data.frame(r)


r <- gather(r, Month, Richness, r1:r3)

df <- data.frame("Site" = env_all$Site_name, "Month" = env_all$Month)
df
df$Richness <- r$Richness
df

rmod <- aov((Richness~Month + Error(Site)), data = df)
summary(rmod)

## Figure 2 in manuscript 

df$Month = factor(df$Month, levels = c("May", "July", "Sept"))

ggplot(df, aes(x = Month, y = Richness)) +
  geom_boxplot(alpha = 0.8, aes(fill = Month)) +
  geom_jitter(position = position_dodge(0.8)) +
  geom_line(aes(group = Site), color = "grey", alpha = 0.8) +
  theme_classic(14) +
  ylab("Richness") +
  guides(shape = guide_legend(title="Month")) +
  guides(fill = guide_legend(title="Month")) +
  theme(legend.position = "right") +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_fill_manual(values = c("steelblue", "orange", "darkred")) 

##ggsave("Fig2_Richness.jpg", plot = last_plot())

### 3.0 GENERALISTS & RARE TAXA ---- 

## generalists described as occurring in over 70% of sites (14 out of 20)

may_gen <- inverts_may_agPA[,colSums(inverts_may_agPA) >= 14]

july_gen <- inverts_july_agPA[,colSums(inverts_july_agPA) >= 14]

sept_gen <- inverts_sept_agPA[,colSums(inverts_sept_agPA) >= 14]

may_gen
july_gen
sept_gen

## calculating percentage of OTUs that are generalists for each month
ncol(may_gen)/ncol(inverts_may_agPA)*100  ## 1.02%
ncol(july_gen)/ncol(inverts_july_agPA)*100 ## 0.62%
ncol(sept_gen)/ncol(inverts_sept_agPA)*100 ## 1.15%



### Calculating the percentage of rare OTUs, only occur in one site

may_rare <- inverts_may_agPA[,colSums(inverts_may_agPA) == 1]

july_rare <- inverts_july_agPA[,colSums(inverts_july_agPA) == 1]

sept_rare <- inverts_sept_agPA[,colSums(inverts_sept_agPA) == 1]

ncol(may_rare)/ncol(inverts_may_agPA)*100 ## 46.74
ncol(july_rare)/ncol(inverts_july_agPA)*100 ## 44.36
ncol(sept_rare)/ncol(inverts_sept_agPA)*100 ## 44.21


### 4.0 PERMANOVAS ----

a1 <- adonis2(inverts_all_agPA ~ env_all$Site_name, permutations = 999, method = "jaccard")
a1

a2 <- adonis2(inverts_all_agPA ~ env_all$Ag_watershed*env_all$Month, permutations = 999, method = "jaccard")
a2


### 5.0 Forward selection of env variables  ----


### Choosing adjustedRsq 

##combined dataset
test_rda <- rda(inverts_all_agPA ~ . + env_all$Month, data = env_all[,4:15], distance = jaccard)
summary(test_rda)
RsquareAdj(test_rda)

## May dataset
test_rda <- rda(inverts_may_agPA ~ ., data = env_may[,4:15], distance = jaccard)
summary(test_rda)
RsquareAdj(test_rda)


##July dataset
test_rda <- rda(inverts_july_agPA ~ ., data = env_july[,4:15], distance = jaccard)
summary(test_rda)
RsquareAdj(test_rda)

##Sept dataset
test_rda <- rda(inverts_sept_agPA ~ ., data = env_sept[,4:15], distance = jaccard)
summary(test_rda)
RsquareAdj(test_rda)



summary(env_all)

## only water quality & local site variables 
wq_may <- env_may[,7:14]
wq_july <- env_july[,7:14]
wq_sept <- env_sept[,7:14]
wq_all <- env_all[,7:14]


## forward selection of local variables
for.all <- forward.sel(inverts_all_agPA, wq_all, nperm = 999, alpha = 0.05)
for.all

for.may <- forward.sel(inverts_may_agPA, wq_may, nperm = 999, alpha = 0.05, adjR2thresh = 0.168)
for.may

for.july <- forward.sel(inverts_july_agPA, wq_july, nperm = 999, alpha = 0.05, adjR2thresh = 0.113)
for.july

for.sept <- forward.sel(inverts_sept_agPA, wq_sept, nperm = 999, alpha = 0.05, adjR2thresh = 0.149)
for.sept

## repeat forward selection for spatial variables

sp_may <- env_may[,4:6]
sp_may

sp_july <- env_july[,4:6]
sp_sept <- env_sept[,4:6]
sp_all <- env_all[,4:6]

for.all2 <- forward.sel(inverts_all_agPA, sp_all, nperm = 999, alpha = 0.05)
for.all2 ## no variables selected

for.may2 <- forward.sel(inverts_may_agPA, sp_may, nperm = 999, alpha = 0.05)
for.may2

for.july2 <- forward.sel(inverts_july_agPA, sp_july, nperm = 999, alpha = 0.05)
for.july2

for.sept2 <- forward.sel(inverts_sept_agPA, sp_sept, nperm = 999, alpha = 0.05)
for.sept2

### 6.0  RDA all data ----

all_rda <- rda(inverts_all_agPA ~ 
                 env_all$Month, distance = jaccard)

all_rda
anova(all_rda, permutations = 999)

ii=summary(all_rda)  #View analysis results
ii

st=as.data.frame(ii$sites[,1:2])*0.5
yz=as.data.frame(ii$biplot[,1:2])
Month=env_all$Month
Site=env_all$Site_name

Month <- factor(Month, levels = c("May", "July", "Sept"))

## Figure 3 in manuscript

ggplot() +
  ##geom_text_repel(data = st,aes(RDA1,RDA2,label=row.names(st)),size=4)+#Show a Square
  geom_point(data = st,aes(RDA1,RDA2, colour = Month, shape = Month),size=4)+
  ##geom_segment(data = yz[1:2,],aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
              ## arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                 ##            type = "closed"),linetype=1, size=0.6,colour = "grey13")+
  stat_ellipse(data = st, aes(x = RDA1, y = RDA2, colour = Month), level = 0.95, size = 1) +
  ##geom_text_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)))+
  scale_colour_manual(values = c("steelblue", "orange", "darkred")) +
  labs(x=paste("RDA 1 (", format(100 *ii$cont[[1]][2,1], digits=3), "%)", sep=""),
       y=paste("RDA 2 (", format(100 *ii$cont[[1]][2,2], digits=3), "%)", sep=""))+
  theme_classic(14) +
  coord_equal()


### 7.0 VARIATION PARTITIONING ----

## dissimilarity matrices
x1 <- vegdist(inverts_may_agPA, method = "jaccard")
x2 <- vegdist(inverts_july_agPA, method = "jaccard")
x3 <- vegdist(inverts_sept_agPA, method = "jaccard")

## add sig variables from forward selection of local and spatial parameters
v.may <- varpart(x1, ~env_may$Ag_watershed, ~ env_may$Conductivity_uS.cm + env_may$Buffer, 
                 ~env_may$Latitude + env_may$Longitude + env_may$StreamOrder)

v.july <- varpart(x2, ~env_july$Ag_watershed, ~env_july$Conductivity_uS.cm + env_july$FDOM_RFU,
                  ~env_july$Latitude + env_july$Longitude + env_july$StreamOrder)

v.sept <- varpart(x3, ~env_sept$Ag_watershed, ~env_sept$ODO_. + env_sept$FDOM_RFU, 
                  ~env_sept$Latitude + env_sept$Longitude + env_sept$StreamOrder)

## rough plots
plot(v.may)
plot(v.july)
plot(v.sept)


## Figure 4 in manuscript

##panel A

plot(v.may, digits = 2, Xnames = c('Regional Agriculture', 'Local Habitat', 'Space'), 
     bg =  c("steelblue", "orange", "darkred"))

##panel B
plot(v.july, digits = 2, Xnames = c('Regional Agriculture', 'Local Habitat', 'Space'), 
     bg =  c("steelblue", "orange", "darkred"))

##panel C
plot(v.sept, digits = 2, Xnames = c('Regional Agriculture', 'Local Habitat', 'Space'),
     bg =  c("steelblue", "orange", "darkred"))




### 8.0 RDA by Month ----

## 8.1 May ----

## Figure 5 in Manuscript

##including regional agriculture plus relevant spatial and local variables
may_rda <- rda(inverts_may_agPA ~ env_may$Ag_watershed + env_may$Conductivity_uS.cm + env_may$Buffer + 
                 env_may$StreamOrder + env_may$Longitude + env_may$Latitude, distance = jaccard)

may_rda
anova(may_rda, permutations = 999)

plot(may_rda)

ii=summary(may_rda)  #View analysis results
ii

sp=as.data.frame(ii$species[,1:2])   #Depending on the drawing result, the drawing data can be enlarged or reduced to a certain extent, as follows
st=as.data.frame(ii$sites[,1:2])*0.5
yz=as.data.frame(ii$biplot[,1:2])

## Panel A (arrow labels cleaned for legibility outside of R, see geom_text_repel removal)
ggplot() +
  geom_text_repel(data = st,aes(RDA1,RDA2,label=row.names(st)),size=4)+#Show a Square
  geom_point(data = st,aes(RDA1,RDA2, colour = grp),size=4, colour = "steelblue")+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, size=0.6,colour = "grey13")+
  ##geom_text_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)))+
  labs(x=paste("RDA 1 (", format(100 *ii$cont[[1]][2,1], digits=4), "%)", sep=""),
       y=paste("RDA 2 (", format(100 *ii$cont[[1]][2,2], digits=4), "%)", sep=""))+
  theme_classic(14) +
  coord_equal()

## 8.2 July ----

tibble(env_july)

july_rda <- rda(inverts_july_agPA ~ env_july$Ag_watershed + env_july$Conductivity_uS.cm + env_july$FDOM_RFU +
                  env_july$StreamOrder + env_july$Longitude + env_july$Latitude, distance = jaccard)

july_rda

anova(july_rda, permutations = 999)

ii=summary(july_rda)  #View analysis results
ii

sp=as.data.frame(ii$species[,1:2])   #Depending on the drawing result, the drawing data can be enlarged or reduced to a certain extent, as follows
st=as.data.frame(ii$sites[,1:2])*0.5
yz=as.data.frame(ii$biplot[,1:2])

## Panel B
ggplot() +
  geom_text_repel(data = st,aes(RDA1,RDA2,label=row.names(st)),size=4)+#Show a Square
  geom_point(data = st,aes(RDA1,RDA2, colour = grp),size=4, colour = "orange")+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, size=0.6,colour = "grey13")+
  ## geom_text_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)))+
  labs(x=paste("RDA 1 (", format(100 *ii$cont[[1]][2,1], digits=4), "%)", sep=""),
       y=paste("RDA 2 (", format(100 *ii$cont[[1]][2,2], digits=4), "%)", sep=""))+
  theme_classic(14) +
  coord_equal()

ggsave("July_RDA.jpg", plot = last_plot())


## 8.3 Sept ----

sept_rda <- rda(inverts_sept_agPA ~ env_sept$Ag_watershed + env_sept$ODO_. + env_sept$FDOM_RFU +
                  env_sept$StreamOrder + env_sept$Longitude + env_sept$Latitude, distance = jaccard)
sept_rda
anova(sept_rda, permutations = 999)

ii=summary(sept_rda)  #View analysis results
ii

sp=as.data.frame(ii$species[,1:2])   #Depending on the drawing result, the drawing data can be enlarged or reduced to a certain extent, as follows
st=as.data.frame(ii$sites[,1:2])*0.5
yz=as.data.frame(ii$biplot[,1:2])

## Panel C
ggplot() +
  geom_text_repel(data = st,aes(RDA1,RDA2,label=row.names(st)),size=4)+#Show a Square
  geom_point(data = st,aes(RDA1,RDA2, colour = grp),size=4, colour = "darkred")+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, size=0.6,colour = "grey13")+
  ##geom_text_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)))+
  labs(x=paste("RDA 1 (", format(100 *ii$cont[[1]][2,1], digits=4), "%)", sep=""),
       y=paste("RDA 2 (", format(100 *ii$cont[[1]][2,2], digits=4), "%)", sep=""))+
  theme_classic(14) +
  coord_equal()


## 9.0 Local contribution to beta diversity -----

BD_may  <- beta.div(inverts_may_agPA, method =  "jaccard", sqrt.D = FALSE, samp = FALSE, nperm = 999)
BD_july  <- beta.div(inverts_july_agPA, method =  "jaccard", sqrt.D = FALSE, samp = FALSE, nperm = 999)
BD_sept <- beta.div(inverts_sept_agPA, method =  "jaccard", sqrt.D = FALSE, samp = FALSE, nperm = 999)

BD_may$LCBD
BD_july$LCBD
BD_sept$LCBD

LBD <- as.data.frame(cbind(BD_may$LCBD, BD_july$LCBD, BD_sept$LCBD))
LBD

LBD <- rename(LBD, c("May"="V1", "July" = "V2", "Sept" = "V3"))

LBD <- LBD %>% gather(Month, LCBD, May:Sept)
LBD

LBD$Month <- factor(LBD$Month, levels = c("May", "July", "Sept"))

LBD$Ag <- env_may$Ag_watershed

## Figure 6 in manuscript

ggplot(LBD, aes(x = Ag, y = LCBD)) +
  geom_point(size = 3, aes(shape = Month, fill = Month)) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE, aes(colour = Month)) +
  labs(x = "Percent Agriculture",
       y = "LCBD") +
  theme_classic(14) +
  theme(legend.position = "right") +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  facet_grid(cols = vars(Month))  +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_colour_manual(values = c("steelblue", "orange", "darkred")) +
  scale_fill_manual(values = c("steelblue", "orange", "darkred")) +
  geom_hline(yintercept = 0.05) +
  guides(colour = guide_legend(title="Month")) +
  guides(shape = guide_legend(title="Month")) +
  guides(fill = guide_legend(title="Month"))


summary(lm(BD_may$LCBD ~ env_may$Ag_watershed))
summary(lm(BD_july$LCBD ~ env_may$Ag_watershed))
summary(lm(BD_sept$LCBD ~ env_may$Ag_watershed))


## END


