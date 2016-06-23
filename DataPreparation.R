setwd("/Users/Mehdi/Desktop/Bioinformatics/FinalProject")

# read MACCS data and set pert_id as rowname
maccs <- read.csv("MACCS_bitmatrix.csv", header = T)
rownames(maccs) <- maccs$pert_id
maccs <- maccs[, -1]

# read LINCS_GENE_Expression data and set pert_id as rowname
gene.data <- read.csv("LINCS_Gene_Experssion_signatures_CD.csv", header = T)
rownames(gene.data) <- gene.data$pert_id
pert_ids <- gene.data$pert_id
gene.data <- gene.data[, -1]

# read FAERS_offsides_PTs data and set pert_id as rowname
side.effects <- read.csv("FAERS_offsides_PTs.csv", header = T)
rownames(side.effects) <- side.effects$pert_id
side.effects <- side.effects[, -1]

# find drugs in common between FEARS_offsides, Gene_Expression and MACCS datasets
# in addition, matches the order of rows in datasets
common.pert_ids <- intersect(pert_ids, rownames(side.effects))
side.effects.com <- side.effects[common.pert_ids, ]
gene.data.com <- gene.data[common.pert_ids, ]
maccs.com <- maccs[common.pert_ids, ]

# convert gene code to gene name
gene.map <- read.csv("meta_Probes_info.csv", header = T, colClasses = c("character", "character"))
rownames(gene.map) <- gene.map$probe
# remove the added 'x' character in gene codes
gene.codes <- substring(colnames(gene.data.com), 2)
gene.names <- gene.map[gene.codes, 2]
colnames(gene.data.com) <- gene.names

#save new datasets
write.csv(maccs.com, "MACCS_681x166.csv")
write.csv(gene.data.com, "Gene_Expression_681x978.csv")
write.csv(side.effects.com, "FAERS_ADR_681x9404.csv")
