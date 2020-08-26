library(sybil)
library(tidyverse)
library(ggplot2)

dir.create("output", showWarnings = FALSE)

##########################################
# Perform knockouts

# Create a list of the locus tags of genes in each deletion
ESS = read.csv("input/Table S2 Essential Gene Analysis.csv")
MGD = read.csv("input/Table S3 Multiple Gene Deletions.csv")
MGD = MGD %>% filter(!is.na(strain))

deletion.data = data.frame()

for (media in c("LB", "MS")) {
  
  #Reload fresh model
  # Reads file "ADP1_react.tsv"
  ADP1 <- readTSVmod(prefix = "input/ADP1")
  locus.tags.of.genes.in.model = allGenes(ADP1)
  
  # This code block shows the active update reactions
  ex <- findExchReact(ADP1)
  upt <- uptReact(ex)
  print(ex[upt])
  
  # ?modelorg to get accessor functions on ADP1
  
  
  
  if (media == "LB") {
    ## Set nutrients being used to Lysogeny Broth (LB-Miller)
    
    ## Set succinate and ammonium to be zero
    ADP1 <- changeBounds(ADP1, "EXF-SUC", lb = c(0), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-AMMONIUM", lb = c(0), ub = c(Inf))
    
    ## Amino acids, nucleobases, some cofactors
    ADP1 <- changeBounds(ADP1, "EXF-D-ALANINE", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-L-ALPHA-ALANINE", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-ARG", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-ASN", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-L-ASPARTATE", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-GLN", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-GLT", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-GLY", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-HOMO-SER", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-LYS", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-MET", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-PHE", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-PRO", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-D-SERINE", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-THR", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-TRP", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-TYR", lb = c(-10), ub = c(Inf))
    
    #nucleobases
    ADP1 <- changeBounds(ADP1, "EXF-URACIL", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-CYTOSINE", lb = c(-10), ub = c(Inf))
    ADP1 <- changeBounds(ADP1, "EXF-PURINE", lb = c(-10), ub = c(Inf))
    
    #cofactors
    ADP1 <- changeBounds(ADP1, "EXF-NICOTINAMIDE_NUCLEOTIDE", lb = c(-10), ub = c(Inf))
    
    
    #nucleosides
    #ADP1 <- changeBounds(ADP1, "EXF-URIDINE", lb = c(-10), ub = c(Inf))
    #ADP1 <- changeBounds(ADP1, "EXF-CYTIDINE", lb = c(-10), ub = c(Inf))  
    #ADP1 <- changeBounds(ADP1, "EXF-THYMIDINE", lb = c(-10), ub = c(Inf)) 
    #ADP1 <- changeBounds(ADP1, "EXF-ADENOSINE", lb = c(-10), ub = c(Inf))  
    #ADP1 <- changeBounds(ADP1, "EXF-GUANOSINE", lb = c(-10), ub = c(Inf)) 
    
    #deoxynucleosides
    #ADP1 <- changeBounds(ADP1, "EXF-DEOXYADENOSINE", lb = c(-10), ub = c(Inf))
    #ADP1 <- changeBounds(ADP1, "EXF-DEOXYCYTIDINE", lb = c(-10), ub = c(Inf))
    #ADP1 <- changeBounds(ADP1, "EXF-DEOXYGUANOSINE", lb = c(-10), ub = c(Inf))
    #ADP1 <- changeBounds(ADP1, "EXF-DEOXYURIDINE", lb = c(-10), ub = c(Inf))
    
  } else {
    ## Set nutrients being used to minimal succinate
    ADP1 <- changeBounds(ADP1, "EXF-SUC", lb = c(-10), ub = c(Inf))
  }
  
  #Run model with no deletions...
  optL <- optimizeProb(ADP1, algorithm = "fba", retOptSol = FALSE)
  
  
  for (row.index in 1:nrow(MGD)) {
    
    this.row=MGD[row.index,]
    
    #One end or entire gene is within deletion range = overlaps
    deleted.genes = ESS %>% filter( (start >= this.row$start ) & (start <= this.row$end ) | 
                                      (end >= this.row$start ) & (end <= this.row$end ) |
                                      (start >= this.row$start ) & (end <= this.row$end ) ) 
    
    locus.tags.of.deleted.genes = as.character(deleted.genes$id)
    
    locus.tags.of.deleted.genes.in.model = locus.tags.of.deleted.genes[locus.tags.of.deleted.genes %in% locus.tags.of.genes.in.model]
    
    if (length(locus.tags.of.deleted.genes.in.model) > 0) {
      optJ = geneDeletion(ADP1,locus.tags.of.deleted.genes.in.model, combinations=length(locus.tags.of.deleted.genes.in.model))
      this.growth.rate = lp_obj(optJ)
      if (this.growth.rate < 1E-6) {
        this.growth.rate = 0
      }
      this.deleted.fluxes = fluxdels(optJ)
    } else {
      this.growth.rate = optL$obj
      this.deleted.fluxes = c()
    }
    
    new.row = data.frame(
      region = this.row$region,
      strain = this.row$strain,
      media = media,
      deletion = this.row$deletion,
      num.genes.deleted = length(locus.tags.of.deleted.genes),
      num.genes.deleted.in.model = length(locus.tags.of.deleted.genes.in.model), 
      num.fluxes.deleted = length(this.deleted.fluxes[[1]]),
      predicted.growth.rate = this.growth.rate,
      fluxes.deleted = paste(this.deleted.fluxes[[1]],collapse=",")
    )
    
    deletion.data = rbind(deletion.data, new.row)
  }
  
}


## Calculate relative growth rates
max.growth.rates.per.media = deletion.data %>% group_by(media) %>% summarize(max.predicted.growth.rate = max(predicted.growth.rate))
deletion.data = deletion.data %>% left_join(max.growth.rates.per.media, by="media")
deletion.data$relative.predicted.growth.rate = deletion.data$predicted.growth.rate / deletion.data$max.predicted.growth.rate

write_tsv(deletion.data, "output/deletion_flux.tsv")


deletion.data$strain = factor(deletion.data$strain, levels=unique(deletion.data$strain))

          
          
#gene to protein rules
#gpr(ADP1)
#metabolites
metabolites = sort(met_name(ADP1))
write.csv(metabolites, "metabolites.csv")
 
ggplot(deletion.data, aes(x=strain, y=relative.predicted.growth.rate, group=media)) + geom_bar(stat="identity")  + facet_grid(rows = vars(media), scales="free_y") + theme_linedraw() + theme(axis.text.x=element_text(angle = 90, hjust = 0, vjust=0.5)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

