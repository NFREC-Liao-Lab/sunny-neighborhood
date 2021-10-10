# Hui-Ling Liao et. al. Heterospecific neighbor plants impact root microbiome diversity and molecular function of root fungi
# Figure 4a, 4b
# Table S2

library(readxl)
library(vegan)
library(RColorBrewer)
library(gplots)

# Dataset 2A.
# Relative aboundance of fungal taxa for an individual root sample.
# The numbers indicate the percent for the ratio of the D1D2 reads that were blasted to an fungal taxa to the total D1D2 reads extracted from an individual root sample.
# Only the data that >1% reads in at least on sample were shown.


# librdata_header <- read.csv("Neighbor_Dataset_2A.csv", nrow=1)
# data <- read.csv("Neighbor_Dataset_2A.csv", skip=1)
data_header <- data.frame(read_excel("Neighbor_Dataset_2A&B.xlsx", n_max=3, col_names=F, skip=1))
data <- data.frame(read_excel("Neighbor_Dataset_2A&B.xlsx", col_names=F, skip=4))

if(any(dim(data)!=c(44,18))){
  stop("Check data import!!!!!!!!!!!!!!!!")
}

taxa_name <- setNames(data[,2:3], data_header[1, 2:3])
taxa_name$group <- gsub("Endophye", "endophyte", taxa_name$group)
taxa_name$group <- factor(gsub("[_ ].+", "", taxa_name$group))
taxa_name$group5 <- factor(gsub("/.+", "", taxa_name$group))
taxa_name2 <- apply(taxa_name, 1, function(x){paste(x, collapse="_")})

h2 <- data_header[2,]
# h2 <- gsub("P. taeda", "Taead", h2)
h2 <- gsub("P.? ?tri", "P\\. tri", h2)
h2 <- na.omit(h2)

data2 <- data[,-(1:3)]
colnames(data2) <- h2
rownames(data2) <- group_name2
data2[is.na(data2)]<- 0

dataExp <- t(apply(data2, 2, function(x){x/sum(x, na.rm=T)}))
apply(dataExp,1,sum, na.rm=T)
site <- data.frame(ID=gsub("\\%read", "", data_header[1,-(1:3)]), label=factor(rownames(dataExp)))


cPal8 <- brewer.pal(8, "Dark2")


pdf("2A_permanova_each_sites.pdf")
all_labels <- levels(site$label)
for(cc in combn(4,2, simplify=FALSE)){
  index <- site$label %in% all_labels[cc]
  dataSub <- dataExp[index, ]
  abDist <- vegdist(dataSub, method="bray")
  result <- adonis(abDist ~ site$label[index] )#, data=site)
  textplot(c(all_labels[cc], capture.output(result)), cex=0.7)
}
dev.off()


all_labels <- levels(site$label)
resultAll <- list()
for(cc in combn(4,2, simplify=FALSE)){
  index <- site$label %in% all_labels[cc]
  dataSub <- dataExp[index, ]
  site_cc <- site$label[index]
  dd <- rep(NA, length=NCOL(dataSub))
  for(i in 1:NCOL(dataSub)){
     try({
       abDist <- vegdist(dataSub[,i], method="bray")
       # zero-adjusted Bray–Curtis
       abDist[is.na(abDist)] <- 0
       result <- adonis(abDist ~ site_cc ) #, data=site)
       # result <- adonis(dataSub[,i] ~ site_cc, method="bray" ) #, data=site)
       dd[i] <-   result$aov.tab[["Pr(>F)"]][1]
     }, silent=T)
  }
  resultAll[[paste(all_labels[cc], collapse=" : ")]] <- dd

}
result_output <- cbind(taxa_name, do.call("cbind", resultAll))
write.csv(result_output, file="2A_permanova_each_ref.csv")

#   result$aov.tab[["Pr(>F)"]][1]
# }


abDist <- vegdist(dataExp, method="bray")
result <- adonis(abDist ~ site$label )#, data=site)
pdf("2A_permanova_all_sites.pdf")
textplot(capture.output(result))
dev.off()


stress <- vector()
for(k in 1:10){
  mm <- metaMDS((dataExp), k=k, distance="bray", trymax=100, autotransform=F, noshare=F, wascores=T)
  stress[k] <- mm$stress
}
formatC(stress)
png("2A_stress_values.png")
plot(stress, xlab="Dimension", main="NMDS Stress values for 2A")
dev.off()

k <- 2

tiff("2A_Site.tiff", res=300, width=2000, height=2000)

cPalSet1 <- brewer.pal(8, "Set1")
customColour <- as.numeric(site$label)
customColour[customColour==1] <- cPalSet1[4]
customColour[customColour==2] <- cPalSet1[3]
customColour[customColour==3] <- cPalSet1[2]
customColour[customColour==4] <- cPalSet1[1]

mm <- metaMDS((dataExp), k=k, distance="bray", trymax=100, autotransform=F, noshare=F, wascores=T)
choice <- c(1, 2)
plot(mm, type="n", display = "sites", cex=2, choice=choice, xlim=c(-1.4, 1.5))
  # main=sprintf("%s   p-value:%.3f", "", result$aov.tab[["Pr(>F)"]][1]), )
# text(mm, display="sites", labels = site$ID, col=cPal8[as.numeric(site$label)], choice=choice )
# points(mm, col=cPal8[as.numeric(site$label)], cex=1.2)
points(mm, pch=as.numeric(site$label)-1, cex=1.2, col=customColour)
with(site, ordiellipse(mm, label, scaling = "symmetric", kind = "ehull", col = cPalSet1[4:1], lwd=2, label=T))
legend("bottomleft", legend=levels(site$label), pch=1:nlevels(site$label)-1, bty="n", col=cPalSet1[4:1])

dev.off()

# png("2A_Site_BW.png")
tiff("2A_Site_BW.tiff", res=300, width=2000, height=2000)

plot(mm, type="n", display = "sites", cex=2, choice=choice, xlim=c(-1.4, 1.5))
points(mm, pch=as.numeric(site$label)-1, cex=1.2)
with(site, ordiellipse(mm, label, scaling = "symmetric", kind = "ehull", col = cPalSet1[4:1], lwd=2, label=T))
legend("bottomleft", legend=levels(site$label), pch=1:nlevels(site$label)-1, bty="n")
dev.off()




tiff("2A_Species.tiff", res=300, width=2000, height=2000)
# dev.off()
# png("2A_Species.png", width=600, height=600)
new_legend <- levels(taxa_name$group5)
new_legend[3] <- "Endophyte or/and EM"
new_legend[4] <- "Pathogen or/and endophyte"
plot(mm, main="", type="n")
# text(mm, display="sites", labels = site$ID, col=cPal8[8])
points(mm, pch=as.numeric(site$label)-1, col=1, cex=1.2)
text(mm, display="species", labels = taxa_name$ref, col=cPal8[taxa_name$group5], cex=0.7)
legend("bottomleft", legend=new_legend, fill=cPal8[1:nlevels(taxa_name$group5)], bty="n")
dev.off()

write.csv(taxa_name, file="2A_taxa.csv")
