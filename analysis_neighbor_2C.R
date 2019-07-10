# taxa_name$group <- gsub("Endophye", "endophyte", taxa_name$group)
library(readxl)
library(vegan)
library(RColorBrewer)
library(gplots)


# librdata_header <- read.csv("Neighbor_Dataset_2C.csv", nrow=1)
# data <- read.csv("Neighbor_Dataset_2C.csv", skip=1)
data_header <- data.frame(read_excel("Neighbor_Dataset_2A&B.xlsx", n_max=2, col_names=F, skip=1, sheet="2C"))
data <- data.frame(read_excel("Neighbor_Dataset_2A&B.xlsx", col_names=F, skip=3, sheet="2C"))
head(data)

if(any(dim(data)!=c(51,18))){
  stop("Check data import!!!!!!!!!!!!!!!!")
}

taxa_name <- setNames(data[,2:3], data_header[1, 2:3])
taxa_name$group <- taxa_name[[data_header[1,2]]]
index <- names(which(table(taxa_name$group) <= 2))
taxa_name$group[taxa_name$group %in% index] <- "Other"
taxa_name$group <- factor(taxa_name$group)
# taxa_name$group <- factor(gsub("[_ ].+", "", taxa_name$group))
# taxa_name$group5 <- factor(gsub("/.+", "", taxa_name$group))
# taxa_name2 <- apply(taxa_name, 1, function(x){paste(x, collapse="_")})

h2 <- data_header[2,]
# h2 <- gsub("P. taeda", "Taead", h2)
h2 <- gsub("Ptri", "P\\. tri", h2)
h2 <- na.omit(h2)

data2 <- data[,-(1:3)]
colnames(data2) <- h2
# rownames(data) <- group_name2
data[is.na(data)]<- 0

dataExp <- t(apply(data2, 2, function(x){x/sum(x, na.rm=T)}))
apply(dataExp,1,sum, na.rm=T)
site <- data.frame(ID=gsub("Sum of(.+)_\\%read", "\\1", data_header[1,-(1:3)]), label=factor(rownames(dataExp)))


cPal8 <- brewer.pal(8, "Dark2")


# site <- rep(c("A","B"), each=4)

# dataTrt <- data.frame(site=site, rep=paste0(site, 1:4))
# cPal8 <- brewer.pal(8, "Dark2")

# for(sheet in sheets){
#   pdf(sprintf("overallResult_%s.pdf", sheet), width=8, height=8)
#   data <- data.frame(read_xlsx("V10_Fig3_statistics.xlsx", skip=1, sheet=sheet))
#   names(data) <- gsub("ecological.function", "group", names(data))
#   summary(data)
#   data[["group"]] <- gsub("Sap .+", "Sap", data[["group"]])
#   data[["group"]] <- gsub("pathogen", "Pathogen", data[["group"]])
#   data[["group"]] <- gsub("Ecto-,.+", "Ecto", data[["group"]])
#   data[["group"]] <- factor(data[["group"]])
#
#   dataExp <- t(data[, grep("[A|B][1-4]", names(data), value=T)])
#   colnames(dataExp) <- data$ref


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
write.csv(result_output, file="2C_permanova_each_ref.csv")


pdf("2C_permanova_each_sites.pdf")
all_labels <- levels(site$label)
for(cc in combn(4,2, simplify=FALSE)){
  index <- site$label %in% all_labels[cc]
  dataSub <- dataExp[index, ]
  abDist <- vegdist(dataSub, method="bray")
  result <- adonis(abDist ~ site$label[index] )#, data=site)
  textplot(c(all_labels[cc], capture.output(result)), cex=0.7)
}
dev.off()


abDist <- vegdist(dataExp, method="bray")
result <- adonis(abDist ~ site$label )#, data=site)
pdf("2C_permanova_all_sites.pdf")
textplot(capture.output(result))
dev.off()

mm <- metaMDS((dataExp), k=2, distance="euclidean", trymax=100, autotransform=F, noshare=F, wascores=T)

stress <- vector()
for(k in 1:10){
  mm <- metaMDS((dataExp), k=k, distance="bray", trymax=100, autotransform=F, noshare=F, wascores=T)
  stress[k] <- mm$stress
}
plot(stress)

k <- 2

png("2C_Site.png", res=300, width=2000, height=2000)
tiff("2C_Site.tiff", res=300, width=2000, height=2000)

cPalSet1 <- brewer.pal(8, "Set1")
customColour <- as.numeric(site$label)
customColour[customColour==1] <- cPalSet1[4]
customColour[customColour==2] <- cPalSet1[3]
customColour[customColour==3] <- cPalSet1[2]
customColour[customColour==4] <- cPalSet1[1]

mm <- metaMDS((dataExp), k=k, distance="bray", trymax=100, autotransform=F, noshare=F, wascores=T)
choice <- c(1, 2)
plot(mm, type="n", display = "sites", cex=2, choice=choice, xlim=c(-1.2, 0.9))
  # main=sprintf("%s   p-value:%.3f", "", result$aov.tab[["Pr(>F)"]][1]), )
# text(mm, display="sites", labels = site$ID, col=cPal8[as.numeric(site$label)], choice=choice )
# points(mm, col=cPal8[as.numeric(site$label)], cex=1.2)
points(mm, pch=as.numeric(site$label)-1, cex=1.2, col=customColour)
with(site, ordiellipse(mm, label, scaling = "symmetric", kind = "ehull", col = cPalSet1[4:1], lwd=2, label=T))
legend("bottomleft", legend=levels(site$label), pch=1:nlevels(site$label)-1, bty="n", col=cPalSet1[4:1])

dev.off()

# png("2C_Site_BW.png", res=300, width=2000, height=2000)
tiff("2C_Site_BW.tiff", res=300, width=2000, height=2000)

plot(mm, type="n", display = "sites", cex=2, choice=choice, xlim=c(-1.2, 0.9))
points(mm, pch=as.numeric(site$label)-1, cex=1.2)
with(site, ordiellipse(mm, label, scaling = "symmetric", kind = "ehull", col = cPalSet1[4:1], lwd=2, label=T))
legend("bottomleft", legend=levels(site$label), pch=1:nlevels(site$label)-1, bty="n")
dev.off()

  # plot(mm, type="p", display = "sites", cex=2, choice=c(2,3),
  #   main=sprintf("%s   p-value:%.3f", "", result$aov.tab[["Pr(>F)"]][1]))
  # text(mm, display="sites", labels = site$ID, col=cPal8[as.numeric(site$label)], choice=c(2,3) )
  # with(site, ordiellipse(mm, label, scaling = "symmetric", kind = "ehull", col = cPal8[1:nlevels(label)], lwd=3, label=T, choice=c(2,3)))




png("2C_Species.png", res=300, width=2000, height=2000)
tiff("2C_Species.tiff", res=300, width=2000, height=2000)
# dev.off()
# png("2C_Species.png", width=600, height=600)
new_legend <- levels(taxa_name$group)
# new_legend[3] <- "Endophyte or/and EM"
# new_legend[4] <- "Pathogen or/and endophyte"
plot(mm, main="", type="n")
# text(mm, display="sites", labels = site$ID, col=cPal8[8])
points(mm, pch=as.numeric(site$label)-1, col=1, cex=1.2)
text(mm, display="species", labels = taxa_name[["Row Labels"]], col=cPal8[taxa_name$group], cex=0.7)
legend("bottomleft", legend=new_legend, fill=cPal8[1:nlevels(taxa_name$group)], bty="n")
dev.off()

write.csv(taxa_name, file="2C_taxa.csv")

  # , col=cPal8[as.numeric(group)], cex=0.7))
  # with(dataExp, text(mm, display="species", labels = as.character(ref), col=cPal8[as.numeric(group)], cex=0.7))
  # #
  # # dev.off()
  #
}

# TODO: Analysis 2C
# TODO: 2A. 44 by 7 (all + pairwise) PERANOVA p-value table.
# TODO: 2A. by groups (EM, AM..etc) by 7
