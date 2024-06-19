## Generic function to make EnrichR analyses.

## ARGUMENTS:
## geneList: vector with a collection of gene names.
## plotBaseName: output base name; "_EnrichR_DB.png" will be appended.
## dbList: list of databases to be interrogated.
## maxEnrichRpVal: maximal adjusted p_value to filter enrichment results.

## OUTPUT:
## The function produces an adj p_value restricted barplot with the top 50 functions 
## and a complete Excel file for each database used. The barplot is colored by 
## enrichment (equation for enrichment: [(b/n)/(B/N)]; where b is the number of
## overlapping genes, n is the number of target genes, B is the number of genes 
## in annotated set and N is the number of genes in the gene universe).
## Enrichment values are colored using a Yellow-Orange-Red palette.


library(RColorBrewer)

doEnrichRanalyses <- function (geneList, plotBaseName, dbList, maxEnrichRpVal) {
  enriched <- enrichr(geneList, dbList)
  for (db in names(enriched)) {
    dataPlot <- enriched[[db]][enriched[[db]]$Adjusted.P.value < maxEnrichRpVal,]
    if (nrow(dataPlot) == 0) {
      next
    }
    dataPlot <- tail(dataPlot[order(dataPlot$Adjusted.P.value, decreasing = T),],50)
    lowerMargin = (-1.29871794871795 * nrow(dataPlot)) + 74.9358974358974
    
    colGradient <- colorRampPalette(brewer.pal(9, "YlOrRd"))
    
    lim <- round(max(abs(c(max(dataPlot$Odds.Ratio, na.rm=T),min(dataPlot$Odds.Ratio, na.rm=T))))) + 1
    colors <- colGradient(50)[as.numeric(cut(dataPlot$Odds.Ratio,breaks = seq(1,lim,length.out = 51),include.lowest = T))]
    
    correction<-min(dataPlot[which(-log10(dataPlot$Adjusted.P.value) != 0),"Adjusted.P.value"])/10
    
    outPlotName <- paste0(plotBaseName, "_EnrichR_", db, ".png")
    png(outPlotName,30,40,units="cm",res=600)
    par(mar=c(lowerMargin,25,4,6))
    barplot(-log10(dataPlot$Adjusted.P.value) + correction
            ,horiz=T
            ,names.arg = dataPlot$Term
            ,xlab="-log10(adj. p-value)"
            ,las=1
            ,cex.names = 0.8
            ,cex.axis = 0.8
            ,cex.lab = 0.8
            ,mgp = c(1.5,0.5,0)
            ,cex.main=1.2
            ,col = colors
    )
    
    Xsta <- max(-log10(dataPlot$Adjusted.P.value + correction))
    Xsta <- Xsta - Xsta*0.1
    Xmax <- Xsta + Xsta/10
    Ymax <- round(nrow(dataPlot)/3)
    color.legend(Xsta,0,Xmax,Ymax,round(c(1,lim),2),colGradient(50),cex=1.5,align="rb",gradient="y")
    par(xpd=TRUE)
    labelx<-(Xsta+Xmax)/2
    labely<-Ymax + Ymax * 0.1
    text(labelx, Ymax + 1, "Enrichment", cex = 1.5)
    par(xpd=FALSE)
    dev.off()
    system(paste0("convert ",outPlotName," -trim -bordercolor White -border 40x20 +repage ",outPlotName),wait=FALSE)
    write.xlsx(enriched[[db]], paste0(plotBaseName, "_EnrichR_", db, ".xlsx"))
  }
  enrichRdataComplete <- do.call("rbind",enriched)
  return (enrichRdataComplete)
}
