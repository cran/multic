#####################################################################
## Function: plot.family.lods
## Purpose:  Summarize Multic Object with calc.fam.liks = TRUE
## Author: Curt Olswold / Eric Lunde / Beth Atkinson
## Created: 4/5/2005
## Parameters:
##  multic.obj: Multic object that contains the fam lod scores
##  npeakfams: Number of top families to consider at peak {5}
##  title: Descriptive 1st line of title {NULL}
##  title.cex: size of title text {.75}
##  type: which of the three plots to present (or all of them)
##        {c("top", "total", "proportion", "all")}
##  legend: Show a legend {T}
##  legend.loc: location of legend on the existing plot, or on new page
##         {c("extra","upper left","upper middle","upper right")}
######################################################################

plotFamilyLods <- function(x, npeakfams = 5, title = NULL,
                             title.cex = .75,legend=T,
                             legend.loc= c("left","middle","right","extra"),
                             type = c("top", "total", "proportion","all"), ...)
{
  ## Since R CMD check requires the first arg be x, rename x to something
  ## more descriptive, multic.obj.
  multic.obj <- x  
  if( ! is.numeric(multic.obj$fam.log.liks) ) {
    stop("\nplot.family.lods cannot be called on a multic object that",
         "does\n",
         "not contain family log likelihoods.\n",
         "To calculate family log likelihoods, rerun multic with the\n",
         "argument \"calc.fam.log.liks = TRUE\".\n")
  }

  if( dim(multic.obj$fam.log.liks)[3] < 2 ) {
    stop("\nplot.family.lods cannot plot the LOD scores of a polygenic ",
         "object.\n")
  }

  type <- ifelse(length(type)==4,"top",
                 c("top", "total", "proportion","all")[pmatch(type,
                                                c("top", "total", "proportion","all"))])

  legend.loc <- ifelse(length(legend.loc)==4,'left',
                       c("left","middle","right","extra")[pmatch(legend.loc,
                                                c("left","middle","right","extra"))])

  multic.obj.string <- deparse(substitute(x))

  ## If the user did not supply a title, use the multic object's name
  if(is.null(title)) {
    title <- multic.obj.string
  }

  ## Depending on title, set multic.obj.string to something useful
  if(title == multic.obj.string) {
    multic.obj.string <- ""
  } else {
    multic.obj.string <- paste("(", multic.obj.string, ")", sep = "")
  }

  ## get family LODs  
  top.families.lod.scores <-
    get.top.n.families(multic.obj, npeakfams)[, "lod.score"]
  peak <- max(multic.obj$log.liks$lod.score, na.rm = TRUE)

  ## Proportion of famLod to totLod
  family.percentage.at.peak <-
    paste(names(top.families.lod.scores), ': ',
          100 * round(top.families.lod.scores / peak, digits = 2), '%',
          sep = '')

  ## keep top npeakfams (familyLod in previous version)
  family.lods <-
    multic.obj$fam.log.liks[names(top.families.lod.scores), 2, -1]

  ## change negative values of family.lods to zero
  family.lods.neg <- family.lods < 0
  family.lods.noneg <- family.lods
  family.lods.noneg[family.lods.neg] <- 0

  ## Top npeakfams Total LODs (pktotLod in previous version)
  family.total.lods <- apply(family.lods.noneg, 2, sum)

  ## keep top npeakfams for proportions (familyProp in previous version)
  family.percentages <-
    apply(family.lods, 1,
          function(lods, total.lod.scores) { lods / total.lod.scores },
          multic.obj$log.liks$lod.score[-1])

  fp.na <- is.na(family.percentages) | is.infinite(family.percentages)
  family.percentages[fp.na] <- 0
  
  ## Get the numeric version of the centimorgan values
  split.positions <- multic.strings.split(row.names(multic.obj$log.liks)[-1],
                                          sep = ".")
  numeric.split.positions <- split.positions[, -c(1,2), drop=F]
  cm <- as.numeric(apply(numeric.split.positions, MARGIN=1, FUN=paste, collapse='.'))


  ## plot top npeakfams
  if(type %in% c("top","all")) {
    matplot(y = t(family.lods.noneg),
            x = cm,
            type = 'l',
            lwd = 3,
            lty = 1:npeakfams + 2,
            col = 1:npeakfams + 2,
            ylim = c(0, max(family.lods) ),
            ylab = 'LOD',
            xlab = 'Position', ...)
    title(paste(title, '\nTop ', npeakfams, ' Family LOD Scores ',
                multic.obj.string, sep = ''),
          cex = title.cex)
  }
  
  ## plot total + top npeakfams
  if(type %in% c("total","all")) {
    y.data <- rbind(multic.obj$log.liks$lod.score[-1],
                    family.total.lods, family.lods.noneg)
    matplot(y = t(y.data),
            x = cm,
            type = 'l',
            lwd = 3,
            lty = 1:(npeakfams + 2),
            col = 1:(npeakfams + 2),
            ylim = c(0, max(y.data)),
            ylab = 'LOD',
            xlab = 'Position', ...)
    title(paste(title, '\nTotal and Top ', npeakfams, ' Family LOD Scores ',
                multic.obj.string, sep = ''),
          cex = title.cex)
  }
  
  ## plot top npeakfams proportion of total Lod
  if(type %in% c("proportion","all")) {
    matplot(y = family.percentages,
            x = cm,
            type = 'l',
            lwd = 3,
            lty = 1:npeakfams + 2,
            col = 1:npeakfams + 2,
            ylab = 'Proportion',
            xlab = 'Position', ...)
    title(paste(title, '\nProportion of Family LOD to Total LOD ',
                multic.obj.string, sep = ''),
          cex = title.cex)
  }
  
  ## CREATE A LEGEND

  ## Information on current plot
  usr <- par("usr")
  
  if(type =="all" & legend==T) {
    plot(x = 0:1, y = 0:1, ylab = '', xlab = '', axes = FALSE, type = 'n')
    legend(x = .25,y = .95,
           legend = c('Total', paste('Fam Total', npeakfams),
             family.percentage.at.peak),
           col = 1:(npeakfams + 2),
           lty = 1:(npeakfams + 2),
           lwd = 3)
    }

  if(type!="all" & legend==T){

    bit <- .005
    xloc <- ifelse(legend.loc=='left', usr[1]*(1+bit),
                   ifelse(legend.loc=='right',usr[2]*(1-bit),
                          (usr[1]+usr[2])/2))
    yloc <- usr[4]*(1-bit)

    if(type=='total') {
      legend(x = xloc, y = yloc,
             legend = c('Total', paste('Fam Total', npeakfams),
               family.percentage.at.peak),
             col = 1:(npeakfams + 2),
             lty = 1:(npeakfams + 2),
             lwd = 3,bty='n')
    }

    if(type!='total') {
      legend(x = xloc, y = yloc,
             legend = family.percentage.at.peak,
             col = 1:npeakfams + 2,
             lty = 1:npeakfams + 2,
             lwd = 3, bty='n')
    }

  }
  

  invisible(family.lods)
}

## Method renamed to plotFamilyLods because it's not S3.
## Keep alias for backward compatibility.
plot.family.lods <- plotFamilyLods
