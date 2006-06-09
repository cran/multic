add.heritability <- function(polygenic, major.gene1, environmental,
                             sibling.sibling,
                             parent.parent,
                             parent.offspring,
                             major.diagonal.indices, inv.exp.sec.der.random,
                             run.alternative.hyps)
{
  trait.count <- length(major.diagonal.indices)
  ##print("add.herit")
  ##browser()
  
  ## This is just for the polygenic model (null hyp), hence the poly = T  
  null.herit.se <- herit.se(polygenic[, "Estimate", "null"],
                            major.gene1[, "Estimate", "null"],
                            environmental[, "Estimate", "null"],
                            major.diagonal.indices,
                            inv.exp.sec.der.random[, , "null"],
                            poly = T)
  
  ## The second dimension of polygenic in these next assignments are indexed
  ## by 5 = heritability, 6 = standard error of heritability, and 7 =
  ## p.value for the standard error of heritability.
  polygenic[major.diagonal.indices, "h^2", 1] <- null.herit.se$herit.poly
  polygenic[major.diagonal.indices, "se.h^2", 1] <-
    null.herit.se$herit.poly.se
  polygenic[major.diagonal.indices, 7, 1] <-
    round(pchisq.mix((null.herit.se$herit.poly / null.herit.se$herit.poly.se)
                     ^ 2,
                     trait.count = trait.count, lod = F),
          digits = 4)
  
  ## For the space that values could not be calculated (because it doesn't
  ## make sense), fill them with NA's
  polygenic[-major.diagonal.indices, "h^2", 1] <- NA
  polygenic[-major.diagonal.indices, "se.h^2", 1] <- NA
  polygenic[-major.diagonal.indices, 7, 1] <- NA

  major.gene1[, , 1] <- NA
  major.gene1[-major.diagonal.indices, "h^2", -1] <- NA
  major.gene1[-major.diagonal.indices, "se.h^2", -1] <- NA
  major.gene1[-major.diagonal.indices, 7, -1] <- NA

  if(run.alternative.hyps) {
    ## This is just for the sporadic model (alt hyp), hence poly = F
    alt.herit.se <- herit.se(polygenic[, "Estimate", -1],
                             major.gene1[, "Estimate", -1],
                             environmental[, "Estimate", -1],
                             major.diagonal.indices,
                             inv.exp.sec.der.random[, , -1],
                             poly = F)
    
    #print('add.herit')
    #browser()
  
    ## The second dimension of polygenic in these next assignments are indexed
    ## by 5 = heritability, 6 = standard error of heritability, and 7 =
    ## p.value for the standard error of heritability.
    polygenic[major.diagonal.indices, "h^2", -1] <- alt.herit.se$herit.poly
    polygenic[major.diagonal.indices, "se.h^2", -1] <-
      alt.herit.se$herit.poly.se
    polygenic[major.diagonal.indices, 7, -1] <-
      round(pchisq.mix((alt.herit.se$herit.poly / alt.herit.se$herit.poly.se)
                       ^ 2,
                       trait.count = trait.count, lod = F),
            digits = 4)
    
    ## For the space that values could not be calculated (becasue it doesn't
    ## make sense), fill them with NA's
    polygenic[-major.diagonal.indices, "h^2", -1] <- NA
    polygenic[-major.diagonal.indices, "se.h^2", -1] <- NA
    polygenic[-major.diagonal.indices, 7, -1] <- NA
    
    major.gene1[major.diagonal.indices, "h^2", -1] <- alt.herit.se$herit.mg
    major.gene1[major.diagonal.indices, "se.h^2", -1] <-
      alt.herit.se$herit.mg.se
    major.gene1[major.diagonal.indices, 7, -1] <-
      round(pchisq.mix((alt.herit.se$herit.mg / alt.herit.se$herit.mg.se) ^ 2,
                       trait.count = trait.count, lod = F),
            digits = 4)
    
    ## For the space that values could not be calculated (becasue it doesn't
    ## make sense), fill them with NA's
    major.gene1[-major.diagonal.indices, "h^2", -1] <- NA
    major.gene1[-major.diagonal.indices, "se.h^2", -1] <- NA
    major.gene1[-major.diagonal.indices, 7, -1] <- NA
  }
  
  return ( list(polygenic = polygenic, major.gene1 = major.gene1) )
}
