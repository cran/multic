##############################################################################
## Title: herit.se.q
## Description: Calculate the heritability and standard error for the
##              polygene and major.gene major diagonal values.
## Author: Eric Lunde, 04/29/2004
##############################################################################
herit.se <- function(polygene, major.gene1, environment,
                     major.diagonal.indices, inv.exp.sec.der.random,
                     poly = T)
{
  ## Calculation of the heritability and its standard error using the delta
  ## method

  ## First case, under the polygenic model
  if(poly) {
    ##print('poly')
    ##browser()
    
    effects.count <- length(major.diagonal.indices)
    
    ## If the number of polygene values does not the same as the number of
    ## major diagonal indices, strip away the non-diagonal values.  We only
    ## want the values from the major diagonal.
    if(length(polygene) != effects.count) {
      polygene <- polygene[major.diagonal.indices]
    }
    
    ## If the number of environment values does not the same as the number of
    ## major diagonal indices, strip away the non-diagonal values.  We only
    ## want the values from the major diagonal.
    if(length(environment) != effects.count) {
      environment <- environment[major.diagonal.indices]
    }
    
    ## If the number of values in the inv.exp.sec.der.random matrix does not
    ## equal the number of major diagonal indices squared, strip away all the
    ## values that do not correspond to the major diagonal values for the
    ## random effects.
    if(length(inv.exp.sec.der.random) != (2 * effects.count) ^ 2) {
      inv.exp.sec.der.random <-
        inv.exp.sec.der.random[c(major.diagonal.indices,
                                 2 * major.diagonal.indices[effects.count]
                                 + major.diagonal.indices),
                               c(major.diagonal.indices,
                                 2 * major.diagonal.indices[effects.count]
                                 + major.diagonal.indices)]
    }
    
    sums.of.major.diagonal <- polygene + environment
    squares.of.sums <- sums.of.major.diagonal ^ 2

    ## The heritability is the major diagonal polygene values divided by the
    ## sum calculated above
    herit.poly <- polygene/sums.of.major.diagonal
    
    ## g' for sigma^2  cbind the values together so that each row is an index
    ## of the major diagaonal.  The number of rows is equal to the number of
    ## major diagonal indices.  
    der.g <- cbind(polygene/squares.of.sums, - environment/squares.of.sums)

    ## D's dimensions are 2 by 2 by effects.count because for each effect we
    ## want to have 4 (2 by 2) values from the null inv.exp.sec.der.random
    ## matrix
    D <- array(0, dim = c(2, 2, effects.count))
    for(i in 1:effects.count) {
      D[, , i] <- 
        c(inv.exp.sec.der.random[                i,                 i],
          inv.exp.sec.der.random[effects.count + i,                 i],
          inv.exp.sec.der.random[                i, effects.count + i],
          inv.exp.sec.der.random[effects.count + i, effects.count + i])
    }

    ## For each index of D, add to herit.poly.se that index times the values
    ## in der.g that correspond to that index of D.  In the loop we have all
    ## the rows of column i (effects.count elements) times all the rows of
    ## column j (effects.cout elements) times the values of the third
    ## dimension at the i,j indices of D (effects.count elements).
    herit.poly.se <- 0
    for (i in 1:2) {
      for (j in 1:2) {
        herit.poly.se <-
          ifelse(is.na(herit.poly), NA,
                 ifelse(is.na(D[i, j, ]), herit.poly.se, 
                        herit.poly.se + der.g[, i] * der.g[, j] * D[i, j, ]))
      }
    }
    herit.poly.se <- sqrt(herit.poly.se)
    
    return (list("herit.poly"= herit.poly, "herit.poly.se"= herit.poly.se))
  }

  ## Second case, under the major gene model
  if(!poly) {    
   ##  print('spor')
   ##  browser()
    
    effects.count <- length(major.diagonal.indices)
    if(is.matrix(polygene)) {
      num.alt.hyps <- ncol(polygene)
    }else {
      num.alt.hyps <-
        ncol(matrix(polygene,
                    nrow = (effects.count + 1) * effects.count / 2))      
    }
    
    polygene <-
      array(polygene,
            dim = c(major.diagonal.indices[effects.count], num.alt.hyps) )
    major.gene1 <-
      array(major.gene1,
            dim = c(major.diagonal.indices[effects.count], num.alt.hyps) )
    environment <-
      array(environment,
            dim = c(major.diagonal.indices[effects.count], num.alt.hyps) )
    inv.exp.sec.der.random <-
      array(inv.exp.sec.der.random,
            dim = c(3 * major.diagonal.indices[effects.count],
              3 * major.diagonal.indices[effects.count],
              num.alt.hyps ) )

    ## If the number of polygene values does not the same as the number of
    ## major diagonal indices, strip away the non-diagonal values.  We only
    ## want the values from the major diagonal.
    if(length(polygene[, 1]) != effects.count) {
      polygene <- array(polygene[major.diagonal.indices, ],
                        dim = c(effects.count, num.alt.hyps) )
    }
    if(length(major.gene1[, 1]) != effects.count) {
      major.gene1 <- array(major.gene1[major.diagonal.indices,],
                           dim = c(effects.count, num.alt.hyps) )
    }
    if(length(environment[, 1]) != effects.count) {
      environment <- array(environment[major.diagonal.indices, ],
                           dim = c(effects.count, num.alt.hyps) )
    }
    if(length(inv.exp.sec.der.random[, , 1]) != (3 * effects.count) ^ 2) {      
      inv.exp.sec.der.random <-
        array(inv.exp.sec.der.random[c(major.diagonal.indices,
                                       major.diagonal.indices[effects.count]
                                       + major.diagonal.indices,
                                       2 * major.diagonal.indices[effects.count]
                                       + major.diagonal.indices), 
                                     c(major.diagonal.indices,
                                       major.diagonal.indices[effects.count]
                                       + major.diagonal.indices,
                                       2 * major.diagonal.indices[effects.count]
                                       + major.diagonal.indices),
                                     ],
              dim = c(3 * effects.count, 3 * effects.count, num.alt.hyps))
    }

    ## We need to have for each loci and each major diagonal value, the sum of
    ## the polygene, major gene, and environment effects.
    ## sums.of.major.diagonal
    sums.of.major.diagonal <- polygene + major.gene1 + environment
    squares.of.sums <- sums.of.major.diagonal ^ 2
    
    ## Calculation of the heritability, which is each polygene (or major gene)
    ## value being divided by the sum of the polygene, major gene, and
    ## environmental effects for that locus and that major diagonal index.
    herit.poly <- ifelse(polygene != 0, polygene/sums.of.major.diagonal, NA)
    herit.mg <-
      ifelse(major.gene1 != 0, major.gene1/sums.of.major.diagonal, NA)

    ## der.g1 and der.g2 have the dimensions 3 by effects.count by
    ## num.alt.hyps.  The first dimension has 3 values because each
    ## polygene or major gene std.err uses 3 derivative values in its
    ## calculation.  The second dimension has effects.count values because
    ## we need a std.err for each major diagonal index.  The third dimension
    ## has length num.alt.hyps which is how I determine how many loci
    ## (alternative hypthoses) there are.
    der.g1 <- array(0, dim = c(3, effects.count, num.alt.hyps))
    der.g2 <- array(0, dim = c(3, effects.count, num.alt.hyps))
    for(i in 1:effects.count) {
      ## g' for polygenic
      der.g1[1, i, ] <- (polygene[i, ] + environment[i, ])/squares.of.sums[i, ]
      der.g1[2, i, ] <- - major.gene1[i, ]/squares.of.sums[i, ]
      der.g1[3, i, ] <- - major.gene1[i, ]/squares.of.sums[i, ]
      
      ## g' for major gene
      der.g2[1, i, ] <- der.g1[2, i, ]
      der.g2[2, i, ] <- der.g1[1, i, ]
      der.g2[3, i, ] <- der.g1[3, i, ]
    }

    ## Calculate the variance of heritability
    herit.poly.se <- array(0, dim = c(effects.count, num.alt.hyps))
    herit.mg.se <- array(0, dim = c(effects.count, num.alt.hyps))

    ##print('herit.se.q')
    ##browser()

    D <- array(0, dim = c(3, 3, effects.count, num.alt.hyps))

    ## This for loop seperates the inv.exp.sec.der.random into 3 by 3 matrices
    ## for each effect count
    for(i in 1:effects.count) {
      D[, , i, ] <- 
        rbind(inv.exp.sec.der.random[                    i,                     i, ],
              inv.exp.sec.der.random[    effects.count + i,                     i, ],
              inv.exp.sec.der.random[2 * effects.count + i,                     i, ],
              inv.exp.sec.der.random[                    i,     effects.count + i, ],
              inv.exp.sec.der.random[    effects.count + i,     effects.count + i, ],
              inv.exp.sec.der.random[2 * effects.count + i,     effects.count + i, ],
              inv.exp.sec.der.random[                    i, 2 * effects.count + i, ],
              inv.exp.sec.der.random[    effects.count + i, 2 * effects.count + i, ],
              inv.exp.sec.der.random[2 * effects.count + i, 2 * effects.count + i, ])
    }

    for (i in 1:3) {
      for (j in 1:3) {
        herit.poly.se <-
          ifelse( is.na(herit.poly), NA,
                 ifelse(is.na(D[i, j, , ]), herit.poly.se, 
                        herit.poly.se
                        + der.g1[i, , ] * der.g1[j, , ] * D[i, j, , ]))
        
        herit.mg.se <-
          ifelse( is.na(herit.mg), NA,
                 ifelse(is.na(D[i, j, , ]), herit.mg.se,
                        herit.mg.se
                        + der.g2[i, , ] * der.g2[j, , ] * D[i, j, , ]))
      }
    }
    
    herit.poly.se <- sqrt(herit.poly.se)
    herit.mg.se <- sqrt(herit.mg.se)
    
    return (list("herit.poly"= herit.poly, "herit.poly.se"= herit.poly.se,
                 "herit.mg"= herit.mg, "herit.mg.se"= herit.mg.se))
  }
}
