pchisq.mix <- function(lod.lrt, trait.count, lod = T) {
  df.mix <- get.df.mix(trait.count)
  p.mix <- get.p.mix(trait.count)

  if(length(df.mix)!=length(p.mix))
    stop("Lengths of df.mix and p.mix do not match")
  if(all(df.mix==0)) stop("Cannot have all values df.mix = 0")
  if(sum(p.mix)!=1) stop("p.mix does not sum to 1")
  if(lod == T) {
    lrt <- lod.lrt * 2 * log(10)
  }else {
    lrt <- lod.lrt
  }

  zed <- df.mix > 0
  df <-df.mix[zed]
  p <- p.mix[zed]
  
  pval <- NULL
  for (i in 1:length(lod.lrt)){
    pval[i] <- sum(p * (1 - pchisq(lrt[i], df)))
  }
  
  return(pval)
}

get.p.mix <- function(trait.count) {
  if(trait.count < 1 | trait.count > 5) {
    stop(paste("\ntrait.count (", trait.count,
               ") must be between 1 and 5 inclusive\n",
               "pchisq.mix.q key 103\n", sep = ""))
  }
  switch(trait.count,
         c(1/2, 1/2),
         c(1/4, 1/2, 1/4),
         c(1/8, 3/8, 3/8, 1/8),
         c(1/16, 4/16, 6/16, 4/16, 1/16),
         c(1/25, 5/25, 10/25, 10/25, 5/25, 1/25))
}

get.df.mix <- function(trait.count) {
  if(trait.count < 1 | trait.count > 5) {
    stop(paste("\ntrait.count (", trait.count,
               ") must be between 1 and 5 inclusive\n",
               "pchisq.mix.q key 97\n", sep = ""))
  }
  switch(trait.count,
         c(0, 1),
         c(0, 1, 3),
         c(0, 1, 3, 6),
         c(0, 1, 3, 6, 10),
         c(0, 1, 3, 6, 10, 15))
}
