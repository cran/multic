## Load the multic library
library(multic)

## Use solar2multic to translate SOLAR's mibds into mloci.out and phi2 into 
## share.out
solar.output.directory <- "multicInput.solar"
solar2multic(phi2 = "solarOutput/phi2.gz",
             pedigree.file = "solarOutput/simulated.ped",
             pedindex.out = "solarOutput/pedindex.out",
             pedindex.cde = "solarOutput/pedindex.cde",
             ibd.directory = "solarOutput/solarMibds",
             output.directory = solar.output.directory)

## Use sw2mloci to translate SimWalk's IBDs into mloci.out (only one
## mloci.out is needed, but I wanted to show examples of both calls)
## NOTE: the mloci.out.gz made by sw2mloci is NOT meant to work with the
## calls to multic below.
sw2mloci("swOutput", "swOutput/c18.map",
         output.directory = "multicInput.simwalk")

## Create a data.frame with the pedigree and phenotype information
ped.file.name <- "solarOutput/simulated.ped"
ped.file <- read.table(ped.file.name, header = TRUE, sep = ",")
phen.file.name <- "solarOutput/simulated.phen"
phen.file <- read.table(phen.file.name, header = TRUE, sep = ",")
ped.phen <- merge(ped.file, phen.file)

## Make sure data is sorted numerically. Merge sorts lexographically.
ped.phen <- ped.phen[order(ped.phen$famid, ped.phen$id), ]

## Call multic with a univariate model and no covariates
trait1 <-
  multic(trait1 ~ 1,
         data = ped.phen,
         famid, id, fa, mo, sex,
         mloci.out = paste(solar.output.directory, "mloci.out", sep = "/"),
         share.out = paste(solar.output.directory, "share.out", sep = "/"))

## Call multic with a different univariate model and two covariates
trait2 <-
  multic(trait2 ~ sex + age,
         data = ped.phen,
         famid, id, fa, mo, sex,
         mloci.out = paste(solar.output.directory, "mloci.out", sep = "/"),
         share.out = paste(solar.output.directory, "share.out", sep = "/"))

## Call multic with a bivariate model and one covariate
trait1.2 <-
  multic(cbind(trait1, trait2) ~ age,
         data = ped.phen,
         famid, id, fa, mo, sex,
         mloci.out = paste(solar.output.directory, "mloci.out", sep = "/"),
         share.out = paste(solar.output.directory, "share.out", sep = "/"))

## See help(multic) for more example multic calls.
