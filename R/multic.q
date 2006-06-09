multic <- function(formula,
                   data = sys.parent(),
                   famid,
                   id,
                   dadid,
                   momid,
                   sex,                   
                   mloci.out = NULL,
                   share.out = "kinship",
                   longitudinal = FALSE,
                   subset = NULL,
                   ascertainment = NULL,
                   control = multic.control(...),
                   ...
                   )
{
  ## ################ ##
  ## Manipulate Input ##
  ## ################ ##
  
  ## Parse the values we want out of the data object
  call <- match.call(expand = FALSE)
  call <- call[match(c("", "formula", "data", "famid", "id", "dadid",
                       "momid", "sex", "ascertainment", "subset"),
                     names(call), nomatch=0)]
  call[[1]] <- as.name("model.frame")
  ## We need to alter the name subset to another name (subset.) in order to
  ## prevent the model.frame evaluation from doing the subset automatically.
  new.subset.name <- "subset."
  names(call)[match("subset", names(call))] <- new.subset.name

  ## as.factor is needed to get character data evaluated in the data object
  ## R only allows integer, logical, and numeric
  if(using.R()) {
    call$famid <- substitute(as.factor(famid))
    call$id <- substitute(as.factor(id))
    call$dadid <- substitute(as.factor(dadid))
    call$momid <- substitute(as.factor(momid))
  }

  ## Use as.factor for sex during S-Plus and R to make searching for invalid
  ## sex codes uniform across platforms
  call$sex <- substitute(as.factor(sex))
  
  ## This is needed to not eliminate NA's
  call <- lappend(call, "na.action", na.pass)

  model.frame <- eval(call, sys.parent())
  ## It would be nice to add code to fix Curt's case of
  ## multic(famid, id, fa, mo, sex,
  ##        formula = cbind(log(lipo.all1), log(lipo.all2))
  ##                   ~ sex + exf.age1 + sex + exf.age2,
  ##        data = rpdn, share.out = "multicInput/share.out",
  ##        longitudinal = T)
  ## Both sex's get combine into one during the model.frame eval call.
  
  ## Extract famid, id, etc from our model.frame
  famid <- as.character(model.extract(model.frame, "famid"))
  id <- as.character(model.extract(model.frame, "id"))
  dadid <- as.character(model.extract(model.frame, "dadid"))
  momid <- as.character(model.extract(model.frame, "momid"))
  sex <- model.extract(model.frame, "sex")
  if( !missing(ascertainment) ){
    ascertainment <- model.extract(model.frame, "ascertainment")
    ## add valid/invalid ascertainment value checks
  }

  ## Check sex for invalid codes
  valid.sex.codes <- c("m", "M", "f", "F", "1", "2")
  invalid.sex.codes <- is.na(match(levels(sex), valid.sex.codes))
  if(any(invalid.sex.codes)) {
    stop("\n", paste(levels(sex)[invalid.sex.codes], collapse = " "),
         " is/are not valid sex codes.",
         "\nmultic.q key 64")
  }
  sex <- as.integer(sex)
  
  ## Print the call for the user to know where the program is
  cat("Call:\n")
  dput(match.call())
  cat("\n")

  ## ########### ##
  ## Error Check ##
  ## ########### ##

  ## Check famid, id, dadid, momid, and sex to see if they have value - 
  ## (not null nor NA)
  check.for.full.data(famid)
  check.for.full.data(id)
  check.for.full.data(dadid)
  check.for.full.data(momid)
  check.for.full.data(sex)

  ## If the user supplied a data argument, validate it.
  if(is.data.frame(data)) {
    check.data.for.missing.values(data, substitute(data))
  }
  
  ## Prepare share.out for reading (if we are using it)
  share.out.orig <- share.out
  if(share.out != "kinship") {
    share.out <- gunzip(share.out, copy = TRUE)
    
    file.size <- unlist(strsplit(trim(multic.system(paste("wc -c", share.out))), " "))[1]
    if(as.integer(file.size) == 0) {
      clean()
      stop(paste("\n", share.out.orig, " is empty\n", sep = ""))
    }
  }
  
  ## Run some tests on mloci.out
  mloci.out.orig <- mloci.out
  if(!is.null(mloci.out)) { 
    ## Prepare mloci.out for reading
    mloci.out <- gunzip(mloci.out, copy = TRUE)
    
    ## Check mloci.out for non-zero size (if we are using it)
    file.size <- unlist(strsplit(trim(multic.system(paste("wc -c", mloci.out))), " "))[1]

    if(as.integer(file.size) == 0) {
      clean()
      stop(paste("\n", mloci.out.orig, " is empty\n", sep = ""))
    }

    ## Check if it is has # lines
    ## n = 1 is the default for Splus, but n = -1 is the default for R
    if("#" != substring(readLines(mloci.out, n = 1), 1, 1)) {
      clean()
      stop(paste("\n", mloci.out.orig, " is not of proper format.\n",
                 "It does not separate ibds by '#'s.\n",
                 "multic.q key 127\n", sep = ""))
    }

    ## check if number of #'s is a factor of the total number of lines
  }
  
  ## ############## ##
  ## Variable Setup ##
  ## ############## ##

  start.time <- date()
  call <- match.call()
  
  ## Extract the values from the control object
  epsilon <- control$epsilon
  max.iterations <- control$max.iterations
  boundary.fix <- control$boundary.fix
  constraints <- control$constraints
  initial.values <- control$initial.values
  save.output.files <- control$save.output.files
  method <- control$method
  calc.fam.log.liks <- control$calc.fam.log.liks
  calc.residuals <- control$calc.residuals
  keep.input <- control$keep.input
  
  method.value <-
    switch(casefold(method),
           leastsq = 1,
           multic = 2,
           maxfun = 4,
           emvc = 5,
           stop(paste("\nThe estimation method '", method,
                      "' is not a valid method.\n",
                      "Please select 'multic', 'leastsq', ",
                      "'maxfun', or 'emvc'.\n",
                      "multic.q key 60\n", sep = ""))
           )
  if(!is.na(match(method.value, c(1, 4, 5)))) {
    stop(paste("\nThe estimation method \"", method,
               "\" has not been implemented.\n",
               "Please contact Dr. de Andrade (mandrade@mayo.edu) ",
               "\nto encourage its implementation",
               "\nmultic.q key 64\n", sep = ""))
  }
  constraints <-
    matrix(constraints, ncol = 1,
           dimnames = list(c("mu", "poly", "mg", "env", "sib.sib",
             "par.par", "par.off"), "constraint"))
  
  use.ascertainment <- !missing(ascertainment)
  
  ## Determine if this multic call is just to calcuate the polygene model
  ## or both polygene model and major gene models. 
  run.alternative.hyps <- !is.null(mloci.out)
  
  traits <- as.matrix(model.extract(model.frame, 'response'))
  trait.count <- ncol(traits)
  ## formula[[2]] returns only the traits/response text and the [-1]
  ## removes the cbind
  if(trait.count > 1) {
    trait.names <- as.character(formula[[2]][-1])
  } else {
    trait.names <- as.character(formula)[2]
  }
  dimnames(traits) <- list(paste(famid, id, sep = "-"), trait.names)
  ## The following code has been replaced by the above code.  I'll hold on
  ## to the code for a while to prove the replacement works
  ## trait.names <- parse.trait.names(names(model.frame)[1], trait.count)
  
  covariates <- model.matrix(terms(formula, data = data),
                             model.frame)[, -1, drop = FALSE]
  if(length(covariates) == 0) {
    covariates <- NULL
  }  
  ## Get the covariate.count and covariate.names
  if(is.null(covariates)) {
    covariate.count <- 0
    covariate.names <- NULL
  } else {
    covariate.count <- ncol(covariates)
    covariate.names <- dimnames(covariates)[[2]]
    dimnames(covariates) <- list(paste(famid, id, sep = "-"),
                                 covariate.names)
  }
  
  ## Apply the subset if it was specified
  ## "subset." is necessary because model.extract calls substitute before
  ## using it in the function.  Ideally we would like to use the previous
  ## variable new.subset.name.
  subset <- model.extract(model.frame, "subset.")
  if(!is.null(subset)) {
    if(length(subset) == length(famid)) {
      cat("Applying subset:", deparse(call$subset), "\n")
      traits[!subset, ] <- NA
      if(!is.null(covariates)) {
        covariates[!subset, ] <- NA
      }
    } else {
      stop("\nThe length of subset, ", length(subset), ", does not match ",
           "the length of the number\nof observations, ", length(famid),
           ".\nmultic.q key 155\n")
    }
  }
  
  ## If longitudinal, verify the correct number of covariates vs traits
  if(longitudinal) {
    if(covariate.count / trait.count
       != as.integer(covariate.count / trait.count)) {
      stop(paste("\nThe total number of covariates (",
                 covariate.count,
                 ") must be a multiple of the number of\n",
                 "traits (or timepoints) (", trait.count, ")\n", sep = ""))
    }
    
    ## The amount of values per random effect depends on only the number of
    ## traits (for multivariate) or the number of time-points (for longitudinal).
    ## Also, if we are doing a longitudinal study, the number of columns for
    ## trait and covariate is not actually the number of traits and covariates
    ## to be used in the study, but rather they are repeat.count times larger.
    ## So we also correct that here.
    repeat.count <- trait.count
    trait.count <- trait.count / repeat.count
    covariate.count <- covariate.count / repeat.count
    random.effects.count <- repeat.count * (repeat.count + 1) / 2
  } else {
    repeat.count <- 1
    random.effects.count <- trait.count * (trait.count + 1) / 2
  }
  
  ## Run the error checks and lm out here.  This needs to be after the subset
  ## because the subset might eliminate someone who can pass this test.
  lm.ready.check(traits, trait.names, covariates, covariate.names,
                 covariate.count)

  ## Calculate the initial.fit once for both calculate.initial.values and 
  ## calculate.coefficients
  if(is.null(covariates)) {
    initial.fit <- lm(traits ~ 1)
  } else {
    ## If longitudinal == TRUE, only the first time-point's covariates go
    ## into this lm otherwise, they all go
    initial.fit <- lm(traits ~ covariates[, 1:covariate.count])
  }
  
  ## Calculate the initial.values, the starting values for multic
  null.initial.values <-
    calculate.initial.values(initial.fit, covariates, trait.count,
                             random.effects.count, constraints,
                             repeat.count, longitudinal, initial.values)

  ## Calculate the coefficients
  null.coefficients <-
    calculate.coefficients(initial.fit, covariates, trait.count, trait.names,
                           covariate.count, covariate.names, longitudinal)

  ## ############################################################# ##
  ## Align the pedigree and phenotype data with the share.out data ##
  ## ############################################################# ##

  ## If we are using kinship for our phi values, reorder the phenotypes
  ## based on mloci.out, unless we are only doing polygenic, then do nothing.
  ## If we are using share.out for our phi values, reorder the phenotypes
  ## based on share.out
  if(share.out == "kinship") {
    if(length(unique(famid)) != get.family.count(famid, id)) {
      stop("\nIf you are running multic with bootstrapped data, you must ",
           "provide a\n",
           "similarly bootstrapped share.out (see expand.multic and ",
           "expand.share).\n",
           "If you are not running multic with bootstrapped data, the ",
           "family\n",
           "members must be congruent in the data set.\n",
           "multic.q key 266")
    }
    if(run.alternative.hyps) {
      share.order <-
        get.share.order(mloci.out, length(id), get.family.sizes(famid, id),
                        mloci.out.orig, using.mloci.out = TRUE);
    } else {
      share.order <- paste(famid, id, sep = "-")
    }
  } else {
    ## Get the order specified in share.out
    share.order <-
      get.share.order(share.out, length(id), get.family.sizes(famid, id),
                      share.out.orig)
  }
  ## Under some odd circumstances this match would die because it said that 
  ## share.order was not a vector.  I can't explain it.  So, that's the
  ## reason it is wrapped in as.vector.  Eric Lunde. 2005-09-13
  fort.reorder <- match(as.vector(share.order), paste(famid, id, sep = "-"))
  
  if(any(is.na(fort.reorder))) {
    stop(paste("\nMatching famid-id against '", share.out.orig, "' failed.",
               "\nThis probably means that '", share.out.orig,
               "'\nand the pedigree data supplied are not describing the ",
               "same data.\n", 
               "multic.q key 262", sep = ""))
  }
  
  ## Reorder the ped/phen data based on the order from share.out
  famid <- famid[fort.reorder]
  id <- id[fort.reorder]
  dadid <- dadid[fort.reorder]
  momid <- momid[fort.reorder]
  sex <- sex[fort.reorder]
  traits <- matrix(traits[fort.reorder, ],
                   ncol = trait.count * repeat.count,
                   nrow = length(traits) / (trait.count * repeat.count),
                   dimnames = dimnames(traits))
  if( !is.null(covariates) ) {
    covariates <- matrix(covariates[fort.reorder, ],
                         ncol = covariate.count * repeat.count,
                         dimnames = dimnames(covariates))
  }
  if( !missing(ascertainment) ) {
    ascertainment <- ascertainment[fort.reorder]
  }

  ## Create the kinship file to replace share.out here, after the reordering
  ## has taken place
  if(share.out == "kinship") {
    make.kinship.file(famid, id, dadid, momid, share.out,
                      run.alternative.hyps, mloci.out)
  }

  ## This is sort of silly to check share.out after its been alligned to data
  ## This may be just repeatative of the 'if(any(is.na(fort.reorder))) {'
  ## from above.  Look into this.
  ## Check share.out and the input data for data consistancy.
  ##validate.share.out(share.out, famid, id, dadid, momid, sex, share.out.orig)
  
  ## If any trait or covariate is NA for an individual, make all of
  ## that individual's traits and covariates NA
  trait.covariate <- cbind(traits, covariates)
  remove <- apply(is.na(trait.covariate), 1, any)
  traits[remove, ] <- NA
  if(!is.null(covariates)) {
    covariates[remove, ] <- NA
  }
  
  ## Prepare the trait and covariate data to be passed to the compiled
  ## portion of multic (replace the NA's with -9, multic.cpp's missing value)
  traits.with.miss.val <- traits
  traits.with.miss.val[is.na(traits.with.miss.val)] <- missing.value()
  covariates.with.miss.val <- covariates
  if( !is.null(covariates.with.miss.val) ) {
    covariates.with.miss.val[is.na(covariates.with.miss.val)] <-
      missing.value()
  }
  
  ## For the null hypothesis, set major gene 1 (and major gene 2 - even
  ## though the user knows nothing about major gene 2) to "Fixed"
  null.constraints <- c(constraints[1:3], "F", constraints[4:7])
  null.constraints[3] <- "F"
  
  family.sizes <- get.family.sizes(famid, id)
  unique.families <- get.unique.families(famid, id)
  family.count <- get.family.count(famid, id)

  ## ############################################### ##
  ## Fit the model without covariates (if necessary) ##
  ## ############################################### ##
  
  if(covariate.count > 0) {
    ## Clean up the directory before running multic
    remove.temp.files()
    
    cat("Fitting traits without covariates...\n")
    .C("multic",
       as.character(famid),
       as.character(id),
       as.character(dadid),
       as.character(momid),
       as.integer(sex),
       as.numeric(as.vector(traits.with.miss.val)),
       as.numeric(NULL),
       as.integer(ascertainment),
       as.integer(length(id)),
       as.integer(use.ascertainment),
       as.integer(2),
       as.numeric(epsilon),
       as.integer(boundary.fix),
       as.integer(method.value),
       as.integer(trait.count),
       covariate.count = as.integer(0),
       repeat.count = as.integer(repeat.count),
       as.character( c(trait.names, "marker1", NULL) ),
       as.numeric(missing.value()),
       initial.values = as.numeric(null.initial.values),
       as.character(null.constraints),
       as.integer(max.iterations),
       coefficients = as.numeric(NULL),
       as.character(share.out),
       print.progress = as.integer(0),
       calculate.residuals = as.integer(0),
       as.integer(family.sizes),
       as.character(unique.families),
       as.integer(family.count),
       PACKAGE = "multic"
       )
    ## I need to save the polygenic values but Mariza and Beth don't know
    ## what should be done for multivariate and longitudinal analysis
    total.variance <- get.variances.from.file(trait.count, 0)
    cat("\n")
  }
  
  ## ################################################## ##
  ## Calculate the null hypothesis - the polygene model ##
  ## ################################################## ##
  
  ## Clean up the directory before running multic    
  remove.temp.files()
  
  ## Write the unique family ids to the family log likelihood for alternative
  ## hypotheses file here, because we do not have access to those values when
  ## we write the numerical values to fam.likA - Eric Lunde 10-08-03
  ##if(calc.fam.log.liks) {
  ##  sink('fam.likA')
  ##  cat(unique(famid))
  ##  cat('\n')
  ##  sink()
  ##}
  
  ## Save the initial.values for the null hypothesis to help generate
  ## multic.par for the original mutltic to compare results.
  write.initial.values(trait.names, covariate.names, null.initial.values,
                       trait.count, random.effects.count, null.coefficients)
  
  ## This run of multic is to calculate the null hypothesis,
  ## hence the run option is hardcoded to 2
  null <-.C("multic",
            as.character(famid),
            as.character(id),
            as.character(dadid),
            as.character(momid),
            as.integer(sex),
            as.numeric(as.vector(traits.with.miss.val)),
            as.numeric(as.vector(covariates.with.miss.val)),
            as.integer(ascertainment),
            as.integer(length(id)),
            as.integer(use.ascertainment),
            as.integer(2),
            as.numeric(epsilon),
            as.integer(boundary.fix),
            as.integer(method.value),
            as.integer(trait.count),
            as.integer(covariate.count),
            as.integer(repeat.count),              
            as.character( c(trait.names, "marker1", covariate.names) ),
            as.numeric(missing.value()),
            initial.values = as.numeric(null.initial.values),
            as.character(null.constraints),
            as.integer(max.iterations),
            coefficients = as.numeric(null.coefficients),
            as.character(share.out),
            as.integer(1),
            as.integer(calc.residuals),
            as.integer(family.sizes),
            as.character(unique.families),
            as.integer(family.count),
            PACKAGE = "multic"
            )
  
  if(covariate.count > 0) {
    null.variance <- get.variances.from.file(trait.count, covariate.count)
    pro.var.due.to.cov <- 1 - (null.variance / total.variance)
  } else {
    pro.var.due.to.cov <-
      c("No covariates were used in this model.",
        "The proportion of variance due to the covariates cannot be calculated.")
  }
  
  ## ######################################################### ##
  ## Calculate the alternative hypotheses - the sporadic model ##
  ## ######################################################### ##
  
  alt.initial.values <- NULL
  alt.coefficients <- NULL
  alt.hyp.count <- 0
  if(run.alternative.hyps) {
    ## Calculate the location of the major gene indices so we can alter the
    ## polygene and major gene starting values for the alternative hyptheses
    polygene.indices <-
      seq(from = trait.count + 1, length = random.effects.count)
    
    major.gene.1.indices <-
      seq(from = trait.count + random.effects.count + 1,
          length = random.effects.count)
    
    ## alt.initial.values are the initial values for the alternative
    ## hypotheses and like wise for the coefficients
    alt.initial.values <- null$initial.values
    if(covariate.count > 0) {
      alt.coefficients <-
        matrix(null$coefficients, nrow = covariate.count,
               dimnames = dimnames(null.coefficients))
    }
  
    ## Divide the final polygenic values in half amount the polygenic and
    ## major gene 1 influence    
    alt.initial.values[polygene.indices] <-
      alt.initial.values[polygene.indices]/2
    alt.initial.values[major.gene.1.indices] <-
      alt.initial.values[polygene.indices]
    alt.initial.values <- matrix(alt.initial.values, ncol = 1,
                                 dimnames = dimnames(null.initial.values))
    
    ## Save the initial.values for the alternative hypothesis to help generate
    ## multic.par for the original mutltic to compare results.
    write.alt.initial.values(null$initial.values, null$coefficients,
                             alt.initial.values,
                             alt.coefficients,
                             trait.count, random.effects.count)
    
    combine.null.and.alt.initial.values()
    
    ## Create 'tempmloci.out' because the function 'extractLoci' reads and
    ## modifies 'tempmloci.out' not 'mloci.out'
    multic.system(paste('cp', mloci.out, "tempmloci.out"))
    
    ## Create the new form of handling mloci.out and loci.outs
    cat(paste("\nSplitting ", mloci.out.orig, "...\n", sep = ""))
    loci.names <- mloci.split("tempmloci.out")
    multic.system('rm -f tempmloci.out')
    multic.system(paste('rm -f', mloci.out))
    
    alt.hyp.count <- length(loci.names)
    cat("\n")
    
    ## Alter constraints by adding in the place holding major gene 2
    alt.constraints <- c(constraints[1:3], "F", constraints[4:7])
    
    ## These runs of multic are to calculate the alternate hypotheses
    ## The run option is hardcoded to 1 because all of these multic calls
    ## are alternate hypotheses.
    sapply(loci.names, run.alternative.multic,
           famid,
           id,
           dadid,
           momid,
           sex,
           traits.with.miss.val,
           covariates.with.miss.val,
           ascertainment,
           epsilon,
           boundary.fix,
           method.value,
           trait.count,
           covariate.count,
           repeat.count,
           c(trait.names, "marker1", covariate.names),
           alt.initial.values,
           alt.constraints,
           max.iterations,
           alt.coefficients,
           share.out,
           calc.residuals,
           family.sizes,
           unique.families,
           family.count
           )
  }
  
  multic.system(paste('rm -f', share.out))
  cat("\n")
  
  iterations <- load.iterations("iterations.log")  
  
  ## ######################### ##
  ## Build the metadata object ##
  ## ######################### ##
  
  ## Setup metadata variable
  metadata <- list()
  
  metadata <- lappend(metadata, "start.time", start.time)
  metadata <- lappend(metadata, "call", call)
  metadata <- lappend(metadata, "epsilon", epsilon)
  metadata <- lappend(metadata, "boundary.fix", boundary.fix)
  metadata <- lappend(metadata, "method", method)
  metadata <- lappend(metadata, "method.value", method.value)
  metadata <- lappend(metadata, "constraints", constraints)
  metadata <- lappend(metadata, "max.iterations", max.iterations)
  metadata <- lappend(metadata, "calc.fam.log.liks", calc.fam.log.liks)
  metadata <- lappend(metadata, "calc.residuals", calc.residuals)
  metadata <- lappend(metadata, "keep.input", keep.input)
  metadata <- lappend(metadata, "subset", subset)
  metadata <- lappend(metadata, "ascertainment", ascertainment)
  metadata <- lappend(metadata, "use.ascertainment", use.ascertainment)
  metadata <- lappend(metadata, "share.out", share.out.orig)
  metadata <- lappend(metadata, "mloci.out", mloci.out.orig)
  metadata <- lappend(metadata, "longitudinal", longitudinal)
  
  if(keep.input) {
    ped <- data.frame(famid, id, dadid, momid, sex)
    metadata <- lappend(metadata, "pedigree", ped)
    metadata <- lappend(metadata, "traits", traits)
    metadata <- lappend(metadata, "covariates", covariates)
  }
  
  metadata <- lappend(metadata, "trait.names", trait.names)
  metadata <- lappend(metadata, "covariate.names", covariate.names)  
  metadata <- lappend(metadata, "trait.count", trait.count)
  metadata <- lappend(metadata, "covariate.count", covariate.count)
  metadata <- lappend(metadata, "repeat.count", repeat.count)
  
  ## Remove mg2* so the user doesn't know about major gene 2
  non.mg2 <- regexpr("mg2..*", dimnames(null.initial.values)[[1]]) < 0
  null.initial.values <- null.initial.values[non.mg2, , drop = FALSE]
  metadata <- lappend(metadata, "null.initial.values", null.initial.values)
  metadata <- lappend(metadata, "null.coefficients", null.coefficients)
  
  if(run.alternative.hyps) {
    alt.initial.values <- alt.initial.values[non.mg2, , drop = FALSE]
  }
  metadata <- lappend(metadata, "alt.initial.values", alt.initial.values)
  metadata <- lappend(metadata, "alt.coefficients", alt.coefficients)
  
  metadata <- lappend(metadata, "alt.hyp.count", alt.hyp.count)
  metadata <- lappend(metadata, "iterations", iterations)
  
  ## ######################## ##
  ## Create the multic object ##
  ## ######################## ##
  
  multic.object <-
    create.multic.object(famid, id, sex, traits, covariates, metadata)
  multic.object <- lappend(multic.object, "R.sq", pro.var.due.to.cov)
  
  ## Add the finish time to the multic metadata
  metadata <- lappend(metadata, "finish.time", date())
  multic.object <- lappend(multic.object, "metadata", metadata )
  
  ## Remove the multic output files if desired
  if(!save.output.files) {
    remove.temp.files()
  }
  
  ## ######################## ##
  ## Return the multic object ##  
  ## ######################## ##
  return (multic.object)
}

create.multic.object <- function(famid, id, sex, traits, covariates,
                                 metadata)
{
  unique.famid <- get.unique.families(famid, id)
  family.sizes <- get.family.sizes(famid, id)
  
  ## Read the multic output files to build the multic object
  ## If traits is not a matrix, make it a matrix, because we want to use the
  ## dimnames property to get the names of the traitss.
  trait.count <- metadata$trait.count
  trait.names <- metadata$trait.names
  
  covariate.count <- metadata$covariate.count
  covariate.names <- metadata$covariate.names  
  repeat.count <- metadata$repeat.count
  alt.hyp.count <- metadata$alt.hyp.count
  
  if(repeat.count == 1) {
    random.effects.count <- trait.count * (trait.count + 1) / 2
  } else if(repeat.count > 1) {
    random.effects.count <- repeat.count * (repeat.count + 1) / 2
  }
  
  run.alternative.hyps <- alt.hyp.count != 0
  
  ## Calculate the amount of fixed effects, the number of covariates plus the
  ## trait mean for each trait
  fixed.effects.count <- trait.count * (1 + covariate.count)
  
  ## Load the fixed and random effects from the multic output file
  ## summary.log and the distances from the distance file (this distance file
  ## does not have a pre-determined name)
  cat("Loading (fixed, random) effects from file...\n")
  effects <- load.effects(fixed.effects.count, random.effects.count,
                          alt.hyp.count + 1, length(id), trait.count,
                          covariate.count, repeat.count,
                          run.alternative.hyps)
  
  fixed.effects <- effects$fixed.effects
  polygenic <- effects$polygenic
  major.gene1 <- effects$major.gene1
  environmental <- effects$environmental
  sibling.sibling <- effects$sibling.sibling
  parent.parent <- effects$parent.parent
  parent.offspring <- effects$parent.offspring
  log.liks <- effects$log.liks
  
  ## Save the titles of the random effect variables
  polygenic.names <- effects$polygenic.names
  major.gene1.names <- effects$major.gene1.names
  environmental.names <- effects$environmental.names
  sibling.sibling.names <- effects$sibling.sibling.names
  parent.parent.names <- effects$parent.parent.names
  parent.offspring.names <- effects$parent.offspring.names
  
  random.effects.names <- c(polygenic.names, major.gene1.names,
                            environmental.names)
  
  ## Save the names used for the fixed effects and the ibd file names for
  ## future use
  fixed.effects.names <- dimnames(fixed.effects)[[1]]
  ibd.names <- dimnames(fixed.effects)[[3]]
  
  ## Load the family log likelihoods from the multic output files fam.lik0
  ## and fam.likA
  if(metadata$calc.fam.log.liks) {
    cat("Loading family log likelihoods from file...\n")
    fam.log.liks <- load.family.log.likelihoods(get.family.count(famid, id),
                                                ##length(unique(famid)),
                                                alt.hyp.count + 1)
  }else {
    fam.log.liks <- c("Family log likelihoods were not calculated during this call of multic.",
                      "To receive family log likelihoods, specify 'calc.fam.log.liks = T'",
                      "in the call to multic.","multic.q key 619")
  }
  
  ## Load the inverse of the expected second derivitive for the fixed effects
  ## from the null hypothesis and all of the alternative hypotheses that are
  ## stored in the multic output file 'invExpSecDerFixed.log' 
  cat("Loading inverse of the expected second derivative from file...\n")
  inv.exp.sec.der.fixed <- load.inv.exp.sec.der.fixed(fixed.effects.count,
                                                      alt.hyp.count + 1,
                                                      fixed.effects.names,
                                                      ibd.names)
  
  ## Load the inverse of the expected second derivitive for the random
  ## effects from the null hypothesis and all of the alternative hypotheses
  ## that are stored in the multic output file 'invExpSecDerRandom.log'
  inv.exp.sec.der.random <- load.inv.exp.sec.der.random(random.effects.count,
                                                        alt.hyp.count,
                                                        random.effects.names,
                                                        ibd.names)
  
  if(trait.count == 1 &&
     (metadata$repeat.count == 1 || metadata$method.value == 1)) {
    var.sandwich.names <- NULL
    if(metadata$constraints[2] != "F") {
      var.sandwich.names <- append(var.sandwich.names, polygenic.names)
    }
    if(metadata$constraints[3] != "F") {
      var.sandwich.names <- append(var.sandwich.names, major.gene1.names)
    }
    if(metadata$constraints[4] != "F") {
      var.sandwich.names <- append(var.sandwich.names, environmental.names)
    }
    if(metadata$constraints[5] != "F") {
      var.sandwich.names <- append(var.sandwich.names, sibling.sibling.names)
    }
    if(metadata$constraints[6] != "F") {
      var.sandwich.names <- append(var.sandwich.names, parent.parent.names)
    }
    if(metadata$constraints[7] != "F") {
      var.sandwich.names <- append(var.sandwich.names,
                                   parent.offspring.names)
    }
    cat("Loading variance-covariance sandwich matrix from file..\n")
    var.sandwich <- load.var.sandwich(random.effects.count,
                                      alt.hyp.count,
                                      var.sandwich.names,
                                      ibd.names, metadata$constraints)
  }else {
    var.sandwich <- c("No robust variance-covariance sandwich matrix",
                      "was calculated for this multic object.")
  }
  
  ## Calculate the indices that will be used when calculating the upcoming
  ## h^2 calculation.  It is the major diagonal of the random effects.
  ## Since the number of major diagonal values depends on either the number
  ## of traits (multivariate) or repeats (longitudinal) exclusively without
  ## dependance on the other and the value that does not determine the number
  ## of random effects is always 1 (trait.count = 1 in longitudinal and
  ## repeat.count = 1 in multivariate), multiplying these values results in
  ## the maximum of the two and thus the correct number of random effects.
  major.diagonal.indices <- (random.effects.count + 1) -
    cumsum(1:(trait.count * repeat.count))[1:(trait.count * repeat.count)]
  ## Reverse the indices so that they are in increasing order
  major.diagonal.indices <-
    major.diagonal.indices[length(major.diagonal.indices):1]
  
  ## Now, since we've loaded the inverse of the expected second derivative for
  ## the random effects, we can calculate the standard error of the
  ## heritability.  Although, we didn't need to, we also saved the calculation
  ## of the heritabily until now. 
  cat("Calculating heritability values...\n")
  effects.plus.herit <- add.heritability(polygenic, major.gene1,
                                         environmental,
                                         sibling.sibling,
                                         parent.parent,
                                         parent.offspring,
                                         major.diagonal.indices,
                                         inv.exp.sec.der.random,
                                         run.alternative.hyps)
  polygenic <- effects.plus.herit$polygenic
  major.gene1 <- effects.plus.herit$major.gene1

  trait.covariate <- cbind(traits, covariates)
  ## not.nas is T when the value is not missing and F when the value is
  ## missing (NA)
  not.nas <- !is.na(trait.covariate)
  ## has.all.values are the individuals who have all of their trait and
  ## covariate values non missing (not NA)
  has.all.values <- apply(not.nas, 1, all)
  
  if(metadata$calc.residuals) {
    ## I first thought that I had deleted code that I shouldn't have, rather
    ## I have keep uneccessary code just because I liked writing it and felt
    ## bad deleteing it.  This .Call and getResidualPlacement.cpp serve no
    ## purpose other than satisfying code writing - Eric Lunde.
    residual.placement1 <- .Call("getResidualPlacement",
                                 as.integer(family.sizes),
                                 as.logical(has.all.values),
                                 as.integer(alt.hyp.count + 1),
                                 PACKAGE = "multic")
    
    cat("Loading V matrices from file...\n")
    v.matrices <- load.v.matrices(unique.famid)
    
    cat("Loading residuals from file...\n")
    y.beta.diffs <- load.y.beta.diffs(unique.famid)
  }else {
    v.matrices <- c("V matrices were not calculated during this call of multic",
                    "To receive V matrices specify 'calc.residuals = T' in the call to multic.",
                    "multic.q key 700")
    y.beta.diffs <- c("Residuals were not calculated during this call of multic.",
                      "To receive residuals specify 'calc.residuals = T' in the call to multic.",
                      "multic.q key 704") 
  }
  
  ## Create cors (correlations) and give it the pearson and spearman
  ## correlations for traits and covariates.  If there is more than one trait,
  ## calculate the other correlations
  cat("Calculating correlation values...\n")
  cors <- calculate.cor.values(traits, covariates, polygenic, environmental,
                               alt.hyp.count)
    
  ## Calculate mean, std dev, kurtosis, skewness, and other such statistics
  subset.cov <- NULL
  if(!is.null(covariates)) {
    subset.cov <- covariates[has.all.values, , drop = FALSE]    
  }
  descriptives <-
    calculate.descriptives(traits[has.all.values, , drop = FALSE],
                           subset.cov, trait.names, covariate.names)  
  
  men.with.no.missing.values <- has.all.values & sex == 1
  women.with.no.missing.values <- has.all.values & sex == 2
  
  ## Count the number of probands (if any)
  proband.count <- ifelse(metadata$use.ascertainment,
                          sum(metadata$ascertainment == 1),
                          0)
  
  ## Combine all of the counting variables into one vector
  counts <- array(c(get.family.count(famid, id),
                    length(id),
                    length(sex[sex==2]),
                    length(sex[sex==1]),
                    proband.count,
                    sum(women.with.no.missing.values)
                    + sum(men.with.no.missing.values),
                    sum(women.with.no.missing.values),
                    sum(men.with.no.missing.values),
                    trait.count,
                    covariate.count,
                    alt.hyp.count),
                  dimnames=list(c("Pedigrees", "People", "Females","Males",
                    "Probands", "People", "Females", "Males", "Traits",
                    "Covariates", "Locations")))
  
  if(run.alternative.hyps) {
    if(multic.strsplit(metadata$mloci.out, sep = '/')[1] != "") {
      metadata$mloci.out <- paste(multic.system('pwd'), metadata$mloci.out, sep = '/')
    }
  }
  if(multic.strsplit(metadata$share.out, sep = '/')[1] != "") {
    metadata$share.out <- paste(multic.system('pwd'), metadata$share.out, sep = '/')
  }
  
  
  ## Combine all of the multic data into a single list object and return it
  multic.object <- list(fam.log.liks = fam.log.liks,
                        fixed.effects = fixed.effects,
                        polygenic = polygenic,
                        major.gene1 = major.gene1,
                        environmental = environmental,                  
                        sibling.sibling = sibling.sibling,
                        parent.parent = parent.parent,
                        parent.offspring = parent.offspring,
                        log.liks = log.liks,
                        var.fixed = inv.exp.sec.der.fixed,
                        var.random = inv.exp.sec.der.random,
                        var.sandwich = var.sandwich,
                        cors = cors,
                        v.matrices = v.matrices,
                        residuals = y.beta.diffs,
                        descriptives = descriptives,
                        counts = counts,
                        call = metadata$call
                        )
  
  ##  multic.object$metadata$mloci.out <- NULL
  oldClass(multic.object) <- 'multic'
  
  return (multic.object)
}

"[.multic" <- function(multic.obj, i, ..., drop = FALSE) {
  locus.name <- dimnames(multic.obj$log.liks)[[1]][i]
  
  locus.name.matches <- multic.obj$fam.log.liks == locus.name
  fam.log.liks <- multic.obj$fam.log.liks[locus.name.matches, ]
  
  if(is.null(multic.obj$var)) {
    var <- NULL
  } else {
    var <- multic.obj$var[, , i]
  }
  
  z <- list(fam.log.liks = fam.log.liks,
            fixed.effects = multic.obj$fixed.effects[, , i],
            polygenic = multic.obj$polygenic[, , i],
            major.gene1 = multic.obj$major.gene1[, , i],
            environmental = multic.obj$environmental[, , i],
            log.liks = multic.obj$log.liks[i, ],
            var.fixed = multic.obj$var.fixed[, , i],
            var.random = multic.obj$var.random[, , i],
            var = var,
            metadata = multic.obj$metadata,
            descriptives = multic.obj$descriptives,
            counts = multic.obj$counts,
            call = multic.obj$call
            )
  
  attributes(z) <- attributes(multic.obj)
  z
}

###############################################################################
## Function Name: get.family.parent.inconsistancies
## Description: Determines if all fathers in a family are indeed male and if all
##              mothers are female
## Input: family.id - scalar value determining the family number we want to
##                    delimit by in famid, id, dadid, momid, and sex
##        famid - vector of each individual's family id number
##        id - vector of each individual's id number, this may be globally
##             unique, in that no other person in any family has the same id or
##             it may be unique within its family.  It matters not to this
##             algorithm.
##        dadid - vector of each individual's father's id number
##        momid - vector of each individual's mother's id number
##        sex - vector of each individual's sex (1 = male, 2 = female)
## Ouput: (unique.father.ids[incorrect.father.ids],
##         unique.mother.ids[incorrect.mother.ids]) - a list of two vectors, the
##                                                    first contains the unique
##                                                    father ids that have been
##                                                    specified as fathers but
##                                                    are not male and the
##                                                    second, contains mother
##                                                    ids that have been
##                                                    specified as mothers, but
##                                                    are not female.
## Side Effects: NONE.
## Author: Eric Lunde, 9-08-03
###############################################################################
## get.family.parent.inconsistancies AND validate.parents ARE WRONG.  THEY
## NEED TO BE FIXED AND USED OR REMOVED FROM THE FILE.  THEY SHOULD NOT STAY
## HERE UNFUNCTIONAL.  (this is only the case for bootstrap'ed data sets.
## Eric Lunde 2005-08-31
get.family.parent.inconsistancies <- function(family.id, famid, id, dadid,
                                              momid, sex)
{
  ## Acquire the T/F vector to specify one family
  family <- family.id == famid
  
  ## Use that vector to narrow our field of ids, sex, dadids and momids to just
  ## that family
  family.ids <- id[family]
  family.sex <- sex[family]
  family.dad.ids <- dadid[family]
  family.mom.ids <- momid[family]
  
  ## Acquire the ids of fathers as specified by children
  father.ids <- family.dad.ids[family.dad.ids != 0]
  unique.father.ids <- unique(father.ids)

  ## Acquire the ids of mothers as specified by children
  mother.ids <- family.mom.ids[family.mom.ids != 0]
  unique.mother.ids <- unique(mother.ids)

  ## Get a T/F vector to see if those fathers are indeed male
  incorrect.father.ids <-
    family.sex[match(unique.father.ids, family.ids)] != 1

  ## Get a T/F vector to see if those mothers are indeed female
  incorrect.mother.ids <-
    family.sex[match(unique.mother.ids, family.ids)] != 2

  ## Return a list of ids that specify father/male and mother/female
  ## inconsistancies
  return (unique.father.ids[incorrect.father.ids],
          unique.mother.ids[incorrect.mother.ids])
}

###############################################################################
# Function Name: validate.parents
# Description: Determine all father/male and mother/female inconsistancies and
#              report them to the user
# Input: famid - vector of each individual's family id number
#        id - vector of each individual's id number, this may be globally
#             unique, in that no other person in any family has the same id or
#             it may be unique within its family.  It matters not to this
#             algorithm.
#        dadid - vector of each individual's father's id number
#        momid - vector of each individual's mother's id number
#        sex - vector of each individual's sex (1 = male, 2 = female)
# Ouput: NONE.
# Side Effects: If an inconsistancy is found, validate.parents prints those
#               errors to the screen and multic is halted via the stop
#               function.
# Author: Eric Lunde, 9-08-03
###############################################################################
## validate.parents AND get.family.parent.inconsistancies ARE WRONG.  THEY
## NEED TO BE FIXED AND USED OR REMOVED FROM THE FILE.  THEY SHOULD NOT STAY
## HERE UNFUNCTIONAL.  (this is only the case for bootstrap'ed data sets.
## Eric Lunde 2005-08-31
validate.parents <- function(famid, id, dadid, momid, sex) {
  ## We just want to define incorrect.ids and incorrect.family.ids
  incorrect.father.ids <- NULL
  incorrect.mother.ids <- NULL
  family.ids.with.incorrect.fathers <- NULL
  family.ids.with.incorrect.mothers <- NULL

  ## For each unique family
  unique.family.ids <- get.unique.families(famid, id)
  for(i in unique.family.ids) {
    ## Store the incorrect father and mother ids of the current family
    incorrect.parent.ids <-
      get.family.parent.inconsistancies(i, famid, id, dadid, momid, sex)
    
    ## current.father/mother.ids are the incorrect ids for this particular
    ## family
    current.father.ids <- incorrect.parent.ids[[1]]
    current.mother.ids <- incorrect.parent.ids[[2]]
    
    ## Append the new father and mother incorrect ids to the lists of all
    ## incorrect father and mother ids
    incorrect.father.ids <- append(incorrect.father.ids, current.father.ids)
    incorrect.mother.ids <- append(incorrect.mother.ids, current.mother.ids)
    
    ## Append the family number of the father ids to the list of problem
    ## families length(incorrect.father.ids) times so that a given index
    ## corresponds to the family id (incorrect.family.ids) and individual id
    ## (incorrect.father.ids)
    family.ids.with.incorrect.fathers <-
      append(family.ids.with.incorrect.fathers,
             rep(i, length(incorrect.father.ids)))

    ## Do the same for mothers
    family.ids.with.incorrect.mothers <-
      append(family.ids.with.incorrect.mothers,
             rep(i, length(incorrect.mother.ids)))
  }

  ## If there are errors, print them and exit program
  if(any(incorrect.father.ids) || any(incorrect.mother.ids)) {
    ## Get the unique family ids so we can specify which family we want to
    ## display the errors
    unique.incorrect.family.ids <-
      sort(unique(c(family.ids.with.incorrect.fathers,
                    family.ids.with.incorrect.mothers)))
    
    ## report.incorrect.family is defined so we can use apply across the unique
    ## family ids instead of for-looping
    report.incorrect.family <- function(family.id,
                                        family.ids.with.incorrect.fathers,
                                        incorrect.father.ids,
                                        family.ids.with.incorrect.mothers,
                                        incorrect.mother.ids)
    {
      ## Get the incorrect father ids of the family specified by
      ## incorrect.family.id
      incorrect.father.ids.by.family <-
        incorrect.father.ids[family.id == family.ids.with.incorrect.fathers]

      ## Get the incorrect mother ids of the family specified by
      ## incorrect.family.id
      incorrect.mother.ids.by.family <-
        incorrect.mother.ids[family.id == family.ids.with.incorrect.mothers]

      ## Print the family id and the incorrect father ids associated with that
      ## family
      if( !is.na(incorrect.father.ids.by.family)
         && length(incorrect.father.ids.by.family) != 0) {
        cat("In family", family.id, ": persons",
            incorrect.father.ids.by.family, "\nare not defined as male when",
            "children define them as fathers.\n\n")
      }
      
      ## Print the family id and the incorrect mother ids associated with that
      ## family
      if( !is.na(incorrect.mother.ids.by.family)
         && length(incorrect.mother.ids.by.family) != 0) {
        cat("In family", family.id, ": persons",
            incorrect.mother.ids.by.family, "\nare not defined as female when",
            "children define them as mothers.\n\n")
      }
    }
    
    apply(array(unique.incorrect.family.ids), 1, report.incorrect.family,
          family.ids.with.incorrect.fathers, incorrect.father.ids,
          family.ids.with.incorrect.mothers, incorrect.mother.ids)
    
    stop('\nExiting multic.s\n')
  }
}

## If we want to use the sequential ids instead of the large ids that SOLAR
## assigned, we must read 'pedindex.cde', locate the field that tells us
## where to find the sequential id in 'pedindex.out', and make all the
## appropriate translations  
apply.sequential.ids <- function(id, dadid, momid) {
  ## Read pedindex.cde for the locations of sequential and original
  ## identifiers
  pedindex.cde <- importData('pedindex.cde', type='ASCII')
  
  ## id.index is the row of pedindex.cde that holds the id "ID" in the
  ## second column
  id.index <- match('ID', pedindex.cde[2])
    
  ## Calculate the sum of adding all the rows of column 1 up to (but not
  ## including) the row that holds the id "ID"
  origIdStart <- sum(pedindex.cde[seq(1, id.index-1), 1])
  
  ## Find how many characters are in the field labeled by "ID"
  origIdLen <- pedindex.cde[match('ID', pedindex.cde[2]), 1]
  
  ## Find how many characters are in the field labeled by "IBDID"
  seqIdLen <- pedindex.cde[1, match('IBDID', pedindex.cde[2])]
  
  ## Read pedindex.out for the sequential and original identifiers
  form <- paste("%", seqIdLen, "f %", origIdStart-seqIdLen, "* %", origIdLen,
                "f", sep='')
  ids <- importData('pedindex.out', type='FASCII', format=form)
  
  ## Create translation from original id to sequential id function
  translateId <-  function(originalId, idTable) {
    if(originalId == 0) {
      return (0)
    }
    return ( idTable[match(originalId, idTable[ ,2]), 1])
  }
  
  ## Sort ids - not necessary, consider deleting, Eric Lunde 9-03
  orderedIds <- ids[order(ids[ ,2]), ]
  
  ## Reassign the id, dadid, and momid to have the sequential identifiers
  ## seqIds <- apply(array(id), 1, translateId, orderedIds)
  ## seqDadIds <- apply(array(dadid), 1, translateId, orderedIds)
  ## seqMomIds <- apply(array(momid), 1, translateId, orderedIds)
  id <- apply(array(id), 1, translateId, orderedIds)
  dadid <- apply(array(dadid), 1, translateId, orderedIds)
  momid <- apply(array(momid), 1, translateId, orderedIds)

  return(id, dadid, momid)
}

## upp.tri.as.vector is a small utility function to extract from a matrix
## object (x) the upper triangle and return the values in a vector.
upp.tri.as.vector <- function(x) {
  x <- as.matrix(x)

  ## If x contains only one element, return it as a vector.
  if(length(x) == 1) {
    return ( as.vector(x) )
  }

  ## Copy the entire first row
  result <- x[1,]

  ## For the rest of the rows, only append the entries from the major diagonal
  ## leftward
  for(i in 2:dim(x)[1]) {
    result <- append(result, x[i, -(1:(i-1))])
  }
  return (result)
}

## get.n returns the integer value of the number of individuals with
## non-missing values in the given x vector (this is used with x being equal
## to the trait or covariate data vectors
get.n <- function(x) {
  return (sum(as.integer(!is.na(x))))
}

#############################################################################
## calculate.initial.values creates a single-column matrix with the initial
## values for the polygene-multic calculation.
##
## If the user specified the initial values via the multic.control object
## or ... parameter, calculate.initial.values populates its result with
## those values.  Basically, it copies the necessary values off the front
## of the control.initial.values argument and then removes those values from
## the control.initial.values argument.  This is a complicated way of simple
## copy or assignment command.  All this extra work is done to make use of
## the way calculate.initial.values creates the dimnames for the matrix.  So,
## if it looks more complicated than it ought to be, that is the reason.
#############################################################################
calculate.initial.values <- function(initial.fit, covariates,
                                     trait.count, random.effects.count,
                                     constraints, repeat.count,
                                     longitudinal, control.initial.values)
{
  ## If the initial.values were provided by the user, make sure they are the
  ## correct length.
  if( !is.null(control.initial.values) ) {
    correct.length = trait.count + random.effects.count * 6
    if(length(control.initial.values)
       != correct.length) {
      stop("\nThe user-specified initial.values does not have the\n",
           "correct length.  It is ", length(control.initial.values),
           " and shoud be ", correct.length, ".\n",
           "multic.q key 1244\n")
    }
  }

  ## Set up the initial.values vector with the mu values
  if( is.null(control.initial.values) ) {
    if(is.null(covariates)) {
      ## mu
      initial.values <- coefficients(initial.fit)
    } else {
      ## R returns a vector if trait.count == 1, S-PLUS return a matrix
      initial.values <- as.matrix(coefficients(initial.fit))[1, ]
    }
    if(longitudinal) {
      initial.values <- initial.values[1]
    }
  } else {
    ## mu
    initial.values <- control.initial.values[1:trait.count]
    control.initial.values <-
      control.initial.values[-(1:(trait.count*repeat.count))]
  }
  
  ## Find out how many non "F"'s there are, not counting mu or mg
  constraint.count <- sum(constraints[c(-1, -3)] != "F")
  if(using.R()) {
    var <- var(residuals(initial.fit), na.rm = TRUE)
  } else {
    var <- var(residuals(initial.fit), na.method = 'omit')
  }
  values <- upp.tri.as.vector(var / constraint.count)
  zeroes <- rep(0, random.effects.count)
  initial.values.names <- paste("mu", 1:trait.count, sep = "")
  var.covar.labels <- get.var.covar.labels(trait.count * repeat.count)

  ## poly
  if( is.null(control.initial.values) ) {
    initial.values <- append(initial.values,
                             ifelse(rep(constraints[2] != "F",
                                        random.effects.count),
                                    values,
                                    zeroes))
  } else {
    initial.values <-
      append(initial.values, control.initial.values[1:random.effects.count])
    control.initial.values <-
      control.initial.values[-(1:random.effects.count)]
  }
  initial.values.names <- append(initial.values.names,
                                 paste("poly", var.covar.labels, sep = ""))
  ## mg
  if( is.null(control.initial.values) ) {
    initial.values <- append(initial.values, zeroes)
  } else {
    initial.values <-
      append(initial.values, control.initial.values[1:random.effects.count])
    control.initial.values <-
      control.initial.values[-(1:random.effects.count)]
  }
  initial.values.names <- append(initial.values.names,
                                 paste("mg", var.covar.labels, sep = ""))
  ## mg2
  initial.values <- append(initial.values, zeroes)
  initial.values.names <- append(initial.values.names,
                                 paste("mg2", var.covar.labels, sep = ""))
  
  ## env
  if( is.null(control.initial.values) ) {
    initial.values <- append(initial.values,
                             ifelse(rep(constraints[4] != "F",
                                        random.effects.count),
                                    values,
                                    zeroes))
  } else {
    initial.values <-
      append(initial.values, control.initial.values[1:random.effects.count])
    control.initial.values <-
      control.initial.values[-(1:random.effects.count)]
  }
  initial.values.names <- append(initial.values.names,
                                 paste("env", var.covar.labels, sep = ""))
  ## sib
  if( is.null(control.initial.values) ) {
    initial.values <- append(initial.values,
                             ifelse(rep(constraints[5] != "F",
                                        random.effects.count),
                                    values,
                                    zeroes))
  } else {
    initial.values <-
      append(initial.values, control.initial.values[1:random.effects.count])
    control.initial.values <-
      control.initial.values[-(1:random.effects.count)]
  }
  initial.values.names <- append(initial.values.names,
                                 paste("sib.sib", var.covar.labels, sep = ""))
  ## pp
  if( is.null(control.initial.values) ) {
    initial.values <- append(initial.values,
                             ifelse(rep(constraints[6] != "F",
                                        random.effects.count),
                                    values,
                                    zeroes))
  } else {
    initial.values <-
      append(initial.values, control.initial.values[1:random.effects.count])
    control.initial.values <-
      control.initial.values[-(1:random.effects.count)]
  }
  initial.values.names <- append(initial.values.names,
                                 paste("par.par", var.covar.labels, sep = ""))
  ## po
  if( is.null(control.initial.values) ) {
    initial.values <- append(initial.values,
                             ifelse(rep(constraints[7] != "F",
                                        random.effects.count),
                                    values,
                                    zeroes))
  } else {
    initial.values <-
      append(initial.values, control.initial.values[1:random.effects.count])
    control.initial.values <-
      control.initial.values[-(1:random.effects.count)]
  }
  initial.values.names <- append(initial.values.names,
                                 paste("par.off", var.covar.labels, sep = ""))

  initial.values <-
    matrix(initial.values, ncol = 1,
           dimnames = list(initial.values.names, "initial.value"))

  return (initial.values)
}

calculate.coefficients <- function(initial.fit, covariates, 
                                   trait.count, trait.names,
                                   covariate.count, covariate.names,
                                   longitudinal) {
  coefficients <- NULL
  
  if(is.null(covariates)) {
    return(NULL)
  } else {
    ## Could these be combined ?
    if(longitudinal) {
      coefficients <- coefficients(initial.fit)[-1, 1]
    } else {
      ## R returns a vector if trait.count == 1, S-PLUS return a matrix
      coefficients <- as.matrix(coefficients(initial.fit))[-1, ]
    }
  }

  ## Check that the coefficients is the correct length
  coefficients.correct.length <- (trait.count * covariate.count)
  if(length(coefficients) != coefficients.correct.length) {
    print(coefficients)
    stop(paste("\nThe variable 'coefficients' is not of the correct length\n",
               "coefficients length = ", length(coefficients), "\n",
               "The correct length is ", coefficients.correct.length, "\n",
               "multic.q key 345\n\n", sep=""))
  }


  coefficients <- ifelse(as.integer(coefficients) >= 10,
                         coefficients /
                         ( 10 ^ (nchar(as.integer(coefficients)) - 1) ),
                         coefficients)

  trait.names <- trait.names[1:trait.count]
  covariate.names <- covariate.names[1:covariate.count]
  
  coefficients <- matrix(coefficients, nrow = covariate.count,
                         dimnames = list(covariate.names, trait.names))
  
  return (coefficients)
}

calculate.descriptives <- function(trait, covariate, trait.names,
                                   covariate.names) {
  ## Generate the trait and covariate information for the print/summary method
  ## First generate the mean, standard deviation, minimum, and maximum for
  ## all traits and covariates
  trait[trait == missing.value()] <- NA
  ## trait.n means the number of individuals with non-missing values for that
  ## trait
  trait.n <- apply(trait, 2,
                   function(trait) { return (sum(!is.na(trait))) } )
  trait.means <- apply(trait, 2, mean, na.rm = TRUE)
  if(using.R()) {
    trait.var <- apply(trait, 2, var, use = "complete.obs")
  } else {
    trait.var <- apply(trait, 2, var, na.method="omit")
  }
  trait.std.dev <- sqrt(trait.var)
  trait.min <- apply(trait, 2, min, na.rm = TRUE)
  trait.max <- apply(trait, 2, max, na.rm = TRUE)
  trait.kurtosis <- apply(trait, 2, multic.kurtosis, na.rm = TRUE)
  trait.skewness <- apply(trait, 2, multic.skewness, na.rm = TRUE)
  trait.summary <- cbind(trait.n, trait.means, trait.std.dev, trait.min,
                         trait.max, trait.kurtosis, trait.skewness)
  ## Build data.frame from the calculated elements
  trait.cov.values <- trait.summary

  if( !is.null(covariate) ) {
    covariate[covariate == missing.value()] <- NA
    ## cov.n means the number of individuals with non-missing values for that
    ## covariate
    cov.n <- apply(covariate, 2,
                   function(covariate) { return (sum(!is.na(covariate))) } )
    cov.means <- apply(covariate, 2, mean, na.rm = TRUE)
    if(using.R()) {
      cov.var <- apply(covariate, 2, var, use = "complete.obs")
    } else {
      cov.var <- apply(covariate, 2, var, na.method="omit")
    }
    cov.std.dev <- sqrt(cov.var)
    cov.min <- apply(covariate, 2, min, na.rm = TRUE)
    cov.max <- apply(covariate, 2, max, na.rm = TRUE)
    cov.kurtosis <- apply(covariate, 2, multic.kurtosis, na.rm = TRUE)
    cov.skewness <- apply(covariate, 2, multic.skewness, na.rm = TRUE)
    covariate.summary <- cbind(cov.n, cov.means, cov.std.dev, cov.min,
                               cov.max, cov.kurtosis, cov.skewness)
    ## Re-Build data.frame from the calculated elements
    trait.cov.values <- rbind(trait.summary, covariate.summary)
  }
  
  trait.cov.matrix <- matrix(trait.cov.values, ncol=7,
                             dimnames=list(
                               c(trait.names, covariate.names),
                               c("n", "Mean", "Std Dev", "Minimum",
                                 "Maximum", "Kurtosis", "Skewness")))
  
  trait.cov.summary <- data.frame(trait.cov.matrix)

  return (trait.cov.summary)
}

## dadid, momid, and sex can be used for other (already written)
## validation functions, or they can be removed
validate.share.out <- function(share.out, famid, id, dadid, momid, sex,
                               share.out.orig) {
  ## share.out.orig is only here to provide the user with the illusion that
  ## the share.out they specified is being validated, when in fact, a copy
  ## named localshare.out is actually being validated.
  cat(paste("Validating ", share.out.orig, "...\n", sep = ""))

  ## Calculate the number of unique family combinations within a family.
  ## In other words, based on the family sizes, how many lines should
  ## share.out have.
  family.sizes <- get.family.sizes(famid, id)
  unique.famid <- get.unique.families(famid, id)
  
  intra.family.combinations <- family.sizes * (family.sizes - 1) / 2
  total.intra.family.combinations <- sum(intra.family.combinations)
  
  ## Calculate how many lines the file share.out has. 
  wc.share.out <- multic.system(paste('wc', share.out))
  ##  lines.of.share.out <- as.integer(substring(wc.share.out, 1, 8))
  lines.of.share.out <-
    as.integer(multic.strsplit(wc.share.out, sep = " ")[1])
  
  ## If these two values are not equal, there is an error.
  if( total.intra.family.combinations != lines.of.share.out ) {
    stop(paste("\nThe number of lines of ", share.out, " (",
               lines.of.share.out, ")\n",
               "must equal the sum total of all the unique family member\n",
               "combinations (", total.intra.family.combinations, ") ",
               "calculated by the sum over all family\nsizes ni, ",
               "ni * (ni-1) / 2.\nmultic.q key 1406", sep = ""))
  }

  ## Calculate how many words the file share.out has.
  ## words.of.share.out <- as.integer(substring(wc.share.out, 9, 16))
  words.of.share.out <-
    as.integer(multic.strsplit(wc.share.out, sep = " ")[2])
  
  ## share.out should have 6 times as many words as it has lines, because
  ## every line should have 6 words.
  if(lines.of.share.out * 6 != words.of.share.out) {
    stop(paste("The number of words of ", share.out, " (",
               words.of.share.out, ") must be 6 times the\n",
               "number of lines of ", share.out, " (", lines.of.share.out,
               ")\n"))
  }

  unique.ids <- paste(famid, id, sep = '-')

  ## Futher test the validity of the data and share.out by comparing the ids
  ## from famid and id to the order listed in share.out
  passed.validation <- .C("validateShareOut",
                          family.sizes = as.integer(family.sizes),
                          family.count = as.integer(length(family.sizes)), 
                          id = as.character(unique.ids),
                          id.length = as.integer(length(unique.ids)),
                          share.out = as.character(share.out),
                          passed.validation = as.integer(integer(1)),
                          PACKAGE = "multic"
                          )$passed.validation

  ## I chose to return a 1/0 (true/false) value from 'validateShareOut' so
  ## that Splus did not crash when it failed validation, but rather just
  ## stopped the function and re-prompt the user.
  if( !as.logical(passed.validation) ) {
    stop(paste("\nThe validation of", share.out, "failed.\n",
               "multic.q key 1293"))
  }
}

write.initial.values <- function(trait.names, covariate.names,
                                 initial.values,
                                 trait.count, random.effects.count,
                                 coefficients) {
  sink("null.initial.values")
  cat("trait.names:    ", trait.names, "\n")
  cat("covariate.names:", covariate.names, "\n\n")
  cat("*** polygene model initial values ***\n")
  cat("mu:             ", initial.values[1:trait.count], "\n")
  cat("polygene:       ",
      initial.values[(trait.count + 1):(trait.count +
                                        random.effects.count)], "\n")
  cat("major gene:     ",
      initial.values[(trait.count
                      + random.effects.count + 1):
                     (trait.count
                      + 2 * random.effects.count)], "\n")
  ##cat("major gene 2:   ",
  ##    initial.values[(trait.count
  ##                    + 2 * random.effects.count + 1):
  ##                   (trait.count
  ##                    + 3 * random.effects.count)], "\n")
  cat("environmental:  ",
      initial.values[(trait.count
                      + 3 * random.effects.count + 1):
                     (trait.count
                      + 4 * random.effects.count)], "\n")
  cat("sibling:        ",
      initial.values[(trait.count
                      + 4 * random.effects.count + 1):
                     (trait.count
                      + 5 * random.effects.count)], "\n")
  cat("parent-parent:  ",
      initial.values[(trait.count
                      + 5 * random.effects.count + 1):
                     (trait.count
                      + 6 * random.effects.count)], "\n")
  cat("par-offspring:  ",
      initial.values[(trait.count
                      + 6 * random.effects.count + 1):
                     (trait.count
                      + 7 * random.effects.count)], "\n")
  cat("betas:          ", coefficients, "\n")
  sink()
}

write.alt.initial.values <- function(null.initial.values, null.coefficients,
                                     alt.initial.values, alt.coefficients,
                                     trait.count, random.effects.count) {
  sink("alt.initial.values")
  cat("\n*** polygene model final values ***\n")
  cat("mu:             ",
      null.initial.values[1:trait.count], "\n")
  cat("polygene:       ",
      null.initial.values[(trait.count + 1):
                          (trait.count
                           + random.effects.count)], "\n")
  cat("major gene:     ",
      null.initial.values[(trait.count
                           + random.effects.count + 1):
                          (trait.count
                           + 2 * random.effects.count)], "\n")
  ##cat("major gene 2:   ",
  ##    null.initial.values[(trait.count
  ##                         + 2 * random.effects.count + 1):
  ##                        (trait.count
  ##                         + 3 * random.effects.count)], "\n")
  cat("environmental:  ",
      null.initial.values[(trait.count
                           + 3 * random.effects.count + 1):
                          (trait.count
                           + 4 * random.effects.count)], "\n")
  cat("sibling:        ",
      null.initial.values[(trait.count
                           + 4 * random.effects.count + 1):
                          (trait.count
                           + 5 * random.effects.count)], "\n")
  cat("parent-parent:  ",
      null.initial.values[(trait.count
                           + 5 * random.effects.count + 1):
                          (trait.count
                           + 6 * random.effects.count)], "\n")
  cat("par-offspring:  ",
      null.initial.values[(trait.count
                           + 6 * random.effects.count + 1):
                          (trait.count
                           + 7 * random.effects.count)], "\n")
  cat("betas:          ", null.coefficients, "\n")
  
  cat("\n*** major gene model initial values ***\n")
  cat("mu:             ",
      alt.initial.values[1:trait.count], "\n")
  cat("polygene:       ",
      alt.initial.values[(trait.count + 1):
                         (trait.count
                          + random.effects.count)], "\n")
  cat("major gene:     ",
      alt.initial.values[(trait.count
                          + random.effects.count + 1):
                         (trait.count
                          + 2 * random.effects.count)], "\n")
  ##cat("major gene 2:   ",
  ##    alt.initial.values[(trait.count
  ##                        + 2 * random.effects.count + 1):
  ##                       (trait.count
  ##                        + 3 * random.effects.count)], "\n")
  cat("environmental:  ",
      alt.initial.values[(trait.count
                          + 3 * random.effects.count + 1):
                         (trait.count
                          + 4 * random.effects.count)], "\n")
  cat("sibling:        ",
      alt.initial.values[(trait.count
                          + 4 * random.effects.count + 1):
                         (trait.count
                          + 5 * random.effects.count)], "\n")
  cat("parent-parent:  ",
      alt.initial.values[(trait.count
                          + 5 * random.effects.count + 1):
                         (trait.count
                          + 6 * random.effects.count)], "\n")
  cat("par-offspring:  ",
      alt.initial.values[(trait.count
                          + 6 * random.effects.count + 1):
                         (trait.count
                          + 7 * random.effects.count)], "\n")
  cat("betas:          ", alt.coefficients, "\n")
  sink()
}

combine.null.and.alt.initial.values <- function() {
  initial.values.text <- multic.system('cat null.initial.values')
  initial.values.text <- paste(initial.values.text, collapse = '\n')
  remove.file('null.initial.values')
  
  alt.initial.values.text <- multic.system('cat alt.initial.values')
  alt.initial.values.text <-
    paste(alt.initial.values.text, collapse = '\n')
  remove.file('alt.initial.values')
  
  all.initial.values.text <-
    paste(initial.values.text, alt.initial.values.text, sep = "\n")
  sink("initial.values")
  cat(all.initial.values.text)
  cat("\n")
  sink()
}

lappend <- function(list.obj, name, object) {
  if(!is.null(object)) {
    new.index <- length(list.obj) + 1
    list.obj[[new.index]] <- object
    names(list.obj)[new.index] <- name
    return(list.obj)
  }
  return (list.obj)
}

mloci.split <- function(mloci.file.name) {
  loci.names <- .Call("splitTempMloci", as.character(mloci.file.name),
                      PACKAGE = "multic")
  return (loci.names)
}

max.multic <- function(..., na.rm = TRUE) {
  ## I have no idea why these ..1 and ... work.  But they do.  Eric Lunde
  ## 2005-08-05
  if(!is.multic(..1)) {
    stop("\n", substitute(...), "is not a multic object\n")
  }

  multic.object <- ..1

  max.lod <- max(as.numeric(multic.object$log.liks$lod.score), na.rm = TRUE)
  index.of.max.lod <- match(max.lod, multic.object$log.liks$lod.score)
  formula <- paste(as.character(multic.object$call$formula)[-1],
                   collapse = " ~ ")
  
  ## formula, chrom #, Position, LOD, p-value,
  output <- data.frame(formula = formula,
                       marker = dimnames(multic.object$log.liks)[[1]][index.of.max.lod],
                       chrom = "chrom",
                       position = "position",
                       max.lod = max.lod,
                       p.value = "p.value")
  return(output)
}

calculate.correlations <- function(effect, trait.count) {
  cors <- .Call("calculateCorrelations", effect, trait.count,
                PACKAGE = "multic")
  return (cors)
}

get.share.order <- function(share.out, length.id, family.sizes,
                            share.out.orig, using.mloci.out = FALSE) {

  expected.line.count <- sum(family.sizes * (family.sizes - 1) / 2)

  actual.line.count <-
    as.integer(multic.system(paste("cat", share.out, "| wc -l")))
  if(using.mloci.out) {
    ## pound.count is needed the actual line count is meant to be the
    ## number of lines in share.out.  We need to divide the number of lines
    ## in mloci by the number of pounds to get a usable value.
    pound.count <-
      as.integer(multic.system(paste("grep '#'", share.out, "| wc -l")))
    ## The "-1" removes the "# " line from consideration.
    actual.line.count <- actual.line.count/pound.count - 1
  }

  if(expected.line.count != actual.line.count) {
    stop(paste("\nBased on the famid provided (probably through the data ",
               "argument),\n",
               share.out.orig, " is not the expected size.\nIt has ",
               actual.line.count, " lines and is expected to have ",
               expected.line.count, ".\n\n",
               "Probable causes are:\n",
               "- the input is not sorted by famid\n",
               "- the ", share.out.orig, " refers to a different\n",
               "  set of families than the data agrument does\n",
               "- if you are bootstrapping, either the data argument or\n",
               "  ", share.out.orig, " have not been properly expanded.\n",
               "multic.q key 1796", sep = ""))
  }

  share.order <- .Call("getShareOutOrder", share.out, length.id,
                       using.mloci.out,
                       PACKAGE = "multic")
  return (share.order)
}

## Create a function so we can sapply across the loci names
run.alternative.multic <- function(loci.name,
                                   famid,
                                   id,
                                   dadid,
                                   momid,
                                   sex,
                                   traits.with.miss.val,
                                   covariates.with.miss.val,
                                   ascertainment,
                                   epsilon,
                                   boundary.fix,
                                   method.value,
                                   trait.count,
                                   covariate.count,
                                   repeat.count,
                                   trait.covariate.names,
                                   alt.initial.values,
                                   constraints,
                                   max.iterations,
                                   alt.coefficients,
                                   share.out,
                                   calc.residuals,
                                   family.sizes,
                                   unique.families,
                                   family.count) {

  ## Calculate how many lines the file share.out has.  This needs to be one
  ## less than each loci.out
  wc.share.out <- multic.system(paste('wc ', share.out, sep = ""))
  lines.of.share.out <-
    as.integer(multic.strsplit(wc.share.out, sep = " ")[1])
  
  if(!file.exists(loci.name)) {
    cat("Cannot find '", loci.name, "'\n",sep = "")
    return(NULL)
  }
  
  ## Copy the file loci.name to loci.out
  multic.system(paste("cp", loci.name, "loci.out"))
  
  ## Calculate how many lines the file loci.out has.  To be of the right
  ## format, loci.out needs to have one more line than share.out
  wc.loci.out <- multic.system('wc loci.out')
  lines.of.loci.out <- as.integer(multic.strsplit(wc.loci.out, sep = " ")[1])
  
  ## If loci.out is empty, we are done processing tempmloci.out and should
  ## end the program
  if(lines.of.loci.out == 0) {
    break
  } else if( (lines.of.loci.out - 1) != lines.of.share.out) {
    cat(paste("\n", share.out, " and loci.out do not have the same number ",
              "of lines.\n", sep = ""))
    cat(paste(share.out, " has ", lines.of.share.out, " and loci.out has ",
              lines.of.loci.out-1, " (not counting the header)\n",
              sep = ""))
    cat("multic.q key 426\n")
    cat(paste("Terminating multic on", date(), "\n"))
    stop()
  }
  
  ## Get the cM value (for mibd._.*) or the mark number for ibd.* so we can
  ## be selective as to which ibds to run alternative hypothesis on.
  ## Eric Lunde, 03/07/2004
  if(FALSE) {
    if(substring(loci.name, 1, 4) == 'ibd.') {
      ibd.file.number <- .C("getIthToken",
                            ibd.file.name = as.character(loci.name),
                            file.count = as.integer(1),
                            i = as.integer(1),
                            delimiter = as.character('.'),
                            PACKAGE = "multic"
                            )$ibd.file.name
    } else if(substring(loci.name, 1, 5) == 'mibd.') {
      ibd.file.number <-
        as.integer(.C("getIthToken",
                      ibd.file.name = as.character(loci.name),
                      file.count = as.integer(1),
                      i = as.integer(2),
                      delimiter = as.character('.'),
                      PACKAGE = "multic"
                      )$ibd.file.name
                   )
    }else {
      stop(paste("\n'", loci.name, "' is an unrecognized file name.\n",
                 "multic.q key 480\n", sep = ""))
    }
  }
  
  ##alt.multic <-
  .C("multic",
     as.character(famid),
     as.character(id),
     as.character(dadid),
     as.character(momid),
     as.integer(sex),
     as.numeric(as.vector(traits.with.miss.val)),
     as.numeric(as.vector(covariates.with.miss.val)),
     as.integer(ascertainment),
     as.integer(length(id)),
     as.integer(!is.null(ascertainment)),
     as.integer(1),
     as.numeric(epsilon),
     as.integer(boundary.fix),
     as.integer(method.value),
     as.integer(trait.count),
     as.integer(covariate.count),
     as.integer(repeat.count),
     as.character( trait.covariate.names),
     as.numeric(missing.value()),
     as.numeric(alt.initial.values),
     as.character(constraints),
     as.integer(max.iterations),
     as.numeric(alt.coefficients),
     as.character(share.out),
     as.integer(1),
     as.integer(calc.residuals),
     as.integer(family.sizes),
     as.character(unique.families),
     as.integer(family.count),
     PACKAGE = "multic"
     )
  
  invisible()
}

calculate.cor.values <- function(trait, covariate, polygenic, environmental,
                                 alt.hyp.count) {
  cors <- list()
  trait.count <- ncol(trait)

  if(using.R()) {
    cors$pearson <- cor(cbind(trait, covariate), use = "complete.obs")
  } else {
    cors$pearson <- cor(cbind(trait, covariate), na.method = "omit")
  }
  
  ## Use rank to get the ranks and then use cor to get the spearman
  ## correlation
  trait[trait == missing.value()] <- NA
  covariate[covariate == missing.value()] <- NA
  if(using.R()) {
    ranked <- apply(cbind(trait, covariate), 2, rank, na.last = "keep")
    cors$spearman <- cor(ranked, use = "complete.obs")
  } else {
    ranked <- apply(cbind(trait, covariate), 2, rank)
    
    ## If NA's were present, they appear after ranking as the highest value
    ## in that column.  We must replace the highest value in each column with
    ## NA if they were there before ranking.
    has.na <- apply(cbind(trait, covariate), 2, function(x) any(is.na(x)) )
    if(any(has.na)) {  
      maxes <- apply(ranked[, has.na, drop = FALSE], 2, max)
      replacement.matrix <- t(apply(ranked[, has.na, drop = FALSE], 1,
                                    function(x, y) x == y, maxes))
      ranked[replacement.matrix] <- NA
    }
    cors$spearman <- cor(ranked, na.method = "omit")
  }

  if(trait.count > 1) {
    poly.effects <- polygenic[, "Estimate", ]
    env.effects <- environmental[, "Estimate", ]
    if(alt.hyp.count == 0) {
      poly.effects <- matrix(poly.effects, ncol = 1,
                             dimnames = list(names(poly.effects), "null"))
      env.effects <- matrix(env.effects, ncol = 1,
                            dimnames = list(names(env.effects), "null"))
    }
    cors$genetic <- calculate.correlations(poly.effects, trait.count)
    cors$genetic[cors$genetic == missing.value()] <- NA
    cors$environment <- calculate.correlations(env.effects, trait.count)
    cors$environment[cors$environment == missing.value()] <- NA

    ## We need to massage the data we already have to use the
    ## calculate.correlations function again.
    genetic <- 1 / polygenic[, "h^2", ]
    genetic[genetic == Inf] <- NA
    if(alt.hyp.count == 0) {
      genetic <- matrix(genetic, ncol = 1,
                        dimnames = list(names(genetic), "null"))
    }
    genetic[get.covariance.indices(trait.count), ] <- cors$genetic
    
    environment <- 1 / (1 - polygenic[, "h^2", ])
    environment[environment == Inf] <- NA
    if(alt.hyp.count == 0) {
      environment <- matrix(environment, ncol = 1,
                            dimnames = list(names(environment), "null"))
    }
    environment[get.covariance.indices(trait.count), ] <- cors$environment
    
    cors$phenotype <- calculate.correlations(genetic, trait.count) +
      calculate.correlations(environment, trait.count)
    cors$phenotype[cors$phenotype == missing.value()] <- NA
    ## rename the dimnames for the phenotype columns
  }else {
    cors$genetic <- "No genetic correlations were calculated for this multic object."
    cors$environment <- "No environmental correlations were calculated for this multic object."
    cors$phenotype <- "No phenotypic correlations were calculated for this multic object."
  }

  return (cors)
}

get.variances.from.file <- function(trait.count, covariate.count) {
  random.effects.count <- (trait.count - 1) * trait.count / 2
  
  ## Read the lines from summary.log
  lines <- readLines("summary.log", n = -1)
  
  ## Trim all leading and trailing white space
  trimmed.lines <- trim(lines)
  
  ## Determine where the polygenic data is
  var.covar.lines <-
    trimmed.lines[seq(from = 1 + trait.count * (covariate.count + 1) + 1,
                      length = 7 * (random.effects.count
                        + trait.count))]
  
  ## Split the line into tokens
  data <- matrix(multic.strsplit(var.covar.lines, sep = " "),
                 byrow = TRUE, ncol = 3)
  
  ## Get and return the second token
  variance <- sum(as.numeric(data[, 2]))
  
  return (variance)
}

get.major.diagonal.indices <- function(trait.count) {
  switch(trait.count,
         c(1),
         c(1, 3),
         c(1, 4, 6),
         c(1, 5, 8, 10),
         c(1, 6, 10, 13, 15)
         )
}

get.covariance.indices <- function(trait.count) {
  switch(trait.count,
         NA,
         c(2),
         c(2, 3, 5),
         c(2, 3, 4, 6, 7, 9),
         c(2, 3, 4, 5, 7, 8, 9, 11, 12, 14)
         )
}

get.var.covar.labels <- function(trait.count) {
  switch(trait.count,
         c("1"),
         c("1", "12", "2"),
         c("1", "12", "13", "2", "23", "3"),
         c("1", "12", "13", "14", "2", "23", "24", "3", "34", "4"),
         c("1", "12", "13", "14", "15", "2", "23", "24", "25", "3", "34",
           "35", "4", "45", "5")
         )  
}

## This function's purpose has been replaced.  It's not deleted only so I
## don't delete it before I should.  aka make a mistake
parse.trait.names <- function(str, trait.count) {
  if(trait.count == 1) {
    return (str)
  }else {
    trait.names <- str
    
    ## Get rid of the "cbind("
    trait.names <- multic.strsplit(trait.names, sep = "(")
    trait.names <- paste(trait.names[-1], collapse = "(")
    
    ## Get rid of the last ")"
    trait.names <- multic.strsplit(trait.names, sep = ")")
    trait.names <- paste(trait.names[-length(trait.names)], collapse = ")")
    
    ## Separate the names by "," and put back together correctly
    tokens <- multic.strsplit(trait.names, sep = ",")
    
    open.parens <- close.parens <- 0
    trait.names <- character(trait.count)
    traits.found <- 0
    begin.index <- 1
    
    for(i in 1:length(tokens)) {
      open.parens <-
        open.parens + length(multic.strsplit(tokens[i], sep = "(")) - 1
      close.parens <-
        close.parens + length(multic.strsplit(tokens[i], sep = ")")) - 1
      
      if(open.parens == close.parens) {
        trait.names[traits.found + 1] <- paste(tokens[begin.index:i],
                                               collapse = ", ")
        traits.found <- traits.found + 1
        begin.index <- i + 1
        open.parens <- close.parens <- 0
      }
    }
    
    return(trait.names)
  }
}

missing.value <- function() {
  return (-9)
}

check.for.full.data <- function(x) {
  xname <- substitute(x)
  if(any(is.null(x) | x == "NA" | is.na(x))) {
    stop(paste("\n", x.name, " has at least one NULL or NA value.  ",
               "multic cannot continue.",
               "\nmultic.q key 38\n", sep = ""))
  }
  if(any(is.na(x))) {
    stop(paste("\n", x.name, " does not have value.  ",
               "multic cannot continue.\n",
               "multic.q key 70\n", sep = ""))
  } 
}

load.family.log.likelihoods <- function(family.count, hyp.count) {
  if( !file.exists('fam.lik') ) {
    cat(paste("\nThe multic output file 'fam.lik0' is not present.\n",
              "The multic object will not contain data concerning the family ",
              "log likelihoods.\n\n",
              sep=""))
    return (NULL)
  }

  result <- .Call("loadFamilyLogLikelihoods",
                  as.integer(family.count),
                  as.integer(hyp.count),
                  PACKAGE = "multic")
  
  return (result)
}

is.multic <- function(obj) {
  return(class(obj) == "multic")
}

load.iterations <- function(file) {
  lines <- multic.system(paste("cat", file))
  line.count <- length(lines)

  marker.names <- lines[seq(from = 1, by = 2, length = line.count/2)]
  iterations <- as.integer(lines[seq(from = 2, by = 2,
                                     length = line.count/2)])

  iterations <- data.frame(iterations, row.names = marker.names)

  return (iterations)
}

formula.multic <- function(x, ...) {
  if(!is.multic(x)) {
    stop("\n", substitute(x), "is not a multic object\n")
  }
  return (x$metadata$call$formula)
}

multic.control <- function(epsilon = 1e-5,
                           max.iterations = 50,
                           boundary.fix = TRUE,
                           constraints = c("E", "E", "E", "E", "F", "F", "F"),
                           initial.values = NULL,
                           save.output.files = FALSE,
                           method = c("multic", "leastsq", "maxfun", "emvc"),
                           calc.fam.log.liks = FALSE,
                           calc.residuals = FALSE,
                           keep.input = calc.residuals)
{
  if(epsilon <= 0.) {
    stop("\nthe value of epsilon supplied is zero or negative")
  }
  if(max.iterations < 1.) {
    stop("\nthe value of max.iterations supplied is zero or negative")
  }
  boundary.fix <- as.logical(boundary.fix)
  if(length(constraints) != 7) {
    stop('\nthe length of constraints supplied is not 7')
  }else {
    ## Check for characters other than "C", "E", and "F"
  }
  save.output.files <- as.logical(save.output.files)
  
  ## If the max.iterations was missing (user didn't specify), I should set it
  ## based on which method they chose.
  method <- match.arg(method)
  
  calc.fam.log.liks <- as.logical(calc.fam.log.liks)
  calc.residuals <- as.logical(calc.residuals)
  keep.input <- as.logical(keep.input)

  control <- list(epsilon = epsilon,
                  max.iterations = max.iterations,
                  boundary.fix = boundary.fix,
                  constraints = constraints,
                  initial.values = initial.values,
                  save.output.files = save.output.files,
                  method = method,
                  calc.fam.log.liks = calc.fam.log.liks,
                  calc.residuals = calc.residuals,
                  keep.input = keep.input)

  return (control)
}

make.block.id.pair <- function(unique.ids, family.sizes) {
  block.id.pair <- .Call("makeBlockIdPair", unique.ids, family.sizes,
                         PACKAGE = "multic")
  return (block.id.pair)
}

make.kinship.file <- function(famid, id, dadid, momid, share.out,
                              run.alternative.hyps, mloci.out) {
  ## if we are using mloci, we can get the id-pairs for our
  ## share.out/kinship file from mloci (read.table, skip = 1, n =
  ## sum(table(famid) * (table(famid) - 1) / 2) ),  if we are not using
  ## mloci, we'll have to generate them ourselves.
  family.sizes <- get.family.sizes(famid, id)
  uuids <- paste(famid, id, sep = "-")
  ## Perhaps quite = TRUE can be eaten by Splus's ...
  if(run.alternative.hyps) {
    if(using.R()) {
      id.pairs <- scan(mloci.out,
                       what = list(character(0), character(0), numeric(0)),
                       n = sum(family.sizes * (family.sizes - 1) / 2) * 3,
                       skip = 1,
                       quiet = TRUE)
    } else {
      id.pairs <- scan(mloci.out,
                       what = list(character(0), character(0), numeric(0)),
                       n = sum(family.sizes * (family.sizes - 1) / 2) * 3,
                       skip = 1)
    }
    
    id.pairs <- data.frame(id.pairs)
    names(id.pairs) <- list("id1", "id2", "phi")
    ## read.table(mloci.out, skip = 1,
    ##                        n = sum(family.sizes * (family.sizes - 1) / 2))
  } else {
    ## Since table sorts the result by famid, we need to reorder
    ## family.sizes to stay consistant with the famid order
    family.sizes <- get.family.sizes(famid, id)

    id.pairs <- make.block.id.pair(uuids, as.integer(family.sizes))
  }
  
  if(using.R()) {
    if(! require(kinship, quietly = TRUE)) {
      stop("kinship could not be found, multic.q key 338")
    }
  } else {
    library(kinship)
  }
  cat("Using kinship to generate phi values...\n")
  bds <- makekinship(famid, uuids,
                     paste(famid, dadid, sep = "-"),
                     paste(famid, momid, sep = "-"))
  
  ## Since bds may not be in the necessary order (it is sorted by
  ## character strings not integers, I need to reorder the matrix to get
  ## the phi values.
  bds.mat <- as.matrix(bds)
  order <- match(uuids, dimnames(bds.mat)[[1]])
  ordered.mat <- bds.mat[order, order]
  
  id1 <- factor(id.pairs[, 1], levels = uuids)
  id2 <- factor(id.pairs[, 2], levels = uuids)
  blocks <- ordered.mat[cbind(as.integer(id1), as.integer(id2))]
  
  ## Multiply by 2 to get 2 * phi (phi2)
  blocks <- 2 * blocks
  
  ## The I()'s are included because R does not want character values in a
  ## data.frame.  The I()'s allow this.
  ped <- data.frame(I(famid), I(id), fa = I(dadid), I(momid))
  
  ## Calculate the sib-sib, parent-parent, and parent-offspring 0/1 values
  siblings <- is.sibling(id.pairs[, 1], id.pairs[, 2], ped)
  spouses <- is.spouse(id.pairs[, 1], id.pairs[, 2], ped)
  parent.offsprings <-
    is.parent.offspring(id.pairs[, 1], id.pairs[, 2], ped)
  
  kinship <- cbind(as.character(id.pairs[, 1]),
                   as.character(id.pairs[, 2]), blocks,
                   siblings, spouses, parent.offsprings)
  
  ## Keep in mind that share.out's value is "kinship"
  if(using.R()) {
    write.table(kinship, share.out, row.names = FALSE,
                col.names = FALSE, sep = "\t", quote = FALSE)
  } else {
    write.table(kinship, share.out, dimnames.write = FALSE, sep = "\t")
  }
}

## family.count, unique.families, and family.sizes were created to provide
## the correct values when using bootstrapped data and files.  For example,
## a family may be repeated, thus unique(famid) would remove all but the
## first instance of that family.  So a family would only be counted once,
## when it should be counted multiple times.  Also if the same family was
## repeated in the data set, unique(famid) does not know where one logical
## family ends and the next begins.  So, we could no longer use unique to get
## information from the data set.
get.family.count <- function(famids, ids) {
  return(length(get.family.sizes(famids, ids)))
}

get.unique.families <- function(famids, ids) {
  return(names(get.family.sizes(famids, ids)))
}

get.family.sizes <- function(famids, ids) {
  family.sizes <- .Call("getFamilySizes",
                        as.character(famids),
                        as.character(ids),
                        PACKAGE = "multic")
  return (family.sizes)
}

lm.ready.check <- function(traits, trait.names, covariates, covariate.names,
                           covariate.count) {
  ## make sure all columns in traits are not all null
  ## make sure all columns in traits when cbind'ed to covariates and NA rows
  ##  removed are not all null
  if(is.null(covariates)) {
    trait.ok <- TRUE
  } else {
    trait.ok <- apply(traits, 2,
                      function(x, covariates) {
                        data <- cbind(x, covariates)
                        observation.count <- nrow(na.omit(data))
                        return (observation.count > 0)
                      }, covariates[, 1:covariate.count])
  }

  if(any(!trait.ok)) {
    stop(paste("\nThe traits: ",
               paste(trait.names[!trait.ok], collapse = ", "),
               " have 0 observations when compared\nwith the covariates: ",
               paste(covariate.names[1:covariate.count], collapse = ", "), 
               "\nmultic.q key 2608", sep = ""))
  }
}

## Compare two multic objects for equality
all.equal.multic <- function(target, current, ..., raw = FALSE) {
  if(!raw) {
    ## Both dget for R and S-PLUS convert all NA vectors to type logical
    ## instead of their original, numeric.  For the known all-NA attributes
    ## of a multic object, convert them.
    
    ## These go away when the object contails more than the polygenic model
    ## major.gene1[,,1]
    mode(target$major.gene1) <- "numeric"
    ## log.liks$lod.score
    mode(target$log.liks$lod.score) <- "numeric"
    ## log.liks$p.value
    mode(target$log.liks$p.value) <- "numeric"
    
    ## After dget'ing the target object, cors shows up as a "named" object
    ## instead of "list"
    class(target$cors) <- "list"

    ## Ignore metadata$start.time and metadata$finish.time differences
    target$metadata$start.time <- NULL
    target$metadata$finish.time <- NULL

    current$metadata$start.time <- NULL
    current$metadata$finish.time <- NULL
  }

  all.equal.results <- all.equal.default(target, current)
  return(all.equal.results)
}

## Validates that a data argument does not contain -9's instead of NA's
check.data.for.missing.values <- function(data, data.name) {
  is.not.valid <- any(missing.value() == data) && all(!is.na(data))
  error.message <-
    paste("The data set '", data.name, "' contains ",
          missing.value(), " without containing any NA values.\n",
          "This is most likely not a valid data set.\n",
          "multic.q key 2591", sep = "")
  if(is.not.valid) {
    cat(error.message, "\n\n", sep = "")
    warning(error.message)
  }  
  invisible()
}

## Given a file name (relative to the current directory, return the path
## to the file from the root.
canonical.path.name <- function(file.name) {
  pwd <- multic.system("pwd")
  full.path.name <- paste(pwd, file.name, sep = "/")
  split.file.name <- multic.strsplit(full.path.name, sep = "/")
  while(sum(".." == split.file.name) > 0) {
    index <- match("..", split.file.name)
    split.file.name <- split.file.name[-c(index-1, index)]
  }
  canonical.path.name <- paste(split.file.name, collapse = "/")
  if("" == canonical.path.name) {
    canonical.path.name <- "/"
  }
  return (canonical.path.name)
}
