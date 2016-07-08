#' Compute averaged beta values and classify CpG
#' \code{DAT} performs the Double Average differencial methylation analysis.
#' @param betaValues A matrix containing beta values for each CpG (row) and sample (col) as can be obtained using the minfi::getBeta() package function.
#' @param expDat A dataframe containing experimental information about the samples as can be obtained from the read.450k() function in the minfi package
#' @param window.width Size of the window during the rolling average computation.
#' @param threads Number of threads to use for parallel computing.
#' @return A list containing \enumerate{
#'      \item Names of samples in each condition
#'      \item Classification of CpG (named character vector)
#'      \item Beta values differences between condition (named numeric vector)
#'      \item Summary statistics
#'      }
#' @examples
#' DAT(betas[1:10000,],expSamples,threads=10)
#' @export
#' @importFrom  dplyr left_join
#' @importFrom dplyr arrange
#' @import magrittr
#' @import minfi
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @importFrom caTools runmean
#' @importFrom Hmisc describe
#' @import foreach
#' @import parallel
#' @import doParallel


DAT <- function(betaValues, expDat, window.width=3, threads = 1){
  if(!is.matrix(betaValues)){stop("betaValues is not a matrix")}
  if(!is.numeric(betaValues)){stop("betaValues is not numeric")}
  if(!all(betaValues>=0 & betaValues <=1,na.rm = T)){stop("The range of beta values should be [0,1]")}
  if(ncol(betaValues) != nrow(expDat)){stop(paste("Number of samples doesn't correspond to experiment!",ncol(betaValues),"!=",nrow(expDat)))}
  if(length(table(expSamples$Sample_Group))!=2){stop(" Bad number of groups, only two conditions accepted")}
  if(any(table(expSamples$Sample_Group) < 3)){warning(" The number of samples in at least one group is lower than 3")}

  # Add annotation from manifest to betavalues dataframe
  betaValues <- data.frame(Name=row.names(betaValues),betaValues,check.names = F)
  nSamp <- ncol(betaValues)-1
  cat("Loading manifest data ... \n")
  anno <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)[,1:5])


  combined <- left_join(betaValues,anno, by = "Name") %>% arrange(chr, pos)
  #row.names(combined) <- combined$Name

  sp <- split(combined,f = combined$chr)
  cat("Calculating rolling average per sample ...\n")
  ## Calculate the average beta values per CpG per chromosome
  sp_res <- lapply(sp, FUN = function(x) {

    namesx <- x[,1]
    x <- as.matrix(x[,2:(nSamp+1)],dimnames=list(namesx, colnames(x[2:(nSamp+1)])))

    x <- apply(x,MARGIN = 2,FUN = function(y){
      caTools::runmean(y, k = 3,alg = "C")
    })
    row.names(x) <- namesx
    return(x)

  })
  betas_AVG <- do.call(rbind, sp_res)


  cat("Classifying CpG ...")

  # Need to identify column/sample group -> don't forget to set the levels of the factor !!!!!
  lev1_grp <- which(colnames(betas_AVG) %in% expDat$Sample_Name[expDat$Sample_Group==levels(expDat$Sample_Group)[1]])
  lev2_grp <- which(colnames(betas_AVG) %in% expDat$Sample_Name[expDat$Sample_Group==levels(expDat$Sample_Group)[2]])
  res_all <- list("meta" = list("controls" = colnames(betas_AVG)[lev1_grp], "treated" = colnames(betas_AVG)[lev2_grp]))
  betas_AVG <- betas_AVG[,c(lev1_grp,lev2_grp)]



  if(parallel::detectCores()==1){
    th <- 1
  } else {
      th <- min(parallel::detectCores()-1,threads)
  }

  grpsp <- split(1:nrow(betas_AVG), ceiling(seq_along(1:nrow(betas_AVG))/(nrow(betaValues)/th)))

  cl <- parallel::makeCluster(length(grpsp), outfile="")
  registerDoParallel(cl)
  cat("Using",length(cl),"threads\n")


  resi <- foreach(y = 1:length(grpsp)) %dopar% {

    res <- list("cat"=NULL,"bvalue"=NULL)
    sec.AVG <- apply(betas_AVG[grpsp[[y]],], MARGIN = 1, FUN = function(x){
      avgx <- mean(x)

      if(all(x[1:length(lev1_grp)] < avgx) & all(x[(length(lev1_grp)+1):length(x)] > avgx)){
        "INCREASED"
      } else {
        if(all(x[1:length(lev1_grp)] > avgx) & all(x[(length(lev1_grp)+1):length(x)] < avgx)){
          "DECREASED"
        } else{
          "UNCHANGED"
        }
      }
    })
    res$cat <- sec.AVG

    cat("Calculating beta values differences ...\n")
    delta_betas <- apply(betas_AVG[grpsp[[y]],],MARGIN = 1, function(x){
      # add paired comparison ?
      mean(x[(length(lev1_grp)+1):length(x)])-mean(x[1:length(lev1_grp)])
    })



    cat("Performing statistics ...\n")
    res$bvalue <- delta_betas
    return(res)
  }
  # --------------------------------------------------
  keys <- unique(unlist(lapply(resi, names)))
  resi2 = setNames(do.call(Map, c(c, resi)),keys)
  summary_res <- tapply(resi2$bvalue,INDEX = resi2$cat,FUN = describe)
  summary_res <- as.list(summary_res)
  summary_res <- do.call(rbind,lapply(summary_res,"[[",4))
  class(summary_res) <- "numeric"
  resi2$summary <- summary_res
  cat("Closing co-workers connexions ...\n")
  stopCluster(cl)
  res_all <-c(res_all,resi2)
  return(res_all)
  cat("The End ...\n")
}


