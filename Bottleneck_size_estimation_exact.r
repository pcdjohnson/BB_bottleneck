# R program to estimate bottleneck size for BRSV NGS project
# This code is adapted from:
# https://github.com/weissmanlab/BB_bottleneck/blob/master/Bottleneck_size_estimation_exact.r
# I also vectorised and parallelised the functions parts of the code, so it runs about 20x faster on 8 cores.
# The input data for this project are a metadata xlsx file with data on each sequence, and
# NGS diversity data output by VSensus (https://github.com/rjorton/VSensus).
# Paul Johnson

rm(list = ls())
library(rmutil)
library(parallel)
library(openxlsx)


# minimum coverage
cov.thresh <- 500

# analyse experiment G or the outbreak samples O?
expt <- "G" 

# mismatch proportions below this threshold are set to zero
err.thresh <- 0.01

# drop sites that differ at consensus level, to ask to what extent
# is any relationship between genetic and transmission dependent on
# consensus differences as opposed to sub-consensus diversity
drop.consensus.differences <- FALSE

start.time <- Sys.time()


# define functions

erfinv <- function(x) qnorm((x + 1)/2)/sqrt(2)

generate_log_likelihood_exact <- function(donor_freqs_observed, recipient_total_reads, 
                                          recipient_var_reads_observed, Nb_val, var_calling_threshold, n_variants) {
  likelihood_vector <- rep(0, n_variants)
  log_likelihood_vector <- rep(0, n_variants)
  
  nu_donor <- donor_freqs_observed[, 1]
  variant_reads <- recipient_var_reads_observed[, 1]
  total_reads <- recipient_total_reads[, 1] 
  vr.gt <- (variant_reads >= var_calling_threshold*total_reads)
  # implement variant calling threshold
  for (k in 0:Nb_val){  
    alpha <- k
    beta <- (Nb_val - k)
    if (alpha == 0)
    { alpha <- 0.00001 }
    if (beta == 0)
    { beta <- 0.00001 }
    m <- alpha/(alpha + beta)
    s <- (alpha + beta)
    likelihood_vector[vr.gt] <- likelihood_vector[vr.gt] + 
      (dbetabinom( variant_reads[vr.gt], total_reads[vr.gt], m, s, log = FALSE)*dbinom(k, size=Nb_val, prob= nu_donor[vr.gt]))
  }
  log_likelihood_vector[vr.gt] = log(likelihood_vector[vr.gt])  
  
  # implement variant calling threshold
  likelihood_vector[!vr.gt] = 0
  log_likelihood_vector[!vr.gt] = 0
  for (k in 0:Nb_val){  
    alpha <- k
    beta <- (Nb_val - k)	
    if (alpha == 0)
    { alpha <- 0.00001 }
    if (beta == 0)
    { beta <- 0.00001 }
    m <- alpha/(alpha + beta)
    s <- (alpha + beta)
    likelihood_vector[!vr.gt] <- likelihood_vector[!vr.gt] + 
      (pbetabinom( floor(var_calling_threshold*total_reads[!vr.gt]), total_reads[!vr.gt], m, s)*dbinom(k, size=Nb_val, prob= nu_donor[!vr.gt]))
  }
  log_likelihood_vector[!vr.gt] = log(likelihood_vector[!vr.gt])
  # Now we sum over log likelihoods of the variants at different loci to get the total log likelihood for each value of Nb
  sum(log_likelihood_vector)
}


generate_log_likelihood_function_exact <- 
  function(donor_freqs_observed, recipient_total_reads, recipient_var_reads_observed, Nb_range, var_calling_threshold, confidence_level = 0.95, n_variants) {
    log_likelihood_function <-
      unlist(mclapply(Nb_range, function(Nb_val, ...) {
        print(Nb_val)
        generate_log_likelihood_exact(donor_freqs_observed, recipient_total_reads, recipient_var_reads_observed, Nb_val, var_calling_threshold, n_variants)
      }, mc.cores = detectCores()))
    return(log_likelihood_function)
  }

restrict_log_likelihood <- function(log_likelihood_function, Nb_min, Nb_max) { # restricts log likelihood to the interval of interest 
  for (h in 1:(Nb_min )){  
    if(h< Nb_min)
    {log_likelihood_function[h] = - 999999999}	      # kludge for ensuring that these values less than Nb_min don't interfere with our search for the max of log likelihood in the interval of Nb_min to Nb_max
  }
  
  return(log_likelihood_function)
}


return_bottleneck_size <- function(log_likelihood_function, Nb_range) { 
  return(Nb_range[log_likelihood_function == max(log_likelihood_function)])
}

return_CI <- 
  function(log_likelihood_function, Nb_range, confidence_level = 0.95) { ## returns lower bound of confidence interval 
    max_log_likelihood = return_bottleneck_size(log_likelihood_function, Nb_range)  ## This is the point on the x-axis (bottleneck size) at which log likelihood is maximized
    max_val =  max(log_likelihood_function)  ## This is the maximum value of the log likelihood function, found when the index is our bottleneck estimate
    CI_height = max_val - erfinv(confidence_level)*sqrt(2)  # This value (  height on y axis) determines the confidence intervals using the likelihood ratio test
    Nb_range_lo <- Nb_range[Nb_range <= max_log_likelihood]
    Nb_range_hi <- Nb_range[Nb_range >= max_log_likelihood]
    CI_lo <-
      Nb_range_lo[order(abs(log_likelihood_function[Nb_range <= max_log_likelihood] - CI_height))][1]
    CI_hi <-
      Nb_range_hi[order(abs(log_likelihood_function[Nb_range >= max_log_likelihood] - CI_height))][1]
    c(CI_lo, CI_hi)
  }

############################################################
############################################################

# load data donor & recipient data 

# path to experiment
expt.path <- 
  paste0("~/OneDrive - University of Glasgow/Projects/BRSVtransmission/NGSdata/", expt, "/")

# list of NGS sample names
file.tab <- file.info(list.files(expt.path, full.names = TRUE))
sample.list <- sapply(strsplit(rownames(file.tab)[file.tab$isdir], "//"), "[[", 2)

##### restrict to Paper 1 sequences #####
metadata.all <-
  read.xlsx("~/Dropbox/projects/BRSVtransmission/BRSV_NGS_metadata/BRSV_NGS_project_metadata_2021-12-01.xlsx",
            sheet = 1)
sample.list <- sample.list[sample.list %in% metadata.all$Sequence[metadata.all$Paper1 %in% "In"]]

# restrict to NGS sequences that have produced diversity data (mpile_div.txt file)
mpile.list <- paste0(expt.path, sample.list, "/", sample.list,"_mpile_div.txt")
sample.list <- sample.list[file.exists(mpile.list)]
mpile.list <- mpile.list[file.exists(mpile.list)]

# make table of all pairs of sequences
all.pairs <- expand.grid(donor = sample.list, recip = sample.list, stringsAsFactors = FALSE)
all.pairs <- all.pairs[all.pairs$donor != all.pairs$recip, ]
dim(all.pairs)

# loop over all pairs
res.tab <- sapply(1:nrow(all.pairs), function(r) {
  print(r)
  print(all.pairs[r, ])
  donor <- all.pairs$donor[r]  # donor sequence
  donor.file <- paste0(expt.path, donor, "/", donor,"_mpile_div.txt")
  if(!file.exists(donor.file)) return(NULL)
  donor.data <- read.delim(donor.file)
  dim(donor.data)
  
  recip <- all.pairs$recip[r] 
  recip.file <- paste0(expt.path, recip, "/", recip,"_mpile_div.txt")
  if(!file.exists(recip.file)) return(NULL)
  recip.data <- read.delim(recip.file)
  dim(recip.data)
  
  print(c(nrow(donor.data), nrow(recip.data)))
  if(nrow(donor.data) != nrow(recip.data)) return(c(bottleneck_size = NA, CIlo = NA, CIhi = NA))
  
  # derive consensus sequences
  actg <- c("A", "C", "G", "T")
  donor.data$consensus <- actg[apply(donor.data[, actg], 1, which.max)]
  recip.data$consensus <- actg[apply(recip.data[, actg], 1, which.max)]
  consensus.diff <- donor.data$consensus != recip.data$consensus
  print(paste(sum(consensus.diff), "consensus differences"))
  
  # optionally, ignore sites where there are consensus-level differences
  if(drop.consensus.differences && sum(consensus.diff) > 0) {
    donor.data <- donor.data[consensus.diff, ]
    recip.data <- recip.data[consensus.diff, ]
  }
  dim(donor.data)
  dim(recip.data)
  
  donor_freqs_recip_freqs_and_reads_observed <-
    data.frame(V1 = donor.data$Mismatch, 
               V2 = recip.data$Mismatch, 
               V3 = recip.data$BaseCov, 
               V4 = recip.data$NonRefN,
               V5 = donor.data$BaseCov, 
               V6 = donor.data$NonRefN)
  
  if(mean(donor_freqs_recip_freqs_and_reads_observed$V3) < cov.thresh |
     quantile(donor_freqs_recip_freqs_and_reads_observed$V3, 0.1) < cov.thresh/4 |
     mean(donor_freqs_recip_freqs_and_reads_observed$V5) < cov.thresh |
     quantile(donor_freqs_recip_freqs_and_reads_observed$V5, 0.1) < cov.thresh/4) {
    return(c(bottleneck_size = NA, CIlo = NA, CIhi = NA))
  }
  
  plot_bool  <- TRUE # args$plot_bool
  var_calling_threshold  <- err.thresh # args$var_calling_threshold
  Nb_step <- 1
  Nb_min <- 1 # args$Nb_min
  Nb_max <- 1000 # args$Nb_max
  Nb_range <- seq(Nb_min, Nb_max, by = Nb_step)
  dim(donor_freqs_recip_freqs_and_reads_observed)
  donor_freqs_recip_freqs_and_reads_observed <- 
    donor_freqs_recip_freqs_and_reads_observed[donor_freqs_recip_freqs_and_reads_observed$V1 >= var_calling_threshold, ]
  dim(donor_freqs_recip_freqs_and_reads_observed)
  apply(donor_freqs_recip_freqs_and_reads_observed[, 1:2], 2, function(x) (which(x > var_calling_threshold)))
  donor_freqs_observed <- as.data.frame(donor_freqs_recip_freqs_and_reads_observed[,1])
  n_variants <- nrow(donor_freqs_recip_freqs_and_reads_observed)
  recipient_total_reads <- as.data.frame(donor_freqs_recip_freqs_and_reads_observed[,3]) #read.table(args[2])
  recipient_var_reads_observed <- as.data.frame(donor_freqs_recip_freqs_and_reads_observed[,4])#read.table(args[3])
  
  
  
  log_likelihood_function <- generate_log_likelihood_function_exact(donor_freqs_observed, recipient_total_reads, recipient_var_reads_observed, Nb_range, var_calling_threshold, confidence_level = 0.95, n_variants)
  
  bottleneck_size <- return_bottleneck_size(log_likelihood_function,  Nb_range)
  CI_index <- return_CI(log_likelihood_function,  Nb_range)
  
  ##########################
  ##############################################################################################  ABOVE THIS LINE DETERMINES PEAK LOG LIKELIHOOD AND CONFIDENCE INTERVALS
  # Now we plot the result
  if(plot_bool == TRUE) {
    pdf(file=paste0(donor, "to", recip, "exact_plot.pdf"))
    plot(Nb_range, log_likelihood_function, type = "b", xlab = expression(N[b]), ylab = "log likelihood", cex = 0.6)
    abline(v = bottleneck_size, col="black", lty = 2)  # Draws a verticle line at Nb value for which log likelihood is maximized
    abline(v = CI_index, col="green", lty = 2) # confidence intervals
    mtext(c(bottleneck_size, CI_index), side = 1, at = c(bottleneck_size, CI_index), line = -1)
    title("Log likelihood of bottleneck size, showing ML estimate with 95% CI")
    mtext(paste("Variant calling threshold:", var_calling_threshold), side = 3, line = 0.5)
    dev.off()
  }
  
  print(Sys.time() - start.time)
  print(paste("Bottleneck size", donor, "to", recip))
  print(bottleneck_size)
  print("confidence interval left bound")
  print(CI_index[1])
  print("confidence interval right bound")
  print(CI_index[2])
  
  c(bottleneck_size = bottleneck_size, CIlo = CI_index[1], CIhi = CI_index[2])
})

print(Sys.time() - start.time)


all.pairs.res <- na.omit(cbind(all.pairs, t(res.tab)))
write.csv(all.pairs.res, 
          file = paste0(expt, ifelse(drop.consensus.differences, 
                                     "_noconsensus", ""), 
                        "_bottleneck_size.csv"),
          row.names = FALSE)

# how similar are bottlenecks in opposite directions?
all.pairs.res$unidirectional <- 
  apply(apply(all.pairs.res[, c("donor", "recip")], 1, sort), 2, paste, collapse = "-")
x <- all.pairs.res$bottleneck_size[!duplicated(all.pairs.res$unidirectional)]
names(x) <- all.pairs.res$unidirectional[!duplicated(all.pairs.res$unidirectional)]
y <- all.pairs.res$bottleneck_size[duplicated(all.pairs.res$unidirectional)]
names(y) <- all.pairs.res$unidirectional[duplicated(all.pairs.res$unidirectional)]
plot(jitter(x), jitter(y[names(x)]), log = "xy")
cor(x, y[names(x)], method = "spearman")
abline(0, 1)
