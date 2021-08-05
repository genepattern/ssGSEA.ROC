suppressMessages(suppressWarnings(install.packages("getopt", repos = "https://cloud.r-project.org/", 
 quiet = TRUE)))
suppressMessages(suppressWarnings(install.packages("optparse", repos = "https://cloud.r-project.org/", 
 quiet = TRUE)))
suppressMessages(suppressWarnings(library("getopt")))
suppressMessages(suppressWarnings(library("optparse")))

suppressMessages(suppressWarnings(install.packages("ROCR", repos = "https://cloud.r-project.org/", quiet = TRUE)))
suppressMessages(suppressWarnings(library("ROCR")))

arguments = commandArgs(trailingOnly = TRUE)

option_list <- list(
make_option("--ssmatrix", dest = "ssmatrix"),
make_option("--clsfile", dest = "clsfile"),
make_option("--reverse", dest = "reverse"),
make_option("--plotnsets", dest = "plotnsets")
)

opt <- parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE, 
 args = arguments)$options

ssmatrix = opt$ssmatrix
clsfile = opt$clsfile
reverse = as.logical(opt$reverse)
plotnsets = as.numeric(opt$plotnsets)

set.seed(147)

ssauc <- function(ssmatrix, clsfile, reverse, nperm = 1000, permutation.type = 0, 
 fraction = 1, replace = F, adjust.FDR.q.val = F) {
 
 filenames <- basename(ssmatrix)
 
 dataset <- read.table(ssmatrix, header = F, stringsAsFactors = FALSE, sep = "\t", 
  fill = TRUE, quote = "", na = NA)
 cls.cont <- readLines(clsfile)
 
 dataset <- dataset[-c(1), ]
 colnames(dataset) <- dataset[c(1), ]
 dataset <- dataset[-c(1), ]
 rownames(dataset) <- dataset[, c(1)]
 dataset_descriptions <- as.data.frame(dataset[, c(1, 2)], stringsAsFactors = FALSE)
 dataset <- dataset[, -c(1, 2)]
 
 num.lines <- length(cls.cont)
 cls.cont[[3]] <- gsub("\\t", " ", cls.cont[[3]])  #Converts any tabs to spaces
 class.list <- unlist(strsplit(cls.cont[[3]], " "))  #Splits CLS on spaces
 
 cls.cont[[2]] <- gsub("\\t", " ", cls.cont[[2]])  #Converts any tabs to spaces
 class.ids <- unlist(strsplit(cls.cont[[2]], " "))  #Splits CLS on spaces
 class.ids <- class.ids[2:length(class.ids)]
 class.ids <- class.ids[class.ids != ""]
 
 s <- length(class.list)
 t <- table(class.list)[c(unique(class.list))]
 l <- length(t)
 
 sigfigs <- nchar(strsplit(as.character(1/nperm), "\\.")[[1]][2])
 
 phen <- vector(length = l, mode = "character")
 phen.label <- vector(length = l, mode = "numeric")
 class.v <- vector(length = s, mode = "numeric")
 for (i in 1:l) {
  phen[i] <- noquote(names(t)[i])
  phen.label[i] <- i - 1
 }
 for (i in 1:s) {
  for (j in 1:l) {
   if (class.list[i] == phen[j]) {
    class.v[i] <- phen.label[j]
   }
  }
 }
 
 if (reverse == FALSE) {
  class.labels <- class.v
 } else if (reverse == TRUE) {
  class.labels <- 1 - class.v
  phen.label <- rev(phen.label)
  t <- rev(t)
  class.ids <- rev(class.ids)
 }
 
 
 names(t) <- class.ids
 
 col.index <- order(class.labels, decreasing = F)
 class.labels <- class.labels[col.index]
 
 class.labels = 1 - class.labels
 
 dataset <- dataset[, c(col.index)]
 dataset_transposed <- as.data.frame(t(dataset), stringsAsFactors = FALSE)
 dataset_calculations <- as.data.frame(cbind(dataset_descriptions, AUC = as.numeric(""), 
  stringsAsFactors = FALSE))
 dataset_calculations <- as.data.frame(cbind(dataset_calculations, `Matthews Correlation (MCC)` = as.numeric(""), 
  stringsAsFactors = FALSE))
 dataset_calculations <- as.data.frame(cbind(dataset_calculations, `cutoff value (Youden Index)` = as.numeric(""), 
  stringsAsFactors = FALSE))
 dataset_calculations <- as.data.frame(cbind(dataset_calculations, Sensitivity = as.numeric(""), 
  stringsAsFactors = FALSE))
 dataset_calculations <- as.data.frame(cbind(dataset_calculations, Specificity = as.numeric(""), 
  stringsAsFactors = FALSE))
 dataset_calculations <- as.data.frame(cbind(dataset_calculations, PPV = as.numeric(""), 
  stringsAsFactors = FALSE))
 dataset_calculations <- as.data.frame(cbind(dataset_calculations, NPV = as.numeric(""), 
  stringsAsFactors = FALSE))
 
 for (i in 1:length(colnames(dataset_transposed))) {
  pred <- prediction(as.numeric(dataset_transposed[, i]), as.integer(class.labels))
  perf <- performance(pred, "tpr", "fpr")
  auc_perf <- performance(pred, measure = "auc")
  auc_perf <- auc_perf@y.values[[1]]
  dataset_calculations[i, c("AUC")] <- auc_perf
  perf2 <- performance(pred, "mat")
  if(all(is.finite(perf2@y.values[[1]]))) {
  dataset_calculations[i, c("Matthews Correlation (MCC)")] <- as.numeric(perf2@y.values[[1]][which.max(abs(perf2@y.values[[1]]))])
  } else {
  dataset_calculations[i, c("Matthews Correlation (MCC)")] <- paste0(unique(perf2@y.values[[1]]),collapse=";")
  }
  if (dataset_calculations[i, c("Matthews Correlation (MCC)")] >= 0) {
   sensis <- performance(pred, measure = "sens")@y.values[[1]]
   specis <- performance(pred, measure = "spec")@y.values[[1]]
  } else {
   sensis <- -performance(pred, measure = "sens")@y.values[[1]]
   specis <- -performance(pred, measure = "spec")@y.values[[1]]
  }
  summing <- sensis + specis
  youden_ind <- which.max(summing)
  cutoffs <- performance(pred, measure = "sens")@x.values[[1]]
  dataset_calculations[i, c("cutoff value (Youden Index)")] <- cutoffs[youden_ind]
  dataset_calculations[i, c("Sensitivity")] <- performance(pred, measure = "sens")@y.values[[1]][youden_ind]
  dataset_calculations[i, c("Specificity")] <- performance(pred, measure = "spec")@y.values[[1]][youden_ind]
  dataset_calculations[i, c("PPV")] <- performance(pred, measure = "ppv")@y.values[[1]][youden_ind]
  dataset_calculations[i, c("NPV")] <- performance(pred, measure = "npv")@y.values[[1]][youden_ind]
 }
 dataset_order <- order(as.numeric(dataset_calculations[, c("Matthews Correlation (MCC)")]), 
  decreasing = TRUE)
 dataset_calculations <- dataset_calculations[dataset_order, ]
 
 rownames(dataset_calculations) <- 1:nrow(dataset_calculations)
 
 
 
 # Permutations
 
 if ((t[[1]] >= 7) & (t[[2]] >= 7)) {
  
  N <- length(dataset[, 1])
  Ns <- length(dataset[1, ])
  
  subset.mask <- matrix(0, nrow = Ns, ncol = nperm)
  reshuffled.class.labels1 <- matrix(0, nrow = Ns, ncol = nperm)
  reshuffled.class.labels2 <- matrix(0, nrow = Ns, ncol = nperm)
  class.labels1 <- matrix(0, nrow = Ns, ncol = nperm)
  class.labels2 <- matrix(0, nrow = Ns, ncol = nperm)
  
  auc.phi.matrix <- matrix(0, nrow = N, ncol = nperm)
  mcc.phi.matrix <- matrix(0, nrow = N, ncol = nperm)
  
  M1 <- matrix(0, nrow = N, ncol = nperm)
  M2 <- matrix(0, nrow = N, ncol = nperm)
  S1 <- matrix(0, nrow = N, ncol = nperm)
  S2 <- matrix(0, nrow = N, ncol = nperm)
  
  gc()
  
  C <- split(class.labels, class.labels)
  class1.size <- length(C[[1]])
  class2.size <- length(C[[2]])
  class1.index <- seq(1, class1.size, 1)
  class2.index <- seq(class1.size + 1, class1.size + class2.size, 1)
  
  for (r in 1:nperm) {
   class1.subset <- sample(class1.index, size = ceiling(class1.size * fraction), 
    replace = replace)
   class2.subset <- sample(class2.index, size = ceiling(class2.size * fraction), 
    replace = replace)
   class1.subset.size <- length(class1.subset)
   class2.subset.size <- length(class2.subset)
   subset.class1 <- rep(0, class1.size)
   for (i in 1:class1.size) {
    if (is.element(class1.index[i], class1.subset)) {
      subset.class1[i] <- 1
    }
   }
   subset.class2 <- rep(0, class2.size)
   for (i in 1:class2.size) {
    if (is.element(class2.index[i], class2.subset)) {
      subset.class2[i] <- 1
    }
   }
   subset.mask[, r] <- as.numeric(c(subset.class1, subset.class2))
   fraction.class1 <- class1.size/Ns
   fraction.class2 <- class2.size/Ns
   
   if (permutation.type == 0) {
    # random (unbalanced) permutation
    full.subset <- c(class1.subset, class2.subset)
    label1.subset <- sample(full.subset, size = Ns * fraction.class1)
    reshuffled.class.labels1[, r] <- rep(0, Ns)
    reshuffled.class.labels2[, r] <- rep(0, Ns)
    class.labels1[, r] <- rep(0, Ns)
    class.labels2[, r] <- rep(0, Ns)
    for (i in 1:Ns) {
      m1 <- sum(!is.na(match(label1.subset, i)))
      m2 <- sum(!is.na(match(full.subset, i)))
      reshuffled.class.labels1[i, r] <- m1
      reshuffled.class.labels2[i, r] <- m2 - m1
      if (i <= class1.size) {
     class.labels1[i, r] <- m2
     class.labels2[i, r] <- 0
      } else {
     class.labels1[i, r] <- 0
     class.labels2[i, r] <- m2
      }
    }
   } else if (permutation.type == 1) {
    # proportional (balanced) permutation
    
    class1.label1.subset <- sample(class1.subset, size = ceiling(class1.subset.size * 
      fraction.class1))
    class2.label1.subset <- sample(class2.subset, size = floor(class2.subset.size * 
      fraction.class1))
    reshuffled.class.labels1[, r] <- rep(0, Ns)
    reshuffled.class.labels2[, r] <- rep(0, Ns)
    class.labels1[, r] <- rep(0, Ns)
    class.labels2[, r] <- rep(0, Ns)
    for (i in 1:Ns) {
      if (i <= class1.size) {
     m1 <- sum(!is.na(match(class1.label1.subset, i)))
     m2 <- sum(!is.na(match(class1.subset, i)))
     reshuffled.class.labels1[i, r] <- m1
     reshuffled.class.labels2[i, r] <- m2 - m1
     class.labels1[i, r] <- m2
     class.labels2[i, r] <- 0
      } else {
     m1 <- sum(!is.na(match(class2.label1.subset, i)))
     m2 <- sum(!is.na(match(class2.subset, i)))
     reshuffled.class.labels1[i, r] <- m1
     reshuffled.class.labels2[i, r] <- m2 - m1
     class.labels1[i, r] <- 0
     class.labels2[i, r] <- m2
      }
    }
   }
  }
  # message('Computing permutations 1-100 of ', nperm, '...')
  
  for (d in 1:nperm) {
   if ((d == 101) || (d == 201) || (d == 301) || (d == 401) || (d == 501) || 
    (d == 601) || (d == 701) || (d == 801) || (d == 901)) {
    # message('Computing permutations ', d, '-', d + 99, ' of ', nperm, '...')
   }
   
   reshuffled.class.labels1 = (1 - reshuffled.class.labels1)
   
   for (i in 1:length(colnames(dataset_transposed))) {
    pred <- prediction(as.numeric(dataset_transposed[, i]), as.integer(reshuffled.class.labels1[, 
      d]))
    auc_perf <- performance(pred, measure = "auc")
    auc.phi.matrix[i, d] <- as.numeric(auc_perf@y.values[[1]])
    mcc_perf <- performance(pred, "mat")
    mcc.phi.matrix[i, d] <- as.numeric(mcc_perf@y.values[[1]][which.max(abs(mcc_perf@y.values[[1]]))])
    
   }
   auc.phi.matrix <- auc.phi.matrix[dataset_order, ]
   mcc.phi.matrix <- mcc.phi.matrix[dataset_order, ]
   
  }
  
  auc.p.vals <- matrix(0, nrow = N, ncol = 2)
  mcc.p.vals <- matrix(0, nrow = N, ncol = 2)
  
  # message('Computing pValues from AUC null distribution...')
  
  for (i in 1:N) {
   pos.phi <- NULL
   neg.phi <- NULL
   for (j in 1:nperm) {
    if (auc.phi.matrix[i, j] >= 0.5) {
      pos.phi <- c(pos.phi, auc.phi.matrix[i, j])
    } else {
      neg.phi <- c(neg.phi, auc.phi.matrix[i, j])
    }
   }
   ES.value <- dataset_calculations[i, c("AUC")]
   if (ES.value >= 0.5) {
    auc.p.vals[i, 1] <- format(round(signif(sum(pos.phi >= ES.value)/length(pos.phi)), 
      digits = sigfigs), nsmall = sigfigs)
   } else {
    auc.p.vals[i, 1] <- format(round(signif(sum(neg.phi <= ES.value)/length(neg.phi)), 
      digits = sigfigs), nsmall = sigfigs)
   }
  }
  
  # message('Computing pValues from MCC null distribution...')
  
  for (i in 1:N) {
   pos.phi <- NULL
   neg.phi <- NULL
   for (j in 1:nperm) {
    if (mcc.phi.matrix[i, j] >= 0) {
      pos.phi <- c(pos.phi, mcc.phi.matrix[i, j])
    } else {
      neg.phi <- c(neg.phi, mcc.phi.matrix[i, j])
    }
   }
   ES.value <- dataset_calculations[i, c("Matthews Correlation (MCC)")]
   if (ES.value >= 0) {
    mcc.p.vals[i, 1] <- format(round(signif(sum(pos.phi >= ES.value)/length(pos.phi)), 
      digits = sigfigs), nsmall = sigfigs)
   } else {
    mcc.p.vals[i, 1] <- format(round(signif(sum(neg.phi <= ES.value)/length(neg.phi)), 
      digits = sigfigs), nsmall = sigfigs)
   }
  }
  
  
  dataset_calculations <- cbind(dataset_calculations, `AUC NOM pValue` = format(auc.p.vals[, 
   1], nsmall = sigfigs))
  dataset_calculations <- cbind(dataset_calculations, `MCC NOM pValue` = format(mcc.p.vals[, 
   1], nsmall = sigfigs))
  
  # message('Computing FDR q-values from MCC distributions...')
  Ng <- dim(dataset_calculations)[1]
  
  NES <- vector(length = Ng, mode = "numeric")
  phi.norm.mean <- vector(length = Ng, mode = "numeric")
  obs.phi.norm.mean <- vector(length = Ng, mode = "numeric")
  phi.norm.median <- vector(length = Ng, mode = "numeric")
  obs.phi.norm.median <- vector(length = Ng, mode = "numeric")
  phi.norm.mean <- vector(length = Ng, mode = "numeric")
  # obs.phi.mean <- vector(length = Ng, mode = 'numeric')
  FDR.mean <- vector(length = Ng, mode = "numeric")
  FDR.median <- vector(length = Ng, mode = "numeric")
  # phi.norm.median.d <- vector(length = Ng, mode = 'numeric')
  # obs.phi.norm.median.d <- vector(length = Ng, mode = 'numeric')
  
  # Obs.ES.index <- dataset_calculations[, c('AUC')] Orig.index <- seq(1, Ng)
  # Orig.index <- Orig.index[Obs.ES.index] Orig.index <- order(Orig.index,
  # decreasing = F) Obs.ES.norm.sorted <- Obs.ES.norm[Obs.ES.index]
  # Obs.ES.norm.sorted <- dataset_calculations[, c('AUC')] #AUC gs.names.sorted <-
  # gs.names[Obs.ES.index]
  
  for (k in 1:Ng) {
   # NES[k] <- Obs.ES.norm.sorted[k]
   NES[k] <- as.numeric(dataset_calculations[, c("Matthews Correlation (MCC)")][k])
   ES.value <- NES[k]
   count.col <- vector(length = nperm, mode = "numeric")
   obs.count.col <- vector(length = nperm, mode = "numeric")
   for (i in 1:nperm) {
    # phi.vec <- phi.norm[, i]
    phi.vec <- mcc.phi.matrix[, i]
    # obs.phi.vec <- obs.phi.norm[, i]
    obs.phi.vec <- as.numeric(dataset_calculations[, c("Matthews Correlation (MCC)")])
    if (ES.value >= 0) {
      count.col.norm <- sum(phi.vec >= 0)
      obs.count.col.norm <- sum(obs.phi.vec >= 0)
      count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec >= ES.value)/count.col.norm, 
     0)
      obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.phi.vec >= 
     ES.value)/obs.count.col.norm, 0)
    } else {
      count.col.norm <- sum(phi.vec < 0)
      obs.count.col.norm <- sum(obs.phi.vec < 0)
      count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec <= ES.value)/count.col.norm, 
     0)
      obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.phi.vec <= 
     ES.value)/obs.count.col.norm, 0)
    }
   }
   phi.norm.mean[k] <- mean(count.col)
   obs.phi.norm.mean[k] <- mean(obs.count.col)
   phi.norm.median[k] <- median(count.col)
   obs.phi.norm.median[k] <- median(obs.count.col)
   FDR.mean[k] <- ifelse(phi.norm.mean[k]/obs.phi.norm.mean[k] < 1, phi.norm.mean[k]/obs.phi.norm.mean[k], 
    1)
   FDR.median[k] <- ifelse(phi.norm.median[k]/obs.phi.norm.median[k] < 1, 
    phi.norm.median[k]/obs.phi.norm.median[k], 1)
  }
  
  # adjust q-values
  
  if (adjust.FDR.q.val == T) {
   pos.nes <- length(NES[NES >= 0])
   min.FDR.mean <- FDR.mean[pos.nes]
   min.FDR.median <- FDR.median[pos.nes]
   for (k in seq(pos.nes - 1, 1, -1)) {
    if (FDR.mean[k] < min.FDR.mean) {
      min.FDR.mean <- FDR.mean[k]
    }
    if (min.FDR.mean < FDR.mean[k]) {
      FDR.mean[k] <- min.FDR.mean
    }
   }
   
   neg.nes <- pos.nes + 1
   min.FDR.mean <- FDR.mean[neg.nes]
   min.FDR.median <- FDR.median[neg.nes]
   for (k in seq(neg.nes + 1, Ng)) {
    if (FDR.mean[k] < min.FDR.mean) {
      min.FDR.mean <- FDR.mean[k]
    }
    if (min.FDR.mean < FDR.mean[k]) {
      FDR.mean[k] <- min.FDR.mean
    }
   }
  }
  
  FDR.mean <- signif(FDR.mean, digits = 5)
  FDR.median <- signif(FDR.median, digits = 5)
  
  dataset_calculations <- cbind(dataset_calculations, `MCC FDR` = FDR.mean)
  # dataset_calculations <- cbind(dataset_calculations, `MCC FDR (Median)` =
  # FDR.median)
 }
 
 dataset_calculations <- as.data.frame(cbind(dataset_calculations, `ssGSEA Score Wilcox pValue` = "", 
  stringsAsFactors = FALSE))
 for (i in 1:dim(dataset_calculations)[1]) {
  dataset_calculations[i, dim(dataset_calculations)[2]] <- as.numeric(wilcox.test(as.numeric(dataset_transposed[1:eval(t[1]), 
   dataset_calculations[i, 1]]), as.numeric(dataset_transposed[eval(1 + 
   t[1]):eval(t[1] + t[2]), dataset_calculations[i, 1]]))[3])
 }
 
 return(list(filenames = filenames, auc = dataset_calculations, dataset = dataset_transposed, 
  sigfigs = sigfigs, class.labels = class.labels, t = t))
}

result <- ssauc(ssmatrix, clsfile, reverse)

filenames <- result$filenames

dataset_calculations <- result$auc

dataset <- result$dataset
sigfigs <- result$sigfigs
class.labels <- result$class.labels
t <- result$t

names <- c(paste0(names(t)[1], " (N=", t[[1]], ")"), paste0(names(t)[2], " (N=", 
 t[[2]], ")"))
names2 <- paste0(names(t)[1], " (N=", t[[1]], ") vs. ", names(t)[2], " (N=", t[[2]], 
 ")")

# plotnsets=dim(dataset)[2]

if (plotnsets > dim(dataset)[2]) {
 # message('More sets selcted than available in input, plotting results for all
 # sets.')
 plotnsets <- dim(dataset)[2]
}
plotsets <- unique(c(head(dataset_calculations[, c(1)], plotnsets), tail(dataset_calculations[, 
 c(1)], plotnsets)))

direction = paste0(names(t)[1], "vs", names(t)[2])

write.table(dataset_calculations, paste0(filenames, ".", direction, ".Results.tsv"), 
 quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")


pdf(file = paste0(filenames, ".", direction, ".Plots.pdf"), width = 14)
for (i in 1:length(plotsets)) {
 
 if ((as.numeric(dataset_calculations[dataset_calculations[c(1)] == plotsets[i], 
  c("AUC")]) >= 0.5) && (as.numeric(dataset_calculations[dataset_calculations[c(1)] == 
  plotsets[i], c("Matthews Correlation (MCC)")]) >= 0)) {
  direction_color = "red"
  legendloc = "bottomright"
  adjloc = c(-0.05, 1)
 } else if ((as.numeric(dataset_calculations[dataset_calculations[c(1)] == plotsets[i], 
  c("AUC")]) < 0.5) && (as.numeric(dataset_calculations[dataset_calculations[c(1)] == 
  plotsets[i], c("Matthews Correlation (MCC)")]) < 0)) {
  direction_color = "blue"
  legendloc = "topleft"
  adjloc = c(1.05, -1)
 } else {
  direction_color = "darkgray"
  legendloc = "bottomright"
  adjloc = c(-0.05, 1)
 }
 pred <- prediction(as.numeric(eval(parse(text = paste0("dataset$", dataset_calculations[dataset_calculations[c(1)] == 
  plotsets[i], 1])))), as.integer(class.labels))
 perf <- performance(pred, "tpr", "fpr")
 par(mfrow = c(1, 2))
 plot(perf, main = paste0(dataset_calculations[dataset_calculations[c(1)] == plotsets[i], 
  1], "\n", names2, " ROC"), type = "l", col = direction_color, cex.main = 1)
 points(1 - as.numeric(dataset_calculations[dataset_calculations[c(1)] == plotsets[i], 
  c("Specificity")]), as.numeric(dataset_calculations[dataset_calculations[c(1)] == 
  plotsets[i], c("Sensitivity")]), col = "slategrey", pch = 20)
 text(1 - as.numeric(dataset_calculations[dataset_calculations[c(1)] == plotsets[i], 
  c("Specificity")]), as.numeric(dataset_calculations[dataset_calculations[c(1)] == 
  plotsets[i], c("Sensitivity")]), col = "slategrey", labels = paste0("cutoff value (Youden Index): ", 
  round(as.numeric(dataset_calculations[dataset_calculations[c(1)] == plotsets[i], 
   c("cutoff value (Youden Index)")]), 5)), adj = adjloc, cex = 0.8)
 
 legend(legendloc, c(paste0("AUC: ", format(round(as.numeric(dataset_calculations[dataset_calculations[c(1)] == 
  plotsets[i], c("AUC")]), 5), nsmall = 5)), paste0("MCC: ", format(round(as.numeric(dataset_calculations[dataset_calculations[c(1)] == 
  plotsets[i], c("Matthews Correlation (MCC)")]), 5), nsmall = 5))), bty = "n", 
  cex = 0.8)
 abline(a = 0, b = 1)
 boxplot(as.numeric(dataset[1:eval(t[1]), dataset_calculations[dataset_calculations[c(1)] == 
  plotsets[i], 1]]), as.numeric(dataset[eval(1 + t[1]):eval(t[1] + t[2]), dataset_calculations[dataset_calculations[c(1)] == 
  plotsets[i], 1]]), names = names, ylab = "Enrichment Score (ssGSEA)", xlab = paste0("Wilcox pValue: ", 
  format(as.numeric(dataset_calculations[dataset_calculations[c(1)] == plotsets[i], 
   c("ssGSEA Score Wilcox pValue")]), format = "e", digits = sigfigs, flag = "#")))
 title(paste0(dataset_calculations[dataset_calculations[c(1)] == plotsets[i], 
  1], "\nssGSEA Scores"), cex.main = 1)
}
suppressMessages(suppressWarnings(dev.off()))
# message(paste0('Plots saved as: ', filenames, '.', direction, '.Plots.pdf'))
