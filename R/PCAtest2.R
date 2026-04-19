#' PCAtest: Statistical Significance of PCA
#'
#'   PCAtest uses random permutations to build null distributions for
#'   several statistics of a PCA analysis: Psi (Vieira 2012), Phi (Gleason and
#'   Staelin 1975), the rank-of-roots (ter Braak 1988), the index of the
#'   loadings (Vieira 2012), and the correlations of the PC with the variables
#'   (Jackson 1991). Comparing these distributions with the observed values of
#'   the statistics, the function tests: (1) the hypothesis that there is more
#'   correlational structure among the observed variables than expected by
#'   random chance, (2) the statistical significance of each PC, and (3) the
#'   contribution of each observed variable to each significant PC. The
#'   function also calculates the sampling variance around mean observed
#'   statistics based on bootstrap replicates.
#'
#' @param x A matrix or dataframe with variables in the columns and the observations in the rows.
#' @param nperm Number of random permutations to build null distributions of the statistics.
#' @param nboot Number of bootstrap replicates to build 95%-confidence intervals of the observed statistics.
#' @param alpha Nominal alpha level for statistical tests.
#' @param indload A logical indicating whether to calculate the index loadings of the variables with the significant PCs.
#' @param varcorr A logical indicating whether to calculate the correlations of the variables with the significant PCs.
#' @param counter A logical specifying whether to show the progress of the random sampling (bootstrap and permutations) on the screen.
#' @param plot A logical specifying whether to plot the null distributions, observed statistics, and 95%-confidence intervals of statistics based on random permutation and bootstrap resampling.
#'
#' @details
#' PCAtest uses the function stats::prcomp to run a PCA using
#'  the arguments scale = TRUE and center = TRUE.
#' PCAtest plots four types of graphs in a single page: (1) a histogram showing the null distribution and the observed value of the Psi statistic, (2) a histogram showing the null distribution and the observed value of the Phi statistic, (3) a bar plot of the percentage of explained variance of each PC1, PC2, ..., etc., showing the sampling variance based on bootstrap replicates and random permutations with 95%-confidence intervals, and (4) a bar plot of the index of the loadings of each observed variable for PC1, showing the sampling variance of bootstrap replicates and random permutations with 95%- confidence intervals. If more than one PC is significant, additional plots for the index of the loadings are shown in as many new pages as necessary given the number of significant PCs. If the PCA is not significant, based on the Psi and Phi testing results, only histograms (1) and (2) are shown.
#'
#' @return
#' An object of class “list” with the following elements:
#' \describe{
#'   \item{psiobs}{The observed Psi statistic.}
#'
#'   \item{phiobs}{The observed Phi statistic.}
#'
#'   \item{psi}{The null distribution of Psi values.}
#'
#'   \item{phi}{The null distribution of Phi values.}
#'
#'   \item{pervarobs}{The percentage of variance explained by each PC based on the observed data.}
#'
#'   \item{pervarboot}{The percentage of variance explained by each PC based on the bootstrapped data.}
#'
#'   \item{pervarbootci}{Confidence intervals of the percentage of variance explained by each PC based on the bootstrapped data.}
#'
#'   \item{pervarperm}{The percentage of variance explained by each PC based on the randomized data.}
#'
#'   \item{pervarpermci}{Confidence intervals of the percentage of variance explained by each PC based on the randomized data.}
#'
#'   \item{indexloadobs}{The index of the loadings of the observed data.}
#'
#'   \item{indexloadboot}{The index of the loadings of the bootstrapped data.}
#'
#'   \item{indexloadbootci}{Confidence intervals of the index of the loadings based on the bootstrapped data.}
#'
#'   \item{indexloadperm}{The index of the loadings based on the randomized data.}
#'
#'   \item{indexloadpermci}{Confidence intervals of the index of the loadings based on the randomized data.}
#'
#'   \item{corobs}{If varcorr=TRUE, the correlations of the observed variables with each significant PC.}
#'
#'   \item{corboot}{If varcorr=TRUE, the correlations of the observed variables with each significant PC based on the bootstrapped data.}
#'
#'   \item{corbootci}{If varcorr=TRUE, the confidence intervals of the correlations of the variables with each significant PC based on the bootstrapped data.}
#'
#'   \item{corperm}{If varcorr=TRUE, the correlations of the observed variables with each significant PC based on randomized data.}
#'
#'   \item{corpermci}{If varcorr=TRUE, the confidence intervals of the correlations of the variables with each significant PC based on randomized data.}
#'
#'}
#'
#' @author
#' Arley Camargo
#'
#' @references
#' \itemize{
#'   \item Gleason, T. C. and Staelin R. (1975) A proposal for handling missing data. Psychometrika, 40, 229–252.
#'   \item Jackson, J. E. (1991) A User’s Guide to Principal Components. John Wiley & Sons, New York, USA.
#'   \item Ringnér, M. (2008) What is principal component analysis? Nature Biotechnology, 26, 303–304.
#'   \item ter Braak, C. F. J. (1990) Update notes: CANOCO (version 3.1). Agricultural Mattematic Group, Report LWA-88-02, Wagningen, Netherlands.
#'   \item Vieira, V. M. N. C. S. (2012) Permutation tests to estimate significances on Principal Components Analysis. Computational Ecology and Software, 2, 103–123.
#'   \item Wong, M. K. L. and Carmona, C. P. (2021) Including intraspecific trait variability to avoid distortion of functional diversity and ecological inference: Lessons from natural assemblages. Methods in Ecology and Evolution. \url{https://doi.org/10.1111/2041- 210X.13568}.
#'}
#' @export
#'
#' @examples
#'#PCA analysis of five uncorrelated (r=0) variables
#'library(MASS)
#'mu <- rep(0,5)
#'Sigma <- matrix(c(rep(c(1,0,0,0,0,0),4),1),5)
#'ex0 <- mvrnorm(100, mu = mu, Sigma = Sigma )
#'result<-PCAtest(ex0, 100, 100, 0.05, varcorr=FALSE, counter=FALSE, plot=TRUE)
#'
#'#PCA analysis of five correlated (r=0.25) variables
#'Sigma <- matrix(c(rep(c(1,0.25,0.25,0.25,0.25,0.25),4),1),5)
#'ex025 <- mvrnorm(100, mu = mu, Sigma = Sigma )
#'result<-PCAtest(ex025, 100, 100, 0.05, varcorr=FALSE, counter=FALSE, plot=TRUE)
#'
#'#PCA analysis of five correlated (r=0.5) variables
#'Sigma <- matrix(c(rep(c(1,0.5,0.5,0.5,0.5,0.5),4),1),5)
#'ex05 <- mvrnorm(100, mu = mu, Sigma = Sigma )
#'result<-PCAtest(ex05, 100, 100, 0.05, varcorr=FALSE, counter=FALSE, plot=TRUE)
#'
#'#PCA analysis of seven morphological variables from 29 ant species (from
#'#Wong and Carmona 2021, https://doi.org/10.1111/2041-210X.13568)
#'data("ants")
#'result<-PCAtest(ants, 100, 100, 0.05, varcorr=FALSE, counter=FALSE, plot=TRUE)

PCAtest2 <- function(
  x,
  nperm = 1000,
  nboot = 1000,
  alpha = 0.05,
  indload = TRUE,
  varcorr = FALSE,
  counter = TRUE,
  plot = TRUE,
  scale = TRUE,
  center = TRUE
) {

  # empirical eigenvalues, loadings, Psi, and Phi

  pcaemp <- stats::prcomp(x, scale = scale, center = center)
  eigenvalues <- pcaemp$sdev^2

  if (nrow(x) < ncol(x)) {
    eigenvalues <- eigenvalues[-length(eigenvalues)]
  }

  eigenobs <- eigenvalues
  pervarobs <- eigenvalues / sum(eigenvalues) * 100

  # empirical index loadings

  indexloadobs <- t(pcaemp$rotation^2 %*% diag(eigenvalues^2))

  # empirical correlations

  if (isTRUE(varcorr)) {
    corobs <- t(pcaemp$rotation %*% diag(sqrt(eigenvalues)))
  }

  # empirical Psi
  Psiobs <- Psi <- sum((eigenvalues - 1)^2)

  # empirical Phi
  Phi <- sum(eigenvalues^2)
  Phi <- sqrt((Phi - ncol(x)) / (ncol(x) * (ncol(x) - 1)))
  Phiobs <- Phi

  # bootstraped data for confidence intervals of empirical eigenvalues, index loadings,
  # and correlations

  pervarboot <- matrix(nrow = nboot, ncol = length(eigenvalues))
  indexloadboot <- array(dim = c(nboot, dim(indexloadobs)))
  corboot <- array(dim = c(nboot, dim(corobs)))

  cat("\nSampling bootstrap replicates... Please wait\n")

  for (i in 1:nboot) {
    if (isTRUE(counter)) {
      cat("\r", i, "of", nboot, "bootstrap replicates\r")
      utils::flush.console()
    }

    bootdata <- x[sample(nrow(x), size = nrow(x), replace = TRUE), ]
    pcaboot <- stats::prcomp(bootdata, scale = scale, center = center)
    eigenvalues <- pcaboot$sdev^2

    if (nrow(x) < ncol(x)) {
      eigenvalues <- eigenvalues[-length(eigenvalues)]
    }

    pervarboot[i, ] <- eigenvalues / sum(eigenvalues) * 100

    if (isTRUE(indload)) {
      indexload <- t(pcaboot$rotation^2 %*% diag(eigenvalues^2))
      indexloadboot[i, , ] <- indexload
    }

    if (isTRUE(varcorr)) {
      corload <- t(pcaemp$rotation %*% diag(sqrt(eigenvalues)))
      corboot[i, , ] <- corload
    }
  }

  cat(
    "\nCalculating confidence intervals of empirical statistics... Please wait\n"
  )

  # confidence intervals of percentage of variation
  confint <- apply(
    pervarboot,
    MARGIN = 2,
    FUN = stats::quantile,
    probs = c(0.025, 0.975)
  )


  if (isTRUE(indload)) {
    # confidence intervals of loadings
    confintindboot <- cbind(as.vector(apply(indexloadboot, c(3, 2), function(x) stats::quantile(x, probs = alpha / 2))),
                            as.vector(apply(indexloadboot, c(3, 2), function(x) stats::quantile(x, probs = 1 - alpha / 2))))
  }
  if (isTRUE(varcorr)) {
    # confidence intervals of correlations
    confintcorboot <- cbind(as.vector(apply(corboot, c(3, 2), function(x) stats::quantile(x, probs = alpha / 2))),
                            as.vector(apply(corboot, c(3, 2), function(x) stats::quantile(x, probs = 1 - alpha / 2))))
  }

  # null distributions based on randomizations of Psi, Phi, eigenvalues, percentage of variation,
  # and index loadings

  Psi <- numeric(nperm)
  Phi <- numeric(nperm)
  eigenrand <- matrix(nrow = nperm, ncol = length(eigenvalues))
  eigenprob <- numeric(length(eigenvalues))
  pervarperm <- matrix(nrow = nperm, ncol = length(eigenvalues))
  indexloadperm <- array(dim = c(nperm, dim(indexloadobs)))
  corperm <- array(dim = c(nperm, dim(corobs)))

  cat("\nSampling random permutations... Please wait\n")

  for (i in 1:nperm) {
    if (isTRUE(counter)) {
      cat("\r", i, "of", nperm, "random permutations                                                 \r")
      utils::flush.console()
    }

    repvalue <- 0
    perm <- apply(x, MARGIN = 2, FUN = sample)
    pcaperm <- stats::prcomp(perm, scale = scale, center = center)
    eigenvalues <- pcaperm$sdev^2 # eigenvalues

    if (nrow(x) < ncol(x)) {
      eigenvalues <- eigenvalues[-length(eigenvalues)]
    }

    pervarperm[i, ] <- eigenvalues / sum(eigenvalues) * 100
    repvalue <- sum((eigenvalues - 1)^2)
    Psi[i] <- repvalue
    Phi[i] <- sqrt((repvalue - ncol(x)) / (ncol(x) * (ncol(x) - 1)))

    eigenrand[i, ] <- eigenvalues

    if (isTRUE(indload)) {
      indexload <- t(pcaboot$rotation^2 %*% diag(eigenvalues^2))
      indexloadperm[i, , ] <- indexload
    }

    if (isTRUE(varcorr)) {
      corload <- t(pcaemp$rotation %*% diag(sqrt(eigenvalues)))
      corperm[i, , ] <- corload
    }
  }

  cat(
    "\nComparing empirical statistics with their null distributions... Please wait\n"
  )

  confintperm <- apply(
    pervarperm,
    MARGIN = 2,
    FUN = stats::quantile,
    probs = c(alpha, 1 - alpha / 2)
  )
  if (isTRUE(indload)) {
    confintind <- cbind(as.vector(apply(indexloadperm, c(3, 2), function(x) stats::quantile(x, probs = alpha / 2))),
                        as.vector(apply(indexloadperm, c(3, 2), function(x) stats::quantile(x, probs = 1 - alpha / 2))))
    meanind <- as.vector(apply(indexloadperm, c(3, 2), mean))
  }

  if (varcorr==T) {
    confintcor <- cbind(as.vector(apply(corperm, c(3, 2), function(x) stats::quantile(x, probs = alpha / 2))),
                        as.vector(apply(corperm, c(3, 2), function(x) stats::quantile(x, probs = 1 - alpha / 2))))
    meancor <- as.vector(apply(corperm, c(3, 2), mean))
  }
  # comparing empirical Psi, Phi and eigenvalues with their null distributions to calculate
  # p-values

  Psiprob <- length(which(Psi > Psiobs)) / nperm
  Phiprob <- length(which(Phi > Phiobs)) / nperm

  for (k in 1:length(eigenvalues)) {
    eigenprob[k] <- length(which(eigenrand[, k] > eigenobs[k])) / nperm
  }

  if (Psiprob < alpha & Phiprob < alpha) {
    # test PC axes if both Psi and Phi are significant

    # find out which PCs are significant

    #	sigaxes <- 0
    #	for (i in 1:length(eigenprob)) {
    #  	if (eigenprob[i] < alpha) {
    #  		sigaxes <- sigaxes + 1
    #  		}
    #	}

    sigaxes <- c()
    for (i in 1:length(eigenprob)) {
      if (eigenprob[i] < alpha) {
        sigaxes[i] <- i
      } else {
        sigaxes[i] <- 0
      }
    }

    # find out which index loadings are significant for each significant axis

    if (isTRUE(indload)) {
      sigload <- list()
      for (j in 1:length(sigaxes)) {
        if (sigaxes[j] != 0) {
          conteo <- c()
          for (i in 1:nperm) {
            conteo <- c(
              conteo,
              which(indexloadperm[[i]][j, ] > indexloadobs[j, ])
            )
          }
          sigload[[j]] <- setdiff(
            c(1:dim(x)[2]),
            names(table(conteo)[table(conteo) / nperm > alpha])
          )
        } else {
          sigload[[j]] <- c("NA")
        }
      }
    }

    # find out which correlations are significant for each significant axis

    if (varcorr == T) {
      sigcor <- list()
      for (j in 1:length(sigaxes)) {
        if (sigaxes[j] != 0) {
          conteo <- c()
          for (i in 1:nperm) {
            conteo <- c(conteo, which(corperm[[i]][j, ] > corobs[j, ]))
          }
          sigcor[[j]] <- setdiff(
            c(1:dim(x)[2]),
            names(table(conteo)[table(conteo) / nperm > alpha])
          )
        } else {
          sigcor[[j]] <- c("NA")
        }
      }
    }
  }

  # screen output

  cat(
    paste(
      "\n",
      "========================================================",
      sep = ""
    ),
    paste(
      "Test of PCA significance: ",
      dim(x)[2],
      " variables, ",
      dim(x)[1],
      " observations",
      sep = ""
    ),
    paste(
      nboot,
      " bootstrap replicates, ",
      nperm,
      " random permutations",
      sep = ""
    ),
    paste("========================================================", sep = ""),
    sep = "\n"
  )

  cat(
    paste(
      "\n",
      "Empirical Psi = ",
      format(round(Psiobs, 4), nsmall = 4),
      ", Max null Psi = ",
      format(round(max(Psi), 4), nsmall = 4),
      ", Min null Psi = ",
      format(round(min(Psi), 4), nsmall = 4),
      ", p-value = ",
      Psiprob,
      sep = ""
    ),
    paste(
      "Empirical Phi = ",
      format(round(Phiobs, 4), nsmall = 4),
      ", Max null Phi = ",
      format(round(max(Phi), 4), nsmall = 4),
      ", Min null Phi = ",
      format(round(min(Phi), 4), nsmall = 4),
      ", p-value = ",
      Phiprob,
      sep = ""
    ),
    sep = "\n"
  )

  if (Psiprob >= alpha & Phiprob >= alpha) {
    # if both Psi and Phi are not significant
    cat(paste("\n", "PCA is not significant!", sep = ""))
  }

  if (Psiprob < alpha & Phiprob < alpha) {
    # test PC axes if both Psi and Phi are significant

    for (i in 1:length(eigenobs)) {
      cat(paste(
        "\n",
        "Empirical eigenvalue #",
        i,
        " = ",
        round(eigenobs[i], 5),
        ", Max null eigenvalue = ",
        round(max(eigenrand[, i]), 5),
        ", p-value = ",
        eigenprob[i],
        sep = ""
      ))
    }

    cat("\n")

    for (i in sigaxes[sigaxes != 0]) {
      cat(paste(
        "\n",
        "PC ",
        i,
        " is significant and accounts for ",
        round(pervarobs[i], digits = 1),
        "% (95%-CI:",
        round(confint[1, i], digits = 1),
        "-",
        round(confint[2, i], digits = 1),
        ")",
        " of the total variation",
        sep = ""
      ))
    }

    if (length(sigaxes[sigaxes != 0]) > 1) {
      cat("\n")
      cat(paste(
        "\n",
        length(sigaxes[sigaxes != 0]),
        " PC axes are significant and account for ",
        round(sum(pervarobs[sigaxes[sigaxes != 0]]), digits = 1),
        "% of the total variation",
        sep = ""
      ))
    }

    if (indload == T) {
      cat("\n")

      for (i in sigaxes[sigaxes != 0]) {
        if (length(sigload[[i]]) == 0) {
          cat(paste(
            "\n",
            "None of the variables have significant loadings on PC ",
            i,
            sep = ""
          ))
        } else {
          if (length(sigload[[i]]) == 1) {
            cat(paste(
              "\n",
              "Variable ",
              sigload[[i]],
              " has a significant loading on PC ",
              i,
              sep = ""
            ))
          } else {
            cat(paste(
              "\n",
              "Variables ",
              paste(sigload[[i]][1:length(sigload[[i]]) - 1], collapse = ", "),
              ", and ",
              paste(sigload[[i]][length(sigload[[i]])]),
              " have significant loadings on PC ",
              i,
              sep = ""
            ))
          }
        }
      }
    }

    if (varcorr == T) {
      cat("\n")

      for (i in sigaxes[sigaxes != 0]) {
        if (length(sigcor[[i]]) == 0) {
          cat(paste(
            "\n",
            "None of the variables have significant correlations with PC ",
            i,
            sep = ""
          ))
        } else {
          if (length(sigcor[[i]]) == 1) {
            cat(paste(
              "\n",
              "Variable ",
              sigcor[[i]],
              " has a significant correlation with PC ",
              i,
              sep = ""
            ))
          } else {
            cat(paste(
              "\n",
              "Variables ",
              paste(sigcor[[i]][1:length(sigcor[[i]]) - 1], collapse = ", "),
              ", and ",
              paste(sigcor[[i]][length(sigcor[[i]])]),
              " have significant correlations with PC ",
              i,
              sep = ""
            ))
          }
        }
      }
    }

    cat("\n\n")
  }

  # plots

  if (plot == T) {
    graphics::par(mfrow = c(2, 2), mar = c(5, 4, 1, 2) + 0.1)

    # plot of empirical and randomized Psi

    h <- graphics::hist(Psi, plot = FALSE)
    h$density = h$counts / sum(h$counts) * 100
    plot(
      h,
      freq = FALSE,
      col = "gray45",
      xlab = "Psi",
      xlim = c(0, max(max(Psi), Psiobs)),
      ylab = "Percentage of permutations",
      main = ""
    )
    graphics::arrows(
      x0 = Psiobs,
      y0 = max(h$density) / 10,
      y1 = 0,
      col = "red",
      lwd = 2,
      length = 0.1
    )
    graphics::legend(
      "topright",
      legend = c("Null distribution", "Empirical value"),
      fill = c("gray45", "red"),
      bty = "n",
      cex = 0.8
    )

    # plot of empirical and randomized Phi

    h <- graphics::hist(Phi, plot = FALSE)
    h$density <- h$counts / sum(h$counts) * 100
    plot(
      h,
      freq = FALSE,
      col = "gray45",
      xlab = "Phi",
      xlim = c(0, max(max(Phi), Phiobs)),
      ylab = "Percentage of permutations",
      main = ""
    )
    graphics::arrows(
      x0 = Phiobs,
      y0 = max(h$density) / 10,
      y1 = 0,
      col = "red",
      lwd = 2,
      length = 0.1
    )

    # plot of bootstrapped and randomized percentage of variation for all PCs

    if (Psiprob < alpha & Phiprob < alpha) {
      # test PC axes if both Psi and Phi are significant

      plot(
        pervarobs,
        ylab = "Percentage of total variation",
        xlab = "PC",
        bty = "n",
        ylim = c(0, max(confint)),
        type = "b",
        pch = 19,
        lty = "dashed",
        col = "red",
        xaxt = "n"
      )
      graphics::axis(1, at = 1:length(eigenobs))
      graphics::lines(
        apply(pervarperm, MARGIN = 2, FUN = mean),
        type = "b",
        pch = 19,
        lty = "dashed",
        col = "gray45"
      )
      suppressWarnings(graphics::arrows(
        x0 = c(1:length(eigenvalues)),
        y0 = confint[2, ],
        y1 = confint[1, ],
        code = 3,
        angle = 90,
        length = 0.05,
        col = "red"
      ))
      suppressWarnings(graphics::arrows(
        x0 = c(1:length(eigenvalues)),
        y0 = confintperm[2, ],
        y1 = confintperm[1, ],
        code = 3,
        angle = 90,
        length = 0.05,
        col = "gray45"
      ))

      # plot of bootstrapped and randomized index loadings for significant PCs only
      if (indload == T) {
        k = 1
        for (i in 1:length(sigaxes[sigaxes != 0])) {
          if (i %in% seq(2, max(2, length(sigaxes[sigaxes != 0])), 4)) {
            # if sigaxes is 2, 6, 10,... make another window to display more plots

            grDevices::dev.new()
            graphics::par(mfrow = c(2, 2), mar = c(5, 4, 1, 2) + 0.1)
            plot(
              indexloadobs[sigaxes[sigaxes != 0][i], ],
              ylab = paste("Index loadings of PC", sigaxes[sigaxes != 0][i]),
              xlab = "Variable",
              bty = "n",
              ylim = c(0, max(confintindboot)),
              pch = 19,
              col = "red",
              xaxt = "n"
            )
            graphics::axis(1, at = 1:dim(x)[2])
            graphics::lines(
              meanind[k:(k - 1 + dim(x)[2])],
              type = "p",
              pch = 19,
              col = "gray45"
            )
            suppressWarnings(graphics::arrows(
              x0 = c(1:dim(x)[2]),
              y0 = confintindboot[k:(k - 1 + dim(x)[2]), 1],
              y1 = confintindboot[k:(k - 1 + dim(x)[2]), 2],
              code = 3,
              angle = 90,
              length = 0.05,
              col = "red"
            ))
            suppressWarnings(graphics::arrows(
              x0 = c(1:dim(x)[2]),
              y0 = confintind[k:(k - 1 + dim(x)[2]), 1],
              y1 = confintind[k:(k - 1 + dim(x)[2]), 2],
              code = 3,
              angle = 90,
              length = 0.05,
              col = "gray45"
            ))
          } else {
            plot(
              indexloadobs[sigaxes[sigaxes != 0][i], ],
              ylab = paste("Index loadings of PC", sigaxes[sigaxes != 0][i]),
              xlab = "Variable",
              bty = "n",
              ylim = c(0, max(confintindboot)),
              pch = 19,
              col = "red",
              xaxt = "n"
            )
            graphics::axis(1, at = 1:dim(x)[2])
            graphics::lines(
              meanind[k:(k - 1 + dim(x)[2])],
              type = "p",
              pch = 19,
              col = "gray45"
            )
            suppressWarnings(graphics::arrows(
              x0 = c(1:dim(x)[2]),
              y0 = confintindboot[k:(k - 1 + dim(x)[2]), 1],
              y1 = confintindboot[k:(k - 1 + dim(x)[2]), 2],
              code = 3,
              angle = 90,
              length = 0.05,
              col = "red"
            ))
            suppressWarnings(graphics::arrows(
              x0 = c(1:dim(x)[2]),
              y0 = confintind[k:(k - 1 + dim(x)[2]), 1],
              y1 = confintind[k:(k - 1 + dim(x)[2]), 2],
              code = 3,
              angle = 90,
              length = 0.05,
              col = "gray45"
            ))
          }
          k = k + dim(x)[2]
        }
      }
    }
  }

  results <- list()

  results[["Empirical Psi"]] <- Psiobs
  results[["Empirical Phi"]] <- Phiobs
  results[["Null Psi"]] <- Psi
  results[["Null Phi"]] <- Phi

  if (Psiprob < alpha & Phiprob < alpha) {
    # test PC axes if both Psi and Phi are significant

    results[["Percentage of variation of empirical PCs"]] <- pervarobs
    results[["Percentage of variation of bootstrapped data"]] <- pervarboot
    results[[
      "Observed confidence intervals of percentage of variation"
    ]] <- confint
    results[["Percentage of variation of randomized data"]] <- pervarperm
    results[[
      "Randomized confidence intervals of percentage of variation"
    ]] <- confintperm

    if (indload == T) {
      results[["Index loadings of empirical PCs"]] <- indexloadobs
      results[["Index loadings with bootstrapped data"]] <- indexloadboot
      results[[
        "Observed confidence intervals of index loadings"
      ]] <- confintindboot
      results[["Index loadings with randomized data"]] <- indexloadperm
      results[[
        "Randomized confidence intervals of index loadings"
      ]] <- confintind
    }

    if (varcorr == T) {
      results[["Correlations of empirical PCs with variables"]] <- corobs
      results[["Correlations in bootstrapped data"]] <- corboot
      results[[
        "Observed confidence intervals of variable correlations"
      ]] <- confintcorboot
      results[["Correlations in randomized data"]] <- corperm
      results[[
        "randomized confidence intervals of variable correlations"
      ]] <- confintcor
    }
  }

  return(results)
}
