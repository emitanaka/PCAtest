
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
    pcaemp$rotation <- pcaemp$rotation[, -ncol(pcaemp$rotation)]
  }

  eigenobs <- eigenvalues
  pervarobs <- eigenvalues / sum(eigenvalues) * 100

  # empirical index loadings

  indexloadobs <- t(pcaemp$rotation^2 %*% diag(eigenvalues^2))

  # empirical correlations

  #if (isTRUE(varcorr)) {
    corobs <- t(pcaemp$rotation %*% diag(sqrt(eigenvalues)))
  #}

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
    discard <- apply(bootdata, MARGIN = 2, function(x) sd(x) == 0)
    while(any(discard)) {
      bootdata <- x[sample(nrow(x), size = nrow(x), replace = TRUE), ]
      discard <- apply(bootdata, MARGIN = 2, function(x) sd(x) == 0)
    }      
    pcaboot <- stats::prcomp(bootdata, scale = scale, center = center)
    eigenvalues <- pcaboot$sdev^2

    if (nrow(x) < ncol(x)) {
      eigenvalues <- eigenvalues[-length(eigenvalues)]
      pcaboot$rotation <- pcaboot$rotation[, -ncol(pcaboot$rotation)]
    }

    pervarboot[i, ] <- eigenvalues / sum(eigenvalues) * 100

    if (isTRUE(indload)) {
      indexload <- t(pcaboot$rotation^2 %*% diag(eigenvalues^2))
      indexloadboot[i, , ] <- indexload
    }

    if (isTRUE(varcorr)) {
      corload <- t(pcaboot$rotation %*% diag(sqrt(eigenvalues)))
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

    perm <- apply(x, MARGIN = 2, FUN = sample)
    pcaperm <- stats::prcomp(perm, scale = scale, center = center)
    eigenvalues <- pcaperm$sdev^2 # eigenvalues

    if (nrow(x) < ncol(x)) {
      eigenvalues <- eigenvalues[-length(eigenvalues)]
      pcaperm$rotation <- pcaperm$rotation[, -ncol(pcaperm$rotation)]
    }

    pervarperm[i, ] <- eigenvalues / sum(eigenvalues) * 100
    
    Psi[i] <- sum((eigenvalues - 1)^2)
    Phi[i] <- sqrt((sum(eigenvalues^2) - ncol(x)) / (ncol(x) * (ncol(x) - 1)))

    eigenrand[i, ] <- eigenvalues

    if (isTRUE(indload)) {
      indexload <- t(pcaperm$rotation^2 %*% diag(eigenvalues^2))
      indexloadperm[i, , ] <- indexload
    }

    if (isTRUE(varcorr)) {
      corload <- t(pcaperm$rotation %*% diag(sqrt(eigenvalues)))
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

  if (isTRUE(varcorr)) {
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

    if (isTRUE(varcorr)) {
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
