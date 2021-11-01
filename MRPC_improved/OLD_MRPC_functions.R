



#=====================================================================================================
> MRPC
function (data, suffStat, GV, FDR = 0.05, alpha = 0.05, indepTest = c("gaussCItest", 
                                                                      "disCItest", "citest"), labels, p, fixedGaps = NULL, 
          fixedEdges = NULL, NAdelete = TRUE, m.max = Inf, u2pd = c("relaxed", 
                                                                    "rand", "retry"), skel.method = c("stable", 
                                                                                                      "original", "stable.fast"), conservative = FALSE, 
          maj.rule = FALSE, solve.confl = FALSE, FDRcontrol = c("LOND", 
                                                                "ADDIS"), tau = 0.5, lambda = 0.25, verbose = FALSE) 
{
  cl <- match.call()
  if (!missing(p)) 
    stopifnot(is.numeric(p), length(p <- as.integer(p)) == 
                1, p >= 2)
  if (missing(labels)) {
    if (missing(p)) 
      stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  }
  else {
    stopifnot(is.character(labels))
    if (missing(p)) {
      p <- length(labels)
    }
    else if (p != length(labels)) 
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else message("No need to specify 'p', when 'labels' is given")
  }
  u2pd <- match.arg(u2pd)
  skel.method <- match.arg(skel.method)
  if (u2pd != "relaxed") {
    if (conservative || maj.rule) 
      stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")
    if (solve.confl) 
      stop("Versions of PC using lists for the orientation rules (and possibly bi-directed edges)\n can only be run with 'u2pd = relaxed'")
  }
  if (conservative && maj.rule) 
    stop("Choose either conservative PC or majority rule PC!")
  if (verbose) 
    cat("test for independence:", indepTest, "\n")
  skel <- ModiSkeleton(data, suffStat, FDR = FDR, alpha = alpha, 
                       indepTest = indepTest, labels = labels, method = skel.method, 
                       fixedGaps = fixedGaps, fixedEdges = fixedEdges, NAdelete = NAdelete, 
                       m.max = m.max, FDRcontrol = FDRcontrol, tau = tau, lambda = lambda, 
                       verbose = verbose)
  skel@call <- cl
  if (!conservative && !maj.rule) {
    switch(u2pd, relaxed = EdgeOrientation(skel, GV = GV, 
                                           suffStat, FDR, alpha, FDRcontrol = FDRcontrol, indepTest, 
                                           tau = tau, lambda = lambda, verbose = verbose))
  }
  else {
    pc. <- pc.cons.intern(skel, suffStat, match.fun(indepTest), 
                          alpha = alpha, version.unf = c(2, 1), maj.rule = maj.rule, 
                          verbose = verbose)
    EdgeOrientation(pc.$sk, FDRcontrol = FDRcontrol, tau = tau, 
                    lambda = lambda, verbose = verbose)
  }
}


#=====================================================================================================


  > ModiSkeleton
function (data, suffStat, FDR, alpha, indepTest = c("gaussCItest", 
                                                    "disCItest", "citest"), labels, p, method = c("stable", 
                                                                                                  "original", "stable.fast"), m.max = Inf, fixedGaps = NULL, 
          fixedEdges = NULL, NAdelete = TRUE, FDRcontrol = c("LOND", 
                                                             "ADDIS"), tau = 0.5, lambda = 0.25, verbose = FALSE) 
{
  cl <- match.call()
  if (!missing(p)) 
    stopifnot(is.numeric(p), length(p <- as.integer(p)) == 
                1, p >= 2)
  if (missing(labels)) {
    if (missing(p)) 
      stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  }
  else {
    stopifnot(is.character(labels))
    if (missing(p)) 
      p <- length(labels)
    else if (p != length(labels)) 
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
  }
  seq_p <- seq_len(p)
  method <- match.arg(method)
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedGaps), c(p, p))) 
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps))) 
    stop("fixedGaps must be symmetric")
  else G <- !fixedGaps
  diag(G) <- FALSE
  if (any(is.null(fixedEdges))) {
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p))) 
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  else if (!identical(fixedEdges, t(fixedEdges))) 
    stop("fixedEdges must be symmetric")
  {
    sepset <- lapply(seq_p, function(.) vector("list", 
                                               p))
    pMax <- matrix(-Inf, nrow = p, ncol = p)
    m <- 0L
    Alpha <- 0L
    K <- 0
    R <- pval <- kappai <- Ci <- Si <- Ci_plus <- numeric(dim(data)[2]^2)
    gammai <- kappai_star <- alphai <- numeric(dim(data)[2]^2)
    normalizer <- 0.4374901658
    exponent <- 1.6
    gammai_sum <- 0
    diag(pMax) <- 1
    done <- FALSE
    ord <- 0L
    n.edgetests <- numeric(1)
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord1 <- ord + 1L] <- 0
      done <- TRUE
      ind <- which(upper.tri(G), arr.ind = TRUE)
      ind <- ind[order(ind[, 1]), ]
      remEdges <- nrow(ind)
      if (verbose) 
        cat("Order=", ord, "; remaining edges:", 
            remEdges, "\n", sep = "")
      if (method == "stable") {
        G.l <- split(G, gl(p, p))
      }
      for (i in 1:remEdges) {
        if (verbose && (verbose >= 2 || i%%100 == 0)) 
          cat("|i=", i, "|iMax=", remEdges, 
              "\n")
        x <- ind[i, 1]
        y <- ind[i, 2]
        if (G[y, x] && !fixedEdges[y, x]) {
          nbrsBool <- if (method == "stable") 
            nbrsBool1 <- G.l[[x]]
          nbrsBool2 <- G.l[[y]]
          nbrsBool1[y] <- FALSE
          nbrs_x <- seq_p[nbrsBool1]
          nbrsBool2[x] <- FALSE
          nbrs_y <- seq_p[nbrsBool2]
          nbrs <- unique(union(nbrs_x, nbrs_y))
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord) 
              done <- FALSE
            S <- seq_len(ord)
            repeat {
              n.edgetests[ord1] <- n.edgetests[ord1] + 
                1
              m <- m + 1
              if (m >= length(R)) {
                R <- c(R, numeric(dim(data)[2]^2))
                pval <- c(pval, numeric(dim(data)[2]^2))
                kappai <- c(kappai, numeric(dim(data)[2]^2))
                kappai_star <- c(kappai_star, numeric(dim(data)[2]^2))
                Ci <- c(Ci, numeric(dim(data)[2]^2))
                Si <- c(Si, numeric(dim(data)[2]^2))
                Ci_plus <- c(Ci_plus, numeric(dim(data)[2]^2))
                gammai <- c(gammai, numeric(dim(data)[2]^2))
                alphai <- c(alphai, numeric(dim(data)[2]^2))
              }
              if (indepTest == "citest") {
                x <- data[, ind[i, 1]]
                y <- data[, ind[i, 2]]
                z <- data[, nbrs[S]]
                if (length(S) == 0) {
                  P <- ci.test(x, y)
                  pval[m] <- P$p.value
                }
                else {
                  P <- ci.test(x, y, z)
                  pval[m] <- P$p.value
                }
                x <- ind[i, 1]
                y <- ind[i, 2]
              }
              if (indepTest == "gaussCItest") {
                pval[m] <- gaussCItest(x, y, nbrs[S], 
                                       suffStat)
              }
              if (indepTest == "disCItest") {
                pval[m] <- disCItest(x, y, nbrs[S], suffStat)
              }
              if (verbose) 
                cat("x=", x, " y=", y, " S=", 
                    nbrs[S], "\n")
              if (is.na(pval[m])) 
                pval[m] <- as.numeric(NAdelete)
              if (pMax[x, y] < pval[m]) 
                pMax[x, y] <- pval[m]
              if (verbose) 
                cat("Test number =", m, "\n")
              if (verbose) 
                cat("pval =", pval[m], "\n")
              if (FDRcontrol == "LOND") {
                alphai[m] <- SeqFDR(m, FDR, a = 2, R)
                Alpha <- alphai[m]
              }
              else if (FDRcontrol == "ADDIS") {
                if (m == 1) {
                  w0 <- tau * lambda * FDR/2
                  Ci_sum <- 0
                  Si_sum <- 0
                  K <- 0
                  gammai[1] <- normalizer/1^exponent
                  alphai[1] <- w0 * gammai[1]
                  Alpha <- alphai[1]
                  R[1] <- pval[1] <= alphai[1]
                }
                else {
                  run_addis <- addis(alpha = FDR, tau = tau, 
                                     lambda = lambda, iter = m, w0 = w0, 
                                     pval = pval, alphai = alphai, gammai = gammai, 
                                     kappai = kappai, kappai_star = kappai_star, 
                                     K = K, Ci = Ci, Si = Si, Ri = R, 
                                     Ci_plus = Ci_plus, Ci_sum = Ci_sum, 
                                     Si_sum = Si_sum, gammai_sum = gammai_sum, 
                                     normalizer = normalizer, exponent = exponent)
                  alphai[[m]] <- run_addis[[1]]
                  gammai[[m]] <- run_addis[[2]]
                  K <- run_addis[[3]]
                  R[[m]] <- run_addis[[5]]
                  Si[[m - 1]] <- run_addis[[6]]
                  Ci[[m - 1]] <- run_addis[[7]]
                  Ci_sum <- run_addis[[9]]
                  Si_sum <- run_addis[[10]]
                  gammai_sum <- run_addis[[12]]
                  if (K != 0) {
                    kappai[[K]] <- run_addis[[4]]
                    kappai_star[[K]] <- run_addis[[11]]
                    Ci_plus[1:K] <- run_addis[[8]]
                  }
                  Alpha <- alphai[[m]]
                }
              }
              else {
                Alpha <- alpha
              }
              if (verbose) 
                cat("Alpha value =", Alpha, "\n")
              if (pval[m] <= Alpha) {
                R[m] <- 1
                if (verbose) 
                  cat("Since pval<Alpha,test is rejected: Nodes are dependent", 
                      "\n")
              }
              else {
                R[m] <- 0
                if (verbose) 
                  cat("Since pval>Alpha,test is accepted:Nodes are independent", 
                      "\n")
              }
              if (pval[m] >= Alpha) {
                G[x, y] <- G[y, x] <- FALSE
                sepset[[x]][[y]] <- nbrs[S]
                break
              }
              else {
                nextSet <- getNextSet(length_nbrs, ord, 
                                      S)
                if (nextSet$wasLast) 
                  break
                S <- nextSet$nextSet
              }
            }
          }
        }
      }
      ord <- ord + 1L
    }
    for (i in 1:(p - 1)) {
      for (j in 2:p) {
        pMax[i, j] <- pMax[j, i] <- max(pMax[i, j], pMax[j, 
                                                         i])
      }
    }
    }
  Gobject <- if (sum(G) == 0) {
    new("graphNEL", nodes = labels)
  }
  else {
    colnames(G) <- rownames(G) <- labels
    as(G, "graphNEL")
  }
  new("MRPCclass", graph = Gobject, call = cl, n = integer(0), 
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests, 
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1), 
      test = m, alpha = Alpha, R = R, K = K, pval = pval, normalizer = normalizer, 
      exponent = exponent, alphai = alphai, kappai = kappai, 
      kappai_star = kappai_star, Ci = Ci, Si = Si, Ci_plus = Ci_plus, 
      gammai = gammai, gammai_sum = gammai_sum)
}

#=====================================================================================================


  > EdgeOrientation
function (gInput, GV, suffStat, FDR, alpha, indepTest, FDRcontrol, 
          tau = 0.5, lambda = 0.25, verbose = FALSE) 
{
  g <- as(gInput@graph, "matrix")
  g1 <- g
  p <- nrow(g)
  tarmat <- matrix(0, nrow(g), ncol(g))
  rownames(tarmat) <- rownames(g)
  colnames(tarmat) <- colnames(g)
  g[lower.tri(g)] <- 0
  edges <- which(g == 1, arr.ind = TRUE)
  if (GV > 0) {
    edgesWithBothVs <- edges[which(edges[, 1] <= GV & edges[, 
                                                            2] <= GV), ]
    edgesWithVFirst <- edges[which(edges[, 1] <= GV & edges[, 
                                                            2] > GV), ]
    edgesWithVSecond <- edges[which(edges[, 1] > GV & edges[, 
                                                            2] <= GV), ]
    if (length(edgesWithBothVs) > 2) {
      tarmat[edgesWithBothVs] <- 1
      tarmat[edgesWithBothVs[, 2:1]] <- 1
    }
    else {
      tarmat[edgesWithBothVs[1], edgesWithBothVs[2]] <- 1
      tarmat[edgesWithBothVs[2], edgesWithBothVs[1]] <- 1
    }
    if (length(edgesWithVFirst) > 2) {
      tarmat[edgesWithVFirst] <- 1
    }
    else {
      tarmat[edgesWithVFirst[1], edgesWithVFirst[2]] <- 1
    }
    if (length(edgesWithVSecond) > 2) {
      tarmat[edgesWithVSecond[, 2:1]] <- 1
    }
    else {
      tarmat[edgesWithVSecond[2], edgesWithVSecond[1]] <- 1
    }
  }
  if (verbose) 
    cat("\n V-structures are as follows :\n")
  m <- gInput@test
  Alpha <- gInput@alpha
  pval <- gInput@pval
  R <- gInput@R
  ind <- which(g1 == 1, arr.ind = TRUE)
  V <- colnames(g)
  w0 <- tau * lambda * FDR/2
  alphai <- gInput@alphai
  kappai <- gInput@kappai
  kappai_star <- gInput@kappai_star
  K <- gInput@K
  Ci <- gInput@Ci
  Ci_sum <- sum(Ci)
  Si <- gInput@Si
  Si_sum <- sum(Si)
  gammai_sum <- gInput@gammai_sum
  Ci_plus <- gInput@Ci_plus
  gammai <- gInput@gammai
  normalizer <- gInput@normalizer
  exponent <- gInput@exponent
  for (i in seq_len(nrow(ind))) {
    x <- ind[i, 1]
    y <- ind[i, 2]
    allZ <- setdiff(which(g1[y, ] == 1), x)
    for (z in allZ) {
      if ((g1[x, z] == 0 & g1[x, y] == 1 & g1[y, z] == 
           1) & !(tarmat[y, x] == 1) & !(tarmat[z, y] == 
                                         1) & !(tarmat[y, z] == 1) & !(y %in% gInput@sepset[[x]][[z]] || 
                                                                       y %in% gInput@sepset[[z]][[x]])) {
        m <- m + 1
        if (m >= length(R)) {
          R <- c(R, numeric(m/2))
          pval <- c(pval, numeric(m/2))
          kappai <- c(kappai, numeric(m/2))
          kappai_star <- c(kappai_star, numeric(m/2))
          Ci <- c(Ci, numeric(m/2))
          Si <- c(Si, numeric(m/2))
          Ci_plus <- c(Ci_plus, numeric(m/2))
          gammai <- c(gammai, numeric(m/2))
          alphai <- c(alphai, numeric(m/2))
        }
        if (indepTest == "gaussCItest") {
          pval[m] <- gaussCItest(x, z, y, suffStat)
        }
        if (indepTest == "disCItest") {
          pval[m] <- disCItest(x, z, y, suffStat)
        }
        if (FDRcontrol == "LOND") {
          alphai[m] <- SeqFDR(m, FDR, a = 2, R)
          Alpha <- alphai[m]
        }
        else if (FDRcontrol == "ADDIS") {
          run_addis <- addis(alpha = FDR, tau = tau, 
                             lambda = lambda, iter = m, w0 = w0, pval = pval, 
                             alphai = alphai, gammai = gammai, kappai = kappai, 
                             kappai_star = kappai_star, K = K, Ci = Ci, 
                             Si = Si, Ri = R, Ci_plus = Ci_plus, Ci_sum = Ci_sum, 
                             Si_sum = Si_sum, gammai_sum = gammai_sum, 
                             normalizer = normalizer, exponent = exponent)
          alphai[[m]] <- run_addis[[1]]
          gammai[[m]] <- run_addis[[2]]
          K <- run_addis[[3]]
          R[[m]] <- run_addis[[5]]
          Si[[m - 1]] <- run_addis[[6]]
          Ci[[m - 1]] <- run_addis[[7]]
          Ci_sum <- run_addis[[9]]
          Si_sum <- run_addis[[10]]
          gammai_sum <- run_addis[[12]]
          if (K != 0) {
            kappai[[K]] <- run_addis[[4]]
            kappai_star[[K]] <- run_addis[[11]]
            Ci_plus[1:K] <- run_addis[[8]]
          }
          Alpha <- alphai[[m]]
        }
        else {
          Alpha <- alpha
        }
        if (verbose) {
          cat("x=", x, " y=", z, " S=", 
              y, "\n")
          cat("Test number =", m, "\n")
          cat("Additional pval value =", pval[m], 
              "\n")
          cat("Alpha value =", Alpha, "\n")
        }
        if (pval[m] <= Alpha) {
          R[m] <- 1
          if (verbose) {
            cat(V[x], "->", V[y], "<-", V[z], 
                "\n")
            cat("Since pval<Alpha,additional test is rejected;", 
                "Nodes", V[x], "and", V[z], 
                "are dependent given", V[y], "\n")
          }
          tarmat[x, y] <- tarmat[z, y] <- 1
        }
        else {
          R[m] <- 0
          if (verbose) 
            cat("Since pval>Alpha,additional test is accepted;", 
                "Nodes", V[x], "and", V[z], 
                "are independent given", V[y], "\n")
        }
      }
    }
  }
  if (any(tarmat == 1) & GV > 0) {
    WW1 <- unique(which(tarmat == 1, arr.ind = T)[, 2])
    WW1 <- setdiff(WW1, 1:GV)
    if (length(WW1) != 0) {
      for (v1 in 1:length(WW1)) {
        WW2 <- unique(which(tarmat[, WW1[v1]] == 1, arr.ind = T))
        WW3 <- unique(which(g1[, WW1[v1]] == 1, arr.ind = T))
        WW3 <- setdiff(WW3, WW2)
        if (length(WW3) != 0) {
          for (v2 in 1:length(WW3)) {
            x <- WW2[1]
            y <- WW1[v1]
            z <- WW3[v2]
            if ((!is.na(x) & !is.na(y) & !is.na(z)) & 
                (x != 0 & y != 0 & z != 0) & (x != y & 
                                              y != z & z != x)) {
              if (g1[y, z] == 1 & tarmat[y, z] != 1 & 
                  tarmat[z, y] != 1 & ((y %in% gInput@sepset[[x]][[z]]) || 
                                       (y %in% gInput@sepset[[z]][[x]]))) {
                tarmat[y, z] <- 1
                tarmat[z, y] <- 0
              }
              if (g1[y, z] == 1 & g1[x, z] != 1 & tarmat[x, 
                                                         y] == 1 & tarmat[y, x] != 1 & tarmat[y, 
                                                                                              z] != 1 & tarmat[z, y] != 1 & !(y %in% 
                                                                                                                              gInput@sepset[[x]][[z]]) & !(y %in% gInput@sepset[[z]][[x]])) {
                m <- m + 1
                if (m >= length(R)) {
                  R <- c(R, numeric(m/2))
                  pval <- c(pval, numeric(m/2))
                  kappai <- c(kappai, numeric(m/2))
                  kappai_star <- c(kappai_star, numeric(m/2))
                  Ci <- c(Ci, numeric(m/2))
                  Si <- c(Si, numeric(m/2))
                  Ci_plus <- c(Ci_plus, numeric(m/2))
                  gammai <- c(gammai, numeric(m/2))
                  alphai <- c(alphai, numeric(m/2))
                }
                if (indepTest == "gaussCItest") {
                  pval[m] <- gaussCItest(x, z, y, suffStat)
                }
                if (indepTest == "disCItest") {
                  pval[m] <- disCItest(x, z, y, suffStat)
                }
                if (FDRcontrol == "LOND") {
                  alphai[m] <- SeqFDR(m, FDR, a = 2, 
                                      R)
                  Alpha <- alphai[m]
                }
                else if (FDRcontrol == "ADDIS") {
                  run_addis <- addis(alpha = FDR, tau = tau, 
                                     lambda = lambda, iter = m, w0 = w0, 
                                     pval = pval, alphai = alphai, gammai = gammai, 
                                     kappai = kappai, kappai_star = kappai_star, 
                                     K = K, Ci = Ci, Si = Si, Ri = R, 
                                     Ci_plus = Ci_plus, Ci_sum = Ci_sum, 
                                     Si_sum = Si_sum, gammai_sum = gammai_sum, 
                                     normalizer = normalizer, exponent = exponent)
                  alphai[[m]] <- run_addis[[1]]
                  gammai[[m]] <- run_addis[[2]]
                  K <- run_addis[[3]]
                  R[[m]] <- run_addis[[5]]
                  Si[[m - 1]] <- run_addis[[6]]
                  Ci[[m - 1]] <- run_addis[[7]]
                  Ci_sum <- run_addis[[9]]
                  Si_sum <- run_addis[[10]]
                  gammai_sum <- run_addis[[12]]
                  if (K != 0) {
                    kappai[[K]] <- run_addis[[4]]
                    kappai_star[[K]] <- run_addis[[11]]
                    Ci_plus[1:K] <- run_addis[[8]]
                  }
                  Alpha <- alphai[[m]]
                }
                else {
                  Alpha <- alpha
                }
                if (verbose) {
                  cat("Additional pval value =", 
                      pval[m], "\n")
                  cat("Alpha value =", Alpha, "\n")
                }
                if (pval[m] <= Alpha) {
                  R[m] <- 1
                  if (verbose) {
                    V <- colnames(g)
                    cat("\n", V[x], "->", 
                        V[y], "<-", V[z], "\n")
                    cat("Since pval<Alpha,additional test is rejected", 
                        "Nodes", V[x], "and", 
                        V[z], "are dependent given", 
                        V[y], "\n")
                  }
                  tarmat[z, y] <- 1
                }
                else {
                  tarmat[y, z] <- 1
                  R[m] <- 0
                }
              }
              if (g1[y, z] == 1 & g1[x, z] == 1 & g1[x, 
                                                     y] == 1 & !(z %in% gInput@sepset[[x]][[y]]) & 
                  !(y %in% gInput@sepset[[x]][[z]])) {
                tarmat[y, z] <- 1
                tarmat[z, y] <- 1
              }
            }
          }
        }
      }
      T1 <- which(tarmat == 1, arr.ind = T)
      G1 <- which(g == 1, arr.ind = T)
      D1 <- setdiff(G1, T1)
      if (length(D1) != 0) {
        for (d1 in 1:length(D1)) {
          Rem1_row <- which(g[, D1[d1]] == 1, arr.ind = T)
          Rem1_col <- which(g[D1[d1], ] == 1, arr.ind = T)
          Rem1 <- c(Rem1_row, Rem1_col)
          Rem2 <- which(tarmat[, Rem1] == 1, arr.ind = T)
          if (length(Rem1) != 0 & length(Rem2) != 0) {
            for (d11 in 1:length(Rem1)) {
              x <- Rem2[d11]
              y <- Rem1[d11]
              z <- D1[d1]
            }
            if ((!is.na(x) & !is.na(y) & !is.na(z)) & 
                (x != 0 & y != 0 & z != 0) & (x != y & 
                                              y != z & z != x)) {
              if (g1[y, z] == 1 & tarmat[y, z] != 1 & 
                  tarmat[z, y] != 1 & ((y %in% gInput@sepset[[x]][[z]]) || 
                                       (y %in% gInput@sepset[[z]][[x]]))) {
                tarmat[y, z] <- 1
                tarmat[z, y] <- 0
              }
              if (g1[y, z] == 1 & g1[x, z] != 1 & tarmat[x, 
                                                         y] == 1 & tarmat[y, x] != 1 & tarmat[y, 
                                                                                              z] != 1 & tarmat[z, y] != 1 & !(y %in% 
                                                                                                                              gInput@sepset[[x]][[z]]) & !(y %in% gInput@sepset[[z]][[x]])) {
                m <- m + 1
                if (m >= length(R)) {
                  R <- c(R, numeric(m/2))
                  pval <- c(pval, numeric(m/2))
                  kappai <- c(kappai, numeric(m/2))
                  kappai_star <- c(kappai_star, numeric(m/2))
                  Ci <- c(Ci, numeric(m/2))
                  Si <- c(Si, numeric(m/2))
                  Ci_plus <- c(Ci_plus, numeric(m/2))
                  gammai <- c(gammai, numeric(m/2))
                  alphai <- c(alphai, numeric(m/2))
                }
                if (indepTest == "gaussCItest") {
                  pval[m] <- gaussCItest(x, z, y, suffStat)
                }
                if (indepTest == "disCItest") {
                  pval[m] <- disCItest(x, z, y, suffStat)
                }
                if (FDRcontrol == "LOND") {
                  alphai[m] <- SeqFDR(m, FDR, a = 2, 
                                      R)
                  Alpha <- alphai[m]
                }
                else if (FDRcontrol == "ADDIS") {
                  run_addis <- addis(alpha = FDR, tau = tau, 
                                     lambda = lambda, iter = m, w0 = w0, 
                                     pval = pval, alphai = alphai, gammai = gammai, 
                                     kappai = kappai, kappai_star = kappai_star, 
                                     K = K, Ci = Ci, Si = Si, Ri = R, 
                                     Ci_plus = Ci_plus, Ci_sum = Ci_sum, 
                                     Si_sum = Si_sum, gammai_sum = gammai_sum, 
                                     normalizer = normalizer, exponent = exponent)
                  alphai[[m]] <- run_addis[[1]]
                  gammai[[m]] <- run_addis[[2]]
                  K <- run_addis[[3]]
                  R[[m]] <- run_addis[[5]]
                  Si[[m - 1]] <- run_addis[[6]]
                  Ci[[m - 1]] <- run_addis[[7]]
                  Ci_sum <- run_addis[[9]]
                  Si_sum <- run_addis[[10]]
                  gammai_sum <- run_addis[[12]]
                  if (K != 0) {
                    kappai[[K]] <- run_addis[[4]]
                    kappai_star[[K]] <- run_addis[[11]]
                    Ci_plus[1:K] <- run_addis[[8]]
                  }
                  Alpha <- alphai[[m]]
                }
                else {
                  Alpha <- alpha
                }
                if (pval[m] <= Alpha) {
                  R[m] <- 1
                  if (verbose) {
                    V <- colnames(g)
                    cat("\n", V[x], "->", 
                        V[y], "<-", V[z], "\n")
                    cat("Since pval<Alpha,additional test is rejected", 
                        "Nodes", V[x], "and", 
                        V[z], "are dependent given", 
                        V[y], "\n")
                  }
                  tarmat[z, y] <- 1
                }
                else {
                  tarmat[y, z] <- 1
                  R[m] <- 0
                }
              }
              if (g1[y, z] == 1 & g1[x, z] == 1 & g1[x, 
                                                     y] == 1 & !(z %in% gInput@sepset[[x]][[y]]) & 
                  !(y %in% gInput@sepset[[x]][[z]])) {
                tarmat[y, z] <- 1
                tarmat[z, y] <- 1
              }
            }
          }
        }
      }
      ind <- which(g1 == 1, arr.ind = TRUE)
      for (i in seq_len(nrow(ind))) {
        x <- ind[i, 1]
        y <- ind[i, 2]
        if (tarmat[x, y] == 0 & tarmat[y, x] == 0) {
          tarmat[x, y] <- 1
          tarmat[y, x] <- 1
        }
      }
    }
  }
  if (any(tarmat == 1) & GV == 0) {
    WW1 <- unique(which(tarmat == 1, arr.ind = T)[, 2])
    if (length(WW1) != 0) {
      for (v1 in 1:length(WW1)) {
        WW2 <- unique(which(tarmat[, WW1[v1]] == 1, arr.ind = T))
        WW3 <- unique(which(g1[, WW1[v1]] == 1, arr.ind = T))
        WW3 <- setdiff(WW3, WW2)
        if (length(WW3) != 0) {
          for (v2 in 1:length(WW3)) {
            x <- WW2[1]
            y <- WW1[v1]
            z <- WW3[v2]
            if ((!is.na(x) & !is.na(y) & !is.na(z)) & 
                (x != 0 & y != 0 & z != 0) & (x != y & 
                                              y != z & z != x)) {
              if (g1[y, z] == 1 & tarmat[y, z] != 1 & 
                  tarmat[z, y] != 1 & ((y %in% gInput@sepset[[x]][[z]]) || 
                                       (y %in% gInput@sepset[[z]][[x]]))) {
                tarmat[y, z] <- 1
                tarmat[z, y] <- 0
              }
              if (g1[y, z] == 1 & g1[x, z] != 1 & tarmat[x, 
                                                         y] == 1 & tarmat[y, x] != 1 & tarmat[y, 
                                                                                              z] != 1 & tarmat[z, y] != 1 & !(y %in% 
                                                                                                                              gInput@sepset[[x]][[z]]) & !(y %in% gInput@sepset[[z]][[x]])) {
                m <- m + 1
                if (m >= length(R)) {
                  R <- c(R, numeric(m/2))
                  pval <- c(pval, numeric(m/2))
                  kappai <- c(kappai, numeric(m/2))
                  kappai_star <- c(kappai_star, numeric(m/2))
                  Ci <- c(Ci, numeric(m/2))
                  Si <- c(Si, numeric(m/2))
                  Ci_plus <- c(Ci_plus, numeric(m/2))
                  gammai <- c(gammai, numeric(m/2))
                  alphai <- c(alphai, numeric(m/2))
                }
                if (indepTest == "gaussCItest") {
                  pval[m] <- gaussCItest(x, z, y, suffStat)
                }
                if (indepTest == "disCItest") {
                  pval[m] <- disCItest(x, z, y, suffStat)
                }
                if (FDRcontrol == "LOND") {
                  alphai[m] <- SeqFDR(m, FDR, a = 2, 
                                      R)
                  Alpha <- alphai[m]
                }
                else if (FDRcontrol == "ADDIS") {
                  run_addis <- addis(alpha = FDR, tau = tau, 
                                     lambda = lambda, iter = m, w0 = w0, 
                                     pval = pval, alphai = alphai, gammai = gammai, 
                                     kappai = kappai, kappai_star = kappai_star, 
                                     K = K, Ci = Ci, Si = Si, Ri = R, 
                                     Ci_plus = Ci_plus, Ci_sum = Ci_sum, 
                                     Si_sum = Si_sum, gammai_sum = gammai_sum, 
                                     normalizer = normalizer, exponent = exponent)
                  alphai[[m]] <- run_addis[[1]]
                  gammai[[m]] <- run_addis[[2]]
                  K <- run_addis[[3]]
                  R[[m]] <- run_addis[[5]]
                  Si[[m - 1]] <- run_addis[[6]]
                  Ci[[m - 1]] <- run_addis[[7]]
                  Ci_sum <- run_addis[[9]]
                  Si_sum <- run_addis[[10]]
                  gammai_sum <- run_addis[[12]]
                  if (K != 0) {
                    kappai[[K]] <- run_addis[[4]]
                    kappai_star[[K]] <- run_addis[[11]]
                    Ci_plus[1:K] <- run_addis[[8]]
                  }
                  Alpha <- alphai[[m]]
                }
                else {
                  Alpha <- alpha
                }
                if (verbose) {
                  cat("Additional pval value =", 
                      pval[m], "\n")
                  cat("Alpha value =", Alpha, "\n")
                }
                if (pval[m] <= Alpha) {
                  R[m] <- 1
                  if (verbose) {
                    V <- colnames(g)
                    cat("\n", V[x], "->", 
                        V[y], "<-", V[z], "\n")
                    cat("Since pval<Alpha,additional test is rejected", 
                        "Nodes", V[x], "and", 
                        V[z], "are dependent given", 
                        V[y], "\n")
                  }
                  tarmat[z, y] <- 1
                }
                else {
                  tarmat[y, z] <- 1
                  R[m] <- 0
                }
              }
              if (g1[y, z] == 1 & g1[x, z] == 1 & g1[x, 
                                                     y] == 1 & !(z %in% gInput@sepset[[x]][[y]]) & 
                  !(y %in% gInput@sepset[[x]][[z]])) {
                tarmat[y, z] <- 1
                tarmat[z, y] <- 1
              }
            }
          }
        }
      }
      T1 <- which(tarmat == 1, arr.ind = T)
      G1 <- which(g == 1, arr.ind = T)
      D1 <- setdiff(G1, T1)
      if (length(D1) != 0) {
        for (d1 in 1:length(D1)) {
          Rem1 <- which(g[, D1[d1]] == 1, arr.ind = T)
          Rem2 <- which(tarmat[, Rem1] == 1, arr.ind = T)
          if (length(Rem1) != 0 & length(Rem2) != 0) {
            for (d11 in 1:length(Rem1)) {
              x <- Rem2[d11]
              y <- Rem1[d11]
              z <- D1[d1]
            }
            if ((!is.na(x) & !is.na(y) & !is.na(z)) & 
                (x != 0 & y != 0 & z != 0) & (x != y & 
                                              y != z & z != x)) {
              if (g1[y, z] == 1 & tarmat[y, z] != 1 & 
                  tarmat[z, y] != 1 & ((y %in% gInput@sepset[[x]][[z]]) || 
                                       (y %in% gInput@sepset[[z]][[x]]))) {
                tarmat[y, z] <- 1
                tarmat[z, y] <- 0
              }
              if (g1[y, z] == 1 & g1[x, z] != 1 & tarmat[x, 
                                                         y] == 1 & tarmat[y, x] != 1 & tarmat[y, 
                                                                                              z] != 1 & tarmat[z, y] != 1 & !(y %in% 
                                                                                                                              gInput@sepset[[x]][[z]]) & !(y %in% gInput@sepset[[z]][[x]])) {
                m <- m + 1
                if (m >= length(R)) {
                  R <- c(R, numeric(m/2))
                  pval <- c(pval, numeric(m/2))
                  kappai <- c(kappai, numeric(m/2))
                  kappai_star <- c(kappai_star, numeric(m/2))
                  Ci <- c(Ci, numeric(m/2))
                  Si <- c(Si, numeric(m/2))
                  Ci_plus <- c(Ci_plus, numeric(m/2))
                  gammai <- c(gammai, numeric(m/2))
                  alphai <- c(alphai, numeric(m/2))
                }
                if (indepTest == "gaussCItest") {
                  pval[m] <- gaussCItest(x, z, y, suffStat)
                }
                if (indepTest == "disCItest") {
                  pval[m] <- disCItest(x, z, y, suffStat)
                }
                if (FDRcontrol == "LOND") {
                  alphai[m] <- SeqFDR(m, FDR, a = 2, 
                                      R)
                  Alpha <- alphai[m]
                }
                else if (FDRcontrol == "ADDIS") {
                  run_addis <- addis(alpha = FDR, tau = tau, 
                                     lambda = lambda, iter = m, w0 = w0, 
                                     pval = pval, alphai = alphai, gammai = gammai, 
                                     kappai = kappai, kappai_star = kappai_star, 
                                     K = K, Ci = Ci, Si = Si, Ri = R, 
                                     Ci_plus = Ci_plus, Ci_sum = Ci_sum, 
                                     Si_sum = Si_sum, gammai_sum = gammai_sum, 
                                     normalizer = normalizer, exponent = exponent)
                  alphai[[m]] <- run_addis[[1]]
                  gammai[[m]] <- run_addis[[2]]
                  K <- run_addis[[3]]
                  R[[m]] <- run_addis[[5]]
                  Si[[m - 1]] <- run_addis[[6]]
                  Ci[[m - 1]] <- run_addis[[7]]
                  Ci_sum <- run_addis[[9]]
                  Si_sum <- run_addis[[10]]
                  gammai_sum <- run_addis[[12]]
                  if (K != 0) {
                    kappai[[K]] <- run_addis[[4]]
                    kappai_star[[K]] <- run_addis[[11]]
                    Ci_plus[1:K] <- run_addis[[8]]
                  }
                  Alpha <- alphai[[m]]
                }
                else {
                  Alpha <- alpha
                }
                if (pval[m] <= Alpha) {
                  R[m] <- 1
                  if (verbose) {
                    V <- colnames(g)
                    cat("\n", V[x], "->", 
                        V[y], "<-", V[z], "\n")
                    cat("Since pval<Alpha,additional test is rejected", 
                        "Nodes", V[x], "and", 
                        V[z], "are dependent given", 
                        V[y], "\n")
                  }
                  tarmat[z, y] <- 1
                }
                else {
                  tarmat[y, z] <- 1
                  R[m] <- 0
                }
              }
              if (g1[y, z] == 1 & g1[x, z] == 1 & g1[x, 
                                                     y] == 1 & !(z %in% gInput@sepset[[x]][[y]]) & 
                  !(y %in% gInput@sepset[[x]][[z]])) {
                tarmat[y, z] <- 1
                tarmat[z, y] <- 1
              }
            }
          }
        }
      }
      ind <- which(g1 == 1, arr.ind = TRUE)
      for (i in seq_len(nrow(ind))) {
        x <- ind[i, 1]
        y <- ind[i, 2]
        if (tarmat[x, y] == 0 & tarmat[y, x] == 0) {
          tarmat[x, y] <- 1
          tarmat[y, x] <- 1
        }
      }
    }
  }
  if (all(tarmat == 0)) {
    tarmat <- g1
  }
  if (GV > 0 & (any(tarmat[1:GV, ] == 1) || any(tarmat[, 1:GV] == 
                                                1)) & all(tarmat[-c(1:GV), -c(1:GV)] == 0)) {
    tarmat1 <- g1
    tarmat1[1:GV, ] <- tarmat[1:GV, ]
    tarmat1[, 1:GV] <- tarmat[, 1:GV]
    tarmat <- tarmat1
  }
  gInput@graph <- as(tarmat, "graphNEL")
  gInput@R <- R[1:m]
  gInput@K <- K
  gInput@pval <- pval[1:m]
  gInput@kappai <- kappai[1:K]
  gInput@kappai_star <- kappai_star[1:K]
  gInput@Ci <- Ci[1:m]
  gInput@Si <- Si[1:m]
  gInput@Ci_plus <- Ci_plus[1:K]
  gInput@gammai <- gammai[1:m]
  gInput@gammai_sum <- gammai_sum
  gInput
}



#=====================================================================================================