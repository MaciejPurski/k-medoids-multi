kmedoidsMulti <- function(df, k, threads=1) {
  # function to calculate euclidean distance
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  
  # function to obtain dissimilarity from the dissimilarity matrix
  diss <- function(i, j) {
    min.elem <- min(i, j)
    max.elem <- max(i, j)
    diss.matrix[[min.elem]][max.elem - min.elem + 1]
  }
  
  # function to compute dissimilarity with closest medoid
  dissToClosestMedoid <- function(j) {
    min(sapply(C, function(n) diss(j, n)))
  }
  
  # function used in BUILD step in order to find out
  # if an object j will be assigned co candidate i or one of the existing medoids
  computeM <- function(i, j) {
    if (i == O[j])
      return(0)
    
    min(M.closest[j], diss(i, O[j]))
  }
  
  updateElements <- function() {
    # calculate distances to medoids
    d <- sapply(1:nrow(O), function(i) sapply(C, function(j) diss(O[i, "id"], j)))
    
    # find closest medoid
    O$med <<- apply(d, 2, function(i) C[which.min(i)])
    O$dissToMed <<- apply(O, 1, function(i) diss(i[["id"]], i[["med"]]))
    
    # find distance to second closest medoid
    O$dissToSecond <<- sapply(1:nrow(O), function(i) min(d[,i][d[,i]!=min(d[,i])]))
  }
  
  # function to compute result of swapping i medoid with h non-medoid
  compute <- function(i, h) {
    # difference when h is not counted anymore to whole cost
    diff <- -h[["dissToMed"]]
    
    # difference when i is no longer a medoid and starts counting to the cost
    d <- sapply(C, function(j) diss(i, j))
    
    # find second closest medoid (first one is 0)
    sec <- min(d[d != min(d)])
    diff <- diff + min(sec, diss(i, h[["id"]]))
    
    diff <- diff + sum(apply(O, 1, function(j) {
      if (j[["id"]] == h[["id"]])
        return(0)
      if (i == j[["med"]]) {
        min(diss(j[["id"]], h[["id"]]), j[["dissToSecond"]]) - j[["dissToMed"]]
      } else {
        min(diss(j[["id"]], h[["id"]]) - j[["dissToMed"]], 0)
      }
    }))
    
    diff
  }
  
  cl <- parallel::makeCluster(threads)
  
  # calculate dissimilarity matrix
  diss.matrix <- parallel::parLapply(cl, seq(to=nrow(df)), function(i) {
                    sapply(seq(from=i, nrow(df)), function(j) euc.dist(df[i,], df[j, ]))})


  # initialize O - not medoids, C - medoids
  O <- 1:nrow(df)
  C <- c()

  # calculate first medoid
  sumOfDistances <- sapply(O, function(i) sum(sapply(O, function(j) {diss(i, j)})))
  C <- append(C, which.min(sumOfDistances))
  O <- O[-C]
  
  # BUILD STEP
  # goes until k medoids are found
  while (length(C) != k) {
    # search for candidates
    M.closest <- parallel::parSapply(cl, O, function(j) dissToClosestMedoid(j))

    M <- parallel::parSapply(cl, O, function(i) sum(sapply(1:length(O), function(j) computeM(i, j))))
    best.candidate <- which.min(M)
    C <- append(C, O[best.candidate])
    O <- O[-best.candidate]
  }
  
  buildObj <- sum(sapply(O, function(i) dissToClosestMedoid(i))) / nrow(df)

  # O becomes data frame consisting of all the informations needed on medoids calculation
  O <- data.frame(id = O, med = 1:length(O), dissToMed = 1:length(O), dissToSecond = 1:length(O))
  
  # SWAP PHASE
  while (TRUE) {
    updateElements()
    
    # T is a matrix which keeps results of swapping i medoid with h non-medoid
    T <- sapply(C, function(i) apply(O, 1, function(h) compute(i, h)))

    T.min <- min(T)

    if (T.min >= 0)
      break

    indices <- which(T == T.min, arr.ind = TRUE)
    
    # swap medoids
    newMedoid <- O[indices[1], "id"]
    O[indices[1], "id"] <- C[indices[2]]
    C[indices[2]] <- newMedoid
  }
  
  parallel::stopCluster(cl)
  
  swapObj <- sum(O$dissToMed) / nrow(df)
  
  res <- list(
          medoids = df[C,],
          id.med = C,
          clustering = setNames(apply(O, 1, function(i) which(C %in% i[["med"]])),
                                apply(O, 1, function(i) row.names(df)[i[["id"]]])),
          objective = setNames(c(buildObj, swapObj),
                               c("build", "swap"))
  )

  return(res)
}