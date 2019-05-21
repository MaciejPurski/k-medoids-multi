library(cluster)

data("USArrests")
df <- scale(USArrests)



args = commandArgs(trailingOnly=TRUE)

k <- args[1]
ex <- args[2]
k <- 4
ex <- 20
df <- df[1:ex,]


pam.res <- pam(df, k, trace.lev = 100,stand=FALSE)

# my pam implementation
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))


# calculate dissimilarities
diss.matrix <- lapply(seq(to=nrow(df)), function(i) {
                if (i + 1 <= nrow(df)) sapply(seq(from=i+1, nrow(df)), function(j) euc.dist(df[i,], df[j, ]))})

diss.matrix[[50]] <- 0L

# function to obtain dissimilarities
diss2 <- function(i, j) {
  if (i == j)
    return(0)
  diss.matrix[[min(i,j)]][max(i,j) - min(i,j)]
}

diss <- function(i, j) {
  return(euc.dist(df[i,], df[j, ]))
}

dissToClosestMedoid <- function(j) {
  min(sapply(C, function(n) diss(j, n)))
}

closestMedoid <- function(j) {
  C[which.min(sapply(C, function(n) diss(j, n)))]
}


dissToSecondClosest <- function(j) {
  t <- sapply(C, function(n) diss(j, n))
  sort(t)

  print(t[2])
  return(t[2])
}

computeM <- function(i, j) {
  if (i == j)
    return(0)
  
  min(dissToClosestMedoid(j), diss(i, j))
}

compute_C2 <- function(i, j, h, d, s) {
  if (j == h)
    return(0)

  if (diss(j,i) > d) {
    return(min(diss(j, h) - d, 0))
  }
  
  if (diss(j, i) == d) {
    return(min(diss(j,h), s) - d)
  }
}

compute_C4 <- function(i, j, h, d, s) {
  if (j == h)
    return(0)
  
  if (diss(j,i) > d) {
    if (diss(j,h) >= d)
      return(0)
    if (diss(j,h) < d)
      return(diss(j,h) - d)
  }
  
  if (diss(j, i) == d) {
    if (diss(j,h) < s)
      return(diss(j,h) - d)
    if (diss(j,h) >= s)
      return(s - d)
  }
}

compute_C5 <- function(i, h) {
  before <- sum(sapply(O, function(j) dissToClosestMedoid(j)))
  print("before")
  print(before)
  O[which(O %in% h)] <- i
  C[which(C %in% i)] <- h

  print("after")
  after <- sum(sapply(O, function(j) dissToClosestMedoid(j)))
  print(after)
  
  diff <- after - before
  
  O[which(O %in% i)] <- h
  C[which(C %in% h)] <- i
  
  return(diff)
}

compute_C6 <- function(i, j, h) {
  if (j == h)
    return(0)
  
  if (diss(j,i) > dissToClosestMedoid(j)) {
    return(min(diss(j, h) - dissToClosestMedoid(j), 0))
  }
  
  if (diss(j, i) == dissToClosestMedoid(j)) {
    return(min(diss(j,h), dissToSecondClosest(j)) - dissToClosestMedoid(j))
  }
}

compute_C6 <- function(i, j, h) {
  if (O[j] == h)
    return(0)
  if (i == meds[j])
    return(min(diss(O[j], h), D2[j]) - diss(O[j], meds[j]))
  else
    return(min(diss(O[j], h) - diss(O[j], meds[j]), 0))
}



# BUILD STEP
O <- 1:nrow(df)
C <- c()

sumOfDistances <- sapply(O, function(i) sum(sapply(O, function(j) {diss(i, j)})))

C <- append(C, which.min(sumOfDistances))
O <- O[-C]

while (length(C) != k) {
  # search for candidates
  M <- sapply(O, function(i) sum(sapply(O, function(j) computeM(i, j))))
  best.candidate <- which.min(M)
  C <- append(C, O[best.candidate])
  O <- O[-best.candidate]
}

while (TRUE) {
  dists <- lapply(O, function(j) sort(sapply(C, function(i) diss(j, i))))
  
  D <- sapply(dists, function(j) j[1])
  S <- sapply(dists, function(j) j[2])
  
  meds <- sapply(O, closestMedoid)
  D2 <- sapply(O, dissToSecondClosest)
  
  sum_diss <- sum(sapply(O, function(j) D[which(O %in% j)]))


  T <- sapply(C, function(i)
        sapply(O, function(h) sum(sapply(O, function(j)
          compute_C2(i, j, h, D[which(O %in% j)], S[which(O %in% j)])))))

  
  print(C)
  print("diss: ")
  print(sum_diss)
  print(T)
  print("min chosen")

  T.min <- min(T)
  print(T.min)
  if (T.min >= 0)
    break
  # powinno byc wymienione 11 z 6 a jest 33 z 8
  
  indices <- which(T == T.min, arr.ind = TRUE)
  
  print(cat("swap ", indices[1], " ", indices[2]))
  print(cat("ind ", O[indices[1]], " ", C[indices[2]]))
  
  newMedoid <- O[indices[1]]
  O[indices[1]] <- C[indices[2]]
  C[indices[2]] <- newMedoid
}

df[C, ]
pam.res$medoids

