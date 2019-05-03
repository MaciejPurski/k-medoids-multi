library(cluster)

data("USArrests")
df <- scale(USArrests)



args = commandArgs(trailingOnly=TRUE)

k <- args[1]
ex <- args[2]
df <- df[1:ex,]


pam.res <- pam(df, k, trace.lev = 100,stand=FALSE)

# my pam implementation
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

diss <- function(i, j) {
  return(euc.dist(df[i,], df[j, ]))
}

updateElements <- function() {
  # calculate distances to medoids
  d <- sapply(1:nrow(O), function(i) sapply(C, function(j) diss(O[i, "id"], j)))
  O$med <- apply(d, 2, function(i) C[which.min(i)])
  O$dissToMed <- apply(O, 1, function(i) diss(i[["id"]], i[["med"]]))
  O$dissToSecond <- sapply(1:nrow(O), function(i) min(d[,i][d[,i]!=min(d[,i])]))
  
  O
}

compute <- function(i, h) {
  # difference when h is not counted anymore to whole cost
  diff <- -h[["dissToMed"]]
  
  # difference when i is no longer a medoid and starts counting to the cost
  d <- dissToMedoids(i)
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

dissToMedoids <- function(i) {
  sapply(C, function(j) diss(i, j))
}

compute2 <- function(i, j, h) {
  if (j == h)
    return(0)
}

dissToClosestMedoid <- function(j) {
  min(sapply(C, function(n) diss(j, n)))
}

computeM <- function(i, j) {
  if (i == j)
    return(0)
  
  min(dissToClosestMedoid(j), diss(i, j))
}

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

print(C)

O <- data.frame(id = O, med = 1:length(O), dissToMed = 1:length(O), dissToSecond = 1:length(O))

while (TRUE) {
  O <- updateElements()
  
  T <- sapply(C, function(i) apply(O, 1, function(h) compute(i, h)))

  print("diss: ")
  print(sum(O$dissToMed))
  print("min chosen")
  
  T.min <- min(T)
  print(T.min)
  if (T.min >= 0)
    break
  # powinno byc wymienione 11 z 6 a jest 33 z 8
  
  indices <- which(T == T.min, arr.ind = TRUE)
  
  print(cat("swap ", indices[1], " ", indices[2]))
  print(cat("ind ", O[indices[1], "id"], " ", C[indices[2]]))
  
  newMedoid <- O[indices[1], "id"]
  O[indices[1], "id"] <- C[indices[2]]
  C[indices[2]] <- newMedoid

}

df[C, ]
pam.res$medoids
