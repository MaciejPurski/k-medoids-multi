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

diss2 <- function(i, j) {
  return(euc.dist(df[i,], df[j, ]))
}

diss.matrix <- lapply(seq(to=nrow(df)), function(i) {
                sapply(seq(from=i, nrow(df)), function(j) euc.dist(df[i,], df[j, ]))})
diss.matrix[[50]] <- 0L
View(diss.matrix)

# function to obtain dissimilarities
diss <- function(i, j) {
  if (i == j)
    return(0)
  diss.matrix[[min(i,j)]][max(i,j) - min(i,j)]
}

updateElements <- function() {
  # calculate distances to medoids
  d <- sapply(1:nrow(O), function(i) sapply(C, function(j) diss(O[i, "id"], j)))
  O$med <- apply(d, 2, function(i) C[which.min(i)])
  O$dissToMed <- apply(O, 1, function(i) diss(i[["id"]], i[["med"]]))
  O$dissToSecond <- sapply(1:nrow(O), function(i) min(d[,i][d[,i]!=min(d[,i])]))
  
  O
}



dissToMedoids <- function(i) {
  sapply(C, function(j) diss(i, j))
}

compute2 <- function(i, j, h) {
  if (j == h)
    return(0)
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
