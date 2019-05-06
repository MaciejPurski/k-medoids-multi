library(cluster)

data("USArrests")
df <- scale(USArrests)

args = commandArgs(trailingOnly=TRUE)

k <- args[1]
k <- 8

pam.res <- pam(df, 8, trace.lev = 100)

# my pam implementation
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))



diss <- lapply(seq(to=nrow(df)), function(i) {
              if (i + 1 > nrow(df)) return()
              sapply(seq(from=i+1, nrow(df)), function(j) euc.dist(df[i,], df[j, ]))}
              )
length(unlist(diss))
# Function to compute dissimilarity of an object x
# to the closest object in S
dissToClosest <- function(x, C, d) {
  min(unlist(lapply(C, function(y) euc.dist(d[x, ], d[y, ]) )))
}

computeC <- function(i, j, C, d) {
  min(dissToClosest(j, C, d), euc.dist(d[i,], d[j,]))
}

distToMedoids <- function(j, C, d) {
  sapply(C, function(n) euc.dist(d[j, ], d[n, ]))
}

compute_C <- function(i, j, h, d, s, df) {
  C_jih <- 0

  if (euc.dist(df[j, ], df[i, ]) > d &
      euc.dist(df[j, ], df[h, ]) > d) {
      C_jih <- 0
  } else if (euc.dist(df[j,], df[i, ]) == d) {
    if (euc.dist(df[j, ], df[h, ]) < s) {
      C_jih <- euc.dist(df[j, ], df[h, ]) - d
    } else {
      C_jih <- s - d
    }
  } else if (euc.dist(df[j,], df[h, ]) < d) {
    C_jih <- euc.dist(df[j,], df[h, ]) - d
  }
  
  C_jih
}

# BUILD STEP
# O is a list of indices of unselected rows
# S is a list of indices of selected rows - clusters
O <- 1:nrow(df)
C <- list()



sumOfDistances <- lapply(O,
                  function(i) Reduce("+", lapply(O,
                  function(j) euc.dist(df[i,], df[j,]))))

# C - selected
C <- append(C, which.min(sumOfDistances))
C
# O - unselected
O <- O[-unlist(C)]

# Compute matrix C
while (length(C) != k) {
  M <- lapply(O, function(i) lapply(O[!O %in% i],
                                    function(j) computeC(i, j, C, df)))

  # find best candidate to add to set S, s is the index in U list

  s <- which.min(lapply(M, function(x) Reduce("+", x)))
  C <- append(C, O[s])
  O <- O[-s]
}

i = 0
while (TRUE) {

dists <- lapply(O, function(j) sort(distToMedoids(j, C, df)))

D <- sapply(dists, function(j) j[1])
S <- sapply(dists, function(j) j[2])


# consider all possible swaps
T <- sapply(C, function(i) sapply(O,
                function(h) sum(sapply(1:length(O),
                  function(j) {
                    if (O[[j]] != h)
                      compute_C(i, O[[j]], h, D[j], S[j], df)
                    else
                      0
                  }))))

T.min <- min(T)

print(T.min)
print(df[unlist(C),])
if (T.min > 0)
  break

indices <- which(T == T.min, arr.ind = TRUE)

newMedoid <- O[[indices[1]]]
O[[indices[1]]] <- C[[indices[2]]]
C[[indices[2]]] <- newMedoid
}

df[unlist(C),]
pam.res$medoids
pam.res$objective
