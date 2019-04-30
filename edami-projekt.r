library(cluster)

data("USArrests")
df <- scale(USArrests)

pam.res <- pam(df, 4)


# my pam implementation

k <- 4
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

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

D <- lapply(dists, function(j) j[1])
S <- lapply(dists, function(j) j[2])


# consider all possible swaps
T <- sapply(C, function(i) sapply(O,
                function(h) sum(sapply(1:length(O),
                  function(j) {
                    if (O[[j]] != i)
                      compute_C(i, O[[j]], h, D[[j]], S[[j]], df)
                    else
                      0
                  }))))

T.min <- min(T)


if (T.min >= 0)
  break

indices <- which(T == T.min, arr.ind = TRUE)

nO <- indices[1, 2]
nC <- indices[1, 1]

O <- append(O, C[[nO]])
C <- append(C, O[[nC]])

O <- O[-nC]
C <- C[-nO]

print(unlist(O))
print(unlist(C))
}

df[unlist(C),]
pam.res$medoids
