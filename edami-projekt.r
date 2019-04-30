library(cluster)

data("USArrests")
df <- scale(USArrests)
head(df, n = 3)

pam.res <- pam(df, 8, do.swap = FALSE, trace.lev=10, stand=FALSE)
print(pam.res)
str(pam.res)
View(df)


# my pam implementation

k <- 8
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

dists <- lapply(O, function(j) sort(distToMedoids(j, C, df)))
View(dists)
