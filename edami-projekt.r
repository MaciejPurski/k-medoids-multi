library(cluster)


data("USArrests")
df <- scale(USArrests)
head(df, n = 3)

pam.res <- pam(df, 3, do.swap = FALSE)
print(pam.res)
str(pam.res)
View(df)


# my pam implementation

k <- 3
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

# Function to compute dissimilarity of an object x
# to the closest object in S
dissToClosest <- function(x, S, d) {
  min(unlist(lapply(S, function(y) { euc.dist(d[x, ], d[y, ]) })))
}

computeC <- function(i, j, S, d) {
  max(dissToClosest(j, S, d) - euc.dist(d[i,], d[j,]), 0)
}

# BUILD STEP
# U is a list of indices of unselected rows
# S is a list of indices of selected rows, e.t rows which are centers
U <- 1:nrow(df)
S <- list()


sumOfDistances <- lapply(U,
                  function(i) Reduce("+", lapply(U,
                  function(j) euc.dist(df[i,], df[j,]))))

# S - selected
S <- append(S, which.min(sumOfDistances))

# U - unselected
U <- U[!U %in% S]

# Compute matrix C
while (length(S) != k) {
  C <- lapply(U, function(i) lapply(U[!U %in% i],
               function(j) computeC(i, j, S, df)))

  # find best candidate to add to set S, s is the index in U list
  s <- which.max(lapply(C, function(x) Reduce("+", x)))
  S <- append(S, U[s])
  U <- U[-s]
}

unlist(S)
df[unlist(S),]
pam.res$medoids
