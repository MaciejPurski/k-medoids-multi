library(cluster)
library(kmedoids)

data("USArrests")
df <- scale(USArrests)

obj.fun = rep(0, 10)

for (i in 1:10) {
  print(paste(i, "/10"))
  pam.res = cluster::pam(df, i)
  res = kmedoids::kmedoidsMulti(df, i)
  
  if (setequal(pam.res$medoids, res$medoids) &
      setequal(pam.res$id.med, res$id.med)) {
    print("PASSED")
  } else {
    print("FAILED")
  }
  
  obj.fun[i] <- res$objective[["swap"]]
}

plot(1:10, obj.fun, "b", xaxt = "n", main = "PAM clustering USArrests, n = 50", xlab = "k", ylab = "objective function")
axis(1, at=1:10, labels=1:10)
abline(v = 4, lty = 2, col = "blue")

res = kmedoidsMulti(df, 4)
pam.res = pam(df, 4)
