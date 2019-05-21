# script which performs testing on multiple threads
library("kmedoids")
library("cluster")
df <- read.csv("HTRU_2.csv")

df$X0 <- NULL

df <- scale(df)

test.run <- function(s, k, reps) {
  sapply(1:8, function(t) mean(sapply(1:reps,
                                      function(i)  {
                                        res <- system.time(kmedoidsMulti(df[sample(nrow(df), s), ], 2, t))
                                        print(c(t, i))
                                        print(res[["elapsed"]])
                                        return(res[["elapsed"]])
                                        
                                      } )))
}

test.pam <- function(s, k, reps) {
  mean(sapply(1:reps, function(i)  {
                         res <- system.time(pam(df[sample(nrow(df), s), ], 2))
                         print(i)
                         print(res[["elapsed"]])
                         return(res[["elapsed"]])
                                        
                        } ))
}

t1 <- test.run(100, 2, 10)

pam1 <- test.pam(100, 2, 10)

t2 <- test.run(200, 2, 10)

pam2 <- test.pam(200, 2, 10)


t3 <- test.run(500, 2, 5)

pam3 <- test.pam(500, 2, 5)

t4 <- test.run(1000, 2, 5)

pam4 <- test.pam(1000, 2, 5)

t5 <- test.run(2000, 2, 5)

pam5 <- test.pam(2000, 2, 5)

t <- list(t1, t2, t3, t4, t5)

t6 <- test.run(5000, 2, 5)

pam6 <- test.pam(5000, 2, 5)

t6 <- test.run(5000, 2, 5)
