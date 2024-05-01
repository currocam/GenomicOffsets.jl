# This script relies on pre-existing libraries to make regression test
using RCall, Serialization

R"""
library(LEA)
data("offset_example")
Y <- offset_example$geno
X <- offset_example$env
Xpred <- offset_example$env.pred
"""
@rget Y; @rget X; @rget Xpred

example_data = (Y=Y, X=X, Xpred=Xpred)
open("example_data.jld", "w") do io
    serialize(io, example_data)
end

# RONA from Gain et al
R"""
rona_fn <- function(Y, X, Xpred) {
    nb_var <- ncol(X)
    n <- nrow(Y)
    mod_lm <- lm(as.matrix(Y) ~ ., data = data.frame(X)) 
    sm <- summary(mod_lm)
    B <- sapply(sm, FUN = function(x) x$coeff[1:(nb_var + 1), 1])
    X <-  cbind(rep(1.0, n), X)
    Xpred <- cbind(rep(1.0, n), Xpred)
    Y.fit <- as.matrix(X) %*% as.matrix(B)
    Y.pred <- as.matrix(Xpred) %*% as.matrix(B)
    allele_frequency_shift <- abs(Y.fit - Y.pred)
    rowMeans(allele_frequency_shift)
    }
rona <- rona_fn(Y, X, Xpred)
"""

@rget rona
open("rona.jld", "w") do io
    serialize(io, rona)
end