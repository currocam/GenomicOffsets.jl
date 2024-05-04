# This script relies on pre-existing libraries to make regression test
using RCall, Serialization

R"""
library(LEA)
data("offset_example")
Y <- offset_example$geno
X <- offset_example$env
Xpred <- offset_example$env.pred
"""
@rget Y;
@rget X;
@rget Xpred;

example_data = (Y=Y, X=X, Xpred=Xpred)
open("example_data.jld", "w") do io
    return serialize(io, example_data)
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
    return serialize(io, rona)
end

# RDA from Gain et al
R"""
rda_fn <- function(Y, X, Xpred){  
  nb_var <- ncol(X)
  n <- nrow(Y)
  mod_lm <- lm(as.matrix(Y) ~ ., data = data.frame(X)) 
  sm <- summary(mod_lm)
  B <- sapply(sm, FUN = function(x) x$coeff[1:(nb_var + 1), 1])
  X <-  as.matrix(cbind(rep(1.0, n), X))
  Xpred <-  as.matrix(cbind(rep(1.0, n), Xpred))
  pc = prcomp(X %*% B)
  proj.x = predict(pc, X %*% B)
  proj.xpred = predict(pc, Xpred %*% B)
  rda.go = rowSums((proj.x - proj.xpred)[,1:ncol(X)]^2)/ncol(Y)
  rda.go
}
rda <- rda_fn(Y, X, Xpred)
"""

@rget rda
open("rda.jld", "w") do io
    return serialize(io, rda)
end

R"""
tw_fn <- function(Y){
    pca.res <- prcomp(Y)
    eigenvalues <- pca.res$sdev^2
    eigenvalues <- eigenvalues[order(-eigenvalues)]
    tracywindom_statistic <- AssocTests::tw(
        eigenvalues =eigenvalues, eigenL = length(eigenvalues)
        )
}
tw <- tw_fn(Y)
"""

@rget tw

open("tracywidom.jld", "w") do io
    return serialize(io, tw)
end

R"""
lfmm2_k2 <- LEA::lfmm2(Y, X, 2, effect.sizes = TRUE)
lfmm2_k2 <- list(
    U = as.matrix(lfmm2_k2@U),
    V = as.matrix(lfmm2_k2@V),
    B = as.matrix(lfmm2_k2@B),
    pvalues = LEA::lfmm2.test(lfmm2_k2, input=Y, env=X, full=TRUE, genomic.control=TRUE)$pvalues
)
lfmm2_k3 <- LEA::lfmm2(Y, X, 3, effect.sizes = TRUE)
lfmm2_k3 <- list(
    U = as.matrix(lfmm2_k3@U),
    V = as.matrix(lfmm2_k3@V),
    B = as.matrix(lfmm2_k3@B),
    pvalues = LEA::lfmm2.test(lfmm2_k3, input=Y, env=X, full=TRUE, genomic.control=TRUE)$pvalues
)
lfmm2_k25 <- LEA::lfmm2(Y, X, 25, effect.sizes = TRUE)
lfmm2_k25 <- list(
    U = as.matrix(lfmm2_k25@U),
    V = as.matrix(lfmm2_k25@V),
    B = as.matrix(lfmm2_k25@B),
    pvalues = LEA::lfmm2.test(lfmm2_k25, input=Y, env=X, full=TRUE, genomic.control=TRUE)$pvalues
)
"""

@rget lfmm2_k2;
@rget lfmm2_k3;
@rget lfmm2_k25;

lfmm2s = (lfmm2_k2=lfmm2_k2, lfmm2_k3=lfmm2_k3, lfmm2_k25=lfmm2_k25)

open("lfmm2.jld", "w") do io
    return serialize(io, lfmm2s)
end

R"""
geometric_k2 <- LEA::genetic.gap(Y, X, X, Xpred, K = 2, scale = TRUE)$offset
geometric_k3 <- LEA::genetic.gap(Y, X, X, Xpred, K = 3, scale = TRUE)$offset
geometric_k25 <- LEA::genetic.gap(Y, X, X, Xpred, K = 25, scale = TRUE)$offset
"""

@rget geometric_k2;
@rget geometric_k3;
@rget geometric_k25;

geometric = (geometric_k2=geometric_k2, geometric_k3=geometric_k3,
             geometric_k25=geometric_k25)

open("geometric.jld", "w") do io
    return serialize(io, geometric)
end

R"""
l2norm <- \(u, v) sqrt(sum((u-v)^2))
gf_fn <- function(Y, X, Xpred){
    Xcols <- paste0("X", 1:ncol(X))
    colnames(X) <- Xcols
    colnames(Xpred) <- Xcols
    Ycols <- paste0("Y", 1:ncol(Y))
    colnames(Y) <- Ycols
    df_gf <- data.frame(Y, X)
    gf <- gradientForest::gradientForest(data=df_gf, predictor.vars=Xcols, response.vars = Ycols, ntree=100)
    currentcumimp <- predict(gf, data.frame(X))
    futurcumimp <- predict(gf, data.frame(Xpred))
    nb_ind <- nrow(X)
    genetic_offset <- c()
    for (i in seq(1,nb_ind)){
        genetic_offset <- c(genetic_offset, l2norm(futurcumimp[i,], currentcumimp[i,]))
    }
    genetic_offset
}
gradient_offset <- gf_fn(Y, X, Xpred)
"""

@rget gradient_offset

open("gradientforest.jld", "w") do io
    return serialize(io, gradient_offset)
end