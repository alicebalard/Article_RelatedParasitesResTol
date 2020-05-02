# contre-example: y/x et x sont correlés, mais x et y ne le sont pas
# e.g. runif for resistance with very small variance
# T/R trade-off but T or R not different between groups

getCorRT <- function(vec){
  n=10
  
  res <- melt(data.frame(sapply(vec, function(x) rnegbin(n, mu=x, theta = 2))))
  # wl <- rep(runif(n, min = 15, max = 20), length(vec))
  wl <- rep(rnorm(n, mean = 15), length(vec))
  
  df <- data.frame(x = res$value,
                   y = wl,
                   gr=rep(LETTERS[1:6], n)[order(rep(LETTERS[1:6], n))])
  
  dfSum <- data.frame(x = as.vector(by(df[,"x"], df[,"gr"], mean)),
                      y = as.vector(by(df[,"y"], df[,"gr"], mean)),
                      gr = unique(df$gr))
  
  mod <-lm(y ~ 0 + x:gr, data=df)
  
  dfRT <- data.frame(meanRes = -as.vector(by(df[,"x"], df[,"gr"], mean)),
                     tol =-mod$coefficients,
                     gr = unique(df$gr))
  
  p1 <- ggplot(dfSum, aes(x, y)) +
    geom_smooth(method = "lm", col = "grey") +
    geom_point(data = df, aes(col=gr), size = 1, pch = 4) +
    geom_point(aes(col=gr, group=gr), size = 4, pch = 1) +
    xlab("parasite load") + ylab("% weight loss")
  
  p2 <- ggplot(dfRT, aes(meanRes, tol)) +
    geom_smooth(method = "lm", col = "grey") +
    geom_point(aes(col=gr, group=gr), size = 4) +
    xlab("resistance") + ylab("tolerance")
  
  plot <- cowplot::plot_grid(p1+theme(legend.position = "none"), p2)
  c1 <- cor.test(dfSum$x, dfSum$y, method = "spearman")
  c2 <- cor.test(dfRT$meanRes, dfRT$tol, method = "spearman")
  
  return(list(c1$p.value, c2$p.value, plot, res, wl))
}

vec = c(1,2,3,1,2,1) # resistance id
test <- getCorRT(vec)
test[[2]]
test[[3]]
  # Run 1000 simulations: how often ResTol correlation is significant?

getCorRT(vec)

?boot
