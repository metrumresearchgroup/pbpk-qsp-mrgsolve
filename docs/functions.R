logbr <- function() {
  x <- 10^seq(-5,5)
  sort(c(x,3*x))
}
library(ggplot2)
noline <- element_blank()
theme_plain <- function(...) {
  theme_bw() + theme(panel.grid.major=noline,panel.grid.minor=noline, 
                     plot.margin=margin(0.5,0.5,1,0.5,unit="cm"),...)
}
rotx <- function(angle=30) theme(axis.text.x = element_text(angle = angle, hjust = 1))
roty <- function(angle=30) theme(axis.text.y = element_text(angle = angle, hjust = 1))
typef <- function(x) factor(x, c(1,2), c("Pitavastatin alone", "Pitavastatin + CsA"))

