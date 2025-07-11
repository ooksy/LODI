######################################## empirical data check ###################

normalized.count = Y
for(i in 1:n){
  for(j in 1:Jsum){
    normalized.count[i, j] = normalized.count[i, j]/exp(hat.ri[M[j],i])
  }
}

par(mfrow=c(4,4))
for(j in 1:(J[1])){
  hist(log(normalized.count[,j]+0.01), main = paste0('Group ', M[j],' OTU ', j), 
       xlab = 'log normalized count', breaks = 20)
}


burn=1:nsamp
library(ggplot2)
library(latex2exp)
library(circlize)
library(ComplexHeatmap)

emp.factor = log(rowMeans(Y))
empirical.cov = cov(log(Y+0.01) - t(emp.factor)[,rep(1:m, J)])
quantile(empirical.cov[row(empirical.cov)!=col(empirical.cov)], c(.025, .975))
summary(empirical.cov[row(empirical.cov)!=col(empirical.cov)])
col_fun = colorRamp2(c(-30, 0, 30), c("green", "white", "red"))
Heatmap(empirical.cov, name = "Cor", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  at = c(-30, 0, 30),
  #  labels =c("-30", "-20", "-10", "0", "10", "20","30"),
  title = "Cov",
  legend_height = unit(15, "cm")
), column_names_gp = gpar(fontsize = 8.5), row_names_gp = gpar(fontsize = 5))

col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))
Heatmap(cor(empirical.cov), name = "Cor", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  at = c(-1, 0, 1),
  title = "Cor",
  legend_height = unit(15, "cm")
), column_names_gp = gpar(fontsize = 8.5), row_names_gp = gpar(fontsize = 5))



############################ check ri ###############################

burn=1:nsamp
library(ggplot2)
library(latex2exp)
library(circlize)
library(ComplexHeatmap)

par(mfrow=c(1,1))
#### traceplot 
plot(burn, ri.st[1,1,burn],'l', xlab = 'iter', ylab = TeX("$r_{11}$")) 


#### compare ri estimates to empirical estimates of ri 
mi=1
df.ri = matrix(NA, nrow = n, ncol = 4)
df.ri[, 2] = sapply(1:n, function(x) mean(ri.st[mi, x, burn]))
df.ri[, 1] = sapply(1:n, function(x) quantile(ri.st[mi, x, burn], probs = 0.025))
df.ri[, 3] = sapply(1:n, function(x) quantile(ri.st[mi, x,burn], probs = 0.975))
df.ri[,4] = hat.ri[mi,]
pos.mean.ri.1 = df.ri[,2]
df.ri = data.frame(df.ri)
colnames(df.ri) <- c('2.5%', 'Pos.Mean', '97.5%', 'True')
ggplot(df.ri, aes(x = True)) + 
  geom_point(aes(y = Pos.Mean, colour = "Pos.Mean"), size = 2) +
  geom_line(aes(x = True, y = True, colour = "True")) + 
  theme_bw() +
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.025, linewidth =0.5,
                colour = "light grey") +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(TeX("$\\hat{r_{i1}}$"))+
  xlab(TeX("$r^{EM}_{i1}$")) + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text=element_text(size=25),
        text=element_text(size=25))


############################ check alpha_ij traceplot ###################################

burn=1:nsamp
library(ggplot2)
library(latex2exp)

mi = 1
new.ls = array(NA, dim = c(n, J[mi], length(burn)))
for(i in burn){
  new.ls[,,i-burn[1]+1] = alphasij.st[S,J.ls[[mi]],i]
}
df.riplusthetaj.mean = rowMeans(new.ls, dims = 2)
df.riplusthetaj.low = apply(new.ls, 1:2, function(x) quantile(x,0.025))
df.riplusthetaj.up = apply(new.ls, 1:2, function(x) quantile(x,0.975))
df.riplusthetaj.median = apply(new.ls, 1:2, function(x) quantile(x,0.5))
Jcol = which(M==mi)
par(mfrow=c(4,4))
for(l in 1:L.alpha){
  plot(burn, alphasij.st[1,l,],'l', xlab = 'iter')
}

############################ check ri+alpha_ij ###################################

burn=1:nsamp
library(ggplot2)
library(latex2exp)

mi = 1
new.ls = array(NA, dim = c(n, J[mi], length(burn)))
for(i in burn){
  new.ls[,,i-burn[1]+1] = alphasij.st[S,J.ls[[mi]],i] + matrix(ri.st[mi,,i], n, J[mi])
}
df.riplusthetaj.mean = rowMeans(new.ls, dims = 2)
df.riplusthetaj.low = apply(new.ls, 1:2, function(x) quantile(x,0.025))
df.riplusthetaj.up = apply(new.ls, 1:2, function(x) quantile(x,0.975))
df.riplusthetaj.median = apply(new.ls, 1:2, function(x) quantile(x,0.5))
df.riplusthetaj.true =  matrix(NA, nrow = n, ncol = J[mi])
Jcol = which(M==mi)
for(i in 1:n){
  for(j in 1:length(Jcol)){
    df.riplusthetaj.true[i, j] = log(Y[i,j]+0.01)
  }
}

df.riplusthetaj = c(df.riplusthetaj.low)
df.riplusthetaj = cbind(df.riplusthetaj, c(df.riplusthetaj.mean))
#df.riplusthetaj = cbind(df.riplusthetaj, c(df.riplusthetaj.median))
df.riplusthetaj = cbind(df.riplusthetaj, c(df.riplusthetaj.up))
df.riplusthetaj = cbind(df.riplusthetaj, c(df.riplusthetaj.true))
df.riplusthetaj = data.frame(df.riplusthetaj)
colnames(df.riplusthetaj) <- c('2.5%', 'Pos.Mean', '97.5%', 'True')
ggplot(df.riplusthetaj, aes(x = True)) +
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.1, linewidth =0.1,
                colour = "light grey") +
  geom_point(aes(y = Pos.Mean, colour = "Pos.Mean"), size = 1) +
  geom_line(aes(x = True, y = True, colour = "True")) + 
  theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(TeX("$\\widehat{r_{i1}+\\alpha_{i1j}}$"))+
  xlab(TeX("$\\log(Y_{i1j}+0.01)$")) + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text=element_text(size=25),
        text=element_text(size=25))



####################### Pos Median Cov  ###############################

burn=1:nsamp
library(ggplot2)
library(latex2exp)

mar.pos.cov.array = array(NA, dim=c(Jsum, Jsum, length(burn)))
for(i in burn){
  mar.pos.cov.array[,,i-burn[1]+1] = tcrossprod(Lambda.st[,,i]) + diag(rep(sig2.st[,i], J))
  # mar.pos.cov.array[,,i-burn[1]+1] = tcrossprod(Lambda.st[,,i]) + diag(rep(sig2.st[i], J))
}
mar.pos.cov = apply(mar.pos.cov.array, c(1,2), median)

empirical.cov = cov(log(Y+0.01) - t(hat.ri)[,rep(1:m, J)])
diff.cov = mar.pos.cov - empirical.cov
df.diff.cov = data.frame(x = diff.cov[upper.tri(diff.cov, diag = F)])
ggplot(df.diff.cov, aes(x=x)) + geom_histogram(color="black", fill="white", size=1.1) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab('Density')+
  #xlim(-2, 2)+
  xlab(TeX("$\\hat{\\Sigma_{jj'}}-\\Sigma_{jj'}^{EM}$")) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=25),
        text=element_text(size=35)) 
#ggsave(paste0("./", seed, '-11.png'))

df.cov = data.frame(x = mar.pos.cov[upper.tri(mar.pos.cov, diag = F)], y = empirical.cov[upper.tri(empirical.cov)])
colnames(df.cov) = c('pos', 'true')
ggplot(df.cov, aes(x = true)) + 
  geom_point(aes(y = pos, colour = "pos"), size = .5, position = position_jitter(w = .01, h = .01)) +
  geom_line(aes(x = true, y = true, colour = "true")) + 
  theme_bw() + 
  theme(legend.position="none") +
  #ylim(-10,10)+
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(TeX("$\\hat{\\Sigma}_{jj'}$"))+
  xlab(TeX("$\\Sigma^{EM}_{jj'}$")) + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text=element_text(size=25),
        text=element_text(size=35)) + 
  theme(axis.title.y = element_text(angle = 0, vjust =0.5))
#ggsave(paste0("./", seed, '-12.png'))


## heatmap of cov, lower is empirical, upper is pso estimate
poscov.vs.trcov = empirical.cov
poscov.vs.trcov[upper.tri(poscov.vs.trcov)] = mar.pos.cov[upper.tri(mar.pos.cov)]
library(circlize)
col_fun = colorRamp2(c(-25, 0, 25), c("green", "white", "red"))
library(ComplexHeatmap)
quantile(poscov.vs.trcov[row(poscov.vs.trcov)!=col(poscov.vs.trcov)], c(.025, .975))
summary(poscov.vs.trcov[row(poscov.vs.trcov)!=col(poscov.vs.trcov)])
Heatmap(poscov.vs.trcov, name = "Cor", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  at = c(-25,-20, -10, 0, 10, 20, 25),
  #  labels =c("-30", "-20", "-10", "0", "10", "20","30"),
  title = "Cov",
  legend_height = unit(15, "cm")
), column_names_gp = gpar(fontsize = 8.5), row_names_gp = gpar(fontsize = 5))


###################### Pos Median Cor ###########################

burn=1:nsamp
library(ggplot2)
library(latex2exp)
mar.pos.cor.array = array(NA, dim=c(Jsum, Jsum, length(burn)))
for(i in burn){
  mar.pos.cor.array[,,i] = cov2cor(tcrossprod(Lambda.st[,,i]) + diag(rep(sig2.st[,i], J)))
}
mar.pos.cor = apply(mar.pos.cor.array, c(1,2), median)
mar.pos.cor.upper = apply(mar.pos.cor.array, c(1,2), function(x){
  quantile(x, 0.975)
})
mar.pos.cor.lower = apply(mar.pos.cor.array, c(1,2), function(x){
  quantile(x, 0.025)
})

true.cor = cov2cor(empirical.cov)
df.diff.cor = data.frame(x = mar.pos.cor[upper.tri(mar.pos.cor, diag = T)] - (true.cor)[upper.tri(true.cor, diag = T)])
RMSE <- function(x) sqrt(mean(x^2))
RMSE(df.diff.cor$x)

ggplot(df.diff.cor, aes(x=x)) + geom_histogram(color="black", fill="white", size=1.1) + theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab('Density')+
  xlim(c(-1, 1))+
  xlab(TeX("$\\hat{\\rho}_{jj'}-\\rho_{jj'}^{tr}$")) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=25),
        text=element_text(size=35)) 

df.cor = data.frame(x = mar.pos.cor[upper.tri(mar.pos.cor, diag = F)], y = true.cor[upper.tri(true.cor, diag = F)])
colnames(df.cor) = c('pos', 'true')
ggplot(df.cor, aes(x = true)) + 
  geom_point(aes(y = pos, colour = "pos"), size = .5, position = position_jitter(w = 0.00015, h = 0)) +
  geom_line(aes(x = true, y = true, colour = "true")) + 
  theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(TeX("$\\hat{\\rho}_{jj'}$"))+
  xlab(TeX("$\\rho^{EM}_{jj'}$")) + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text=element_text(size=25),
        text=element_text(size=35)) + 
  theme(axis.title.y = element_text(angle = 0, vjust =0.5))


## heatmap of cor, lower is empirical, upper is pso estimate
poscor.vs.trcor = true.cor
poscor.vs.trcor[upper.tri(poscor.vs.trcor)] = mar.pos.cor[upper.tri(mar.pos.cor)]
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))
library(ComplexHeatmap)
#png(paste0("./", seed, '-16.png'), width = 5383/10, height = 4479/10, units = "px")
Heatmap(poscor.vs.trcor, name = "Cor", col = col_fun, cluster_rows = F,  cluster_columns = F, heatmap_legend_param = list(
  title = "Cor",
  legend_height = unit(15, "cm")
), column_names_gp = gpar(fontsize = 8.5), row_names_gp = gpar(fontsize = 5))


