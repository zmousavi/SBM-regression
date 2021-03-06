rm (list=ls())
require("Matrix")

#Idea: 'm' SBM graphs with same block matrix. 
#Simulate block matrix entries using Minh's var(p_hat-p) theorem.
#Note: we do not generate SBM graphs here.
#compute an estimate for var(p_hat-p)
#use the estimate to adjust measurement error 

pi_k= 0.5
pi_l = 0.5 

graph_size = 100

beta_0 = 1
beta_1=1
beta_2=1
beta_3 =1
beta = matrix(c(beta_0, beta_1, beta_2, beta_3))

beta_true_se = c()
beta_naive_se = c()
beta_adj_se = c()


m = 100 # number of graphs 
set.seed(83)
p11_true = runif(n=m,  0.5, 0.6) 
set.seed(79)
p12_true = runif(n=m,  0.5, 0.6) 
set.seed(439)
p22_true = runif(n=m, 0.5, 0.6) 

P_true = cbind(p11_true, p12_true, p22_true) #B matrix entries 

#theoeretical var(p_hat-p) according to minh
var_p11_minh =   (2*(1-p11_true)*p11_true)/(pi_k^2*(graph_size)^2)
var_p12_minh =  (1*(1-p12_true)*p12_true)/(pi_k*pi_l*(graph_size)^2)
var_p12_minh =  var_p22_minh * var_p11_minh + var_p22_minh
var_p22_minh =  (2*(1-p22_true)*p22_true)/(pi_l^2*(graph_size)^2)
U = cbind(var_p11_minh, var_p12_minh, var_p22_minh)

beta_true = list()
beta_naive = list()
beta_adj=list()
beta_true_se = list()
beta_naive_se = list()
beta_adj_se=list()



beta_all_true = c()
beta_all_naive = c()
beta_all_adj=c()
beta_all_true_se =c()
beta_all_naive_se =c()
beta_all_adj_se =c()


mc_runs = 100
for (z in 1:mc_runs){

  if (z%%100 ==0){
    print (z)
  } 
  
  
P_naive = cbind(rnorm(m, mean=P_true[,1], sd=sqrt(var_p11_minh)),
                rnorm(m, mean=P_true[,2], sd=sqrt(var_p12_minh)),
                rnorm(m, mean=P_true[,3], sd=sqrt(var_p22_minh)))
  
X_true=P_true
X_naive = P_naive
y = beta_0 + X_true%*%beta[2:4] + rnorm(length(m), sd = sqrt(1e-2*mean(1e-4))) 


#estimated var(p_hat-p)
var_p11_minh_hat =   (2*(1-P_naive[,1])*P_naive[,1])/(pi_k^2*(graph_size)^2)
var_p12_minh_hat =  (1*(1-P_naive[,2])*P_naive[,2])/(pi_k*pi_l*(graph_size)^2)
var_p22_minh_hat =  (2*(1-P_naive[,3])*P_naive[,3])/(pi_l^2*(graph_size)^2)
U_hat = cbind(var_p11_minh_hat, var_p12_minh_hat, var_p22_minh_hat)
sd_U = sqrt(U_hat)

Sigma_U = diag(diag(crossprod(sd_U)))
solve(t(X_true)%*%X_true)%*%t(X_true)%*% y
X_true_big = cbind(rep(1, m), X_true)
beta_true[[z]] = solve(t(X_true_big)%*%X_true_big)%*%t(X_true_big)%*% y

X_naive_big = cbind(rep(1, m), X_naive)

Sigma_U_big = bdiag (0, Sigma_U)
beta_naive[[z]] = solve(t(X_naive_big)%*%X_naive_big)%*%t(X_naive_big)%*% y
beta_adj[[z]] = solve(t(X_naive_big)%*%X_naive_big-Sigma_U_big )%*%t(X_naive_big)%*% y


beta_true_se[[z]] = (beta_true[[z]]-beta)^2
beta_naive_se[[z]]  = (beta_naive[[z]]-beta)^2
beta_adj_se[[z]] = (beta_adj[[z]]-beta)^2

}

b0_true = c()
b1_true = c()
b2_true = c()
b3_true = c()

b0_naive = c()
b1_naive = c()
b2_naive = c()
b3_naive = c()

b0_adj = c()
b1_adj = c()
b2_adj = c()
b3_adj = c()

b0_true_se = c()
b1_true_se = c()
b2_true_se = c()
b3_true_se = c()

b0_naive_se = c()
b1_naive_se = c()
b2_naive_se = c()
b3_naive_se = c()

b0_adj_se = c()
b1_adj_se = c()
b2_adj_se = c()
b3_adj_se = c()



for (z in 1:length(beta_true)){
  
  if (z%%100 ==0){
    print (z)
  } 
 
b0_true[z] =  beta_true[[z]][1]
b1_true[z] =  beta_true[[z]][2]
b2_true[z] = beta_true[[z]][3]
b3_true[z] = beta_true[[z]][4]

b0_true_se[z] =  beta_true_se[[z]][1]
b1_true_se[z] =  beta_true_se[[z]][2]
b2_true_se[z] = beta_true_se[[z]][3]
b3_true_se[z] = beta_true_se[[z]][4]

b0_naive[z] = beta_naive[[z]][1]
b1_naive[z] = beta_naive[[z]][2]
b2_naive[z] = beta_naive[[z]][3]
b3_naive[z] = beta_naive[[z]][4]

b0_naive_se[z] = beta_naive_se[[z]][1]
b1_naive_se[z] = beta_naive_se[[z]][2]
b2_naive_se[z] = beta_naive_se[[z]][3]
b3_naive_se[z] = beta_naive_se[[z]][4]

b0_adj_se[z] = beta_adj_se[[z]][1]
b1_adj_se[z] = beta_adj_se[[z]][2]
b2_adj_se[z] = beta_adj_se[[z]][3]
b3_adj_se[z] = beta_adj_se[[z]][4]

b0_adj[z] = beta_adj[[z]][1]
b1_adj[z] = beta_adj[[z]][2]
b2_adj[z] = beta_adj[[z]][3]
b3_adj[z] = beta_adj[[z]][4]


}

###########
###########
###########

setwd("~/Desktop/redata/Feb6")
dev.new()
pdf("b0_estimate_sbm_sim.pdf", 7, 5)
#par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b0_true, b0_naive,  b0_adj,  notch=TRUE, 
        #main=bquote(paste("b1 Estimate","\n graph_size:", n , "\n mc_runs:",length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)
dev.off()

dev.new()
pdf("b11_estimate_sim.pdf", 7, 5)
boxplot(b1_true, b1_naive,  b1_adj,  notch=TRUE, 
        #main=bquote(paste("b1 Estimate","\n graph_size:", n , "\n mc_runs:",length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)
dev.off()

dev.new()
pdf("b12_estimate_sim.pdf", 7, 5)
boxplot(b2_true, b2_naive,  b2_adj,  notch=TRUE, 
        main=bquote(paste("b2 Estimate", "\n graph_size:", n ,"\n mc_runs:", length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted" ))
abline(h = b2_true, col="red", lty=2, lwd=0.5)
dev.off()

dev.new()
pdf("b22_estimate_sim.pdf", 7, 5)
boxplot(b3_true, b3_naive,  b3_adj,  notch=TRUE, 
        main=bquote(paste("b3 Estimate", "\n graph_size:", n ,"\n mc_runs:", length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = b3_true, col="red", lty=2, lwd=0.5)
dev.off()

dev.new()
pdf("b0_MSE_sbm_sim.pdf", 7, 5)
boxplot(b0_true_se, b0_naive_se,  b0_adj_se, notch=TRUE, 
        #   main=bquote(paste("b1 Square Error", "\n graph_size:", n ,"\n mc_runs:", length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = 0, col="red", lty=2, lwd=0.5)
dev.off()

dev.new()
pdf("b11_MSE_sim.pdf", 7, 5)
boxplot(b1_true_se, b1_naive_se,  b1_adj_se, notch=TRUE, 
        #   main=bquote(paste("b1 Square Error", "\n graph_size:", n ,"\n mc_runs:", length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = 0, col="red", lty=2, lwd=0.5)
dev.off()

dev.new()
pdf("b12_MSE_sim.pdf", 7, 5)
boxplot(b2_true_se, b2_naive_se,  b2_adj_se, notch=TRUE, 
       # main=bquote(paste("b2 Square Error", "\n graph_size:", n ,"\n mc_runs:", length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = 0, col="red", lty=2, lwd=0.5)
dev.off()

dev.new()
pdf("b22_MSE_sim.pdf", 7, 5)
par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b3_true_se, b3_naive_se,  b3_adj_se, notch=TRUE, 
      #  main=bquote(paste("b3 Square Error", "\n graph_size:", n ,"\n mc_runs:", length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = 0, col="red", lty=2, lwd=0.5)
dev.off()




#CI 
t.test(b0_true)$conf.int
t.test(b0_naive)$conf.int
t.test(b0_adj)$conf.int

t.test(b1_true)$conf.int
t.test(b1_naive)$conf.int
t.test(b1_adj)$conf.int

t.test(b2_true)$conf.int
t.test(b2_naive)$conf.int
t.test(b2_adj)$conf.int

t.test(b3_true)$conf.int
t.test(b3_naive)$conf.int
t.test(b3_adj)$conf.int

#SE (sign test )
N = sum((b0_adj_se - b0_naive_se) !=0)
k = sum(b0_adj_se < b0_naive_se)
1 - pbinom(k, size=N, prob=0.5)

N = sum((b1_adj_se - b1_naive_se) !=0)
k = sum(b1_adj_se < b1_naive_se)
1 - pbinom(k, size=N, prob=0.5)

N = sum((b2_adj_se - b2_naive_se) !=0)
k = sum(b2_adj_se < b2_naive_se)
1 - pbinom(k, size=N, prob=0.5)

N = sum((b3_adj_se - b3_naive_se) !=0)
k = sum(b3_adj_se < b3_naive_se)
1 - pbinom(k, size=N, prob=0.5)






