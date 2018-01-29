rm (list=ls())


#Idea: 'm' SBM graphs with same block matrix, B, with entries [p11, p12; p21, p22] 
#Simulate var(p_hat-p) per Minh's theorem in order to generate p_hat.
#Note: we do not generate SBM graphs here.

#set up graph regression with p_hat 
#Using Minh's theorem construct an estimate for var(p_hat-p) 
#use the estimate to adjust measurement error 

pi_k= 0.5
pi_l = 0.5 

graph_size = 50

beta_0 = 1
beta_1=1
beta_2=1
beta_3 =1
beta = matrix(c(beta_0, beta_1, beta_2, beta_3))

beta_true_se = c()
beta_naive_se = c()
beta_adj_se = c()


m = 100 # number of graphs 
p11_true = runif(n=m,  0.5, 0.6) 
p12_true = runif(n=m,  0.5, 0.6) 
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
z=1
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
par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b0_true, b0_naive,  b0_adj,notch=TRUE, 
        main=bquote(paste("b0 Estimate", "\n graph_size:", graph_size , "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = beta_0, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b1_true, b1_naive,  b1_adj,notch=TRUE, 
        main=bquote(paste("b1 Estimate", "\n graph_size:", graph_size , "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b2_true, b2_naive,  b2_adj, notch=TRUE, 
        main=bquote(paste("b2 Estimate", "\n graph_size:", graph_size , "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted "))
abline(h = b2_true, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b3_true, b3_naive,  b3_adj,notch=TRUE, 
        main=bquote(paste("b3 Estimate", "\n graph_size:", graph_size , "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = b3_true, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b0_true_se, b0_naive_se,  b0_adj_se,notch=TRUE, 
        main=bquote(paste("b0 Square Error", "\n graph_size:", graph_size , "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = 0, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b1_true_se, b1_naive_se,  b1_adj_se,notch=TRUE, 
        main=bquote(paste("b1 Square Error", "\n graph_size:", graph_size , "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = 0, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b2_true_se, b2_naive_se,  b2_adj_se, notch=TRUE, 
        main=bquote(paste("b2 Square Error", "\n graph_size:", graph_size , "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = 0, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b3_true_se, b3_naive_se,  b3_adj_se,notch=TRUE, 
        main=bquote(paste("b3 Square Error", "\n graph_size:", graph_size , "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = 0, col="red", lty=2, lwd=0.5)

########
########
########


# par(mar=c(7.1,4.1,4.1,2.1))
# boxplot(beta_all_true_se, beta_all_naive_se,  beta_all_adj_se,notch=TRUE, 
#         main=bquote(paste("beta_all Square Error", "\n graph_size:", graph_size , "\n mc_runs:", mc_runs, "\n m:", m)), 
#         cex.main=0.5,
#         las=2, names=c("True", "Naive", "Adjusted"))
# abline(h = b0_true, col="red", lty=2, lwd=0.5)


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










###########
###########
###########
###########
###########
###########
###########
###########
###########
par(mar=c(7.1,4.1,4.1,2.1))
 boxplot(b1_true, b1_naive,  b1_adj[abs(b1_adj)<mean(b1_adj)+3*sd(b1_adj)],notch=TRUE, 
         main=bquote(paste("b1 Estimate", "\n mc_runs:", mc_runs, "\n m:", m)), 
         cex.main=0.5,
         las=2, names=c("True", "Naive", "|Adjusted| < \n 3sd(Adjusted)"))
 abline(h = b1_true, col="red", lty=2, lwd=0.5)
 
 par(mar=c(7.1,4.1,4.1,2.1))
 boxplot(b2_true, b2_naive,  b2_adj[abs(b2_adj)<mean(b2_adj)+3*sd(b2_adj)], notch=TRUE, 
         main=bquote(paste("b2 Estimate", "\n mc_runs:", mc_runs, "\n m:", m)), 
         cex.main=0.5,
         las=2, names=c("True", "Naive", "|Adjusted| < \n 3sd(Adjusted)"))
 abline(h = b2_true, col="red", lty=2, lwd=0.5)
 
 par(mar=c(7.1,4.1,4.1,2.1))
 boxplot(b3_true, b3_naive,  b3_adj[abs(b3_adj)<mean(b3_adj)+3*sd(b3_adj)],notch=TRUE, 
         main=bquote(paste("b3 Estimate", "\n mc_runs:", mc_runs, "\n m:", m)), 
         cex.main=0.5,
         las=2, names=c("True", "Naive", "|Adjusted| < \n 3sd(Adjusted)"))
 abline(h = b3_true, col="red", lty=2, lwd=0.5)
 
 par(mar=c(7.1,4.1,4.1,2.1))
 boxplot(b1_true_se, b1_naive_se,  b1_adj_se[abs(b1_adj_se)<mean(b1_adj_se)+3*sd(b1_adj_se)],notch=TRUE, 
         main=bquote(paste("b1 Square Error", "\n mc_runs:", mc_runs, "\n m:", m)), 
         cex.main=0.5,
         las=2, names=c("True", "Naive", "|Adjusted| < \n 3sd(Adjusted)"))
 abline(h = 0, col="red", lty=2, lwd=0.5)
 
 par(mar=c(7.1,4.1,4.1,2.1))
 boxplot(b2_true_se, b2_naive_se,  b2_adj_se[abs(b2_adj_se)<mean(b2_adj_se)+3*sd(b2_adj_se)], notch=TRUE, 
         main=bquote(paste("b2 Square Error", "\n mc_runs:", mc_runs, "\n m:", m)), 
         cex.main=0.5,
         las=2, names=c("True", "Naive", "|Adjusted| < \n 3sd(Adjusted)"))
 abline(h = 0, col="red", lty=2, lwd=0.5)
 
 par(mar=c(7.1,4.1,4.1,2.1))
 boxplot(b3_true_se, b3_naive_se,  b3_adj_se[abs(b3_adj_se)<mean(b3_adj_se)+3*sd(b1_adj_se)],notch=TRUE, 
         main=bquote(paste("b3 Square Error", "\n mc_runs:", mc_runs, "\n m:", m)), 
         cex.main=0.5,
         las=2, names=c("True", "Naive", "Adjusted < \n 3sd(Adjusted)"))
 abline(h = 0, col="red", lty=2, lwd=0.5)
 
 det(crossprod(P_true))
 
 mean(b3_true) 
 mean(b2_true) 
 mean(b1_true) 
 
 mean(b3_naive) 
 mean(b2_naive) 
 mean(b1_naive) 
 
 mean(b1_adj)
 mean(b2_adj)
 mean(b3_adj)
 
 mean(b1_adj_ideal)
 mean(b2_adj_ideal)
 mean(b3_adj_ideal)
 
 mean(b1_adj_cheat)
 mean(b2_adj_cheat)
 mean(b3_adj_cheat)
 
 b1_adj_trim = b1_adj[abs(b1_adj)<3*sd(b1_adj)] 
 b2_adj_trim = b2_adj[abs(b2_adj)<3*sd(b2_adj)] 
 b3_adj_trim = b3_adj[abs(b3_adj)<3*sd(b3_adj)] 
 
 mean( b1_adj_trim)
 mean(b2_adj_trim)
 mean(b3_adj_trim)
 
 b1_naive_trim = b1_naive[abs(b1_adj)<3*sd(b1_adj)] 
 b2_naive_trim = b2_naive[abs(b2_adj)<3*sd(b2_adj)] 
 b3_naive_trim = b3_naive[abs(b3_adj)<3*sd(b3_adj)] 
 
 b1_true_trim = b1_true[abs(b1_adj)<3*sd(b1_adj)] 
 b2_true_trim = b2_true[abs(b2_adj)<3*sd(b2_adj)] 
 b3_true_trim = b3_true[abs(b3_adj)<3*sd(b3_adj)] 
 
 t.test(b1_true)$conf.int
 t.test(b1_naive)$conf.int
 t.test(b1_adj)$conf.int
 t.test(b1_adj_trim)$conf.int
 
t.test(b1_true, b1_naive, paired=TRUE, alternative="greater" ) 
t.test(b1_adj, b1_naive, paired=TRUE, alternative="greater" ) 
t.test(b1_true, b1_adj, paired=TRUE, alternative="greater" ) 
t.test(b1_true_trim, b1_adj_trim, paired=TRUE, alternative="greater" ) 



b1_adj_se_trim = b1_adj_se[abs(b1_adj_se)<mean(b1_adj_se)+3*sd(b1_adj_se)] 
b2_adj_se_trim = b2_adj_se[abs(b2_adj_se)<mean(b2_adj_se)+3*sd(b2_adj_se)] 
b3_adj_se_trim = b3_adj_se[abs(b3_adj_se)<mean(b3_adj_se)+3*sd(b3_adj_se)] 

b1_true_se_trim = b1_true_se[abs(b1_adj_se)<mean(b1_adj_se)+3*sd(b1_adj_se)] 
b2_true_se_trim = b2_true_se[abs(b2_adj_se)<mean(b2_adj_se)+3*sd(b2_adj_se)] 
b3_true_se_trim = b3_true_se[abs(b3_adj_se)<mean(b3_adj_se)+3*sd(b3_adj_se)] 

b1_naive_se_trim = b1_naive_se[abs(b1_adj_se)<mean(b1_adj_se)+3*sd(b1_adj_se)] 
b2_naive_se_trim = b2_naive_se[abs(b2_adj_se)<mean(b2_adj_se)+3*sd(b2_adj_se)] 
b3_naive_se_trim = b3_naive_se[abs(b3_adj_se)<mean(b3_adj_se)+3*sd(b3_adj_se)] 

N = sum((b1_adj_se - b1_naive_se) !=0)
k = sum(b1_adj_se < b1_naive_se)
1 - pbinom(k, size=N, prob=0.5)

N = sum(b1_adj_se_trim - b1_naive_se_trim !=0)
k = sum(b1_adj_se_trim < b1_naive_se_trim)
1 - pbinom(k, size=N, prob=0.5)

t.test(b2_true)$conf.int
t.test(b2_naive)$conf.int
t.test(b2_adj)$conf.int
t.test(b2_adj_trim)$conf.int

t.test(b1_adj, b1_naive, paired=TRUE, alternative="greater" ) 

t.test(b2_true, b2_naive, paired=TRUE, alternative="greater" ) 
t.test(b2_adj, b2_naive, paired=TRUE, alternative="greater" ) 
t.test(b2_adj_trim, b2_naive_trim, paired=TRUE, alternative="greater" ) 
t.test(b2_true, b2_adj, paired=TRUE, alternative="greater" ) 
t.test(b2_true_trim, b2_adj_trim, paired=TRUE, alternative="greater" ) 

N = sum((b2_adj_se - b2_naive_se) !=0)
k = sum(b2_adj_se < b2_naive_se)
1 - pbinom(k, size=N, prob=0.5)

N = sum((b2_adj_se_trim - b2_naive_se_trim) !=0)
k = sum(b2_adj_se_trim < b2_naive_se_trim)
1 - pbinom(k, size=N, prob=0.5)

t.test(b3_true)$conf.int
t.test(b3_naive)$conf.int
t.test(b3_adj)$conf.int
t.test(b3_adj_trim)$conf.int

t.test(b3_true, b3_naive, paired=TRUE, alternative="greater" ) 
t.test(b3_adj, b3_naive, paired=TRUE, alternative="greater" ) 
t.test(b3_adj_trim, b3_naive_trim, paired=TRUE, alternative="greater" ) 

t.test(b3_true, b3_adj, paired=TRUE, alternative="greater" ) 
t.test(b3_true_trim, b3_adj_trim, paired=TRUE, alternative="greater" )

#sign test 
N = sum((b3_adj_se - b3_naive_se) !=0)
k = sum(b3_adj_se < b3_naive_se)
1 - pbinom(k, size=N, prob=0.5)

N = sum((b3_adj_se_trim - b3_naive_se_trim) !=0)
k = sum(b3_adj_se_trim < b3_naive_se_trim)
1 - pbinom(k, size=N, prob=0.5)


###########
###########
###########



par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b1_true, b1_naive,  b1_adj[abs(b1_adj)<2],  b1_adj_ideal[abs(b1_adj_ideal)<2], notch=TRUE, 
        main=bquote(paste("b1 Estimate", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "|Adjusted| < 2", "Adj_ideal"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b2_true, b2_naive,   b2_adj[abs(b2_adj)<2],  b2_adj_ideal[abs(b2_adj_ideal)<2], notch=TRUE, 
        main=bquote(paste("b2 Estimate", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "|Adjusted| < 2", "Adj_ideal"))
abline(h = b2_true, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b3_true, b3_naive,  b3_adj[abs(b3_adj)<2], b3_adj_ideal[abs(b3_adj_ideal)<2],notch=TRUE, 
        main=bquote(paste("b3 Estimate", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "|Adjusted| < 2", "Adj_ideal"))
abline(h = b3_true, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b1_true_se, b1_naive_se, b1_adj_se[ abs(b1_adj_se)<2],  b1_adj_ideal_se[ abs(b1_adj_ideal_se)<2] ,notch=TRUE, 
        main=bquote(paste("b1 Square Error", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "|Adjusted| < 2", "Adj_ideal"))
abline(h = 0, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b2_true_se, b2_naive_se, b2_adj_se[  abs(b2_adj_se)<2], b2_adj_ideal_se[ abs(b2_adj_ideal_se)<2], notch=TRUE, 
        main=bquote(paste("b2 Square Error", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "|Adjusted| < 2", "Adj_ideal"))
abline(h = 0, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b3_true_se, b3_naive_se,  b3_adj_se[ abs(b3_adj_se)<2], b3_adj_ideal_se[ abs(b3_adj_ideal_se)<2],notch=TRUE, 
        main=bquote(paste("b3 Square Error", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "|Adjusted| < 2", "Adj_ideal"))
abline(h = 0, col="red", lty=2, lwd=0.5)


###########
###########
###########

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b1_true, b1_naive,  b1_adj, b1_adj_ideal ,notch=TRUE, 
        main=bquote(paste("b1 Estimate", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        ylim=c(-10, 10),
        las=2, names=c("True", "Naive", "Adjusted", "Adj_ideal"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b2_true, b2_naive,  b2_adj, b2_adj_ideal, notch=TRUE, 
        main=bquote(paste("b2 Estimate", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted" , "Adj_ideal"))
abline(h = b2_true, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b3_true, b3_naive,  b3_adj, b3_adj_ideal, notch=TRUE, 
        main=bquote(paste("b3 Estimate", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted", "Adj_ideal"))
abline(h = b3_true, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b1_true_se, b1_naive_se,  b1_adj_se, b1_adj_ideal_se, notch=TRUE, 
        main=bquote(paste("b1 Square Error", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted", "Adj_ideal"))
abline(h = 0, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b2_true_se, b2_naive_se,  b2_adj_se, b2_adj_ideal_se, notch=TRUE, 
        main=bquote(paste("b2 Square Error", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted", "Adj_ideal"))
abline(h = 0, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b3_true_se, b3_naive_se,  b3_adj_se, b3_adj_ideal_se, notch=TRUE, 
        main=bquote(paste("b3 Square Error", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        ylim = c(0, 10),
        las=2, names=c("True", "Naive", "Adjusted", "Adj_ideal"))
abline(h = 0, col="red", lty=2, lwd=0.5)


######
######
######



par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b1_true, b1_naive,  b1_adj,  notch=TRUE, 
        main=bquote(paste("b1 Estimate","\n graph_size:", n , "\n mc_runs:",length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b2_true, b2_naive,  b2_adj,  notch=TRUE, 
        main=bquote(paste("b2 Estimate", "\n graph_size:", n ,"\n mc_runs:", length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted" ))
abline(h = b2_true, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b3_true, b3_naive,  b3_adj,  notch=TRUE, 
        main=bquote(paste("b3 Estimate", "\n graph_size:", n ,"\n mc_runs:", length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = b3_true, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b1_true_se, b1_naive_se,  b1_adj_se, notch=TRUE, 
        main=bquote(paste("b1 Square Error", "\n graph_size:", n ,"\n mc_runs:", length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = 0, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b2_true_se, b2_naive_se,  b2_adj_se, notch=TRUE, 
        main=bquote(paste("b2 Square Error", "\n graph_size:", n ,"\n mc_runs:", length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = 0, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b3_true_se, b3_naive_se,  b3_adj_se, notch=TRUE, 
        main=bquote(paste("b3 Square Error", "\n graph_size:", n ,"\n mc_runs:", length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = 0, col="red", lty=2, lwd=0.5)


###########
###########
###########


par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b1_true, b1_naive,  b1_adj, b1_adj_cheat ,notch=TRUE, 
        main=bquote(paste("b1 Estimate", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted", "Adj_cheat"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b2_true, b2_naive,  b2_adj, b2_adj_cheat, notch=TRUE, 
        main=bquote(paste("b2 Estimate", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted" , "Adj_cheat"))
abline(h = b2_true, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b3_true, b3_naive,  b3_adj, b3_adj_cheat, notch=TRUE, 
        main=bquote(paste("b3 Estimate", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted", "Adj_cheat"))
abline(h = b3_true, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b1_true_se, b1_naive_se,  b1_adj_se, b1_adj_cheat_se, notch=TRUE, 
        main=bquote(paste("b1 Square Error", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted", "Adj_cheat"))
abline(h = 0, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b2_true_se, b2_naive_se,  b2_adj_se, b2_adj_cheat_se, notch=TRUE, 
        main=bquote(paste("b2 Square Error", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted", "Adj_cheat"))
abline(h = 0, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b3_true_se, b3_naive_se,  b3_adj_se, b3_adj_cheat_se, notch=TRUE, 
        main=bquote(paste("b3 Square Error", "\n mc_runs:", mc_runs, "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted", "Adj_cheat"))
abline(h = 0, col="red", lty=2, lwd=0.5)


mean(b1_adj[abs(b1_adj)<10])
mean(b1_adj_cheat[abs(b1_adj_cheat)<10])
mean(b1_adj_ideal[abs(b1_adj_ideal)<10])
mean(b1_naive[abs(b1_naive)<10])
mean(b1_true)

mean(b2_adj[abs(b2_adj)<10])
mean(b2_adj_cheat[abs(b2_adj_cheat)<10])
mean(b2_adj_ideal[abs(b2_adj_ideal)<10])
mean(b2_naive[abs(b2_naive)<10])

mean(b3_adj[abs(b3_adj)<10])
mean(b3_adj_cheat[abs(b3_adj_cheat)<10])
mean(b3_adj_ideal[abs(b3_adj_ideal)<10])
mean(b3_naive[abs(b3_naive)<10])
