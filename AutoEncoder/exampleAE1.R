

#import the watermelon dataset
wm=read.csv(file = "C:/Users/Bruin/Desktop/GS Academia/SEM 2/STAT 519/PROJ 1/watermelon/simulated datasets wawamelon/WaterMelonMaster.csv")
wm$species.names=as.factor(wm$species.names)
Y=wm[, 2]
X=wm[,-c(1)]
X$Flesh.color=as.factor(X$Flesh.color)
mod=lm(Lglb~., data = X)

model.mat=model.matrix(mod)

beta=solve(t(model.mat)%*%model.mat)%*%t(model.mat)%*%Y


#using SGD

beta_hat <- matrix(0.1, nrow=ncol(model.mat))
lr=0.0001
N=200
# Repeat below for N-iterations
for (j in 1:N){
  
  # Calculate the cost/error (y_guess - y_truth)
  residual <- (model.mat %*% beta_hat) - Y
  # Calculate the gradient at that point
  delta <- (t(model.mat) %*% residual) * (1/nrow(model.mat))
  # Move guess in opposite direction of gradient
  beta_hat <- beta_hat - (lr*delta)
  print(cbind(beta_hat,delta))
  
}
