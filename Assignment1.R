setwd("D:/Assignment 1")

data = read.csv("gene_data.csv", header = F, fileEncoding = "UTF-8-BOM")

#Task 1

Y1=data[,1]
Y2=data[,2]
Y3=data[,3]
Y4=data[,4]
Y5=data[,5]
Y6=data[,6]
print(Y1)
print(Y2)
print(Y3)
print(Y4)
print(Y5)
print(Y6)

#1.1: Time series plot
plot(Y1,Y2,type='l', xlab ='Time-Series',ylab = 'Data Frame', main ='Time Series Plot of each Gene', col='red' )
plot(Y1,Y3,type='l', xlab ='Time-Series',ylab = 'Data Frame', main ='Time Series Plot of each Gene for V1 and V3', col='green')
plot(Y1,Y4,type='l', xlab ='Time-Series',ylab = 'Data Frame', main ='Time Series Plot of each Gene for V1 and V4', col='blue')
plot(Y1,Y5,type='l', xlab ='Time-Series',ylab = 'Data Frame', main ='Time Series Plot of each Gene for V1 and V5', col='yellow')
plot(Y1,Y6,type='l', xlab ='Time-Series',ylab = 'Data Frame', main ='Time Series Plot of each Gene for V1 and V6', col='orange')

plot(Y1,Y2,type='l', xlab ='Time-Series',ylab = 'Data Frame', main ='Time Series Plot of each Gene', col='red' )
lines(Y1,Y3,type='l', xlab ='Time-Series',ylab = 'Data Frame',  col='green')
lines(Y1,Y4,type='l', xlab ='Time-Series',ylab = 'Data Frame', col='blue')
lines(Y1,Y5,type='l', xlab ='Time-Series',ylab = 'Data Frame', col='yellow')
lines(Y1,Y6,type='l', xlab ='Time-Series',ylab = 'Data Frame', col='orange')
legend('bottomright', c("Y2", "Y3", "Y4", "Y5", "Y6"), lty = 1, col = 2:6)

#1.2: Distribution for each gene
hist(Y2, xlab ='Time-Series',ylab = 'Data Frame', main ='Distribution of each Gene in time-series', col='red' )
hist(Y3, xlab ='Time-Series',ylab = 'Data Frame', main ='Distribution of each Gene in time-series', col='green')
hist(Y4, xlab ='Time-Series',ylab = 'Data Frame', main ='Distribution of each Gene in time-series', col='blue')
hist(Y5, xlab ='Time-Series',ylab = 'Data Frame', main ='Distribution of each Gene in time-series', col='yellow')
hist(Y6, xlab ='Time-Series',ylab = 'Data Frame', main ='Distribution of each Gene in time-series', col='orange')

#1.3: Correlation and Scatter plots
cor(Y2,Y3)
cor(Y2,Y4)
cor(Y2,Y5)
cor(Y2,Y6)


cor(Y3,Y4)
cor(Y3,Y5)
cor(Y3,Y6)

cor(Y4,Y5)
cor(Y4,Y6)

cor(Y5,Y6)

plot(data[2:6], main= 'Scatter Plot for each Gene')

#-----------------------------------------------------------------------------------

#Task 2: Regression - modelling the relationship between gene expression 

#model 1

#model 1.2.1

b1 = matrix(1, length(Y1), 1)
X1 = cbind(Y5,Y4^2, b1)
print(X1)
t1Hat = solve(t(X1) %*% X1) %*% t(X1) %*% Y3
print(t1Hat)

#model 1.2.2

y1Hat = X1 %*% t1Hat
print(X1)
error1 = Y3 - y1Hat
print(X1)
RSS1 = sum((Y3-y1Hat)^2)
print(RSS1)

#model 1.2.3

n1 = NROW(Y1)
n1
sigma1=RSS1/(n1-1)
print(X1)
p1=-n1/2*log(2*pi)-n1/2*log(sigma1)-(1/(2*sigma1))*RSS1
print(p1)

#model 1.2.4

k1 = 3
AIC1 = 2*k1-2*p1
BIC1 = k1*log(n1)-2*p1
print(AIC1)
print(BIC1)

#model 1.2.5 

M1 = as.matrix(error1)
qqnorm(M1, main ='Normal Q-Q Plot for model 1' )
qqline(M1)

#-------------------------------------------------------------------------------

#model 2

#model 2.2.1

b2 = matrix(1, length(Y1), 1)
X2 = cbind(Y5,Y4^2,Y6,b2)
t2Hat = solve(t(X2) %*% X2) %*% t(X2) %*% Y3
print(t2Hat)

#model 2.2.2

y2Hat = X2 %*% t2Hat
error2 = Y3 - y2Hat

RSS2 = sum((Y3-y2Hat)^2)
print(RSS2)

#model 2.2.3

n2 = NROW(Y1)
sigma2=RSS2/(n2-1)
p2=-n2/2*log(2*pi)-n2/2*log(sigma2)-(1/(2*sigma2))*RSS2
print(p2)

#model 2.2.4

k2 = 4
AIC2 = 2*k2-2*p2
BIC2 = k2*log(n2)-2*p2
print(AIC2)
print(BIC2)

#model 2.2.5 

M2 = as.matrix(error2)
qqnorm(M2,, main ='Normal Q-Q Plot for model 2')
qqline(M2)

#-------------------------------------------------------------------------------

#model 3

#model 3.2.1

X3 = cbind(Y4,Y5,Y6^3)
t3Hat = solve(t(X3) %*% X3) %*% t(X3) %*% Y3
print(t3Hat)

#model 3.2.2

y3Hat = X3 %*% t3Hat
error3 = Y3 - y3Hat

RSS3 = sum((Y3-y3Hat)^2)
print(RSS3)

#model 3.2.3

n3 = NROW(Y1)
sigma3=RSS3/(n3-1)
p3=-n3/2*log(2*pi)-n3/2*log(sigma3)-(1/(2*sigma3))*RSS3
print(p3)

#model 3.2.4

k3 = 3
AIC3 = 2*k3-2*p3
BIC3 = k3*log(n3)-2*p3
print(AIC3)
print(BIC3)

#model 3.2.5 

M3 = as.matrix(error3)
qqnorm(M3, main ='Normal Q-Q Plot for model 3')
qqline(M3)

#-----------------------------------------------------------------------------

#model 4

#model 4.2.1

b4 = matrix(1, length(Y1), 1)
X4 = cbind(Y5,Y4^2,Y6^3,b4)
t4Hat = solve(t(X4) %*% X4) %*% t(X4) %*% Y3
print(t4Hat)

#model 4.2.2

y4Hat = X4 %*% t4Hat
error4 = Y3 - y4Hat

RSS4 = sum((Y3-y4Hat)^2)
print(RSS4)

#model 4.2.3

n4 = NROW(Y1)
sigma4=RSS4/(n4-1)
p4=-n4/2*log(2*pi)-n4/2*log(sigma4)-(1/(2*sigma4))*RSS4
print(p4)

#model 4.2.4

k4 = 4
AIC4 = 2*k4-2*p4
BIC4 = k4*log(n4)-2*p4
print(AIC4)
print(BIC4)

#model 4.2.5 

M4 = as.matrix(error4)
qqnorm(M4, main ='Normal Q-Q Plot for model 4')
qqline(M4)

#-----------------------------------------------------------------------------

#model 5

#model 5.2.1

b5 = matrix(1, length(Y1), 1)
X5 = cbind(Y5,Y2^2,Y4^2,b5)
t5Hat = solve(t(X5) %*% X5) %*% t(X5) %*% Y3
print(t5Hat)

#model 5.2.2

y5Hat = X5 %*% t5Hat
error5 = Y3 - y5Hat

RSS5 = sum((Y3-y5Hat)^2)
print(RSS5)

#model 5.2.3

n5 = length(Y1)
n5
sigma5=RSS5/(n5-1)
p5=-n5/2*log(2*pi)-n5/2*log(sigma5)-(1/(2*sigma5))*RSS5
print(p5)

#model 5.2.4

k5 = 4
AIC5 = 2*k5-2*p5
BIC5 = k5*log(n5)-2*p5
print(AIC5)
print(BIC5)

#model 5.2.5

M5 = as.matrix(error5)
qqnorm(M5, main ='Normal Q-Q Plot for model 5')
qqline(M5)

#----------------------------------------------------------------------------------

#model 4.2.6

#model 4.2.7
data_split=nrow(data)
split_data1=floor(data_split*.7)
train=data[1:split_data1,]
split_data2=split_data1+1
test=data[split_data2:301,]


#model 4.2.7.1
train_Y2=train[,2]
train_Y3=train[,3]
train_Y4=train[,4]
train_Y5=train[,5]
train_Y6=train[,6]

b6 = matrix(1, length(train_Y2), 1)
train_X = cbind(train_Y5,train_Y4^2,train_Y6^3,b6)
train_Hat = solve(t(train_X) %*% train_X) %*% t(train_X) %*% train_Y3
print(train_Hat)

#model 4.2.7.2
test_y1=test[,1]
test_Y2=test[,2]
test_Y3=test[,3]
test_Y4=test[,4]
test_Y5=test[,5]
test_Y6=test[,6]

test_b = matrix(1, length(test_Y2), 1)
test_X= cbind(test_Y5,test_Y4^2,test_Y6^3,test_b)

test_Y_Hat= test_X %*% train_Hat
test_error = test_Y3-test_Y_Hat
print(test_error)
print(test_error)

test_RSS = sum((test_Y3-test_Y_Hat)^2)
print(test_RSS)

test_n = length(test_Y3)
test_n
test_sigma=test_RSS/(test_n-1)
p5=-test_n/2*log(2*pi)-test_n/2*log(test_sigma)-(1/(2*test_sigma))*test_RSS
print(p5)

#model 4.2.7.3
y_var_Hat = matrix(0, test_n, 1)

cov_thetaHat = test_sigma * (solve(t(test_X) %*% test_X))

for( i in 1:test_n){
  X_i = matrix( test_X[i,] , 1 , 4 ) 
  y_var_Hat[i,1] = X_i %*% cov_thetaHat %*% t(X_i)
}

CI = 2 * sqrt(y_var_Hat)
print(CI)

plot(test_y1 , test_Y_Hat, type = 'l')
segments(test_y1, test_Y_Hat-CI, test_y1 , test_Y_Hat+CI)

#------------------------------------------------------------------------------
#Task 3
#Task 3.1
absolute = abs(t4Hat)
theta_main1=max(absolute)
theta_main2=max(absolute[absolute != max(absolute)])
theta_main3=min(absolute)
theta_main4=min(absolute[absolute != min(absolute)])
theta_main1
theta_main2
theta_main3
theta_main4



#Task 3.2

theta1lower=theta_main1-.5*(theta_main1)
theta1upper=theta_main1+.5*(theta_main1)
theta2lower=theta_main2-.5*(theta_main2)
theta2upper=theta_main2+.5*(theta_main2)


#Task 3.3
accept=NULL
for(j in 1:10000){
  sample1=runif(1,theta1lower,theta1upper)
  sample2=runif(1,theta2lower,theta2upper)
  t4Hat=c(sample1,theta_main3,theta_main4,sample2)
  y4hat=X4%*%t4Hat
  RSS4 = sum((Y3-y4Hat)^2)
  if(RSS4<4){
    accept=rbind(accept,t4Hat)
  }
}
accept

# Task 3.4
hist(accept[,1],probability = TRUE,main = "Theta 1 Histogram", xlab="Theta1")
lines(density((accept[,1])))

hist(accept[,4],probability = TRUE,main = "Theta 2 Histogram", xlab="Theta2")
lines(density(accept[,4]))
plot(accept[,1],accept[,4],main = "Joint Plot", xlab="Theta1",ylab="Theta2")

