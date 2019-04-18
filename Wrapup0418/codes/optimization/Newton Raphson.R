## Newton Raphson to find the max
## 2019-2-28

f = function(x){(x[1]+1)^2 + (x[2]-1)^2 + 3}

eps = 0.00001 

### the first derivative
f_1 = function(X){
  p = length(X)
  y = c()
  for(i in 1:p){
    X1 = X
    X1[i] = X1[i] + eps
    y = c(y, (f(X1) - f(X))/eps)
  }
  #print(f(X1))
  return(y)
}

f1 = function(x){
  res = (f(x + eps) - f(x))/eps
  return(res)
}

xold = 100

for(i in 1:1000){
  print(i)
  xnew = xold - f(xold)/f1(xold)
  if(abs(xnew - xold)<eps){print(xnew); res = xnew; break}
  xold = xnew
  print(xold)
}

why = c()
for(i in seq(50,51,0.0001)){
  print(i)
  why = c(why, f(i))
}



f2 = function(x){
  return((f1(x + eps) - f1(x))/eps)
}

### the partial derivative
f_1_part = function(X,i){
  X1 = X
  X1[i] = X[i] + eps
  y = (f(X1) - f(X))/eps
 # print(f(X1))
  return(y)
}

### the second derivative
f_2 = function(X){
  p = length(X)
  y = matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      X1 = X
      X1[j] = X[j] + eps
      f1 = f_1_part(X1,i)  
      f2 = f_1_part(X,i)
      y[i,j] = (f1 - f2)/eps
    }
  }
  return(y)
}

### Newton method
begin = Sys.time()
#x_old = c(1,1)
x_old  = 0
value = c(); x_value = c();x_value2 = c()
set.seed(123)
for(times in 1:100){
  print('********')
  print(times)
  Fx = f_1(x_old)
  value = c(value,Fx)
  x_value = c(x_value, x_old[1])
  x_value2 = c(x_value2, x_old[2])
  F2 = f_2(x_old)
  x_new = x_old - solve(F2) %*% Fx
  if(sum(abs(x_old - x_new)) < 10e-7){
    print('coveraged!')
    break
  }
  x_old = x_new
}
end = Sys.time()
x_old

# 100 times not coverage. Took about 3 hours. 

plot(seq(50,50.1280, 0.0001),why, type = 'b',cex = 0.1)