#elevation = cost
#W will be the variable name for the parameters of the function, they will be initialized randomly and learning rate (a) will also be selected
#To calculate the gradient using all the data
#w = w-a*gradient

import numpy as np

def f(w1,w2,x):
    yhat = w1 + w2*x
    return yhat

def dx_w1(w1,w2,x,y):
    yhat = f(w1,w2,x)
    gradient = 2*(yhat - y)
    return gradient
def dx_w2(w1,w2,x,y):
    yhat = f(w1,w2,x)
    gradient = 2*x*(yhat - y)
    return gradient
def gradient_w1(w1,w2,xs,ys):
    N = len(ys)
    total = 0
    for x,y in zip(xs,ys):
        total= total + dx_w1(w1,w2,x,y)
    gradient = total/N
    return gradient

def gradient_w2(w1,w2,xs,ys):
    N = len(ys)
    total = 0
    for x,y in zip(xs,ys):
        total= total + dx_w2(w1,w2,x,y)
    gradient = total/N
    return gradient



def gradient_descent(xs, ys, learnign_rate=0.01, max_num_iteration = 1000):
    w1 = np.random.uniform(0,1,1)
    w2 = np.random.uniform(0,1,1)

    for i in range(max_num_iteration):
        w1 = w1 - learnign_rate*gradient_w1(w1,w2,xs,ys)
        w2 = w2 - learnign_rate*gradient_w2(w1,w2,xs,ys)
        
        if i % 100 == 0:
            print (f"Iteration {i}")
            print (f"W1 = {w1}")
            print (f"W2 = {w2}")
    return (w1,w2)
            

xs= [1,2,3,4,5,6,7]
ys= [1,2,3,4,5,6,7]
(w1,w2) = gradient_descent(xs,ys)
print(w1,w2)
