import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

np.random.seed(1)

N = 100
x = np.random.rand(N,1)*5
y = 2.3 + 5.1*x
y_obs = y + 2*np.random.randn(N,1)

plt.scatter(x,y_obs,label='Observations')
plt.plot(x,y,c='r',label='True function')
plt.legend()
plt.show()

def f(w,b):
    return w*x+b

def loss_function(e):
    L = np.sum(np.square(e))/N
    return L

def dL_dw(e,w,b):
    return -2*np.sum(e*df_dw(w,b))/N

def df_dw(w,b):
    return x

def dL_db(e,w,b):
    return -2*np.sum(e*df_db(w,b))/N

def df_db(w,b):
    return np.ones(x.shape)

def gradient_descent(iter=100,gamma=0.1):
    # get starting conditions
    w = 10*np.random.randn()
    b = 10*np.random.randn()
    
    params = []
    loss = np.zeros((iter,1))
    for i in range(iter):

        params.append([w,b])
        e = y_obs - f(w,b)
        loss[i] = loss_function(e)


        w_new = w - gamma*dL_dw(e,w,b)
        b_new = b - gamma*dL_db(e,w,b)
        w = w_new
        b = b_new
        
    return params, loss
        
params, loss = gradient_descent()

iter=100
gamma = 0.1
w = 10*np.random.randn()
b = 10*np.random.randn()

params = []
loss = np.zeros((iter,1))
for i in range(iter):

    params.append([w,b])
    e = y_obs - f(w,b)
    loss[i] = loss_function(e)

    w_new = w - gamma*dL_dw(e,w,b)
    b_new = b - gamma*dL_db(e,w,b)
    w = w_new
    b = b_new

dL_dw(e,w,b)


plt.plot(loss)

params = np.array(params)
plt.plot(params[:,0],params[:,1])
plt.title('Gradient descent')
plt.xlabel('w')
plt.ylabel('b')
plt.show()


N = 1000
D = 5
X = 5*np.random.randn(N,D)
w = np.random.randn(D,1)
y = X.dot(w)
y_obs = y + np.random.randn(N,1)

(X*w.T).shape

def f(w):
    return X.dot(w)

def loss_function(e):
    L = e.T.dot(e)/N
    return L

def dL_dw(e,w):
    return -2*X.T.dot(e)/N 

def gradient_descent(iter=100,gamma=1e-3):

    w = np.random.randn(D,1)
    params = []
    loss = np.zeros((iter,1))
    for i in range(iter):
        params.append(w)
        e = y_obs - f(w)
        loss[i] = loss_function(e)

        w = w - gamma*dL_dw(e,w)
        
    return params, loss
        
params, loss = gradient_descent()

plt.plot(loss)

params[-1]

model = LinearRegression(fit_intercept=False)
model.fit(X,y)
model.coef_.T

np.hstack([params[-1],model.coef_.T])

def dL_dw(X,e,w):
    return -2*X.T.dot(e)/len(X)

def gradient_descent(gamma=1e-3, n_epochs=100, batch_size=20, decay=0.9):
    epoch_run = int(len(X)/batch_size)
    

    w = np.random.randn(D,1)
    params = []
    loss = np.zeros((n_epochs,1))
    for i in range(n_epochs):
        params.append(w)
        
        for j in range(epoch_run):
            idx = np.random.choice(len(X),batch_size,replace=False)
            e = y_obs[idx] - X[idx].dot(w)

            w = w - gamma*dL_dw(X[idx],e,w)
        loss[i] = e.T.dot(e)/len(e)    
        gamma = gamma*decay 
        
    return params, loss
        
params, loss = gradient_descent()

plt.plot(loss)

np.hstack([params[-1],model.coef_.T])