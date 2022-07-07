import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras

# init model
model = keras.Sequential([keras.layers.Dense(units=1, input_shape=[1])])
model.compile(optimizer='sgd', loss='mean_squared_error')

# read and store data
f = pd.read_csv('data.csv')
xs = [float(i) for i in list(f['x'])]
ys = [float(i) for i in list(f['y'])]

# fit model
model.fit(xs, ys, epochs=500)

# output
print("value for y when x is 10.0: {}".format(model.predict([10.0])))
