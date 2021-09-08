


import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow import keras
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn import model_selection
from keras.models import Sequential
from keras.layers import Dense, Activation
from keras.optimizers import SGD
import matplotlib.figure as figure
from sklearn.model_selection import KFold
from sklearn import metrics
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
import optuna
import keras.backend as K
from keras.layers import Convolution2D, Input, Dense, GlobalAveragePooling2D
from keras.models import Model
from keras.utils import to_categorical
import time




data = np.loadtxt('inputfile_gene.csv', delimiter=',', dtype=float, skiprows=1)

labels = data[:, 0:1] 
features = data[:, 1:] 
X = features
Y = labels.ravel()

autoscaledX = (X - X.mean(axis=0)) / X.std(axis=0, ddof=1)
autoscaledY = (Y - Y.mean()) / Y.std(ddof=1)

X_train, X_test, y_train, y_test = train_test_split(autoscaledX, autoscaledY, test_size=0.3, random_state=0)






def create_model(num_layer, activation, num_filters):
    inputs = Input(shape=(autoscaledX.shape[1],))
    x = Dense(units=num_filters[0], activation=activation)(inputs)
    for i in range(1,num_layer):
        x = Dense(units=num_filters[i], activation=activation)(x)
    predictions = Dense(units=1)(x)
    model = Model(input=inputs, output=predictions)
    return model


def objective(trial):
   
    K.clear_session()


    num_layer = trial.suggest_int("num_layer", 1, 3)

    
    num_filters = [int(trial.suggest_discrete_uniform("num_filter_"+str(i), 100, 500, 100)) for i in range(num_layer)]

    
    activation = trial.suggest_categorical("activation", ["relu", "sigmoid"])

    
    optimizer = trial.suggest_categorical("optimizer", ["rmsprop"])

    model = create_model(num_layer, activation, num_filters)
    model.compile(optimizer=optimizer,
          loss="mse",
          metrics=["mae"])

    model.fit(X_train, y_train, verbose=0, epochs=20, batch_size=16, validation_split=0.1)
    tuna_pred_test = model.predict(X_test)
    
    return (1.0 - (r2_score(y_test, tuna_pred_test)))




starttime = time.time()
study = optuna.create_study()
study.optimize(objective, n_trials=100)




elapsedtime = time.time() - starttime
print ("Elapsed time in hyperparameter optimization: {0} [sec]".format(elapsedtime))




print(study.best_params)
print(study.best_value)
print(study.best_trial)



1-study.best_value






