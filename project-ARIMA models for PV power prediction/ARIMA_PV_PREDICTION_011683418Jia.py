#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 22:44:32 2020
ARIMA Model for course 570 final project
@author: linli 011686418
"""

###import data
import math
import pandas as pd
import numpy as np
from statsmodels.tsa.arima_model import ARIMA
from pandas import DataFrame
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error

# process data, drop lines with'nan'
data = pd.read_excel('solar_data-2017-15min.xlsx',header=0)
data['TIMESTAMP']=data['TIMESTAMP'].str.slice(0,19)
data['TIMESTAMP'] =pd.to_datetime(data['TIMESTAMP'])
data.index= data['TIMESTAMP']
data=data.drop(['TIMESTAMP'],axis=1)
data[data<0]= math.nan
data_pro= data.dropna() 

### parameters configuration
import warnings
# evaluate an ARIMA model for a given order (p,d,q)
def evaluate_arima_model(X, arima_order):
	# prepare training dataset
	train_size = int(len(X) * 0.005)
	train, test = X[0:train_size], X[train_size:int(len(X) * 0.01)]
	history = [x for x in train]
	# make predictions
	predictions = list()    
	for t in range(len(test)):
		model = ARIMA(history, order=arima_order)
		model_fit = model.fit(disp=0)
		yhat = model_fit.forecast()[0]
		predictions.append(yhat)
		history.append(test[t])
	# calculate out of sample error
	error = mean_squared_error(test, predictions)
	return error
 
# evaluate combinations of p, d and q values for an ARIMA model
def evaluate_models(dataset, p_values, d_values, q_values):
	dataset = dataset.astype('float32')
	best_score, best_cfg = float("inf"), None
	for p in p_values:
		for d in d_values:
			for q in q_values:
				order = (p,d,q)
				try:
					mse = evaluate_arima_model(dataset, order)
					if mse < best_score:
						best_score, best_cfg = mse, order
					print('ARIMA%s MSE=%.3f' % (order,mse))
				except:
					continue
	print('Best ARIMA%s MSE=%.3f' % (best_cfg, best_score))
### fit model
# evaluate parameters
p_values = [3,4,5]
d_values = range(0, 2)
q_values = range(0, 2)
warnings.filterwarnings("ignore")
evaluate_models(data_pro.PowerOutput_kW_Avg, p_values, d_values, q_values)

#divide data
X =np.array(data_pro.PowerOutput_kW_Avg) 
EXOG = np.array(data_pro[['SolarIrradiation_Wm2','AmbTemp_C','WindSpeed_ms','WindDir_deg']])
X_list = X.tolist()
EXOG_list = EXOG.tolist()
size = int(len(X) * 0.6)
train, test = X_list[0:size], X_list[size:size+96]
EXOG_train,EXOG_test = EXOG_list[0:size], EXOG_list[size:int(len(X) * 0.8)]

## Model 1-univariate ARIMA model
model = ARIMA(train, order=(3,0,1))
model_fit = model.fit(disp=0)
print(model_fit.summary())
# plot residual errors
residuals = DataFrame(model_fit.resid)
residuals.plot()
plt.title('ARIMA Fit Residual Error Line Plot')
plt.savefig('ARIMA Fit Residual Error Line Plot')
plt.show()
residuals.plot(kind='kde')
plt.title('ARIMA Fit Residual Error Density Plot')
plt.savefig('ARIMA Fit Residual Error Density Plot')
plt.show()
print(residuals.describe())
# predict
EXOG_train,EXOG_test = EXOG_list[0:size], EXOG_list[size:size+96]
start_index = size
end_index = size+95
forecast =model_fit.predict(start=start_index, end= end_index) 

## Model 2-multivariate ARIMA model
# fit model-with exogent
model = ARIMA(train,exog=EXOG_train, order=(3,0,1))
model_fit = model.fit(disp=0)
print(model_fit.summary())
# plot residual errors
residuals = DataFrame(model_fit.resid)
residuals.plot()
plt.title('ARIMA Fit Residual Error Line Plot')
plt.savefig('ARIMA Fit Residual Error Line Plot')
plt.show()
residuals.plot(kind='kde')
plt.title('ARIMA Fit Residual Error Density Plot')
plt.savefig('ARIMA Fit Residual Error Density Plot')
plt.show()
print(residuals.describe())
# predict
EXOG_train,EXOG_test = EXOG_list[0:size], EXOG_list[size:size+96]
start_index = size
end_index = size+95
forecast =model_fit.predict(start=start_index, end= end_index,exog=EXOG_test) 

## Model 3-ARIMA Rolling Forecast model

X = data_pro.PowerOutput_kW_Avg
size = int(len(X) * 0.6)
train, test = X[0:size], X[size:size+96]
history = [x for x in train]
predictions = list()
for t in range(len(test)):
	model = ARIMA(history, order=(3,0,1))
	model_fit = model.fit(disp=0)
	output = model_fit.forecast()
	yhat = output[0]
	predictions.append(yhat)
	obs = test[t]
	history.append(obs)
	print('predicted=%f, expected=%f' % (yhat, obs))
error = mean_squared_error(test, predictions)
print('Test MSE: %.3f' % error)
# plot
plt.plot(range(0,96,1),test[0:96],color = 'blue', label="test")
plt.plot(range(0,96,1),predictions[0:96], color='red', label="predict")
plt.legend(loc = 'upper right')
plt.title('ARIMA Rolling Forecast Line Plot for One Day')
plt.savefig('ARIMA Rolling Forecast Line Plot for One Day')
plt.show()

### Model evaluation
#MAE
from sklearn.metrics import mean_absolute_error
MAE_data = mean_absolute_error(test, predictions)
MAE_data
# MAPE
def mean_absolute_percentage_error(y_true, y_pred): 
    Error = np.sum(np.abs(np.subtract(np.array(y_true), np.array(y_pred))))
    Average = np.sum(y_true)
    return Error/Average
MAPE_data = mean_absolute_percentage_error(test, predictions)
MAPE_data
#RMSE
from sklearn.metrics import mean_squared_error
from math import sqrt
#calculate RMSE
sqrt(mean_squared_error(test, predictions))
#SDE
err = np.subtract(np.array(test), np.array(predictions))
SDE_data =np.std(err)
SDE_data

