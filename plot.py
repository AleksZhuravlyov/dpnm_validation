import pandas as pd

import matplotlib.pyplot as plt

import numpy as np

import json

with open('params.json', "r") as f:
    params = json.load(f)

U_release_av = np.load('U_release_av.npy')
U_outlet_av = np.load('U_outlet_av.npy')
times = np.load('times.npy')

data = pd.DataFrame()

data['time'] = times
data['U_release_av'] = U_release_av
data['U_outlet_av'] = U_outlet_av

data.index = data['time']
del data['time']
data.index.name = 'time'


def calculate_U_av_p(h, mu, dP, dx):
    return h ** 2 / 12 / mu * dP / dx


U_av_p = calculate_U_av_p(2 * params['y_coord_in'],
                          params['viscosity'] * params['density'],
                          (params['p_inlet'] - params['p_outlet']) * params[
                              'density'], params['len_x'])

U_outlet_p = 2 * U_release_av + U_av_p

data['U_outlet_p'] = U_outlet_p

data = data.drop(0.005)

# data.plot()
# plt.show()

error = 2 * (U_outlet_av - U_outlet_p) / (U_outlet_av + U_outlet_p)

mean = np.mean(error)

print(mean)
