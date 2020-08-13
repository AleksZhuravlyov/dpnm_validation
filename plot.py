# MIT License
# 
# Copyright (c) 2020 Aleksandr Zhuravlyov and Zakhar Lanets
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


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
