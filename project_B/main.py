import numpy as np
import yfinance as yf
import pandas as pd

def get_norm_return(pool, start_date, end_date, dropna=True):
    data = yf.download(" ".join(pool), start=start_date, end=end_date)
    adj_close = data['Adj Close']
    if dropna:
        return (adj_close / adj_close.iloc[0]).dropna(axis=1)
    else:
        return adj_close / adj_close.iloc[0]

def random_portfolio(pool, k='rand', k_low=2, k_up=50, w='uniform'):
    if k == 'rand': 
        k_stocks = np.random.randint(k_low, k_up)
    else:
        k_stocks = k
    if w == 'uniform':
        weights = np.ones(k_stocks) / k_stocks
    elif w == 'rand':
        weights = np.random.rand(k_stocks)
        weights = weights / weights.sum()
    else:
        weights = w
    return (np.random.choice(pool, k_stocks, replace=False), weights)