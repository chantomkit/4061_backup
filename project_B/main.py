import numpy as np
import yfinance as yf
import pandas as pd

def stocks_norm_return(pool, start_date, end_date, dropna=True):
    data = yf.download(" ".join(pool), start=start_date, end=end_date)
    adj_close = data['Adj Close']
    if dropna:
        return (adj_close / adj_close.iloc[0]).dropna(axis=1)
    else:
        return adj_close / adj_close.iloc[0]

def init_stock_choice(pool, k='rand', k_low=2, k_up=50):
    if k == 'rand': 
        k_stocks = np.random.randint(k_low, k_up)
    else:
        k_stocks = k
    return np.random.choice(pool, k_stocks, replace=False)

def init_weight_choice(k_stocks, w='uniform'):
    if w == 'uniform':
        weights = np.ones(k_stocks) / k_stocks
    elif w == 'rand':
        weights = np.random.rand(k_stocks)
        weights = weights / weights.sum()
    return weights

def build_portfolio(norm_return, capital, stocks, weights):
    portfolio = norm_return[stocks] * weights * capital
    portfolio['TotalPos'] = portfolio.sum(axis=1)
    portfolio['PercentagePos'] = portfolio.TotalPos / portfolio.TotalPos.iloc[0]
    portfolio['DailyPercentageReturn'] = portfolio['TotalPos'].pct_change(1)
    return portfolio

def sharpe(ret, vol, kval=1, risk_free_ret=0.04):
    return (ret - risk_free_ret) * kval / vol

class MC_portfolio():
    def __init__(self, norm_return, capital):
        self.norm_return = norm_return
        self.capital = capital
        self.trade_periods = {'D': 252, 'W': 52, 'M': 12}

    def mc_stock_choice(self, iter=1000, k='rand', k_low=2, k_up=50, w='uniform', freq='D', history=False):
        self.best_stocks, self.best_sharpe = None, -np.inf
        if history: self.all_stocks, self.all_sharpe = [], []
        for i in range(iter):
            stocks = init_stock_choice(self.norm_return.columns, k=k, k_low=k_low, k_up=k_up)
            weights = init_weight_choice(len(stocks), w=w)

            portfolio = build_portfolio(self.norm_return, self.capital, stocks, weights)
            mreturn, stdreturn = portfolio.DailyPercentageReturn.mean(), portfolio.DailyPercentageReturn.std()
            sp = sharpe(mreturn, stdreturn, kval=np.sqrt(self.trade_periods[freq]), risk_free_ret=0.04/self.trade_periods[freq])
            
            if history:
                self.all_sharpe.append(sp)
                self.all_stocks.append(stocks)

            if sp > self.best_sharpe:
                self.best_sharpe = sp
                self.best_stocks = stocks
        return

    def mc_weight_choice(self, stocks, iter=1000, freq='D'):
        self.best_weights, self.best_sharpe = None, -np.inf
        for i in range(iter):
            weights = init_weight_choice(len(stocks), w='rand')
            norm_return_tmp = self.norm_return[stocks]

            log_return = np.log(norm_return_tmp/norm_return_tmp.shift(1))
            exp_return = np.sum(log_return.mean() * weights) * self.trade_periods[freq]
            exp_volatility = np.sqrt(np.dot(weights.T, np.dot(log_return.cov() * self.trade_periods[freq], weights)))

            sp = sharpe(exp_return, exp_volatility)
            if sp > self.best_sharpe:
                self.best_sharpe = sp
                self.best_weights = weights
        return