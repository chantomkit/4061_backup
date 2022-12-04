import pandas as pd
import matplotlib.pyplot as plt

for i in range(5, 16):
    name = "p_"+str(round(i * 0.1, 2))+"00000.csv"
    df = pd.read_csv(name, header=None)
    x, y = df[0].values, df[1].values
    plt.scatter(x, y, marker='.', label="p="+str(round(i * 0.1, 2)))
plt.legend()
plt.savefig("q2.png")