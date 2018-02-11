import numpy as np

data = [[0.28, 1.31, -6.2, 0.011, 1.03, -0.21, 1.36, 2.17, 0.14],
        [0.07, 0.58, -0.78, 1.27, 1.28, 0.08, 1.41, 1.45, -0.38],
        [1.54, 2.01, -1.63, 0.13, 3.12, 0.16, 1.22, 0.99, 0.69],
        [-0.44, 1.18, -4.32, -0.21, 1.23, -0.11, 2.46, 2.19, 1.31],
        [-0.81, 0.21, 5.73, -2.18, 1.39, -0.19, 0.68, 0.79, 0.87],
        [1.52, 3.16, 2.77, 0.34, 1.96, -0.16, 2.51, 3.22, 1.35],
        [2.2, 2.42, -0.19, -1.38, 0.94, 0.45, 0.6, 2.44, 0.92],
        [0.91, 1.94, 6.21, -0.12, 0.82, 0.17, 0.64, 0.13, 0.97],
        [0.65, 1.93, 4.38, -1.44, 2.31, 0.14, 0.85, 0.58, 0.99],
        [-0.26, 0.82, -0.96, 0.26, 1.94, 0.08, 0.66, 0.51, 0.88]]
data = np.array(data)

def parzen(x, train, h):
    n, d = train.size()
    pnx = 0
    for i in range(n):
        yi = x - train[i, :].T
        pnx += (2*np.pi)**(-float(d)/2)*h**(-d)*np.exp(np.dot(-yi.T, yi)/(2*h**2))
    return pnx/n

w1 = data[:, 1:4]
w2 = data[:, 4:7]
w3 = data[:, 7:]
n = 10
h = 1
tp = [[0.5, 1, 0], [0.31, 1.51, -0.5], [-0.3, 0.44, -0.1]]
for t in range(3):
    p1 = parzen(tp[t,:].T, w1, h)
    p2 = parzen(tp[t, :].T, w2, h)
    p3 = parzen(tp[t, :].T, w3, h)
    maxp, maxi = max([p1, p2, p3])
    print(maxi)


