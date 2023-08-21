#!/usr/bin/env python3

import numpy as np

mu = np.array([[0.1,0.2,0.3],[0.9,0.8,0.7]])
pi = np.array([0.5,0.5])
n = 10

k, d = mu.shape


print(f"#@N:{n}")
print(f"#@D:{d}")

for i in range(n):
    for j in range(d):
        meth = np.random.choice(range(k), p=pi)
        print(f"{i}\tr{i}\t{j}\t{j}\t{meth}")
