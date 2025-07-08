#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(0, 2*np.pi, 1000)
x = 16 * np.sin(t)**3
y = 13 * np.cos(t) - 5 * np.cos(2*t) - 2 * np.cos(3*t) - np.cos(4*t)

plt.figure(figsize=(8, 6))
plt.plot(x, y, color='red', linewidth=3)
plt.fill(x, y, 'pink')
plt.title('I like you', fontsize=16)
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()
