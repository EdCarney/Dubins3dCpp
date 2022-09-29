import sys
import matplotlib.pyplot as plt
import numpy as np

if (len(sys.argv) > 1):
    fileName = sys.argv[1]
else:
    fileName = "path.txt"

with open("path.txt", "r") as f:
    lines = f.readlines()

points = [line.split() for line in lines]

path = np.array([[],[],[]])
for point in points:
        path = np.array([
            np.append(path[0], float(point[0])),
            np.append(path[1], float(point[1])),
            np.append(path[2], float(point[2]))
            ])

fig = plt.figure()
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax = fig.add_subplot(projection='3d')

ax.plot3D(path[0], path[1], path[2], "g-")

plt.show()