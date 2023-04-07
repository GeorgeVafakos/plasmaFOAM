# import fluidfoam
from fluidfoam import readmesh
from fluidfoam import readvector
# from matplotlib.patches import Circle
# import import
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker

# Fonts Configuration
font_title = {'family':'sans-serif', 'weight':'bold', 'size':14}
font_labels = {'family':'sans-serif', 'weight':'normal', 'size':14}

sol = '/home/cfdlab/OpenFOAM/cfdlab-7/opt/tutorials/incompressible/icoFoam/cavity/cavity'

# coordinates
x, y, z = readmesh(sol, structured=True)

nx, ny, nz = x.shape
print("Nx = ", nx, "Ny = ", ny, "Nz = ", nz)

timename = '0.5'
vel = readvector(sol, timename, 'U', structured=True)

# plt.figure()
# levels = np.arange(-0.5, 1.1, 0.1)
# plt.contourf(x[:, :, nz//2], y[:, :, nz//2], vel[0, :, :, nz//2],
#              levels=levels)

# # Setting axis labels
# plt.xlabel('x (m)')
# plt.ylabel('y (m)')

# print(len(x[:, :, nz//2]))
# print(len(y[:, :, nz//2]))
# print(len(vel[1,:, :, nz//2]))
# print(len(vel[0,:, :, nz//2]))

XX = x[:, :, nz//2]
YY = y[:, :, nz//2]
Ux = vel[0,:, :, nz//2]
Uy = vel[1,:, :, nz//2]


fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1,  aspect='equal')
plt.rcParams['contour.negative_linestyle'] = 'solid'
CF =  ax1.contour(x[:, :, nz//2], y[:, :, nz//2], vel[0, :, :, nz//2], 14, colors='k',linewidths=0.8)
PC =  ax1.contourf(x[:, :, nz//2], y[:, :, nz//2], vel[0, :, :, nz//2], 14, cmap='jet',antialiased=True)
# QV =  ax1.quiver(x[:, :, nz//2], y[:, :, nz//2], vel[0, :, :, nz//2], vel[1, :, :, nz//2])
# SL =  ax1.streamplot(XX, YY, Ux, Uy)
plt.title('Re=1000', **font_title)
plt.xlabel('x/L', **font_labels)
plt.ylabel('y/L', **font_labels)
divider = make_axes_locatable(ax1)
cax1 = divider.append_axes("right", size="2%", pad=0.2)
cbar = plt.colorbar(PC, cax = cax1)
fig.tight_layout()



# plt.streamplot(YY, XX, Uy, Ux)

plt.show()
# # velocity field
# vel = readvector(sol, '1200', 'U', False)

# # flip them around
# X = x.reshape(240, 920)
# Y = y.reshape(240, 920)
# Y = Y[::-1, :]

# # magnitude of the velocity field
# magU = np.sqrt( vel[0, :]**2 + vel[1, :]**2 )
# magU.resize(240, 920)
# magU = magU[::-1, :]

# # circle represents the turbine
# circle = Circle((0,0), 5.0, color='black', fill=False)

# # plot the contours
# fig = plt.figure(figsize=(6, 1.5))
# plt.rc('font', size=18)
# plt.rcParams["font.family"] = "monospace"
# plt.title('VAWT wake @ TSR=6.0')
# plt.xlabel('x')
# plt.ylabel('y')
# ax = fig.add_subplot(1,1,1)
# img = ax.imshow(magU,
#                 interpolation='bilinear',
#                 cmap='jet',
#                 extent=[X.min(), X.max(), Y.min(), Y.max()])
# cbar = fig.colorbar(ax=ax,
#                     mappable=img,
#                     orientation='horizontal',
#                     pad=0.3,
#                     label='m/s')
# ax.add_artist(circle)
# plt.show()
