import numpy as np
import matplotlib.pyplot as plt

nxdis = 3
nydis = 3
nzdis = 1

xsiz = 10.
ysiz = 10.
zsiz = 10.

xdb= np.zeros (nxdis*nydis*nzdis)
ydb= np.zeros (nxdis*nydis*nzdis)
zdb= np.zeros (nxdis*nydis*nzdis)

xdis = xsiz  / nxdis
ydis = ysiz  / nydis
zdis = zsiz  / nzdis
i    = -1
xloc = -0.5*(xsiz+xdis)
for ix  in range (nxdis):
    xloc = xloc + xdis
    yloc = -0.5*(ysiz+ydis)
    for iy  in range (nydis):
        yloc = yloc + ydis
        zloc = -0.5*(zsiz+zdis)
        for iz  in range (nzdis):
            zloc = zloc + zdis
            i = i+1
            xdb[i] = xloc + 0.5*xsiz
            ydb[i] = yloc + 0.5*ysiz
            zdb[i] = zloc + 0.5*zsiz



plt.plot (xdb, ydb, 'o')
plt.show()

print xdb
print ydb
print zdb
