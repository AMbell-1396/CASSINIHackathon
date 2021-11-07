from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import mpl_toolkits
#from mpl_toolkits.basemap import Basemap

#CONSTANTS
LL2M = 111139

nc_file = 'FMI-ARC-SEAICE_THICK-CS2SMOS_NRT_OBS_1636193060855.nc'
data = Dataset(nc_file, modem='r')


t = 0
lon = data.variables['yc'][:]
#lon = (lon[t])
print("lon ",lon)
lat = data.variables['xc'][:]
#lat = (lat[t])
print("lat ",lat)

T = data.variables['analysis_sea_ice_thickness'][:]
print(type(T[t]))

#print("time = ", time)

#COORDINATES
# latitude = x axis
# longitude = y axis
#a = np.array([-156.2263, 73.0534])
#b = np.array([-128.7232, 75.1035])
#POINTS A AND B GIVEN AS A PROJECTION
a = np.array([-1000, 400])
b = np.array([600, 600])
#a = lonlan2m(a)
print("a",len(a))
#b = lonlan2m(b)
m = (b[1] - a[1]) / (b[0] - a[0])
print(m)

#def line(a,b,x):
#    m = (b[1] - a[1]) / (b[0] - a[0])
#    c = b[1] - m*b[0]
#    y = m * x + c
#    return y

#find nearest point in data to the actual initial point
near_xa_m = abs(lat - a[0])
near_xa = min(near_xa_m)
near_xb_m = abs(lat - b[0])
near_xb = min(near_xb_m)
near_ya_m = abs(lon - a[1])
print(near_ya_m)
near_ya = min(near_ya_m)
near_yb_m = abs(lon - b[1])
print(near_yb_m)
near_yb = min(near_yb_m)

#index of points in the array
i_near_xa = ma.where(near_xa_m == near_xa)[0][0]
i_near_xb = ma.where(near_xb_m == near_xb)[0][0]
i_near_ya = ma.where(near_ya_m == near_ya)[0][0]
i_near_yb = ma.where(near_yb_m == near_yb)[0][0]
print("NEAR AX", i_near_xa)
print("NEAR AY", i_near_ya)
print("NEAR BX", i_near_xb)
print("NEAR BY", i_near_yb)


def coords_mask(initial_x_index,initial_y_index,final_x_index,final_y_index):
    m = (final_y_index - initial_y_index)/(final_x_index - initial_x_index)
    #initial_x_index = i_near_xa
    x = initial_x_index
    #print("initial x",x)
    array_x_coords = np.array([])
    array_y_coords = np.array([])
    while x < final_x_index:
        #print("inside while")
        y = (m * x) + initial_y_index
        #print("y",y)
        #near_y_m = abs(lat - y)
        #print(near_y_m)
        #near_y = min(near_y_m)
        #print("neary",near_y)
        #print("some", ma.where(near_y_m == 483.25)[0][0])
        #i_near_y = ma.where(near_y_m == near_y)[0]
        i_near_y = round(y)
        array_y_coords = np.append(array_y_coords, i_near_y)
        array_x_coords = np.append(array_x_coords, x)
        x = x+1
    return [array_x_coords.astype(int), array_y_coords.astype(int)]


def matrix_mask(x_coords, y_coords):
    mask = np.zeros((len(lat), len(lon)))
    #diagonal
    for i in range(len(x_coords)):
        mask[x_coords[i]][y_coords[i]] = 1
        i = i+1
    return mask




XY = coords_mask(i_near_xa,i_near_ya,i_near_xb,i_near_yb)
print("XY",XY)
MASK = matrix_mask(XY[0],XY[1]).transpose()
#print(MASK)
plt.subplot(1,2,1)
plt.imshow(MASK, interpolation='None')
print(len(T[t][0]))
print(len(MASK[0]))
print(T[t]*MASK)
#masked_data = ma.masked_array(T[t],MASK)
masked_data = T[t]*MASK
print(masked_data)
plt.subplot(1,2,2)
plt.imshow(masked_data, interpolation='None')
plt.show()

def average_data(data, masked):
    d1_array_averages = np.array([])
    for i in range(len(lon)):
        for j in range(len(lat)):
            if(masked[i][j] == 1):
                left = data[i-1][j]
                right = data[i+1][j]
                up = data[i][j+1]
                down = data[i][j-1]
                center = data[i][j]
                average = (center+left+right+up+down)/5
                d1_array_averages = np.append(d1_array_averages,average)
            j = j+1
        i = i+1
    return d1_array_averages

print(average_data(T[t],MASK))
print(T[t][97][85])
print(T[t][96][85])
print(T[t][98][85])
print(T[t][97][84])
print(T[t][97][86])



