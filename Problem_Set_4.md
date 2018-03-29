```python
from aide_design.play import*
import CoastalFunctions as CF
```

Erica Marroquin
CEE 4350: Coastal Engineering
HW 4, Problem 3

### Part i
```python
g = pc.gravity.magnitude

def h_x(delta, x):
  bathy = s*(x + delta * x_0 * np.exp(-(x-x_0)**2/(wavelength)))
  return bathy

s = 1/50
T = 10 #s
x_0 = -1000 #m
wavelength = 5000 #m
delta_x = 50 #m taking a measurement every 50 m to cut down number of iterations that need to occur
array_elements = 4000/delta_x
array_elements_loop = int(array_elements - 1)
x_array = np.linspace(-4000, 0, num = 80)

#plotting bathymetry in the on-offshore direction (i)
#for delta = 0 condition
plt.plot(x_array, h_x(0, x_array))
plt.xlabel('x, the distance in on-offshore direction')
plt.ylabel('h(x), bathymetry')
plt.title('Bathymetry as a function of on-offshore distance where delta = 0')
plt.show()

#for delta = 0.5 condition
plt.plot(x_array, h_x(0.5, x_array))
plt.xlabel('x, the distance in on-offshore direction')
plt.ylabel('h(x), bathymetry')
plt.title('Bathymetry as a function of on-offshore distance where delta = 0.5')
plt.show()

#for delta = -0.5 condition
plt.plot(x_array, h_x(-0.5, x_array))
plt.xlabel('x, the distance in on-offshore direction')
plt.ylabel('h(x), bathymetry')
plt.title('Bathymetry as a function of on-offshore distance where delta = -0.5')
plt.show()
```

### Part ii
```python
#making different arrays to fill with a loop for numerical integration
#if the arrays were expanded each time the loop iterated, the run time would
#be SO LONG. like in MATLAB, predefined arrays are best.
delta_y = np.ones(array_elements_loop)
delta_y[0] = 0
y = np.ones(array_elements_loop)
y[0] = 0
alpha = np.zeros(array_elements_loop)
alpha[0] = 30 #given angle of incidence
h = np.ones(array_elements_loop)
h[0] = abs(h_x(0, x_array[0]))
k = np.ones(array_elements_loop)
k[0] = CF.wavenumber(T, h[0])
Cg = np.zeros(array_elements_loop)
Cg[0] = np.sqrt(g*h[0])
h_0 = h[0]
Cg_0 = Cg[0]
b = np.zeros(array_elements_loop)
b_0 = np.cos(alpha[0])
a = np.zeros(array_elements_loop)
a[0] = 1
a_0 = a[0]

KAPPA = k[0]*np.sin(a_0) #a constant

def wave_rays(lil_delta):
  for i in range(0, array_elements_loop):
    k[i] = CF.wavenumber(T, abs(h_x(lil_delta, x_array[i]))) #finding a new k for each height
    delta_y[i] = KAPPA/(np.sqrt(k[i]**2 + KAPPA**2))*delta_x #finding the new step change
    y[i] = y[i-1] + delta_y[i] #moving over by previous delta y to get a new height
    alpha[i] = np.arctan(delta_y[i]/delta_x) #finding a new angle of incidence for the new delta y value
    Cg[i] = np.sqrt(g*abs(h_x(lil_delta, x_array[i]))) #calculates new group velocity
    b[i] = np.cos(alpha[i])
    a[i] = a_0*np.sqrt(Cg_0/Cg[i])*np.sqrt(b_0/b[i])  #new amplitude
  return k, delta_y, y, alpha, Cg, b, a

```

### Graphs for Part ii and iii
```python
x_graphing = np.linspace(-4000, 1, num = 79) #to make the vectors the same length (delta_y is 79 elements in the array_elements_loop)

#delta = 0
plt.plot(x_graphing, wave_rays(-0.5)[6]) #wave_rays returns an array, and y values are the second index of the arrays
plt.xlabel('Distance to shoreline in on-offshore direction (m) for delta = -0.5')
plt.ylabel('Wave amplitude of wave ray (m)')
plt.title('Wave amplitude when delta = -0.5')
plt.show()
```
