```python
from aide_design.play import*
import CoastalFunctions as CF
```

Erica Marroquin

Lab 3: Wave Shoaling and Breaking


```python
#given variables
frequency = np.array([1.25, 0.85, 0.5, 0.50]) #Hz
period = 1/frequency #seconds
slope = 1/10
h_cm = 19.7*u.cm #SWL in cm
h = (h_cm.to(u.m)).magnitude #SWL in m
x = np.array([-1, -0.6, -0.2, 0, 0.1, 0.2]) #distances along x-axis in m

#vertical height of each x-value measured
h_wrt_slope = x*slope #m

#measured variables
H0_mm = np.array([60, 74, 34, 25])*(u.mm) #initial height in mm
H0 = (H0_mm.to(u.m)).magnitude #initial height in m
a0 = H0/2 #initial amplitude in m

W1_a_mm = np.array([45, 15, 12, 7.5, 2.5,1])*(u.mm)
W1_a = (W1_a_mm.to(u.m)).magnitude #m
W2_a_mm = np.array([40, 30, 26, 14, 11, 9])*(u.mm)
W2_a = (W2_a_mm.to(u.m)).magnitude #m
W3_a_mm = np.array([21, 19, 20, 9, 5, 3])*(u.mm)
W3_a = (W3_a_mm.to(u.m)).magnitude #m
W4_a_mm = np.array([16, 15, 13.5, 7, 4, 2])*(u.mm)
W4_a = (W4_a_mm.to(u.m)).magnitude #m

# nondimensionalizing water height and amplitude
h_nd = h_wrt_slope/h

W1_a_nd = W1_a/a0[0]
W2_a_nd = W2_a/a0[1]
W3_a_nd = W3_a/a0[2]
W4_a_nd = W4_a/a0[3]

#graphing wave non-dimensionalized amplitude from experimental data
plt.plot(h_nd, W1_a_nd, 'ro', label = "Wave Case 1")
plt.plot(h_nd, W2_a_nd, 'yo', label = "Wave Case 2")
plt.plot(h_nd, W3_a_nd, 'go', label = "Wave Case 3")
plt.plot(h_nd, W4_a_nd, 'bo', label = "Wave Case 4")
plt.xlabel('Non-dimensionalized Water Depth')
plt.ylabel('Non-dimensionalized Amplitude')
plt.legend(loc = 'best')
plt.title('Non-dimensionalized Amplitude from Experimental Data vs Water Depth')
plt.savefig('ND_A_vs_h.png')
plt.show()
```

```python
#calculations for theoretical values
g = (pc.gravity).magnitude

#finding each initial wavenumber
k_theo_init = np.array([CF.wavenumber(period[0],h),
                        CF.wavenumber(period[1],h),
                        CF.wavenumber(period[2],h),
                        CF.wavenumber(period[3],h)])

#finding initial wave phase speed (wave solarity)
sigma = (2*np.pi/(period))
C_p_init = sigma/(k_theo_init)

#finding initial wavelength
wavelen_init = C_p_init*period
wavelen_init

#finding initial n values
n_theo_init = np.array([CF.n(k_theo_init[0],h),
                        CF.n(k_theo_init[1],h),
                        CF.n(k_theo_init[2],h),
                        CF.n(k_theo_init[3],h)])

#initial group velocity
C_g_init = n_theo_init*C_p_init

surf_sim_theo = np.ones(len(wavelen_init))
#calculating Iribarren number
for i in range(0,len(wavelen_init)):
  surf_sim_theo[i] = slope/(np.sqrt(H0[i]/(wavelen_init[i])))
  if surf_sim_theo[i] < 0.5:
    print('This is a spilling breaker with Iribarren number of', surf_sim_theo[i],'.')
  elif surf_sim_theo[i] > 3.3:
    print('This is a surging breaker with Iribarren number of', surf_sim_theo[i],'.')
  else:
    print('This is a plunging breaker with Iribarren number of', surf_sim_theo[i],'.')

#calculating k for each wave case
h_theo = abs(np.linspace(-1, 0.001, num=50))*slope

k1_theo = np.ones(len(h_theo))
k2_theo = np.ones(len(h_theo))
k3_theo = np.ones(len(h_theo))
k4_theo = np.ones(len(h_theo))

for i in range(0,len(h_theo)):
  k1_theo[i] = CF.wavenumber(period[0], h_theo[i])
  k2_theo[i] = CF.wavenumber(period[1], h_theo[i])
  k3_theo[i] = CF.wavenumber(period[2], h_theo[i])
  k4_theo[i] = CF.wavenumber(period[3], h_theo[i])

#calculting n values
n1_theo = CF.n(k1_theo,h_theo)
n2_theo = CF.n(k2_theo,h_theo)
n3_theo = CF.n(k3_theo,h_theo)
n4_theo = CF.n(k4_theo,h_theo)

#non-dimensionalized wave height using theoretical values
h_nd_theo = h_theo/h

#non-dimensionalized wave amplitude using theoretical values
W1_a_nd_theo = np.sqrt(C_g_init[0]/(n1_theo*(np.sqrt(g*h_theo))))
W2_a_nd_theo = np.sqrt(C_g_init[1]/(n2_theo*(np.sqrt(g*h_theo))))
W3_a_nd_theo = np.sqrt(C_g_init[2]/(n3_theo*(np.sqrt(g*h_theo))))
W4_a_nd_theo = np.sqrt(C_g_init[3]/(n4_theo*(np.sqrt(g*h_theo))))
```

```python
#graphing wave non-dimensionalized amplitude using theoretical data
plt.plot(-h_nd_theo, W1_a_nd_theo, label = "Wave Case 1")
plt.plot(-h_nd_theo, W2_a_nd_theo, label = "Wave Case 2")
plt.plot(-h_nd_theo, W3_a_nd_theo, label = "Wave Case 3")
plt.plot(-h_nd_theo, W4_a_nd_theo, label = "Wave Case 4")
plt.xlabel('Non-dimensionalized Water Depth')
plt.ylabel('Non-dimensionalized Wave Amplitude')
plt.legend(loc = 'best')
plt.title('Non-dimensionalized Amplitude from Theoretical Values vs Water Depth')
plt.savefig('ND_A_vs_h_theo')
plt.show()
```

```python
#comparing theoretical and experimental data
plt.plot(-h_nd_theo, W1_a_nd_theo, label = "Theoretical Wave Case 1")
plt.plot(-h_nd_theo, W2_a_nd_theo, label = "Theoretical Wave Case 2")
plt.plot(-h_nd_theo, W3_a_nd_theo, label = "Theoretical Wave Case 3")
plt.plot(-h_nd_theo, W4_a_nd_theo, label = "Theoretical Wave Case 4")
plt.plot(h_nd[0:4], W1_a_nd[0:4], 'ro', label = "Experimental Wave Case 1")
plt.plot(h_nd[0:4], W2_a_nd[0:4], 'yo', label = "Experimental Wave Case 2")
plt.plot(h_nd[0:4], W3_a_nd[0:4], 'go', label = "Experimental Wave Case 3")
plt.plot(h_nd[0:4], W4_a_nd[0:4], 'bo', label = "Experimental Wave Case 4")
plt.xlabel('Non-dimensionalized Water Depth')
plt.ylabel('Non-dimensionalized Wave Amplitude')
plt.legend(loc = 'best')
plt.title('Non-dimensionalized Amplitude vs Water Depth')
plt.savefig('ND_A_vs_h_comparison')
plt.show()

```
