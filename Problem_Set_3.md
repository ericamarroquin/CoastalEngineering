## Erica Marroquin
## CEE 4350: Coastal Engineering
## Problem Set #3

```python
from aide_design.play import*
import CoastalFunctions as CF #a file i made to keep track of all functions
import pypandoc
```

### Question 1a
##### Givens
```python
T = 12*(u.s)
h2 = 2*(u.m)
a_1 = 0.5*(u.m)
sigma = 2*np.pi/(T)
```
The minimum clearance needed to maintain a dry deck would be the max amplitude of the waves against the first row of piles. To get the amplitude of those waves, energy conservation can be used. The amount of energy is the same at both the piles of the dock and the deep water where the buoy is.

The energy conservation equation is as follows

$$\overline{{E_1}^+} {C_g}_1 = \overline{{E_2}^+} {C_g}_2$$


where  $E = \frac{1}{2}\rho g a^{2}$ and

(1) refers to the position at sea and
(2) refers to the position at the piles.

The energy conservation equation can be solved for:
$$a_2 = a_1 \sqrt{\frac{{C_g}_1}{{C_g}_2}}$$

$$C_g = C_p n$$
$$C_p = \frac{\sigma}{k}$$
 and $n = \frac{1}{2}$ in deep water; $n = 1$ in shallow water.

In order to find the wave number in the deep water, the assumption that $tanh(kh) \rightarrow 1$ was made. Thus, the dispersion relationship becomes:
$$\sigma^2 = g k$$

```python
k_1 = (sigma**2)/(CF.g)
C_p1 = sigma/(k_1)
n1 = 0.5 #in deep water
C_g1 = C_p1*n1
```
The assumption was made that the waves hitting the piles are shallow water waves.
```python
C_p2 = np.sqrt(CF.g*h2)
n2 = 1 #in shallow water
C_g2 = C_p2*n2

a_2 = a_1*(np.sqrt(C_g1/C_g2))
print(a_2.magnitude)
```
The amplitude of the shallow water waves at the piles is 0.727 m. Assuming that “clearance” starts at the front row of piles (at $h_2$) the maximum clearance needed to keep the deck dry would be 0.727 m.

### Question 1b
```python
k_2 = CF.wavenumber(T.magnitude,h2.magnitude)
wave_length2 = (2*np.pi)/(k_2)

cond = h2.magnitude/wave_length2

if cond < 0.05: #h/lamba < 1/20
  print('This is a shallow water wave.')
elif cond > 0.5: #h/lamba > 1/2
  print('This is a deep water wave.')
else:
  print('This is a transitional wave.')
```
The waves at the pile location are shallow water waves. According to the water limits if $\frac{h}{\lambda} < \frac{1}{20} = 0.05$ the wave is considered shallow. In this case $\frac{h}{\lambda} = 0.0379$.

### Question 1c
The formula for vertical particle velocity in a standing wave is
$$u = a \sigma \frac{cosh(k(h+z))}{sinh(kh)} cos(kx-\sigma t)$$


The following assumptions can be made:
- At the maximum water height $k(h+z) \rightarrow 0$. The $cosh$ term then goes to 1.
- The maximum values for $u$ happen at the wave crests and troughs. In a $cos$ wave, this occurs at $0, \pi, 2 \pi, ..,$ etc. If $(kx- \sigma t) = 0$ then the $cos$ term goes to 1.
- In a shallow water wave $sinh(kh) \rightarrow kh$.

The formula for vertical particle velocity in shallow water then becomes
$$u = \frac{a \sigma}{kh}$$

```python
u_shallow = (a_2*sigma.magnitude)/(k_2*h2)
print(u_shallow)
```
The vertical particle velocity in shallow water is 1.597 m/s.

The formula for horizontal particle velocity in a standing wave is
$$u = a \sigma \frac{sinh(k(h+z))}{sinh(kh)} sin(kx-\sigma t)$$

The following assumptions can be made:
- Maximum $w$ values occurs at the mean water level. This is $z = 0$. Thus, $sinh(k(h+z)) \rightarrow kh$.
- The maximum values for $w$ happens when the $sin$ wave is at its maximum at $\frac{\pi}{2}, \frac{3\pi}{2},...,$ etc. At these values, the $sin$ term goes to 1.
- In a shallow water wave $sinh(kh) \rightarrow kh$.

The formula for horizontal particle velocity in shallow water then becomes
$$w = a \sigma$$

```python
w_shallow = a_2*sigma.magnitude
print(w_shallow)
```

The horizontal particle velocity in shallow water is 0.3807 m/s.

### Question 2
To get the corresponding wave amplitude at the section where the width is $B_2$, the following assumptions are made:
- The waves are linear
- Group velocity is the same throughout the channel, ${C_g}_1 = {C_g}_2$

To get amplitude based on width instead of height, the energy conservation equation (per unit width) needs to be multiplied by the width of the channel.

$$\overline{{E_1}^+} {C_g}_1 B_1= \overline{{E_2}^+} {C_g}_2 B_2$$
$$\frac{1}{2} \rho g {A_1}^2 B_1 = \frac{1}{2} \rho g {A_2}^2 B_2$$
$$\frac{{A_1}^2}{{A_2}^2} = \frac{B_2}{B_1}$$
The wave amplitude, $A_2$, where the channel width is $B_2$ is
$$A_2 = A_1 \sqrt{\frac{B_2}{B_1}}$$
