#A .py file containing useful functions for performing coastal engineering calculations.
#A huge thank you to Zoe, https://github.com/zam7, for letting me add her functions to this file.

from aide_design.play import*
g = pc.gravity

def wavenumber(T, h):
  """Return the wavenumber of wave using period and water height from bed."""
  k = 10  # this is a guess to find what k is
  diff = (((2*np.pi)/T)**2)-(g.magnitude * k * np.tanh(k*h))
  while diff<0:
      LHS = ((2*np.pi)/T)**2
      RHS = g.magnitude * k * np.tanh(k*h)
      diff = LHS - RHS
      k = k - 0.0001
  return k
