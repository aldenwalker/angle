import sage.all as SAGE 

def r(x):
  """round x to the nearest multiple of 1/2; 
  it is an error if x is an integer"""
  if x.is_integral():
    raise ValueError("x can't be an integer")
  return x.floor() + 1/2

def r_set(x, which):
  """round x to the nearest 1/2 + 2*integer if which==g
  or 1/2 + 2*integer + 1 if which==f; this function is *undefined* if 
  x = 1/2 + 2n and we call it to round to 1/2 + 2n + 1 or vice versa"""
  if which=='g':
    return (x-1/2).round('even') + 1/2
  else:
    return (x-1/2).round('odd') + 1/2


def in_domain(x, which, extended=False):
  """decide whether a point lies in the domain given --- if extended==True, 
  this accepts an endpoint of the domain"""
  if not extended and x.is_integral():
    return False
  if extended and x.is_integral():
    return True
  if which == 'g':
    return x.floor()%2 == 0
  else:
    return x.floor()%2 == 1

def which_domain(x):
  """return which domain the point is in -- it shouldn't be an integer 
  for this to be well-defined"""
  if x.is_integral():
    return None
  else:
    return ('g' if x.floor()%2==0 else 'f')

def dist_to_end_of_domain(x, which):
  """return the distance to end of the domain that the point is in 
  (just ceiling(x) - x unless x is an integer)"""
  if not in_domain(x, which, extended=True):
    raise ValueError("Can't get distance to end of domain point isn't in")
  if x.is_integral():
    return (1 if ((which=='g') != (x%2==1)) else 0) #xor 
  return x.ceiling() - x

def E(x, theta, lam, which=None):
  """actually apply the map to the point"""
  if which != None:
    if not in_domain(x, which, extended=True):
      raise ValueError("Point isn't in the desired domain")
    rounded = r_set(x, which)
    return lam*(x-rounded) + rounded + theta
  else:
    rounded = r(x)
    return lam*(x-rounded) + rounded + theta
  

class AugmentedTrajectory:
  def __init__(self, x, T, theta, lam):
    """we can accept either x as a number to start with
    or a list of all the images right away; if we get a 
    list, we don't check the sanity"""
    self.theta = theta
    self.lam = lam
    self.T = T
    if isinstance(x, list):
      self.X = [ell for ell in x] #duplicate the list
    else:
      self.X = [x]
      while len(self.X) < len(self.T):
        self.X.append( E(x, self.theta, self.lam, T[len(self.X)-1]) )
  
  def __repr__(self):
    return "AugmentedTrajectory(" + str(self.X) + "," + str(self.T) + ',' + str(self.theta) + ',' + str(self.lam) + ')'
  
  def __str__(self):
    return repr(self)
  
  def step_to_theta(self):
    distances = [dist_to_end_of_domain(self.X[i], self.T[i]) for i in xrange(len(self.X))]
    constraints = [-1]
    factor = 1
    n = 0
    while n+1 < len(self.X):
      constraints.append( distances[n+1] / factor )
      n += 1
      factor += self.theta^n
    return distances, constraints
  
  

































