import sage.all as SAGE 

def r(x):
  """round x to the nearest multiple of 1/2; 
  it is an error if x is an integer"""
  if x.is_integral():
    raise ValueError("x can't be an integer")
  return x.floor() + 1/2

def r_set(x, which):
  """round x to the nearest 1/2 + 2*integer if which==g
  or 1/2 + 2*integer + 1 if which==f"""
  pass

def trajectory(x, theta, lam):
  pass

def in_domain(x, which, extended=False):
  if not extended and x.is_integral():
    return False
  if extended and x.is_integral():
    return True
  if which == 'g':
    return x.floor()%2 == 0
  else:
    return x.floor()%2 == 1

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
  def __init__(self, x, T):
    """we can accept either x as a number to start with
    or a list of all the images right away; if we get a 
    list, we don't check the sanity"""
    if isinstance(x, list):
      self.X = [ell for ell in x] #duplicate the list
    else:
      
  

  