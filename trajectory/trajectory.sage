
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
  return x.ceil() - x

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
  def __init__(self, x, theta, lam, depth, starting_map=None, blank=False):
    """An augmented trajectory which records the points and the letters in 
    the trajectory.  If x is an integer, starting_map needs to be specified; 
    if it sees an integer along the way, it'll choose the *smallest* trajectory"""
    self.theta = theta
    self.lam = lam    
    self.depth = depth
    if blank:
      return
    if x.is_integral() and starting_map==None:
      raise ValueError("If x is an integer, the starting map needs to be specified")
    if starting_map != None and not in_domain(x,starting_map,extended=True):
      raise ValueError("x needs to be in the domain of the starting map")
    self.T = [(starting_map if starting_map != None else which_domain(x))]
    self.X = [x]
    while len(self.T) < self.depth:
      self.X.append( E(self.X[-1], self.theta, self.lam, self.T[-1]) )
      if self.X[-1].is_integral():
        self.T.append( ('f' if self.X[-1]%2==0 else 'g') )
      else:
        self.T.append( which_domain(self.X[-1]) )
  
  def __repr__(self):
    return "AugmentedTrajectory(" + str(self.X) + "," + str(self.T) + ',' + str(self.theta) + ',' + str(self.lam) + ')'
  
  def __str__(self):
    return repr(self)
  
  def theta_dist_cons(self):
    """return how far each entry in X is from a breakpoint, and then how 
    much each constrains theta (the minimum of constraints is exactly 
    how far the next breakpoint is"""
    distances = [dist_to_end_of_domain(self.X[i], self.T[i]) for i in xrange(len(self.X))]
    constraints = [-1]
    factor = 1
    n = 0
    while n+1 < len(self.X):
      constraints.append( distances[n+1] / factor )
      n += 1
      factor += self.lam^n
    return distances, constraints
  
  def next_theta_trajectory(self):
    """Returns an augmented trajectory of the same depth 
    which is the smallest possible next trajectory.  That is, get the 
    constraints on theta from theta_dist_cons, advance theta to that value 
    and swap the map on the active constraint, and then generate the 
    rest of the trajectory.  
    
    This is simple if the remaining trajectory doesn't hit an integer. 
    If it *does*, we should return the smallest possible; that is, 
    we should return a trajectory whose distance constraint 
    will be zero."""
    dists,cons = self.theta_dist_cons()
    
    #get the *first* minimum
    min_cons_ind = None
    min_cons = None
    for i in xrange(1,self.depth):
      if min_cons == None or cons[i] < min_cons:
        min_cons = cons[i]
        min_cons_i = i
    
    new_at = AugmentedTrajectory(self.X[0], self.theta, self.lam, self.depth, blank=True)
    
    #get the new theta
    new_at.theta = self.theta + min_cons
    
    #generate the new trajectory
    new_at.X = self.depth*[None]
    new_at.T = self.depth*[None]
    new_at.X[0] = self.X[0]
    new_at.T[0] = self.T[0]
    for i in xrange(0,self.depth-1):
      #get the next x value
      x = E(new_at.X[i], new_at.theta, new_at.lam, new_at.T[i])
      
      #if x is integral, we need to decide what domain it's in
      if x.is_integral(): 
        if i < min_cons_i-1:
          if (x%2==0) != (self.T[i+1]=='f'): #xor
            t = self.T[i+1]
          else:
            raise ValueError("x shouldn't be integral and restricted before the smallest constraint")
        elif i == min_cons_i-1:
          t = ('f' if self.T[min_cons_i] == 'g' else 'g')
        else:
          #set the letter to be the smaller domain
          t = ('f' if x%2==0 else 'g')
      else:
        if i==min_cons_i-1:
          raise ValueError("x should be integral for the smallest constraint")
        t = which_domain(x)
      
      #append to the lists
      new_at.X[i+1] = x
      new_at.T[i+1] = t
    
    return new_at
  
  def constant_trajectory_theta_func(self):
    """Return a function which produces theta given lambda such that 
    the endpoint of the trajectory remains fixed; here the "endpoint"
    of the trajectory is the first place it hits an integer"""
    var('l j');
    n = min([i for i in xrange(1,self.depth) if self.X[i].is_integral()])
    f = self.X[n] - l^n*(self.X[0] - r_set(self.X[0],self.T[0]))
    #print self
    #print l^n*(self.X[0] - r_set(self.X[0],self.T[0]))
    #print f
    for i in xrange(0, n-1):
      f -= l^(n-i-1) * (r_set(self.X[i],self.T[i]) - r_set(self.X[i+1],self.T[i+1]))
      #print self.X[i], self.T[i], r_set(self.X[i], self.T[i]), l^(n-i-1) * (r_set(self.X[i],self.T[i]) - r_set(self.X[i+1],self.T[i+1]))
    f -= r_set(self.X[n-1], self.T[n-1])
    #print r_set(self.X[n-1], self.T[n-1])
    f *= 1/sum(l^j,j,0,n-1)
    #print f
    return (n,f)
    


def constant_trajectory_thetas(x, starting_map, depth):
  """return a list of rational functions of lambda which plot the 
  thetas which hold a given trajectory constant"""
  L = []
  at = AugmentedTrajectory(x, 0, 2, depth, starting_map=starting_map)
  at = at.next_theta_trajectory()
  while True:
    #while 0 in at.theta_dist_cons()[0][1:]:
    #  at = at.next_theta_trajectory()
    if at.theta > 2:
      break
    L.append( (at, at.constant_trajectory_theta_func()) )
    at = at.next_theta_trajectory()
  
  print "Got initial list"
  
  #get the bounds over which the functions are valid
  #first we must remove duplicates from the list
  i = 0
  while i < len(L)-1:
    if (L[i][1][1]-L[i+1][1][1]).is_zero():
      del L[i+1]
    else:
      i += 1
  
  print "Deduped list"
  #print L
  
  #for each depth (increasing), find the functions of that depth
  #and scan through for everything they cut
  L = [[f,en,at,2,1] for at,(en,f) in L]
  for N in xrange(1,depth):
    #find the functions of this depth
    L_inds_this_depth = [i for i in xrange(len(L)) if L[i][1] == N]
    print "Trimming with ", len(L_inds_this_depth), " curves of depth ", N
    for i in L_inds_this_depth:
      f1,en1,at1,high1,low1 = L[i]
      #scan backwards, cutting off (and being cut off by) the functions
      #"seen" means that we've *we* have been cut off (or the intersection is None or 1)
      j=i-1
      min_seen_depth = None
      while j >= 0 and (min_seen_depth == None or min_seen_depth > 1):
        f2,en2,at2,high2,low2 = L[j]
        if min_seen_depth != None and en2 >= min_seen_depth:
          j -= 1
          continue
        #compute the intersection
        try:
          s = find_root(f1-f2,max(low1,low2),2)
        except RuntimeError:
          s = None
        if s == None or abs(s-1) < 1e-8:
            min_seen_depth = en2
        else:
          if en1 < en2: # we are cutting them off
            L[j][-1] = s
          else:         # we are being cut off
            L[i][-1] = s
            low1 = s
            min_seen_depth = en2
        j -= 1
        continue
      #now do the same scan, except forwards
      j = i+1
      min_seen_depth = None
      while j < len(L) and (min_seen_depth == None or min_seen_depth > 1):
        f2,en2,at2,high2,low2 = L[j]
        if min_seen_depth != None and en2 >= min_seen_depth:
          j += 1
          continue
        #compute the intersection
        try:
          s = find_root(f1-f2,max(low1,low2),2)
        except RuntimeError:
          s = None
        if s == None or abs(s-1) < 1e-8:
            min_seen_depth = en2
        else:
          if en1 < en2: # we are cutting them off
            L[j][-1] = s
          else:         # we are being cut off
            L[i][-1] = s
            low1 = s
            min_seen_depth = en2
        j += 1
        continue
  
  return L
            
    




def theta_breakpoint_grid(x, starting_map, lam_step, depth):
  dat = []
  for lam in xsrange(1, 2+lam_step, lam_step):
    at = AugmentedTrajectory(x, 0, lam, depth, starting_map=starting_map)
    while True:
      at = at.next_theta_trajectory()
      if at.theta > 2:
        break
      dists,cons = at.theta_dist_cons()
      if 0 not in dists[1:]:
        dat.append( (at.theta, lam) )
  return dat
      






























