# This file was *autogenerated* from the file trajectory.sage.
from sage.all_cmdline import *   # import sage library
_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_1en8 = RealNumber('1e-8')#the domain of g is [0,1]; the domain of f is [1,2]


def r(x):
  """round x to the nearest multiple of 1/2; 
  it is an error if x is an integer"""
  if x.is_integral():
    raise ValueError("x can't be an integer")
  return x.floor() + _sage_const_1 /_sage_const_2 

def r_set(x, which):
  """round x to the nearest 1/2 + 2*integer if which==g
  or 1/2 + 2*integer + 1 if which==f; this function is *undefined* if 
  x = 1/2 + 2n and we call it to round to 1/2 + 2n + 1 or vice versa"""
  if which=='g':
    return (x-_sage_const_1 /_sage_const_2 ).round('even') + _sage_const_1 /_sage_const_2 
  else:
    return (x-_sage_const_1 /_sage_const_2 ).round('odd') + _sage_const_1 /_sage_const_2 


def in_domain(x, which, extended=False):
  """decide whether a point lies in the domain given --- if extended==True, 
  this accepts an endpoint of the domain"""
  if not extended and x.is_integral():
    return False
  if extended and x.is_integral():
    return True
  if which == 'g':
    return x.floor()%_sage_const_2  == _sage_const_0 
  else:
    return x.floor()%_sage_const_2  == _sage_const_1 

def which_domain(x):
  """return which domain the point is in -- it shouldn't be an integer 
  for this to be well-defined"""
  if x.is_integral():
    return None
  else:
    return ('g' if x.floor()%_sage_const_2 ==_sage_const_0  else 'f')

def dist_to_end_of_domain(x, which):
  """return the distance to end of the domain that the point is in 
  (just ceiling(x) - x unless x is an integer)"""
  if not in_domain(x, which, extended=True):
    raise ValueError("Can't get distance to end of domain point isn't in")
  if x.is_integral():
    return (_sage_const_1  if ((which=='g') != (x%_sage_const_2 ==_sage_const_1 )) else _sage_const_0 ) #xor 
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
      self.X.append( E(self.X[-_sage_const_1 ], self.theta, self.lam, self.T[-_sage_const_1 ]) )
      if self.X[-_sage_const_1 ].is_integral():
        self.T.append( ('f' if self.X[-_sage_const_1 ]%_sage_const_2 ==_sage_const_0  else 'g') )
      else:
        self.T.append( which_domain(self.X[-_sage_const_1 ]) )
  
  @classmethod
  def from_full_trajectory(cls, x, theta, lam, T):
    depth = len(T)
    X = [-_sage_const_1  for ell in T]
    X[_sage_const_0 ] = x
    for i in xrange(depth-_sage_const_1 ):
      if not in_domain(X[i], T[i], extended=True):
        raise ValueError("x must be in domain")
      X[i+_sage_const_1 ] = E(X[i], theta, lam, T[i])
    ans = cls(x, theta, lam, depth, blank=True)
    ans.X = [p for p in X]
    ans.T = [t for t in T]
    return ans
  
  def __repr__(self):
    return "AugmentedTrajectory(" + str(self.X) + "," + str(self.T) + ',' + str(self.theta) + ',' + str(self.lam) + ')'
  
  def __str__(self):
    return repr(self)
  
  def theta_dist_cons(self):
    """return how far each entry in X is from a breakpoint, and then how 
    much each constrains theta (the minimum of constraints is exactly 
    how far the next breakpoint is"""
    distances = [dist_to_end_of_domain(self.X[i], self.T[i]) for i in xrange(len(self.X))]
    constraints = [-_sage_const_1 ]
    factor = _sage_const_1 
    n = _sage_const_0 
    while n+_sage_const_1  < len(self.X):
      constraints.append( distances[n+_sage_const_1 ] / factor )
      n += _sage_const_1 
      factor += self.lam**n
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
    for i in xrange(_sage_const_1 ,self.depth):
      if min_cons == None or cons[i] < min_cons:
        min_cons = cons[i]
        min_cons_i = i
    
    new_at = AugmentedTrajectory(self.X[_sage_const_0 ], self.theta, self.lam, self.depth, blank=True)
    
    #get the new theta
    new_at.theta = self.theta + min_cons
    
    #generate the new trajectory
    new_at.X = self.depth*[None]
    new_at.T = self.depth*[None]
    new_at.X[_sage_const_0 ] = self.X[_sage_const_0 ]
    new_at.T[_sage_const_0 ] = self.T[_sage_const_0 ]
    for i in xrange(_sage_const_0 ,self.depth-_sage_const_1 ):
      #get the next x value
      x = E(new_at.X[i], new_at.theta, new_at.lam, new_at.T[i])
      
      #if x is integral, we need to decide what domain it's in
      if x.is_integral(): 
        if i < min_cons_i-_sage_const_1 :
          if (x%_sage_const_2 ==_sage_const_0 ) != (self.T[i+_sage_const_1 ]=='f'): #xor
            t = self.T[i+_sage_const_1 ]
          else:
            raise ValueError("x shouldn't be integral and restricted before the smallest constraint")
        elif i == min_cons_i-_sage_const_1 :
          t = ('f' if self.T[min_cons_i] == 'g' else 'g')
        else:
          #set the letter to be the smaller domain
          t = ('f' if x%_sage_const_2 ==_sage_const_0  else 'g')
      else:
        if i==min_cons_i-_sage_const_1 :
          raise ValueError("x should be integral for the smallest constraint")
        t = which_domain(x)
      
      #append to the lists
      new_at.X[i+_sage_const_1 ] = x
      new_at.T[i+_sage_const_1 ] = t
    
    return new_at
  
  def constant_trajectory_theta_func(self):
    """Return a function which produces theta given lambda such that 
    the endpoint of the trajectory remains fixed; here the "endpoint"
    of the trajectory is the first place it hits an integer"""
    var('l j');
    n = min([i for i in xrange(_sage_const_1 ,self.depth) if self.X[i].is_integral()])
    f = self.X[n] - l**n*(self.X[_sage_const_0 ] - r_set(self.X[_sage_const_0 ],self.T[_sage_const_0 ]))
    #print self
    #print l^n*(self.X[0] - r_set(self.X[0],self.T[0]))
    #print f
    for i in xrange(_sage_const_0 , n-_sage_const_1 ):
      f -= l**(n-i-_sage_const_1 ) * (r_set(self.X[i],self.T[i]) - r_set(self.X[i+_sage_const_1 ],self.T[i+_sage_const_1 ]))
      #print self.X[i], self.T[i], r_set(self.X[i], self.T[i]), l^(n-i-1) * (r_set(self.X[i],self.T[i]) - r_set(self.X[i+1],self.T[i+1]))
    f -= r_set(self.X[n-_sage_const_1 ], self.T[n-_sage_const_1 ])
    #print r_set(self.X[n-1], self.T[n-1])
    f *= _sage_const_1 /sum(l**j,j,_sage_const_0 ,n-_sage_const_1 )
    #print f
    return (n,f)
    


def constant_trajectory_thetas(x, starting_map, depth):
  """return a list of rational functions of lambda which plot the 
  thetas which hold a given trajectory constant"""
  L = []
  at = AugmentedTrajectory(x, _sage_const_0 , _sage_const_2 , depth, starting_map=starting_map)
  at = at.next_theta_trajectory()
  while True:
    #while 0 in at.theta_dist_cons()[0][1:]:
    #  at = at.next_theta_trajectory()
    if at.theta > _sage_const_2 :
      break
    L.append( (at, at.constant_trajectory_theta_func()) )
    at = at.next_theta_trajectory()
  
  print "Got initial list"
  
  #get the bounds over which the functions are valid
  #first we must remove duplicates from the list
  i = _sage_const_0 
  while i < len(L)-_sage_const_1 :
    if (L[i][_sage_const_1 ][_sage_const_1 ]-L[i+_sage_const_1 ][_sage_const_1 ][_sage_const_1 ]).is_zero():
      del L[i+_sage_const_1 ]
    else:
      i += _sage_const_1 
  
  print "Deduped list"
  #print L
  
  #for each depth (increasing), find the functions of that depth
  #and scan through for everything they cut
  L = [[f,en,at,_sage_const_2 ,_sage_const_1 ] for at,(en,f) in L]
  for N in xrange(_sage_const_1 ,depth):
    #find the functions of this depth
    L_inds_this_depth = [i for i in xrange(len(L)) if L[i][_sage_const_1 ] == N]
    print "Trimming with ", len(L_inds_this_depth), " curves of depth ", N
    for i in L_inds_this_depth:
      f1,en1,at1,high1,low1 = L[i]
      #scan backwards, cutting off (and being cut off by) the functions
      #"seen" means that we've *we* have been cut off
      j=i-_sage_const_1 
      min_seen_depth = None
      while j >= _sage_const_0  and (min_seen_depth == None or min_seen_depth > en1):
        f2,en2,at2,high2,low2 = L[j]
        if min_seen_depth != None and en2 >= min_seen_depth:
          j -= _sage_const_1 
          continue
        #compute the intersection
        try:
          s = find_root(f1-f2,max(low1,low2)+_sage_const_1en8 ,_sage_const_2 )
        except RuntimeError:
          s = None
        if s == None or s < low1 + _sage_const_1en8 :
            min_seen_depth = en2
        else:
          if en1 < en2: # we are cutting them off
            L[j][-_sage_const_1 ] = s
          else:         # we are being cut off
            L[i][-_sage_const_1 ] = s
            low1 = s
            min_seen_depth = en2
        j -= _sage_const_1 
        continue
      #now do the same scan, except forwards
      j = i+_sage_const_1 
      min_seen_depth = None
      while j < len(L) and (min_seen_depth == None or min_seen_depth > en1):
        f2,en2,at2,high2,low2 = L[j]
        if min_seen_depth != None and en2 >= min_seen_depth:
          j += _sage_const_1 
          continue
        #compute the intersection
        try:
          s = find_root(f1-f2,max(low1,low2)+_sage_const_1en8 ,_sage_const_2 )
        except RuntimeError:
          s = None
        if s == None or s < low1 + _sage_const_1en8 :
            min_seen_depth = en2
        else:
          if en1 < en2: # we are cutting them off
            L[j][-_sage_const_1 ] = s
          else:         # we are being cut off
            L[i][-_sage_const_1 ] = s
            low1 = s
            min_seen_depth = en2
        j += _sage_const_1 
        continue
  
  return L
            
    




def theta_breakpoint_grid(x, starting_map, lam_step, depth):
  dat = []
  for lam in xsrange(_sage_const_1 , _sage_const_2 +lam_step, lam_step):
    at = AugmentedTrajectory(x, _sage_const_0 , lam, depth, starting_map=starting_map)
    while True:
      at = at.next_theta_trajectory()
      if at.theta > _sage_const_2 :
        break
      dists,cons = at.theta_dist_cons()
      if _sage_const_0  not in dists[_sage_const_1 :]:
        dat.append( (at.theta, lam) )
  return dat
      



def E_one_point_theta_interval(x, theta_int, lam, which):
  """give the output interval applying all theta values to the 
  point x"""
  #we can just apply it to the first point and the last point and that's it
  return [E(x,theta_int[_sage_const_0 ], lam, which=which), E(x,theta_int[_sage_const_1 ], lam, which=which)]
  

def E_point_interval_theta_interval(point_int, theta_int, lam, which):
  """give the output interval of applying the theta range to the 
  point range (the first theta to the first point, the last theta 
  to the last point, etc)"""
  return [E(point_int[_sage_const_0 ], theta_int[_sage_const_0 ], lam, which), \
          E(point_int[_sage_const_1 ], theta_int[_sage_const_1 ], lam, which)]

def last_point_in_domain(x, which):
  if not in_domain(x, which, extended=True):
    raise ValueError("x needs to be in the domain")
  #print "Getting last point in domain of ", x, " in ", which, ", integral: ", x.is_integral()
  if x.is_integral():
    if (x%_sage_const_2 ==_sage_const_0 ) != (which=='g'):
      return x
    else:
      return x+_sage_const_1 
  else:
    return ceil(x)

def split_interval(theta_int, point_int, which):
  """split the interval point_int into regions which 
  fall into the (extended) domain of the map "which".  Also split the 
  (linearly corresponding) theta interval, and return a list of 
  pairs of theta_int, point_int"""
  if in_domain(point_int[_sage_const_0 ], which):
    start_point = point_int[_sage_const_0 ]
  else:
    next_ok_point = ceil(point_int[_sage_const_0 ])
    if next_ok_point > point_int[_sage_const_1 ]:
      return []
    start_point = next_ok_point
  
  ans = []
  while True:
    #start point has the first point which lies in the domain
    #find the minimum of the endpoint of the interval and the endpoint 
    #of the domain, and chop off that interval, and repeat
    last_domain_point = last_point_in_domain(start_point, which)
    #print "Working on start point ", start_point, " with map ", which, ", last_domain: ", last_domain_point
    if last_domain_point >= point_int[_sage_const_1 ]:
      ans.append([start_point, point_int[_sage_const_1 ]])
      break
    else:
      ans.append([start_point, last_domain_point])
    #find the next interval
    next_domain_beginning = last_domain_point + _sage_const_1 
    if point_int[_sage_const_1 ] < next_domain_beginning:
      break
    start_point = next_domain_beginning
  
  #now scale to get the appropriate theta intervals
  if point_int[_sage_const_1 ] == point_int[_sage_const_0 ]:
    scaling_factor = _sage_const_0 
  else:
    scaling_factor = (theta_int[_sage_const_1 ]-theta_int[_sage_const_0 ])/(point_int[_sage_const_1 ]-point_int[_sage_const_0 ])
  theta_ans = [[theta_int[_sage_const_0 ] + scaling_factor*(a1-point_int[_sage_const_0 ]), \
                theta_int[_sage_const_0 ] + scaling_factor*(a2-point_int[_sage_const_0 ])] for a1,a2 in ans]
  
  return zip(theta_ans, ans)


def intervals_overlap(i1, i2):
  """return true iff the two intervals intersect"""
  if i1[_sage_const_1 ] < i2[_sage_const_0 ] or i2[_sage_const_1 ] < i1[_sage_const_0 ]:
    return False
  else:
    return True



def theta_intervals_for_trajectory(x, lam, desired_trajectory, extended=True):
  """return a set of intervals of theta for which the trajectory 
  will be admissible for parameters theta, lam with input point x.
  If extended=True, then it will make sure that the leftmost and 
  rightmost intervals stop at the real place they should, rather than 0 and 2"""
  if len(desired_trajectory) == _sage_const_0 :
    return [[_sage_const_0 ,_sage_const_2 ]]
  if not in_domain(x, desired_trajectory[_sage_const_0 ], extended=True):
    return []
  if len(desired_trajectory) == _sage_const_1 :
    return [[_sage_const_0 ,_sage_const_2 ]]
  if not extended:
    #compute the interval of points which are accessible from 
    #the starting point with any theta
    initial_interval = E_one_point_theta_interval(x, [_sage_const_0 ,_sage_const_2 ], lam, desired_trajectory[_sage_const_0 ])
    #split the interval into ones which agree with our next letter
    current_intervals = split_interval([_sage_const_0 ,_sage_const_2 ], initial_interval, desired_trajectory[_sage_const_1 ])
  else:
    initial_interval = E_one_point_theta_interval(x, [-_sage_const_1 ,_sage_const_3 ], lam, desired_trajectory[_sage_const_0 ])
    current_intervals = split_interval([-_sage_const_1 ,_sage_const_3 ], initial_interval, desired_trajectory[_sage_const_1 ])
    #go through the intervals and discard all those that don't overlap [0,2]
    current_intervals = [(theta_int, point_int) for (theta_int, point_int) in current_intervals if intervals_overlap(theta_int, [_sage_const_0 ,_sage_const_2 ])]
  i = _sage_const_1 
  while i < len(desired_trajectory)-_sage_const_1 :
    #print "Current intervals: ", current_intervals
    #go through all of our current intervals, act by them, and 
    #pick out the bits which agree with the desired trajectory
    new_intervals = []
    for theta_int, point_int in current_intervals:
      #print "Doing interval ", theta_int, point_int 
      #get the output interval
      oi = E_point_interval_theta_interval(point_int, theta_int, lam, desired_trajectory[i])
      #print "Got output ", oi
      #split it
      si = split_interval(theta_int, oi, desired_trajectory[i+_sage_const_1 ])
      #print "Split to ", si
      new_intervals.extend(si)
        
    current_intervals = new_intervals
    if extended:
      # we need to discard any outside intervals
      current_intervals = [(theta_int, point_int) for (theta_int, point_int) in current_intervals if intervals_overlap(theta_int, [_sage_const_0 ,_sage_const_2 ])]
    i += _sage_const_1 
  
  #the remaining intervals are the intervals we want
  return current_intervals



def find_feasible_regions(x1, T1, x2, T2, disregard_trivial_ints=True):
  """find the feasible region in the theta,lam parameter space
  such that the points x1 and x2 have trajectories T1 and T2"""
  
  #get the theta intervals
  theta_ints1 = theta_intervals_for_trajectory(x1, _sage_const_2 , T1, extended=True)
  theta_ints2 = theta_intervals_for_trajectory(x2, _sage_const_2 , T2, extended=True)
  
  print "Found theta intervals for lambda=2:"
  print theta_ints1
  print theta_ints2
  
  if disregard_trivial_ints:
    theta_ints1 = [(ti,pi) for (ti,pi) in theta_ints1 if ti[_sage_const_0 ] != ti[_sage_const_1 ]]
    theta_ints2 = [(ti,pi) for (ti,pi) in theta_ints2 if ti[_sage_const_0 ] != ti[_sage_const_1 ]]
  
  #build the augmented trajectories and theta as a function of lambda
  #for each theta interval
  theta_func_ints1 = []
  for ti,pi in theta_ints1:
    traj1 = AugmentedTrajectory.from_full_trajectory(x1, ti[_sage_const_0 ], _sage_const_2 , T1)
    traj2 = AugmentedTrajectory.from_full_trajectory(x1, ti[_sage_const_1 ], _sage_const_2 , T1)
    traj1_func = traj1.constant_trajectory_theta_func()
    traj2_func = traj2.constant_trajectory_theta_func()
    theta_func_ints1.append( [ti, traj1_func, traj2_func] )
  theta_func_ints2 = []
  for ti,pi in theta_ints2:
    traj1 = AugmentedTrajectory.from_full_trajectory(x2, ti[_sage_const_0 ], _sage_const_2 , T2)
    traj2 = AugmentedTrajectory.from_full_trajectory(x2, ti[_sage_const_1 ], _sage_const_2 , T2)
    traj1_func = traj1.constant_trajectory_theta_func()
    traj2_func = traj2.constant_trajectory_theta_func()
    theta_func_ints2.append( [ti, traj1_func, traj2_func] )
  
  #for each func pair, find where the functions intersect to get 
  # a minimal lambda value
  for L in theta_func_ints1, theta_func_ints2:
    for FP in L:
      ti, (en1,f1), (en2, f2) = FP
      if (f1-f2).is_zero():
        s = _sage_const_2 
      else:
        try:
          s = find_root(f1-f2,_sage_const_1 +_sage_const_1en8 ,_sage_const_2 )
        except RuntimeError:
          s = None
      if s == None or s < _sage_const_1en8 :
        min_lam = _sage_const_1 
      else:
        min_lam = s
      FP.append(min_lam)
  
  intersecting_pairs = []
  
  #now go through and find pairs of func ints which overlap
  for tfi1 in theta_func_ints1:
    ti1, (en11, f11), (en12, f12), lam1 = tfi1
    for tfi2 in theta_func_ints2:
      ti2, (en21, f21), (en22, f22), lam2 = tfi2
      if intervals_overlap(ti1, ti2):
        if ti1[_sage_const_1 ] <= ti2[_sage_const_1 ]:
          t,ell = ti1[_sage_const_1 ], _sage_const_2 
        else:
          t,ell = ti2[_sage_const_0 ], _sage_const_2 
        intersecting_pairs.append( (t,ell, tfi1, tfi2) )
        continue
      elif ti1[_sage_const_1 ] < ti2[_sage_const_0 ]:
        try:
          s = find_root(f12-f21, max(lam1, lam2), _sage_const_2 )
          t,ell = f12.substitute(l=s), s
        except RuntimeError:
          s = None
      else:
        try:
          s = find_root(f22-f11, max(lam1, lam2), _sage_const_2 )
          t,ell = f22.substitute(l=s), s
        except RuntimeError:
          s = None
      if s != None:
        intersecting_pairs.append( (t,ell, tfi1, tfi2) )
        
        
  
  return theta_func_ints1, theta_func_ints2, intersecting_pairs
      
        













