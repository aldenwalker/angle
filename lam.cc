#include <cmath>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <set>

#include "lam.h"

// every function here assumes that all angles are normalized to be 
// in [0,2pi)
double normalize_angle(double x) {
  double n = fmod(x, 2*PI);
  if (n<0) {
    return n+2*PI;
  } else {
    return n;
  }
}


//return the difference between the angles which is smaller than pi
double angle_diff(double a, double b) {
  double d = b-a;
  double r = fmod(d, 2*PI);
  if (r < 0) {
    r += 2*PI;
  }
  if (r > PI) {
    r = 2*PI - r;
  }
  return r;
}
 
//return the distance from a to b
double angle_dist(double a, double b) {
  if (a < b) {
    return b-a;
  } else {
    return (b+2*PI) - a;
  }
}
 
 //return the angle in the middle of the two angles (closer than PI/2)
 double average_angle(double x, double y) {
  double A = (y>x ? y : x);
  double a = (y>x ? x : y);
  double av = (a+A)/2.0;
  if (A-a > PI) {
    return normalize_angle(av + PI);
  }
  return av;
}
   
bool cyclically_ordered(double x, double y, double z) {
  if (y < x) {
    y += 2*PI;
  }
  if (z < x) {
    z += 2*PI;
  }
  return y < z;
}

bool weak_cyclically_ordered(double x, double y, double z) {
  if (fabs(x-y) < 1e-8 || fabs(y-z) < 1e-8 || fabs(x-z) < 1e-8) {
    return true;
  }
  if (y < x) {
    y += 2*PI;
  }
  if (z < x) {
    z += 2*PI;
  }
  if (fabs(x-y) < 1e-8 || fabs(y-z) < 1e-8 || fabs(x-z) < 1e-8) {
    return true;
  }
  return y <= z;
}

bool intervals_intersect(double a1, double a2, double b1, double b2) {
  //if either b endpoint lives in the a interval, they intersect
  if (cyclically_ordered(a1, b1, a2) || cyclically_ordered(a1, b2, a2)) {
    return true;
  }
  //if not, then the cyclic order of b1, b2 is all that matters
  return cyclically_ordered(b1, a1, b2);
}



void points_on_circle(std::vector<Point2d<float> >& points, int num_per_radian, 
                      double cx, double cy, double r, double a, double a_extent) {
  int num_points = int(fabs(num_per_radian*a_extent*(int(r)+1))) + 2;
  double step = a_extent / (double)(num_points-1);
  points.resize(num_points);
  for (int i=0; i<num_points; ++i) {
    double ang = a + step * (double)i;
    points[i] = Point2d<float>( cx + r*cos(ang), cy + r*sin(ang) );
  }
}


/*************************************************************************/
Graph::Graph() {
  num_edges = num_verts = 0;
  edges.resize(0);
}

Graph Graph::from_edges(const std::vector< std::vector<int> >& e) {
  Graph G;
  G.edges.resize(0);
  G.num_verts = e.size();
  G.edges_starting_with.resize(e.size());
  G.edges_ending_with.resize(e.size());
  for (int i=0; i<(int)e.size(); ++i) {
    for (int j=0; j<(int)e[i].size(); ++j) {
      G.edges_starting_with[i].push_back(G.edges.size());
      G.edges_ending_with[e[i][j]].push_back(G.edges.size());
      G.edges.push_back( std::make_pair(i, e[i][j]) );
    }
  }
  G.num_edges = G.edges.size();
  return G;   
}


Graph Graph::from_edges_with_weights(const std::vector< std::vector<int> >& e,
                                     const std::vector<int>& weights) {
  Graph G = Graph::from_edges(e);
  int old_num_edges = G.num_edges;
  //std::cout << "Before multiplying edges: " << G.mathematica_string() << "\n";
  for (int i=0; i<old_num_edges; ++i) {
    G.multiply_edge(i, weights[G.edges[i].first]);
  }
  return G;
}


void Graph::multiply_edge(int ei, int m) {
  if (m == 1) return;
  //remove the end from where it ends
  std::vector<int>::iterator it;
  int dv = edges[ei].second;
  it = std::find( edges_ending_with[dv].begin(), edges_ending_with[dv].end(), ei );
  if (it == edges_ending_with[dv].end()) {
    std::cout << "Couldn't find edge in dest vert?\n";
  }
  edges_ending_with[dv].erase(it);
  
  int v_to_add = m-1;
  int last_edge_added = ei;
  for (int i=0; i<v_to_add; ++i) {
    int this_vi = edges_starting_with.size();
    int this_ei = edges.size();
    edges[last_edge_added].second = this_vi;
    edges_ending_with.push_back(std::vector<int>(1,last_edge_added));
    edges_starting_with.push_back(std::vector<int>(1,this_ei));
    last_edge_added = this_ei;
    edges.push_back(std::make_pair(this_vi, -1));
  }
  edges[last_edge_added].second = dv;
  edges_ending_with[dv].push_back(last_edge_added);
  num_edges = edges.size();
  num_verts = edges_starting_with.size();
}

double Graph::approximate_leading_eigenvalue(int iters) const {
  std::vector<double> weights(num_verts, 1);
  std::vector<double> old_weights(num_verts);
  if (num_verts == 0) return 0;
  for (int i=0; i<iters+1; ++i) {
    old_weights = weights;
    for (int j=0; j<num_verts; ++j) weights[j] = 0;
    for (int j=0; j<num_edges; ++j) {
      weights[edges[j].second] += old_weights[edges[j].first];
    }
    if (i==iters) {
      //std::cout << "[";
      //for (int j=0; j<num_verts; ++j) {
      //  std::cout << "(" << weights[j] << "," << old_weights[j] << ")";
      //}
      //std::cout << "]\n";
      for (int j=0; j<num_verts; ++j) {
        if (fabs(weights[j])+1e-5 > 1.0/double(num_verts)) {
          return weights[j] / old_weights[j];
        }
      }
    }
    double csum = 0;
    for (int j=0; j<num_verts; ++j) {
      csum += weights[j];
    }
    if (fabs(csum) < 1e-8) return 0;
    for (int j=0; j<num_verts; ++j) {
      weights[j] /= csum;
    }
  }
  std::cout << "Didn't find eigenvalue?\n";
  return -1;
}


void Graph::find_cycle_components(std::vector<std::vector<std::deque<int> > >& comps) {

}


int Graph::rank_if_undirected() const {
  int nv = num_verts;
  if (nv==0 || num_edges == 0) {
    return 0;
  }
  return num_edges - num_verts + 1;
}

std::vector<Graph> Graph::connected_components() const {
  return std::vector<Graph>();
}

std::string Graph::mathematica_string() const {
  std::stringstream ans;
  ans.str("");
  ans << "Graph[{";
  if (edges.size() > 0) {
    ans << "DirectedEdge[" << edges[0].first << "," << edges[0].second << "]";
  }
  for (int i=1; i<(int)edges.size(); ++i) {
    ans << ", DirectedEdge[" << edges[i].first << "," << edges[i].second << "]";
  }
  ans << "}]";
  return ans.str();
}


std::ostream& operator<<(std::ostream& os, const Graph& G) {
  os << "Vertices (" << G.edges_starting_with.size() << ") :\n";
  for (int i=0; i<(int)G.edges_starting_with.size(); ++i) {
    os << i << ": [";
    for (int j=0; j<(int)G.edges_starting_with[i].size(); ++j) {
      os << " " << G.edges[G.edges_starting_with[i][j]].second;
    }
    os << "]\n";
  }
  return os;
}

/*************************************************************************/


Leaf::Leaf() {
  x = y = 0;
  word.resize(0);
}

Leaf::Leaf(double x, double y) {
  this->x = x;
  this->y = y;
  word.resize(0);
}

void Leaf::get_circle_data(double& centerx, double& centery, double& radius, double& a1, double& a_extent) const {
  double alpha = angle_diff(x, y);
  double beta = average_angle(x, y);
  double gamma = PI/2.0 - alpha/2.0;
  if (fabs(alpha - PI) < 1e-8) { 
    a1 = x;
    a_extent = y;
    radius = -1;
    return;
  }
  double to_center = 1.0/cos(alpha/2.0);
  centerx = to_center*cos(beta);
  centery = to_center*sin(beta);
  radius = to_center*sin(alpha/2.0);
  if (radius > 40) {
    a1 = x;
    a_extent = y;
    radius = -1;
    return;
  }
  double middle_angle = beta - PI;
  a1 = middle_angle - gamma;
  a1 = normalize_angle(a1);
  a_extent = 2*gamma;
}


void Leaf::polygon_points(std::vector<Point2d<float> >& points, int num_per_radian) const {
  points.resize(0);
  double cx, cy, r, a, a_ext;
  //std::cout << "Finding points along the leaf: " << *this << "\n";
  get_circle_data(cx, cy, r, a, a_ext);
  //std::cout << "Got the circle data: (" << cx << "," << cy << ") " << r << " " << a << " " << a_ext << "\n";
  if (r < 0) {
    //std::cout << "Drawing a straight leaf\n";
    points.push_back( Point2d<float>(cos(x), sin(x)) );
    points.push_back( Point2d<float>(cos(y), sin(y)) );
  } else if ( fabs((cx + r*cos(a)) - cos(x)) < 1e-4 && fabs((cy + r*sin(a)) - sin(x)) < 1e-4 ) {
    //std::cout << "The leaf goes positive angles\n";
    points_on_circle(points, num_per_radian, cx, cy, r, a, a_ext);
  } else {
    //std::cout << "The leaf goes negative angles\n";
    points_on_circle(points, num_per_radian, cx, cy, r, a + a_ext, -a_ext);
  }
  //std::cout << "Got points: \n";
  //for (int i=0; i<(int)points.size(); ++i) {
  //  std::cout << i << ": " << points[i] << "\n";
  //}
}


bool Leaf::operator<(const Leaf& other) const {
  return angle_diff(x,y) < angle_diff(other.x, other.y);
}

std::ostream& operator<<(std::ostream& os, const std::deque<char>& v) {
  for (int i=0; i<(int)v.size(); ++i) {
    os << (int)v[i];
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, const Leaf& ell) {
  return os << "Leaf(" << ell.x << "," << ell.y << "," << ell.word << ")"; 
}


/**************************************************************************/
ThickLeaf::ThickLeaf(double p1x, double p1y, double p2x, double p2y) {
  this->p1x = p1x;
  this->p1y = p1y;
  this->p2x = p2x;
  this->p2y = p2y;
  this->word.resize(0);
}

ThickLeaf::ThickLeaf(double p1x, double p1y, double p2x, double p2y, const std::deque<char>& w) {
  this->p1x = p1x;
  this->p1y = p1y;
  this->p2x = p2x;
  this->p2y = p2y;
  this->word = w;
}


void ThickLeaf::polygon_points(std::vector<Point2d<float> >& points, int num_per_radian) const {
  Leaf ell1(p1x, p1y);
  Leaf ell2(p2x, p2y);
  
  double gap1 = angle_dist(p2y, p1x);
  double gap2 = angle_dist(p1y, p2x);
  
  //std::cout << "Finding points along: " << *this << "\n";
  
  points.resize(0);
  std::vector<Point2d<float> > new_points;
  ell1.polygon_points(new_points, num_per_radian);
  points.insert(points.end(), new_points.begin(), new_points.end());
  
  points_on_circle(new_points, num_per_radian, 0, 0, 1, p1y, gap2);
  points.insert(points.end(), new_points.begin(), new_points.end());
  
  ell2.polygon_points(new_points, num_per_radian);
  points.insert(points.end(), new_points.begin(), new_points.end());
  
  points_on_circle(new_points,  num_per_radian, 0, 0, 1, p2y, gap1);
  points.insert(points.end(), new_points.begin(), new_points.end());
}



bool ThickLeaf::is_disjoint(const Leaf& ell) const {
  return (cyclically_ordered(p1x, ell.x, p1y) && cyclically_ordered(p1x, ell.y, p1y)) ||
         (cyclically_ordered(p2x, ell.x, p2y) && cyclically_ordered(p2x, ell.y, p2y));
}


void ThickLeaf::split_by_leaves(std::vector<ThickLeaf>& result, const Leaf& ell) const {
  result.resize(0);
  std::deque<char> this_w = word;
  //there are several possibilities: there are 4 cyclic positions for 
  //ell.x, and then there are several for y in each case.  This involves 
  //lots of cyclic order checks.
  
  // (1) ell crosses the entire leaf:
  if ( (cyclically_ordered(p1x, ell.x, p1y) && cyclically_ordered(p2x, ell.y, p2y)) ||
       (cyclically_ordered(p1x, ell.y, p1y) && cyclically_ordered(p2x, ell.x, p2y)) ) {
    return;
  }
  // (2) ell is disjoint from the leaf
  if ( (cyclically_ordered(p1x, ell.x, p1y) && cyclically_ordered(p1x, ell.y, p1y)) ||
       (cyclically_ordered(p2x, ell.x, p2y) && cyclically_ordered(p2x, ell.y, p2y)) ) {
    result.push_back(*this);
    return;
  }
  // (3) ell lies between p1y and p2x (there are two possibilities on the order of ell)
  if (cyclically_ordered(p1y, ell.x, p2x) && cyclically_ordered(p1y, ell.y, p2x) ) {
    if (cyclically_ordered(p1y, ell.x, ell.y)) {
      result.push_back( ThickLeaf(p1x, p1y, ell.x, p2y, this_w) );
      result.push_back( ThickLeaf(p1x, ell.y, p2x, p2y, this_w) );
    } else {
      result.push_back( ThickLeaf(p1x, p1y, ell.y, p2y, this_w) );
      result.push_back( ThickLeaf(p1x, ell.x, p2x, p2y, this_w) );
    }
    return;
  }
  // (4) ell lies between p2y and p1x (again, two possibilities on the order of ell)
  if (cyclically_ordered(p2y, ell.x, p1x) && cyclically_ordered(p2y, ell.y, p1x)) {
    if (cyclically_ordered(p2y, ell.x, ell.y)) {
      result.push_back( ThickLeaf(p2x, p2y, ell.x, p1y, this_w) );
      result.push_back( ThickLeaf(p2x, ell.y, p1x, p1y, this_w) );
    } else {
      result.push_back( ThickLeaf(p2x, p2y, ell.y, p1y, this_w) );
      result.push_back( ThickLeaf(p2x, ell.x, p1x, p1y, this_w) );
    }
    return;
  }
  //(5) ell lies inside the leaf (two directions)
  if (cyclically_ordered(p2y, ell.x, p1x) && cyclically_ordered(p1y, ell.y, p2x)) {
    result.push_back( ThickLeaf(p1x, p1y, ell.y, ell.x, this_w) );
    result.push_back( ThickLeaf(p2x, p2y, ell.x, ell.y, this_w) );
    return;
  }
  if (cyclically_ordered(p2y, ell.y, p1x) && cyclically_ordered(p1y, ell.x, p2x)) {
    result.push_back( ThickLeaf(p1x, p1y, ell.x, ell.y, this_w) );
    result.push_back( ThickLeaf(p2x, p2y, ell.y, ell.x, this_w) );
    return;
  }
  //(6) ell has one endpoint inside and one outside (four possibilities, 
  //    and for each possibility, there are two directions).  For this section,
  //    we use the fact that we've eliminated the above possibilities.  We'll
  //    branch on the options for ell.x
  if (cyclically_ordered(p2y, ell.x, p1x)) {
    if (cyclically_ordered(p1x, ell.y, p1y)) {
      result.push_back( ThickLeaf(ell.x, p1y, p2x, p2y, this_w) );
    } else {
      result.push_back( ThickLeaf(p1x, p1y, p2x, ell.x, this_w) );
    }
    return;
  }
  if (cyclically_ordered(p1x, ell.x, p1y)) {
    if (cyclically_ordered(p2y, ell.y, p1x)) {
      result.push_back( ThickLeaf(ell.y, p1y, p2x, p2y, this_w) );
    } else {
      result.push_back( ThickLeaf(p1x, ell.y, p2x, p2y, this_w) );
    }
    return;
  }
  if (cyclically_ordered(p1y, ell.x, p2x)) {
    if (cyclically_ordered(p1x, ell.y, p1y)) {
      result.push_back( ThickLeaf(p1x, ell.x, p2x, p2y, this_w) );
    } else {
      result.push_back( ThickLeaf(p1x, p1y, ell.x, p2y, this_w) );
    }
    return;
  }
  if (cyclically_ordered(p2x, ell.x, p2y)) {
    if (cyclically_ordered(p1y, ell.y, p2x)) {
      result.push_back( ThickLeaf(p1x, p1y, ell.y, p2y, this_w) );
    } else {
      result.push_back( ThickLeaf(p1x, p1y, p2x, ell.y, this_w) );
    }
    return;
  }
  // if we get here, something is wrong
  std::cout << "Couldn't figure out how to cut leaf?\n";
}

void ThickLeaf::split_by_leaves(std::vector<ThickLeaf>& result, const std::vector<Leaf>& leaves) const {
  result.resize(0);
  result.push_back(*this);
  for (int i=0; i<(int)leaves.size(); ++i) {
    int j=0;
    while (j < (int)result.size()) {
      if (!result[j].is_disjoint(leaves[i])) {
        std::vector<ThickLeaf> split_leaves;
        result[j].split_by_leaves(split_leaves, leaves[i]);
        result.erase(result.begin() + j);
        result.insert(result.begin()+j, split_leaves.begin(), split_leaves.end());
        j = j + split_leaves.size();
      } else {
        ++j;
      }
    }
  }
}



bool ThickLeaf::contains(const ThickLeaf& other) const {
  return (weak_cyclically_ordered(p2y, other.p2y, p1x) &&
          weak_cyclically_ordered(p2y, other.p1x, p1x) &&
          weak_cyclically_ordered(p1y, other.p1y, p2x) &&
          weak_cyclically_ordered(p1y, other.p2x, p2x))   ||
         (weak_cyclically_ordered(p2y, other.p1y, p1x) &&
          weak_cyclically_ordered(p2y, other.p2x, p1x) &&
          weak_cyclically_ordered(p1y, other.p2y, p2x) &&
          weak_cyclically_ordered(p1y, other.p1x, p2x));
}

bool ThickLeaf::contains(const Leaf& other) const {
  return (weak_cyclically_ordered(p2y, other.x, p1x) && 
          weak_cyclically_ordered(p1y, other.y, p2x)) ||
         (weak_cyclically_ordered(p2y, other.y, p1x) &&
          weak_cyclically_ordered(p1y, other.x, p2x));
}

bool ThickLeaf::intersects(const ThickLeaf& other) const {
  //to check if they intersect, it suffices to check whether the ending 
  //intervals intersect on both sides
  //std::cout << "Checking if there is an intersection:\n";
  //std::cout << *this << "\n" << other << "\n";
  bool ans =  (intervals_intersect( p2y, p1x, other.p2y, other.p1x ) && 
          intervals_intersect( p1y, p2x, other.p1y, other.p2x))    ||
         (intervals_intersect( p2y, p1x, other.p1y, other.p2x ) && 
          intervals_intersect( p1y, p2x, other.p2y, other.p1x)) ;
  //std::cout << (ans ? "Yes" : "No") << "\n";
  return ans;
}


bool ThickLeaf::operator<(const ThickLeaf& other) const {
  return angle_dist(p2y, p1x) + angle_dist(p1y, p2x) <
         angle_dist(other.p2y, other.p1x) + angle_dist(other.p1y, other.p2x);
}

void ThickLeaf::find_inclusions(std::vector<std::vector<std::pair<bool,int> > >& inclusions, 
                                const std::vector<ThickLeaf>& A,
                                const std::vector<ThickLeaf>& B,
                                int verbose) {
  if (A.size() != B.size()) {
    std::cout << "inclusions lists aren't the same size?\n";
    return;
  }
  inclusions.resize(A.size());
  for (int i=0; i<(int)A.size(); ++i) {
    inclusions[i].resize(0);
  }
  if (verbose>0) {
    std::cout << "Finding inclusions\n";
  }
  for (int i=0; i<(int)A.size(); ++i) {
    for (int j=0; j<(int)B.size(); ++j) {
      if (A[i].contains(B[j])) {
        if (verbose>0) {
          std::cout << "Contains:\n";
          std::cout << A[i] << "\n" << B[j] << "\n";
        }
        inclusions[i].push_back( std::make_pair(true, j) );
      } else if (A[i].intersects(B[j])) {
        if (verbose>0) {
          std::cout << "Intersects:\n";
          std::cout << A[i] << "\n" << B[j] << "\n";
        }
        inclusions[i].push_back( std::make_pair(false, j) );
      }
    }
  }
}


std::ostream& operator<<(std::ostream& os, const ThickLeaf& ell) {
  return os << "ThickLeaf(" << ell.p1x << "," << ell.p1y << "," << ell.p2x << "," << ell.p2y << "," << ell.word << ")"; 
}



/***************************************************************************/

PartiallyDefinedMap::PartiallyDefinedMap() {
  image_x = image_y = domain_x = domain_y = 0;
}

PartiallyDefinedMap::PartiallyDefinedMap(double x1, double y1, double x2, double y2) {
  domain_x = x1;
  domain_y = y1;
  double domain_size = y1 - x1;
  if (domain_size < 0) {
    domain_size = 2*PI + domain_size;
  }
  //std::cout << "Computed domain size of " << domain_size << " from " << domain_x << " " << domain_y << "\n";
  image_x = x2;
  image_y = y2;
  double image_size = y2-x2;
  if (image_size < 0) {
    image_size = 2*PI + image_size;
  }
  //std::cout << "Computed image size of " << image_size << " from " << image_x << " " << image_y << "\n";
  scale_factor = image_size/domain_size;
  //std::cout << "Got scale factor of: " << scale_factor << "\n";
}

PartiallyDefinedMap::PartiallyDefinedMap(double lambda, double theta, bool f) {
  if (f) { 
    domain_x = normalize_angle( (PI+theta) - (lambda*PI/2.0) );
    domain_y = normalize_angle( (PI+theta) + (lambda*PI/2.0) );
    image_x = PI/2.0;
    image_y = 3*PI/2.0;
  } else {
    domain_x = normalize_angle( theta - (lambda*PI/2.0) );
    domain_y = normalize_angle( theta + (lambda*PI/2.0) );
    image_x = 3*PI/2.0;
    image_y = PI/2.0;
  }
  scale_factor = 1.0/lambda;
}

PartiallyDefinedMap PartiallyDefinedMap::inverse() const {
  return PartiallyDefinedMap(image_x, image_y, domain_x, domain_y);
}

bool PartiallyDefinedMap::leaf_in_domain(const Leaf& ell) const {
  return weak_cyclically_ordered(domain_x, ell.x, domain_y) && 
         weak_cyclically_ordered(domain_x, ell.y, domain_y);
}

bool PartiallyDefinedMap::leaf_in_domain(const ThickLeaf& ell) const {
  return weak_cyclically_ordered(domain_x, ell.p1x, domain_y) && 
         weak_cyclically_ordered(domain_x, ell.p1y, domain_y) && 
         weak_cyclically_ordered(domain_x, ell.p2x, domain_y) && 
         weak_cyclically_ordered(domain_x, ell.p2y, domain_y);
}


double PartiallyDefinedMap::act_on_point(double x) const {
  //std::cout << "Acting by " << *this << " on " << x << "\n";
  if (x < domain_x - 1e-3) {
    x += 2*PI;
  }
  //std::cout << "Difference = " << x-domain_x << "\n";
  double y = (x-domain_x)*scale_factor;
  //std::cout << "Scaled by " << scale_factor << " to " << y << "\n";
  double ans = normalize_angle(image_x + y);
  //std::cout << "Got " << ans << "\n";
  return ans;
}

Leaf PartiallyDefinedMap::act_on_leaf(const Leaf& ell) const {
  //std::cout << "Acting on " << ell << "\n";
  Leaf ans(act_on_point(ell.x), act_on_point(ell.y));
  //std::cout << "Got " << ans << "\n";
  return ans;
}

ThickLeaf PartiallyDefinedMap::act_on_leaf(const ThickLeaf& ell) const {
  ThickLeaf ans(act_on_point(ell.p1x), act_on_point(ell.p1y), 
                act_on_point(ell.p2x), act_on_point(ell.p2y));
  //std::cout << "Acting on leaf: " << ell << "\n";
  //std::cout << "Got: " << ans << "\n";
  return ans;
}

std::ostream& operator<<(std::ostream& os, const PartiallyDefinedMap& f) {
  return os << "PDM(" << f.domain_x << " " << f.domain_y << " -> " << f.image_x << " " << f.image_y << ")";
}


/****************************************************************************/


Lamination::Lamination() {
  lambda = theta = 0;
  f = PartiallyDefinedMap();
  g = PartiallyDefinedMap();
}

Lamination::Lamination(double lambda, double theta) {
  this->lambda = lambda;
  this->theta = theta;
  f = PartiallyDefinedMap(lambda, theta, true);
  g = PartiallyDefinedMap(lambda, theta, false);
}
  
  
bool Lamination::check_can_act(const Leaf& ell, int which_map) {
  return (which_map==0 ? f.leaf_in_domain(ell) : g.leaf_in_domain(ell));
}

Leaf Lamination::act_on_leaf(const Leaf& ell, int which_map) {
  Leaf result = (which_map==0 ? f.act_on_leaf(ell) : g.act_on_leaf(ell));
  result.word = ell.word;
  result.word.push_front(which_map);
  return result;
}

void Lamination::compute_lam_type(int depth, LamType& ellt, int& difficulty) {
  if (lambda < 1 || 2 < lambda || theta < 0 || TWOPI < theta) {
    ellt = NONE;
    difficulty = 0;
    return;
  }
  difficulty = 0;
  std::vector<Leaf> stack(0);
  stack.push_back( Leaf(PI/2.0, 3*PI/2.0) );
  int done_depth = 0;
  double epsilon = (2*PI - (PI*lambda))*0.99999999;
  
  bool can_either_act = check_can_act(stack.back(), 0);
  if (!can_either_act) {
    ellt = NONE;
    return;
  }  
  
  while (stack.size() > 0) {
    Leaf ell = stack.back();
    stack.pop_back();
    ++difficulty;
    bool can_f_act = check_can_act(ell, 0);
    bool can_g_act = check_can_act(ell, 1);
    if (angle_diff(ell.x, ell.y) < epsilon) {
      continue;
    }
    if ((int)ell.word.size() == depth) {
      break;
    }
    if (can_f_act) stack.push_back(act_on_leaf(ell, 0));
    if (can_g_act) stack.push_back(act_on_leaf(ell, 1));
    if ((int)stack.back().word.size() > done_depth) done_depth = stack.back().word.size();
  }
  if (done_depth < depth) ellt = PROPER;
  else ellt = CUT_POINT;
}
    
   
void Lamination::compute_lam_type_with_words(int depth, 
                                             LamType& ellt, 
                                             int& difficulty, 
                                             std::vector<Leaf>& viable_leaves) {
  if (lambda < 1 || 2 < lambda || theta < 0 || TWOPI < theta) {
    ellt = NONE;
    difficulty = 0;
    return;
  }
  difficulty = 0;
  viable_leaves.resize(0);
  std::vector<Leaf> stack(0);
  stack.push_back( Leaf(PI/2.0, 3*PI/2.0) );
  double epsilon = (2*PI - (PI*lambda))*0.99999999;
  
  bool can_either_act = check_can_act(stack.back(), 0);
  if (!can_either_act) {
    ellt = NONE;
    return;
  }  
  
  while (stack.size() > 0) {
    Leaf ell = stack.back();
    stack.pop_back();
    ++difficulty;
    if (angle_diff(ell.x, ell.y) < epsilon) continue;
    if ((int)ell.word.size() == depth) {
      viable_leaves.push_back(ell);
      continue;
    }
    bool can_f_act = check_can_act(ell, 0);
    bool can_g_act = check_can_act(ell, 1);
    if (can_f_act) stack.push_back(act_on_leaf(ell, 0));
    if (can_g_act) stack.push_back(act_on_leaf(ell, 1));
  }
  if (viable_leaves.size() == 0) ellt = PROPER;
  else ellt = CUT_POINT;
}
  

void Lamination::compute_backwards(int depth, 
                                   ThickLeaf& initial, 
                                   std::vector<ThickLeaf>& viable_initials,
                                   std::vector<ThickLeaf>& images, 
                                   int verbose) {
  viable_initials.resize(0);
  images.resize(0);
  initial = ThickLeaf(PI/2.0, 3.0*PI/2.0, 3.0*PI/2.0, PI/2.0);

  PartiallyDefinedMap fi = f.inverse();
  PartiallyDefinedMap gi = g.inverse();
  
  if (verbose > 0) std::cout << "Inverse maps:\n" << fi << "\n" << gi << "\n";
  
  Leaf ell(PI/2.0, 3.0*PI/2.0);
  
  if (verbose > 0)  std::cout << "Demo: l: " << ell << "\nf(l): " << f.act_on_leaf(ell) << "\nfi(l): " << 
                         fi.act_on_leaf(ell) << "\nfi(f(ell)): " << fi.act_on_leaf(f.act_on_leaf(ell)) << 
                        "\nf(fi(ell)): " << f.act_on_leaf(fi.act_on_leaf(ell)) << "\n";
  
  std::vector<Leaf> fg_ell(2);
  fg_ell[0] = f.act_on_leaf(ell);
  fg_ell[1] = g.act_on_leaf(ell);
  initial = (cyclically_ordered(ell.x, fg_ell[0].x, fg_ell[0].y) ? 
                    ThickLeaf( fg_ell[0].x, fg_ell[0].y, ell.y, ell.x ) :
                    ThickLeaf( fg_ell[0].y, fg_ell[0].x, ell.y, ell.x )   );
  std::vector<ThickLeaf> current_leaves(0);
  std::vector<ThickLeaf> next_leaves(0);
  current_leaves.push_back(initial);
  for (int i=0; i<depth; ++i) {
    
    if (verbose > 0) {
      std::cout << "Current leaves:\n";
      for (int j=0; j<(int)current_leaves.size(); ++j) {
        std::cout << j << ": " <<current_leaves[j] << "\n";
      }
    }
    
    next_leaves.resize(0);
  
    //act on each of the current leaves by the only map which can act on them
    //then cut the leaves by the images fell or gell as appropriate
    for (int j=0; j<(int)current_leaves.size(); ++j) {
      ThickLeaf image;
      if (verbose > 1) std::cout << "Acting on leaf: " << current_leaves[j] << "\n";
      if (fi.leaf_in_domain(current_leaves[j])) {
        if (verbose > 1) std::cout << "By fi\n";
        image = fi.act_on_leaf(current_leaves[j]);
        image.word = current_leaves[j].word;
        image.word.push_front(0);
      } else if (gi.leaf_in_domain(current_leaves[j])) {
        if (verbose > 1) std::cout << "By gi\n";
        image = gi.act_on_leaf(current_leaves[j]); 
        image.word = current_leaves[j].word;
        image.word.push_front(1);
      } else {
        std::cout << "Couldn't apply either map?\n";
        return;
      }
      if (verbose > 1) std::cout << "Got image leaf: " << image << "\n";
      std::vector<ThickLeaf> temp_next_leaves;
      image.split_by_leaves(temp_next_leaves, fg_ell);
      if (verbose > 1) {
        std::cout << "Split to the " << temp_next_leaves.size() << " leaves:\n";
        for (int k=0; k<(int)temp_next_leaves.size(); ++k) {
          std::cout << k << ": " << temp_next_leaves[k] << "\n";
        }
      }
      next_leaves.insert(next_leaves.end(), temp_next_leaves.begin(), temp_next_leaves.end());
    }
    current_leaves = next_leaves;
  }
  images = current_leaves;
  if (verbose > 0) {
    std::cout << "Got the image leaves:\n";
    for (int i=0; i<(int)current_leaves.size(); ++i) {
      std::cout << i << ": " << current_leaves[i] << "\n";
    }
  }
  if (verbose > 0) std::cout << "Getting viable initials:\n";
  for (int i=0; i<(int)current_leaves.size(); ++i) {
    ThickLeaf image = current_leaves[i];
    if (verbose> 0) std::cout << "From: " << image << "\n";
    for (int j=0; j<(int)current_leaves[i].word.size(); ++j) {
      image = (current_leaves[i].word[j] == 0 ? f.act_on_leaf(image) : g.act_on_leaf(image));
      if (verbose > 0) std::cout << "Applying " << (current_leaves[i].word[j]==0 ? "f" : "g") << " gives: " << image << "\n";
    }
    viable_initials.push_back(image);
  }
  
}



void Lamination::find_limit_leaves(std::vector<ThickLeaf>& initials,
                                   std::vector<ThickLeaf>& images, 
                                   Graph& inclusion_graph,
                                   bool& complete,
                                   int depth, int verbose) {
  std::vector<ThickLeaf> current_initials(0);
  std::vector<ThickLeaf> current_images(0);
  std::vector<bool> marked_for_action(0);
  std::vector<std::vector<std::pair<bool, int> > > inclusions;
  std::vector<std::vector< std::deque<char> > > inclusion_loops;
  bool went_past_depth = false;
  
  PartiallyDefinedMap fi = f.inverse();
  PartiallyDefinedMap gi = g.inverse();

  Leaf ell(PI/2.0, 3.0*PI/2.0);

  std::vector<Leaf> fg_ell(2);
  fg_ell[0] = f.act_on_leaf(ell);
  fg_ell[1] = g.act_on_leaf(ell);
  ThickLeaf initialf,initialg;
  initialf = (cyclically_ordered(ell.x, fg_ell[0].x, fg_ell[0].y) ? 
                    ThickLeaf( fg_ell[0].x, fg_ell[0].y, ell.y, ell.x ) :
                    ThickLeaf( fg_ell[0].y, fg_ell[0].x, ell.y, ell.x )   );
  initialg = (cyclically_ordered(ell.y, fg_ell[1].x, fg_ell[1].y) ? 
                    ThickLeaf( fg_ell[1].x, fg_ell[1].y, ell.x, ell.y ) :
                    ThickLeaf( fg_ell[1].y, fg_ell[1].x, ell.x, ell.y )   );
  
  current_initials.push_back(initialf);
  current_initials.push_back(initialg);
  current_images.push_back(initialf);
  current_images.push_back(initialg);
  marked_for_action.push_back(true);
  marked_for_action.push_back(true);
  
  while (true) {
    
    if (verbose>0) {
      std::cout << "Current initials and images:\n";
      for (int i=0; i<(int)current_initials.size(); ++i) {
        std::cout << i << ": " << (marked_for_action[i] ? "*\t" : "\t") << current_initials[i] << "\n\t" << current_images[i] << "\n";
      }
    }
    
    std::vector<ThickLeaf> act_stack(0);
    //everything that is marked for action, we need to act on it
    //until everything comes back to inside the initial one
    //
    //first, remove everything to be acted on, and push it on the stack
    for (int j=(int)current_images.size()-1; j>=0; --j) {
      if (marked_for_action[j]) {
        act_stack.push_back( current_images[j] );
        current_images.erase( current_images.begin() + j );
        current_initials.erase( current_initials.begin() + j );
      }
    }
    
    //now act on the stack until things go inside the initial, or, 
    //if they are below the allowed depth, forget them and record that 
    //we aren't rigorous
    std::vector<ThickLeaf> new_images(0);
    while (act_stack.size() > 0) {
      ThickLeaf temp = act_stack.back();
      
      if (verbose>0) { 
        std::cout << "Action stack:\n";
        for (int i=0; i<(int)act_stack.size(); ++i) {
          std::cout << i << ": " << act_stack[i] << "\n";
        }
      }
      
      act_stack.pop_back();
      ThickLeaf temp_image;
      
      if (verbose>0) std::cout << "Acting on " << temp << "\n";
      
      std::vector<ThickLeaf> temp_images(0);
      if (fi.leaf_in_domain(temp)) {
        temp_image = fi.act_on_leaf(temp);
        temp_image.word = temp.word;
        temp_image.word.push_front(0);
      } else if (gi.leaf_in_domain(temp)) {
        temp_image = gi.act_on_leaf(temp);
        temp_image.word = temp.word;
        temp_image.word.push_front(1);
      } else {
        std::cout << "Couldn't apply either map?\n";
        std::cout << lambda << " " << theta << "\n";
        std::cout.precision(15);
        std::cout << temp << "\n";
        return;
      }
      
      if (verbose>0) std::cout << "Got: " << temp_image << "\n";
      
      temp_image.split_by_leaves(temp_images, fg_ell);
      
      if (verbose>0) {
        std::cout << "Split to the leaves:\n";
        for (int i=0; i<(int)temp_images.size(); ++i) {
          std::cout << i << ": " << temp_images[i] << "\n";
        }
      }
      
      for (int j=0; j<(int)temp_images.size(); ++j) {
        if (initialf.contains(temp_images[j]) || initialg.contains(temp_images[j])) {
          if (verbose>0) std::cout << j << " is good\n";
          new_images.push_back(temp_images[j]);
        } else if ((int)temp_images[j].word.size() >= depth) {
          if (verbose>0) std::cout << j << " is too deep\n";
          went_past_depth = true;
        } else {
          if (verbose>0) std::cout << j << " isn't good, so onto the stack\n";
          act_stack.push_back(temp_images[j]);
        }
      }
    }
    
    //push the new images on to the real list of images, and record the initials
    for (int j=0; j<(int)new_images.size(); ++j) {
      ThickLeaf temp_initial = new_images[j];
      for (int k=0; k<(int)new_images[j].word.size(); ++k) { 
        temp_initial = (new_images[j].word[k] == 0 
                         ? f.act_on_leaf(temp_initial) 
                         : g.act_on_leaf(temp_initial));
      }
      current_images.push_back(new_images[j]);
      current_initials.push_back(temp_initial);
    }
    marked_for_action.resize(current_images.size());
    for (int j=0; j<(int)current_images.size(); ++j) {
      marked_for_action[j] = false;
    }
    
    if (verbose>0) {
      std::cout << "New initials and images:\n";
      for (int i=0; i<(int)current_initials.size(); ++i) {
        std::cout << i << ": " << (marked_for_action[i] ? "*\t" : "\t") << current_initials[i] << "\n\t" << current_images[i] << "\n";
      }
    }
    
    //now we have a new list of initials and images.  Get the inclusion
    //lists and see if they are conclusive
    bool all_good = true;
    ThickLeaf::find_inclusions(inclusions, current_images, current_initials,verbose);
    
    if (verbose>0) {
      std::cout << "Inclusions:\n";
      for (int i=0; i<(int)inclusions.size(); ++i) {
        std::cout << i << ": [";
        for (int j=0; j<(int)inclusions[i].size(); ++j) {
          std::cout << " (" << inclusions[i][j].first << "," << inclusions[i][j].second << ")" ;
        }
        std::cout << "]\n";
      }
    }
    
    for (int i=0; i<(int)current_images.size(); ++i) {
      for (int j=0; j<(int)inclusions[i].size(); ++j) {
        if (!inclusions[i][j].first) {
          if ((int)current_images[i].word.size() >= depth) {
            went_past_depth = true;
            continue;
          }
          all_good = false;
          marked_for_action[ inclusions[i][j].second ] = true;
        }
      }
    }
    
    if (all_good) break;
  }
  
  //now we have inclusions where everything is determined; get rid of 
  //the bool decoration and get the graph
  std::vector<std::vector<int> > compact_inclusions(inclusions.size(), std::vector<int>());
  for (int i=0; i<(int)inclusions.size(); ++i) {
    compact_inclusions[i].resize(inclusions[i].size());
    for (int j=0; j<(int)inclusions[i].size(); ++j) {
      compact_inclusions[i][j] = inclusions[i][j].second;
    }
  }
  
  std::vector<int> inclusion_weights(current_images.size());
  for (int i=0; i<(int)current_images.size(); ++i) {
    inclusion_weights[i] = current_images[i].word.size();
  }
  
  if (current_images.size() > 10) {
    inclusion_graph = Graph::from_edges(compact_inclusions);
  } else {
    inclusion_graph = Graph::from_edges_with_weights(compact_inclusions, inclusion_weights);
  }
  
  if (verbose>0) {
    std::cout << "Images (" << (went_past_depth ? "incomplete" : "complete") << "):\n";
    for (int i=0; i<(int)current_initials.size(); ++i) {
      std::cout << i << ":\t" << current_initials[i] << "\n\t" << current_images[i] << "\n";
    }
  
    std::cout << "Inclusion graph:\n";
    std::cout << inclusion_graph;
  }
  
  initials = current_initials;
  images = current_images;
  complete = !went_past_depth;
}
  












