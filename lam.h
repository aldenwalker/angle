#ifndef __LAM__
#define __LAM__

#include <vector>
#include <deque>
#include <string>

#include "point.h"

#define PI 3.14159265358979323846
#define TWOPI 6.2831853071795864769

enum LamType {CUT_POINT, PROPER, NONE};

double angle_diff(double a, double b);
double angle_dist(double a, double b);

struct Graph {
  std::vector< std::vector<int> > edges_starting_with;
  std::vector< std::vector<int> > edges_ending_with;
  std::vector< std::pair<int,int> > edges;
  int num_verts;
  int num_edges;
  Graph();
  
  void multiply_edge(int ei, int m);
  
  static Graph from_edges(const std::vector< std::vector<int> >& e);
  static Graph from_edges_with_weights(const std::vector< std::vector<int> >& e,
                                       const std::vector<int>& weights);
  void find_cycle_components(std::vector<std::vector<std::deque<int> > >& comps);
  int rank_if_undirected() const;
  double approximate_leading_eigenvalue(int iters) const;
  std::vector<Graph> connected_components() const;
  
  std::string mathematica_string() const;
};

std::ostream& operator<<(std::ostream& os, const Graph& G);


struct Leaf {
  double x;
  double y;
  std::deque<char> word;
  Leaf();
  Leaf(double x, double y);
  void polygon_points(std::vector<Point2d<float> >& points, int num_per_radian=5) const;
  void get_circle_data(double& centerx, double& centery, double& radius, double& a1, double& a_extent) const;
  bool operator<(const Leaf& other) const;
};

std::ostream& operator<<(std::ostream& os, const Leaf& ell);

struct ThickLeaf {
  double p1x;
  double p1y;
  double p2x;
  double p2y;
  std::deque<char> word;
  
  ThickLeaf() {}
  ThickLeaf(double p1x, double p1y, double p2x, double p2y);
  ThickLeaf(double p1x, double p1y, double p2x, double p2y, const std::deque<char>& w);
  
  void polygon_points(std::vector<Point2d<float> >& points, int num_per_radian=5) const;
  bool is_disjoint(const Leaf& ell) const;
  void split_by_leaves(std::vector<ThickLeaf>& result, const Leaf& ell) const;
  void split_by_leaves(std::vector<ThickLeaf>& result, const std::vector<Leaf>& leaves) const;
  bool contains(const ThickLeaf& other) const;
  bool contains(const Leaf& other) const;
  bool intersects(const ThickLeaf& other) const;
  
  bool operator<(const ThickLeaf& other) const;
  
  static void find_inclusions(std::vector<std::vector<std::pair<bool, int> > >& inclusions, 
                              const std::vector<ThickLeaf>& A,
                              const std::vector<ThickLeaf>& B,
                              int verbose=0);
  
};

std::ostream& operator<<(std::ostream& os, const ThickLeaf& ell);



struct PartiallyDefinedMap {
  double domain_x;
  double domain_y;
  double image_x;
  double image_y;
  double scale_factor;
  
  PartiallyDefinedMap();
  PartiallyDefinedMap(double lambda, double theta, bool f);
  PartiallyDefinedMap(double x1, double y1, double x2, double y2);
  //PartiallyDefinedMap(double lambda1, double theta1, double lambda2, double theta2);
  PartiallyDefinedMap inverse() const;
  
  bool leaf_in_domain(const Leaf& ell) const;  
  bool leaf_in_domain(const ThickLeaf& ell) const;
  double act_on_point(double x) const;
  Leaf act_on_leaf(const Leaf& ell) const;
  ThickLeaf act_on_leaf(const ThickLeaf& ell) const;
};

std::ostream& operator<<(std::ostream& os, const PartiallyDefinedMap& f);

struct Lamination {
  double lambda;
  double theta;
  PartiallyDefinedMap f;
  PartiallyDefinedMap g;
  
  Lamination();
  Lamination(double lambda, double theta);
  
  bool check_can_act(const Leaf& ell, int which_map);
  Leaf act_on_leaf(const Leaf& ell, int which_map);
  void compute_lam_type(int depth, LamType& ellt, int& difficulty);
  void compute_lam_type_with_words(int depth, LamType& ellt, int& difficulty, std::vector<Leaf>& viable_leaves);
  void compute_backwards(int depth, ThickLeaf& initial, 
                                    std::vector<ThickLeaf>& viable_initials,
                                    std::vector<ThickLeaf>& images,
                                    int verbose=0);
  void find_limit_leaves(std::vector<ThickLeaf>& initials,
                         std::vector<ThickLeaf>& images, 
                         Graph& inclusion_graph,
                         bool& complete,
                         int depth, int verbose=0);
  
};


#endif
