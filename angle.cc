#include <cstdlib>
#include <cstring>
#include <iostream>

#include "angle_gui.h"
#include "lam.h"

int main(int argc, char* argv[]) {
  std::cout.precision(10);
  
  if (argc > 1) {
    if (argv[1][0] == '-' && argv[1][1] == 'L') {
      int d;
      double lam,the;
      bool mathematica = (strlen(argv[1]) > 2 && argv[1][2] == 'm');
      if (argc < 5) {
        std::cout << "Please input a depth, lambda, and theta\n";
        return 1;
      }
      d = atoi(argv[2]);
      lam = atof(argv[3]);
      the = atof(argv[4]);
      Lamination L(lam, the);
      std::vector<ThickLeaf> initials, images;
      Graph inclusion_graph;
      bool complete;
      L.find_limit_leaves(initials, images, inclusion_graph, complete, d, 2);
      std::cout << (complete ? "complete\n" : "incomplete\n");
      if (mathematica) {
        std::cout << inclusion_graph.mathematica_string() << "\n";
      } else {
        std::cout << inclusion_graph;
      }
      return 0;
    }
  }
  
  AngleGui a;
  a.launch();
  return 0;
}
