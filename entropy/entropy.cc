#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>

struct Interval {
  double x;
  double y;
  double length;
  bool empty;
  Interval() {
    empty = true;
  }
  Interval(double x, double y) {
    this->x = x;
    this->y = y;
    length = y-x;
    empty = (length == 0);
  }
  Interval intersection(const Interval& I) const {
    if (y < I.x) {
      return Interval();
    } else if (x < I.x) {
      if (y < I.y) {
        return Interval(I.x, y);
      } else {
        return Interval(I.x, I.y);
      }
    } else if (y < I.y) {
      return Interval(x,y);
    } else if (x < I.y) {
      return Interval(x, I.y);
    } else {
      return Interval();
    }
  }
  Interval operator+(double t) const {
    return Interval(x+t, y+t);
  }
  Interval operator-(double t) const {
    return Interval(x-t,y-t);
  }
  Interval mod_one() const {
    double d1, d2;
    d1 = x;
    d2 = y;
    if (d1 < 0.0) {
      while (d1 < 0.0) {
        d1 += 1.0;
        d2 += 1.0;
      }
      return Interval(d1, d2);
    }
    if (d2 > 1.0) {
      while (d2 > 1.0) {
        d2 -= 1.0;
        d1 -= 1.0;
      }
      return Interval(d1, d2);
    }
    return Interval(d1, d2);
  }
};

std::ostream& operator<<(std::ostream& os, const Interval& I) {
  if (I.empty) return os << "()";
  return os << "(" << I.x << "," << I.y << ")";
}


struct LIM {
  Interval dom;
  Interval im;
  double factor;
  LIM() {}
  LIM(const Interval& d, const Interval& i) {
    dom = d;
    im = i;
    factor = (dom.length > 0 ? im.length / dom.length : 0);
  }
  double image(double x) const {
    return (x - dom.x)*factor + im.x;
  }
  double preimage(double x) const {
    if (factor==0) return dom.x;
    return (x - im.x)*(1.0/factor) + dom.x;
  }
  Interval image(const Interval& I) const {
    Interval in_dom = I.intersection(dom);
    if (in_dom.empty) return in_dom;
    return Interval( image(in_dom.x), image(in_dom.y) );
  }
  Interval preimage(const Interval& I) const {
    Interval in_im = I.intersection(im);
    if (in_im.empty) return in_im;
    return Interval( preimage(in_im.x), preimage(in_im.y) );
  }
};

std::ostream& operator<<(std::ostream& os, const LIM& L) {
  return os << L.dom << " -> " << L.im << " (x" << L.factor << ")";
}

struct LIET {
  std::vector<LIM> comps;
  LIET() {
    comps.resize(0);
  }
  LIET(const std::vector<LIM>& comps) {
    this->comps = comps;
  }
  
  static LIET double_plus_theta(double theta) {
    LIET ans(std::vector<LIM>(0));
    //f map
    double f_breakpoint;
    LIM f1, f2;
    if (theta < 0.25) {
      f_breakpoint = 0.125 - theta/2.0;
      f1 = LIM( Interval(0.0, f_breakpoint), Interval(0.75 + theta, 1.0) );
      f2 = LIM( Interval(f_breakpoint, 0.5), Interval(0.0, 0.75 + theta) );
    } else {
      f_breakpoint = 0.625 - theta / 2.0;
      f1 = LIM( Interval(0.0, f_breakpoint), Interval(theta - 0.25, 1.0) );
      f2 = LIM( Interval(f_breakpoint, 0.5), Interval(0.0, theta - 0.25) );
    }
    ans.comps.push_back(f1);
    ans.comps.push_back(f2);
    //g map
    double g_breakpoint;
    LIM g1, g2;
    if (theta < 0.75) {
      g_breakpoint = 0.875 - theta/2.0;
      g1 = LIM( Interval(0.5, g_breakpoint), Interval(0.25 + theta, 1.0) );
      g2 = LIM( Interval(g_breakpoint, 1.0), Interval(0.0, 0.25 + theta) );
    } else {
      g_breakpoint = 1.375 - theta/2.0;
      g1 = LIM( Interval(0.5, g_breakpoint), Interval(theta - 0.75, 1.0) );
      g2 = LIM( Interval(g_breakpoint, 1.0), Interval(0.0, theta - 0.75) );
    }
    ans.comps.push_back(g1);
    ans.comps.push_back(g2);
    return ans;
  }  
  
  static LIET lambda_plus_theta(double theta, double lambda) {
    LIET ans(std::vector<LIM>(0));
    Interval I1(0.0,0.5);
    Interval I1_offset = I1 + 1.0;
    Interval I1_offset_2 = I1 + 2.0;
    Interval I2(0.5,1.0);
    Interval I2_offset = I2 + 1.0;
    Interval I2_offset_down = I2 - 1.0;
    //f map
    double zero_image = theta + 0.25*(1.0-lambda);
    double half_image = theta + 0.25*(1.0+lambda);
    Interval fI(zero_image, half_image);
    std::vector<Interval> images(5);
    images[0] = fI.intersection( I2_offset_down );
    images[1] = fI.intersection( I1 );
    images[2] = fI.intersection( I2 );
    images[3] = fI.intersection( I1_offset );
    images[4] = fI.intersection( I2_offset );
    int i=0;
    while (images[i].empty) ++i;
    double bp1 = (images[i].y-theta-0.25)/lambda + 0.25;
    if (i<3 && !images[i+2].empty) {
      double bp2 = (images[i+1].y-theta-0.25)/lambda + 0.25;
      ans.comps.push_back( LIM( Interval(0.0,bp1), images[i].mod_one() ) );
      ans.comps.push_back( LIM( Interval(bp1,bp2), images[i+1].mod_one() ) );
      ans.comps.push_back( LIM( Interval(bp2,0.5), images[i+2].mod_one() ) );
    } else {
      ans.comps.push_back( LIM( Interval(0.0,bp1), images[i].mod_one() ) );
      ans.comps.push_back( LIM( Interval(bp1,0.5), images[i+1].mod_one() ) );
    }
    ///////////////
    half_image = theta + 0.75 - lambda*0.25;
    double one_image = theta + 0.75 +lambda*0.25;
    Interval gI(half_image, one_image);
    images[0] = gI.intersection( I1 );
    images[1] = gI.intersection( I2 );
    images[2] = gI.intersection( I1_offset );
    images[3] = gI.intersection( I2_offset );
    images[4] = gI.intersection( I1_offset_2 );
    i=0;
    while (images[i].empty) ++i;
    bp1 = (images[i].y-theta-0.75)/lambda + 0.75;
    if (i<3 && !images[i+2].empty) {
      double bp2 = (images[i+1].y-theta-0.75)/lambda + 0.75;
      ans.comps.push_back( LIM( Interval(0.5,bp1), images[i].mod_one() ) );
      ans.comps.push_back( LIM( Interval(bp1,bp2), images[i+1].mod_one() ) );
      ans.comps.push_back( LIM( Interval(bp2,1.0), images[i+2].mod_one() ) );
    } else {
      ans.comps.push_back( LIM( Interval(0.5,bp1), images[i].mod_one() ) );
      ans.comps.push_back( LIM( Interval(bp1,1.0), images[i+1].mod_one() ) );
    }
    return ans;
  }
  
  double entropy(int depth) {
    //std::vector<Interval> current_partition(comps.size());
    //for (int i=0; i<(int)comps.size(); ++i) {
    //  current_partition[i] = comps[i].dom;
    //}
    std::vector<Interval> current_partition(2);
    current_partition[0] = Interval(0.0,0.5);
    current_partition[1] = Interval(0.5,1.0);
    std::vector<Interval> new_partition;
    for (int i=0; i<depth; ++i) {
      new_partition.resize(0);
      //std::cout << "Current partition: \n";
      for (int j=0; j<(int)current_partition.size(); ++j) {
        //std::cout << j << ": " << current_partition[j] << " -> ";
        for (int k=0; k<(int)comps.size(); ++k) {
          Interval I = comps[k].preimage(current_partition[j]);
          if (!I.empty) new_partition.push_back(I);
          //std::cout << I << " ";
        }
        //std::cout << "\n";
      }
      current_partition = new_partition;
    }
    std::cout << current_partition.size() << "\n";
    double H = 0;
    for (int i=0; i<(int)current_partition.size(); ++i) {
      double ell = current_partition[i].length;
      H += -ell*log(ell);
    }
    return H / (double)depth;
  }
  
  double entropy_tree(int depth) {
    std::vector<std::pair<int, Interval> > stack(0);
    stack.push_back( std::make_pair(0, Interval(0.0,0.5)) );
    stack.push_back( std::make_pair(0, Interval(0.5,1.0)) );
    double H = 0;
    int seen = 0;
    while (stack.size() > 0) {
      int d = stack.back().first;
      Interval I = stack.back().second;
      stack.pop_back();
      for (int i=0; i<(int)comps.size(); ++i) {
        Interval II = comps[i].preimage(I);
        if (II.empty) continue;
        if (d+1 == depth) {
          seen++;
          H += -II.length*log(II.length);
          continue;
        }
        stack.push_back( std::make_pair( d+1, II ) );
      }
    }
    std::cout << seen << "\n";
    return H / (double)depth;
  }
    
};

std::ostream& operator<<(std::ostream& os, const LIET& L) {
  for (int i=0; i<(int)L.comps.size(); ++i) {
    os << i << ": " << L.comps[i] << "\n";
  }
  return os;
}



int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "Usage: ./entropy [-p] depth <theta or step>\n";
    std::cout << "\"./entropy depth theta\" computes the entropy at theta to the given depth\n";
    std::cout << "\"./entropy -p depth step\" computes for 0<theta<1 with the given step size\n";
    return 0;
  }
  int depth;
  bool plot = false;
  bool with_lambda = false;
  double lambda;
  double step;
  double theta;
  if (argv[1][0] == '-') {
    plot = true;
    if (argv[1][2] != '\0') {
      with_lambda = true;
    }
    depth = atoi(argv[2]);
    step = atof(argv[3]);
  } else {
    depth = atoi(argv[1]);
    theta = atof(argv[2]);
    if (argc > 3) {
      with_lambda = true;
      lambda = atof(argv[3]);
    }
  }
  
  if (plot) {
    if (with_lambda) { 
      for (double t = 0.0; t<1.0; t+=step) {
        for (double l = 1.0; l<2.0; l+=step) {
          LIET L = LIET::lambda_plus_theta(t,l);
          double e = L.entropy(depth);
          std::cout << "{" << t << "," << l << "," << e << "},";
        }
      }
      std::cout << "\n";
    } else {
      for (double t = 0.0; t<1.0; t+=step) {
        LIET L = LIET::double_plus_theta(t);
        double e = L.entropy(depth);
        std::cout << "{" << t << "," << e << "},";
      }
      std::cout << "\n";
    }
  } else {
    LIET L = (with_lambda ? LIET:: lambda_plus_theta(theta, lambda) : LIET::double_plus_theta(theta));
    std::cout << L;
    //double e = L.entropy(depth);
    double e2 = L.entropy_tree(depth);
    if (with_lambda) {
      //std::cout << "Entropy of theta=" << theta << ", lambda=" << lambda << " at depth " << depth << ": " << e << "\n";
      std::cout << "Entropy of theta=" << theta << ", lambda=" << lambda << " at depth " << depth << ": " << e2 << "\n";
    } else {
      //std::cout << "Entropy of theta=" << theta << " at depth " << depth << ": " << e << "\n";
      std::cout << "Entropy of theta=" << theta << " at depth " << depth << ": " << e2 << "\n";
    }
  }
  return 0;
}










































