#ifndef __ANGLE_GUI__
#define __ANGLE_GUI__

extern "C" {
#include <X11/Xlib.h>
#include <X11/Xutil.h>
}

#include "point.h"
#include "lam.h"

struct AngleGui;

struct Widget {
  Point2d<int> ul; //the upper left hand corner
  int height;
  int width;
  AngleGui* ag; //the gui, so we can call the member functions
  Pixmap p;
  GC gc;
  void (AngleGui::*event_signal)(void*);
  
  virtual void initial_draw() = 0;
  virtual void redraw() = 0;
  void clear();
  bool contains_pixel(int x, int y);
  bool intersects_rectangle(const Point2d<int>& p, int w, int h);
  virtual void run_event_signal(XEvent* e) {
    (ag->*event_signal)((void*)e);
  }
  Widget() {}
};

struct WidgetDraw : Widget {
  WidgetDraw() {}
  WidgetDraw(AngleGui* i, int w, int h, void (AngleGui::*f)(void*));
  void redraw();
  void initial_draw();
};

struct WidgetText : Widget {
  std::string text;
  Point2d<int> text_data;
  
  WidgetText() {}
  WidgetText(AngleGui* ag, int w, const std::string& text);
  
  void initial_draw();
  void redraw();
  void update_text(const std::string& new_text);
};


struct WidgetIntSelectorEvent {
  enum {UP, DOWN} type;
};

struct WidgetIntSelector : Widget {
  std::string label;
  std::string value_string;
  
  Point2d<int> leftarrow_point; //these are so we know where we clicked
  Point2d<int> rightarrow_point; 
  
  WidgetIntSelector() {}
  WidgetIntSelector(AngleGui* ag, int w, const std::string& label, int start_val, void (AngleGui::*f)(void*));
  void run_event_signal(XEvent* e);
  void update_value(int new_val);
  void initial_draw();
  void redraw();
};

struct WidgetButton : Widget {
  std::string text;
  Point2d<int> text_position;
  
  WidgetButton() {}
  WidgetButton(AngleGui* i, int w, int h, const std::string& t, void (AngleGui::*f)(void*));
  void initial_draw();
  void redraw();
};


struct WidgetCheck : Widget {
  std::string text;
  bool checked;
  Point2d<int> text_position;
  
  WidgetCheck() {}
  WidgetCheck(AngleGui* i, int w, const std::string& t, bool c, void (AngleGui::*f)(void*));
  void redraw();
  void initial_draw();
};






struct AngleGui {
  
  //graphics stuff
  Display* display;
  int screen;
  Window main_window;
  Colormap col_map;
  int get_rgb_color(double r, double g, double b);
  bool main_window_initialized;
  int main_window_height;
  int main_window_width;
  
  std::vector<Widget*> widgets;
  
  WidgetDraw W_param_plot;
  WidgetButton W_param_recenter;
  WidgetButton W_param_zoom_in;
  WidgetButton W_param_zoom_out;
  WidgetIntSelector W_param_depth;
  WidgetIntSelector W_param_mesh;
  WidgetCheck W_param_words;
  WidgetCheck W_param_inclusion;
  WidgetCheck W_param_plot_points;
  WidgetCheck W_param_endpoints;
  WidgetIntSelector W_param_endpoints_depth;
  WidgetButton W_param_write_trajectories;
  WidgetText W_param_mouse_lambda;
  WidgetText W_param_mouse_theta;
  
  WidgetDraw W_lam_plot;
  WidgetText W_lam_lambda_label;
  WidgetText W_lam_theta_label;
  WidgetText W_lam_lam_type;
  WidgetText W_lam_trajectories;
  WidgetIntSelector W_lam_depth;
  WidgetIntSelector W_lam_backward_depth;
  WidgetIntSelector W_lam_limit_leaves_depth;
  WidgetCheck W_lam_highlight;
  
  void S_param_recenter(void* e);
  void S_param_zoom_in(void* e);
  void S_param_zoom_out(void* e);
  void S_param_depth(void* e);
  void S_param_mesh(void* e);
  void S_param_words(void* e);
  void S_param_inclusion(void* e);
  void S_param_endpoints(void* e);
  void S_param_endpoints_depth(void* e);
  void S_param_plot_points(void* e);
  void S_param_draw(void* e);
  void S_param_write_trajectories(void* e);
  
  void S_lam_depth(void* e);
  void S_lam_backward_depth(void* e);
  void S_lam_limit_leaves_depth(void* e);
  void S_lam_highlight(void* e);
  void S_lam_draw(void* e);
  
  void draw_param();
  void draw_lam();
  void recompute_lam_data();
  void reset_param_grid();
  void reset_param_window();
  int compute_color_from_grid(const Point4d<int>& p);
  void reset_highlighted_point(double L, double T);
  void draw_highlighted_point();
  void draw_param_point(double L, double T, int col);
  void draw_param_disk(double L, double T, double r, int col);
  
  double param_theta_l;
  double param_theta_u;
  double param_lambda_l;
  double param_lambda_u;
  int param_depth;
  int param_pixel_width;
  int param_pixel_height;
  int param_num_pixel_groups_width;
  int param_num_pixel_groups_height;
  int param_pixel_group_pixels;
  double param_pixel_group_lambda_size;
  double param_pixel_group_theta_size;
  double param_pixel_lambda_size;
  double param_pixel_theta_size;
  std::vector<std::vector<Point4d<int> > > param_grid;
  bool param_words;
  bool param_inclusion;
  bool param_endpoints;
  int param_endpoints_depth;
  bool param_plot_points;
  
  void param_LT_to_pixel_group(double lambda, double theta, int& x, int& y);
  void param_pixel_group_to_LT(int x, int y, double& lambda, double& theta);
  void param_pixel_to_LT(int x, int y, double& lambda, double& theta);
  void param_LT_to_pixel(double lambda, double theta, int& x, int& y);
  void lam_coords_to_pixel(double dx, double dy, int& x, int& y);
  void lam_coords_to_pixel(double dx, double dy, short& x, short& y);
  void lam_coords_to_pixel(const std::vector<Point2d<float> >& floats, 
                           std::vector<XPoint>& pixels);
  void lam_size_to_pixel(double dm, int& m);
  
  int lam_depth;
  int lam_backward_depth;
  int lam_limit_leaves_depth;
  double lam_lambda;
  double lam_theta;
  int lam_pixel_size;
  bool lam_highlight;
  
  void detach_widget(Widget* w);
  void pack_widget_upper_right(const Widget* w1, Widget* w2);
  void pack_window();
  void launch();
  void main_loop();
};





#endif
