#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "angle_gui.h"

/***************************************************************************
 * widget functions
 ***************************************************************************/
//return the width (.x) and height offset (.y) of the text string
Point2d<int> compute_text_width_and_offset(Display* d, GC* gc, const std::string& text) {
  XFontStruct* font = XLoadQueryFont(d, "fixed");
  XSetFont(d, *gc, font->fid);  
  XCharStruct cs;
  int dir, descent, ascent;
  XTextExtents(font, text.c_str(), text.size(), &dir, &ascent, &descent, &cs);
  return Point2d<int>(cs.rbearing - cs.lbearing, (ascent-descent)/2);
}

bool Widget::contains_pixel(int x, int y) {
  return (ul.x <= x) && (x < ul.x + width) && (ul.y <= y) && (y < ul.y + height);
}

bool Widget::intersects_rectangle(const Point2d<int>& p, int w, int h) {
  return !(ul.x > p.x + w || ul.y > p.y + h || p.x > ul.x + width || p.y > ul.y + height);
}

void Widget::clear() {
  XSetForeground(ag->display, gc, WhitePixel(ag->display, ag->screen));
  XFillRectangle(ag->display, ag->main_window, gc, ul.x, ul.y, width, height);
}


WidgetDraw::WidgetDraw(AngleGui* i, int w, int h, void (AngleGui::*f)(void*)) {
  width = w;
  height = h;
  ag = i;
  event_signal = f;
  p = XCreatePixmap(ag->display, ag->main_window,
                    width, height, DefaultDepth(ag->display, ag->screen));
  gc = XCreateGC(ag->display, RootWindow(ag->display, ag->screen), 0, NULL);
  XSetForeground(ag->display, gc, WhitePixel(ag->display, ag->screen));
  XSetBackground(ag->display, gc, WhitePixel(ag->display, ag->screen));
  XFillRectangle(ag->display, p, gc, 0, 0, width, height);
  XSetForeground(ag->display, gc, BlackPixel(ag->display, ag->screen));
  XDrawRectangle(ag->display, p, gc, 0, 0, width-1, height-1);
}

void WidgetDraw::redraw() {
  initial_draw();
}

void WidgetDraw::initial_draw() {
  XCopyArea(ag->display, p, ag->main_window, gc, 0, 0, width, height, ul.x, ul.y);
}


WidgetText::WidgetText(AngleGui* ag, int w, const std::string& text) {
  this->text = text;
  height = 20;
  this->ag = ag;
  event_signal = NULL;
  gc = XCreateGC(ag->display, RootWindow(ag->display, ag->screen), 0, NULL);
  text_data = compute_text_width_and_offset(ag->display, &gc, text);
  if (w < 0) {
    width = text_data.x + 10;
  } else {
    width = w;
  }
  p = XCreatePixmap(ag->display, ag->main_window,
                    width, height, DefaultDepth(ag->display, ag->screen));
}
  
void WidgetText::initial_draw() {
  XSetForeground(ag->display, gc, WhitePixel(ag->display, ag->screen));
  XFillRectangle(ag->display, p, gc, 0, 0, width, height);
  XSetForeground(ag->display, gc, BlackPixel(ag->display, ag->screen));
  XDrawString(ag->display, p, gc, 5, height/2 + text_data.y, text.c_str(), text.size()); 
  XCopyArea(ag->display, p, ag->main_window, gc, 0, 0, width, height, ul.x, ul.y);
}

void WidgetText::redraw() {
  initial_draw();
}

void WidgetText::update_text(const std::string& new_text) {
  text = new_text;
  text_data = compute_text_width_and_offset(ag->display, &gc, text);
  redraw();
}










WidgetIntSelector::WidgetIntSelector(AngleGui* ag, int w, const std::string& label, int start_val, void (AngleGui::*f)(void*)) {
  std::stringstream T;
  T.str("");
  T << start_val;
  value_string = T.str();
  this->label = label;
  
  this->ag = ag;
  height = 20;
  event_signal = f;
  gc = XCreateGC(ag->display, RootWindow(ag->display, ag->screen), 0, NULL);
  
  Point2d<int> label_data = compute_text_width_and_offset(ag->display, &gc, label);
  Point2d<int> value_data = compute_text_width_and_offset(ag->display, &gc, value_string);
  
  int desired_width = label_data.x + value_data.x + 5*5 + 2*10;
  if (w > 0) {
    width = w; 
  } else {
    width = desired_width;
  }
  
  p = XCreatePixmap(ag->display, ag->main_window,
                    width, height, DefaultDepth(ag->display, ag->screen));
  
}

void WidgetIntSelector::run_event_signal(XEvent* e) {
  if (e->type != ButtonPress) return;
  WidgetIntSelectorEvent se;
  if (ul.x+leftarrow_point.x < e->xbutton.x && e->xbutton.x < ul.x+leftarrow_point.x+10) {
    se.type = WidgetIntSelectorEvent::DOWN;
  } else if (ul.x+rightarrow_point.x-10 < e->xbutton.x && e->xbutton.x < ul.x+rightarrow_point.x) {
    se.type = WidgetIntSelectorEvent::UP;
  } else {
    return;
  }
  (ag->*event_signal)((void*)(&se));
}

void WidgetIntSelector::update_value(int new_val) {
  std::stringstream T;
  T.str("");
  T << new_val;
  value_string = T.str();
  redraw();
}


void WidgetIntSelector::initial_draw() {
  Point2d<int> label_data = compute_text_width_and_offset(ag->display, &gc, label);
  Point2d<int> label_position(5, height/2 + label_data.y);
  
  leftarrow_point = Point2d<int>(5+label_data.x+5, height/2);
  
  Point2d<int> value_data = compute_text_width_and_offset(ag->display, &gc, value_string);
  Point2d<int> value_position(leftarrow_point.x+10+5, height/2 + value_data.y);
  
  rightarrow_point = Point2d<int>(value_position.x + value_data.x + 5 + 10, height/2);
  
  XSetForeground(ag->display, gc, WhitePixel(ag->display, ag->screen));
  XFillRectangle(ag->display, p, gc, 0, 0, width, height);
  XSetForeground(ag->display, gc, BlackPixel(ag->display, ag->screen));
  
  XDrawString(ag->display, p, gc, label_position.x, label_position.y, label.c_str(), label.size()); 
  
  XPoint LA[3];
  LA[0].x = leftarrow_point.x;    LA[0].y = leftarrow_point.y;
  LA[1].x = leftarrow_point.x+10; LA[1].y = leftarrow_point.y+5;
  LA[2].x = leftarrow_point.x+10; LA[2].y = leftarrow_point.y-5;
  XFillPolygon(ag->display, p, gc, LA, 3, Convex, CoordModeOrigin);
  
  XDrawString(ag->display, p, gc, value_position.x, value_position.y, value_string.c_str(), value_string.size()); 
  
  XPoint RA[3];
  RA[0].x = rightarrow_point.x;    RA[0].y = rightarrow_point.y;
  RA[1].x = rightarrow_point.x-10; RA[1].y = rightarrow_point.y-5;
  RA[2].x = rightarrow_point.x-10; RA[2].y = rightarrow_point.y+5;
  XFillPolygon(ag->display, p, gc, RA, 3, Convex, CoordModeOrigin);
  
  XCopyArea(ag->display, p, ag->main_window, gc, 0, 0, width, height, ul.x, ul.y);
}

void WidgetIntSelector::redraw() { 
  initial_draw(); 
}





int AngleGui::get_rgb_color(double r, double g, double b) {
  XColor temp;
  temp.flags = DoRed | DoGreen | DoBlue;
  temp.red = (int)(r*65535);
  temp.green = (int)(g*65535);
  temp.blue = (int)(b*65535);
  if (XAllocColor(display, DefaultColormap(display, screen), &temp) == 0) {
    std::cout << "Color not found?\n";
  }
  return temp.pixel;
}


WidgetButton::WidgetButton(AngleGui* i, int w, int h, const std::string& t, void (AngleGui::*f)(void*)) {
  ag = i;
  text = t;
  event_signal = f;
  
  gc = XCreateGC(ag->display, RootWindow(ag->display, ag->screen), 0, NULL);
  
  Point2d<int> text_data = compute_text_width_and_offset(ag->display, &gc, text);
  
  if (h > 0) {
    height = h;
  } else {
    height = 20;
  }
  if (w > 0) {
    width = w; 
  } else {
    width = text_data.x + 10;
  }
  text_position = Point2d<int>(5, height/2 + text_data.y); 

  p = XCreatePixmap(ag->display, ag->main_window,
                    width, height, DefaultDepth(ag->display, ag->screen));
  
  //clear the pixmap
  XSetForeground(ag->display, gc, WhitePixel(ag->display, ag->screen));
  XSetBackground(ag->display, gc, WhitePixel(ag->display, ag->screen));
  XFillRectangle(ag->display, p, gc, 0, 0, width, height);
  
  //set the real colors
  XSetForeground(ag->display, gc, BlackPixel(ag->display, ag->screen));
  XSetBackground(ag->display, gc, WhitePixel(ag->display, ag->screen));
  
  //draw the square
  XSetLineAttributes(ag->display, gc, 0, LineSolid, CapButt, JoinMiter);
  XDrawLine(ag->display, p, gc, 1, 1, 1, height-2);
  XDrawLine(ag->display, p, gc, 1, height-2, width-2, height-2);
  XDrawLine(ag->display, p, gc, width-2, height-2, width-2, 1);
  XDrawLine(ag->display, p, gc, width-2, 1, 1, 1);
  //draw the label
  XDrawString(ag->display, p, gc, text_position.x, text_position.y, text.c_str(), text.size()); 
}

void WidgetButton::initial_draw() {
  XCopyArea(ag->display, p, ag->main_window, gc, 0, 0, width, height, ul.x, ul.y);
}

void WidgetButton::redraw() {
  initial_draw();
}


WidgetCheck::WidgetCheck(AngleGui* i, int w, const std::string& t, bool c, void (AngleGui::*f)(void*)) {

  ag = i;
  text = t;
  checked = c;
  event_signal = f;
  
  gc = XCreateGC(ag->display, RootWindow(ag->display, ag->screen), 0, NULL);
  
  Point2d<int> text_data = compute_text_width_and_offset(ag->display, &gc, text);
  
  height = 20;
  if (w > 0) {
    width = w; 
  } else {
    width = text_data.x + 25;
  }
  text_position = Point2d<int>(20, height/2 + text_data.y); 
  
  p = XCreatePixmap(ag->display, ag->main_window,
                    width, height, DefaultDepth(ag->display, ag->screen));
  
  initial_draw();
} 

void WidgetCheck::redraw() {
  XSetForeground(ag->display, gc, WhitePixel(ag->display, ag->screen));
  XFillRectangle(ag->display, p, gc, 0, 0, width, height);
  XSetForeground(ag->display, gc, BlackPixel(ag->display, ag->screen));
  if (checked) {
    XFillRectangle(ag->display, p, gc, 5, 5, 11, 11);
  } else {
    XDrawRectangle(ag->display, p, gc, 5, 5, 10, 10);
  }
  XDrawString(ag->display, p, gc, text_position.x, text_position.y, text.c_str(), text.size());  
  XCopyArea(ag->display, p, ag->main_window, gc, 0, 0, width, height, ul.x, ul.y);
}
  

void WidgetCheck::initial_draw() {
  redraw();
}





/****************************************************************************
 * signals 
 ****************************************************************************/
void AngleGui::S_param_draw(void* ee) {
    XEvent* e = (XEvent*)ee;
    if (e->type == KeyPress) return;
  
  //the following is run if the button is pressed or if there is
  //motion where the button is down
  if ( (e->type == ButtonPress && e->xbutton.button == Button1) ||
       (e->type == MotionNotify && ((e->xmotion.state >> 8)&1)) ) {
    int widget_x = e->xbutton.x - W_param_plot.ul.x;
    int widget_y = e->xbutton.y - W_param_plot.ul.y;
    double L,T;
    param_pixel_to_LT(widget_x, widget_y, L, T);
    reset_highlighted_point(L, T);
  }
  /*
  //it's a right mouse click -- zoom in
  } else if (e->type == ButtonPress && e->xbutton.button == Button3) { 
    int widget_x = e->xbutton.x - W_mand_plot.ul.x;
    int widget_y = e->xbutton.y - W_mand_plot.ul.y;
    cpx c = mand_pixel_to_cpx(Point2d<int>(widget_x, widget_y));
    IFS.set_params(c,c);
    mand_zoom(0.5);
    recompute_point_data();
  }
  */
  
  //additionally, if the mouse is moved, we need to update the 
  //text
  if (e->type == MotionNotify) {
    int widget_x = e->xbutton.x - W_param_plot.ul.x;
    int widget_y = e->xbutton.y - W_param_plot.ul.y;
    double lambda,theta;
    param_pixel_to_LT(widget_x, widget_y, lambda, theta);
    std::stringstream T; T.str("");
    T.precision(15);
    T << "Mouse lambda: " << lambda;
    W_param_mouse_lambda.update_text(T.str());
    T.str("");
    T << "Mouse theta: " << theta;
    W_param_mouse_theta.update_text(T.str());
  }
  
}


void AngleGui::S_param_recenter(void* e) {
  XEvent* ee = (XEvent*)e;
  if (ee->type != ButtonPress) return;
  double lam_radius = (param_lambda_u - param_lambda_l)/2.0;
  double theta_radius = (param_theta_u - param_theta_l)/2.0;
  param_lambda_u = lam_lambda + lam_radius;
  param_lambda_l = lam_lambda - lam_radius;
  param_theta_u = lam_theta + theta_radius;
  param_theta_l = lam_theta - theta_radius;
  reset_param_window();
}

void AngleGui::S_param_zoom_in(void* e) {
  XEvent* ee = (XEvent*)e;
  if (ee->type != ButtonPress) return;
  double lam_radius = (param_lambda_u - param_lambda_l)/2.0;
  double theta_radius = (param_theta_u - param_theta_l)/2.0;
  lam_radius *= 0.75;
  theta_radius *= 0.75;
  param_lambda_u = lam_lambda + lam_radius;
  param_lambda_l = lam_lambda - lam_radius;
  param_theta_u = lam_theta + theta_radius;
  param_theta_l = lam_theta - theta_radius;
  reset_param_window();
}

void AngleGui::S_param_zoom_out(void* e) {
  XEvent* ee = (XEvent*)e;
  if (ee->type != ButtonPress) return;
  double lam_radius = (param_lambda_u - param_lambda_l)/2.0;
  double theta_radius = (param_theta_u - param_theta_l)/2.0;
  lam_radius *= 1.0/0.75;
  theta_radius *= 1.0/0.75;
  param_lambda_u = lam_lambda + lam_radius;
  param_lambda_l = lam_lambda - lam_radius;
  param_theta_u = lam_theta + theta_radius;
  param_theta_l = lam_theta - theta_radius;
  reset_param_window();
}

void AngleGui::S_param_depth(void* e) {
  WidgetIntSelectorEvent* ee = (WidgetIntSelectorEvent*)e;
  if (ee->type == WidgetIntSelectorEvent::DOWN) {
    if (param_depth == 0) return;
    param_depth -= 1;
  } else {
    param_depth += 1;
  }
  W_param_depth.update_value(param_depth);
  draw_param();
}

void AngleGui::S_param_mesh(void* e) {
  WidgetIntSelectorEvent* ee = (WidgetIntSelectorEvent*)e;
  if (ee->type == WidgetIntSelectorEvent::DOWN) {
    if (param_pixel_group_pixels == 1) return;
    param_pixel_group_pixels /= 2;
  } else {
    param_pixel_group_pixels *= 2;
  }
  W_param_mesh.update_value(param_pixel_group_pixels);
  param_num_pixel_groups_width = param_pixel_width/param_pixel_group_pixels + 1;
  param_num_pixel_groups_height = param_pixel_height/param_pixel_group_pixels + 1;
  param_pixel_group_lambda_size = (param_lambda_u - param_lambda_l)/(double)param_num_pixel_groups_height;
  param_pixel_group_theta_size = (param_theta_u - param_theta_l)/(double)param_num_pixel_groups_width;
  reset_param_grid();
  draw_param();
}

void AngleGui::S_param_words(void* e) {
  XEvent* ee = (XEvent*)e;
  if (ee->type != ButtonPress) return;
  param_words = !param_words;
  W_param_words.checked = param_words;
  W_param_words.redraw();
  draw_param();
  draw_lam();
}

void AngleGui::S_param_inclusion(void* e) {
  XEvent* ee = (XEvent*)e;
  if (ee->type != ButtonPress) return;
  param_inclusion = !param_inclusion;
  W_param_inclusion.checked = param_inclusion;
  W_param_inclusion.redraw();
  draw_param();
  draw_lam();
}


void AngleGui::S_param_plot_points(void* e) {
  XEvent* ee = (XEvent*)e;
  if (ee->type != ButtonPress) return;
  param_plot_points = !param_plot_points;
  W_param_plot_points.checked = param_plot_points;
  W_param_plot_points.redraw();
  draw_param();
}

void AngleGui::S_lam_draw(void* ee) {
  XEvent* e = (XEvent*)ee;
  if (e->type == KeyPress) return;
}

void AngleGui::S_lam_depth(void* e) {
  WidgetIntSelectorEvent* ee = (WidgetIntSelectorEvent*)e;
  if (ee->type == WidgetIntSelectorEvent::DOWN) {
    if (lam_depth == 0) return;
    lam_depth -= 1;
  } else {
    lam_depth += 1;
  }
  W_lam_depth.update_value(lam_depth);
  draw_lam();
}

void AngleGui::S_lam_backward_depth(void* e) {
  WidgetIntSelectorEvent* ee = (WidgetIntSelectorEvent*)e;
  if (ee->type == WidgetIntSelectorEvent::DOWN) {
    if (lam_backward_depth == 0) return;
    lam_backward_depth -= 1;
  } else {
    lam_backward_depth += 1;
  }
  W_lam_backward_depth.update_value(lam_backward_depth);
  if (lam_backward_depth>0) {
    lam_limit_leaves_depth = 0;
    W_lam_limit_leaves_depth.update_value(lam_limit_leaves_depth);
  }
  draw_lam();
}


void AngleGui::S_lam_limit_leaves_depth(void* e) {
  WidgetIntSelectorEvent* ee = (WidgetIntSelectorEvent*)e;
  if (ee->type == WidgetIntSelectorEvent::DOWN) {
    if (lam_limit_leaves_depth == 0) return;
    lam_limit_leaves_depth -= 1;
  } else {
    lam_limit_leaves_depth += 1;
  }
  W_lam_limit_leaves_depth.update_value(lam_limit_leaves_depth);
  if (lam_limit_leaves_depth>0) {
    lam_backward_depth = 0;
    W_lam_backward_depth.update_value(lam_backward_depth);
  }
  draw_lam();
}


/****************************************************************************
 * translating points
 ****************************************************************************/
void AngleGui::param_LT_to_pixel_group(double lambda, double theta, int& x, int& y) {
  x = (theta - param_theta_l) / param_pixel_group_theta_size;
  y = param_num_pixel_groups_height - ((lambda - param_lambda_l) / param_pixel_group_lambda_size);
}

void AngleGui::param_pixel_group_to_LT(int x, int y, double& lambda, double& theta) {
  theta = param_theta_l + (x+0.5)*param_pixel_group_theta_size;
  lambda = param_lambda_u - ((y+0.5)*param_pixel_group_lambda_size);
}

void AngleGui::param_pixel_to_LT(int x, int y, double& lambda, double& theta) {
  theta = param_theta_l + (x+0.5)*( (param_theta_u-param_theta_l)/param_pixel_width );
  lambda = param_lambda_u - ((y+0.5)*( (param_lambda_u-param_lambda_l)/param_pixel_height ));
}

void AngleGui::param_LT_to_pixel(double lambda, double theta, int& x, int& y) {
  x = (theta - param_theta_l) / param_pixel_theta_size;
  y = param_pixel_height - ((lambda - param_lambda_l) / param_pixel_lambda_size);
}

void AngleGui::lam_coords_to_pixel(double dx, double dy, int& x, int& y) {
  double amount_per_pixel = 2.1/(double)lam_pixel_size;
  x = int( (dx - (-1.05)) / amount_per_pixel );
  y = int( lam_pixel_size - ((dy - (-1.05)) / amount_per_pixel) );
}

void AngleGui::lam_coords_to_pixel(double dx, double dy, short& x, short& y) {
  double amount_per_pixel = 2.1/(double)lam_pixel_size;
  x = short( (dx - (-1.05)) / amount_per_pixel );
  y = short( lam_pixel_size - ((dy - (-1.05)) / amount_per_pixel) );
}


void AngleGui::lam_coords_to_pixel(const std::vector<Point2d<float> >& floats, 
                                   std::vector<XPoint>& pixels) {
  pixels.resize(floats.size());
  for (int i=0; i<(int)floats.size(); ++i) {
    lam_coords_to_pixel(floats[i].x, floats[i].y, pixels[i].x, pixels[i].y);
  }
}

void AngleGui::lam_size_to_pixel(double dm, int& m) {
  double amount_per_pixel = 2.1/(double)lam_pixel_size;
  m = int( dm / amount_per_pixel );
}
  


/****************************************************************************
 * drawing 
 ****************************************************************************/
void AngleGui::draw_lam() {
  Lamination L(lam_lambda, lam_theta);
  Leaf initial_leaf(PI/2.0, 3.0*PI/2.0);
  std::vector<Leaf> leaf_stack(0);
  leaf_stack.push_back(initial_leaf);
  
  //std::cout << "Redrawing lamination at " << lam_lambda << " " << lam_theta << "\n\n\n";
  
  Widget& LP = W_lam_plot;
  XSetForeground(display, LP.gc, WhitePixel(display, screen));
  XFillRectangle(display, LP.p, LP.gc, 0, 0, LP.width, LP.height);
  XSetForeground(display, LP.gc, BlackPixel(display, screen));
  XDrawRectangle(display, LP.p, LP.gc, 0, 0, LP.width-1, LP.height-1);
  XSetFillStyle(display, LP.gc, FillSolid);
  
  //////////////////////////////////////// backward drawing
  //we do this first because we want the leaves to go over the filled area
  if (lam_backward_depth > 0) {
    ThickLeaf initial;
    std::vector<ThickLeaf> viable_initials;
    std::vector<ThickLeaf> images;
    L.compute_backwards(lam_backward_depth, initial, viable_initials, images);
    std::vector<Point2d<float> > polygon_points;
    std::vector<XPoint> draw_points;
    //display the initial leaf
    initial.polygon_points(polygon_points, 10);
    lam_coords_to_pixel(polygon_points, draw_points);
    XSetForeground(display, LP.gc, get_rgb_color(0.9375, 0.875, 0));
    XFillPolygon(display, LP.p, LP.gc, &draw_points[0], draw_points.size(), Nonconvex, CoordModeOrigin);
    
    //display the images
    for (int i=0; i<(int)images.size(); ++i) {
      images[i].polygon_points(polygon_points, 10);
      lam_coords_to_pixel(polygon_points, draw_points);
      XSetForeground(display, LP.gc, get_rgb_color(0.5 + (0.5/(double)(i+1)), 0, 0.5 + (0.5/(double)(i+1))) );
      XFillPolygon(display, LP.p, LP.gc, &draw_points[0], draw_points.size(), Nonconvex, CoordModeOrigin);
    }
    
    //display the viable initial leaves
    for (int i=0; i<(int)viable_initials.size(); ++i) {
      viable_initials[i].polygon_points(polygon_points, 10);
      lam_coords_to_pixel(polygon_points, draw_points);
      XSetForeground(display, LP.gc, get_rgb_color(0, 0.5 + (0.5/(double)(i+1)), 0.5 + (0.5/(double)(i+1))) );
      XFillPolygon(display, LP.p, LP.gc, &draw_points[0], draw_points.size(), Nonconvex, CoordModeOrigin);
    }

  } else if (lam_limit_leaves_depth > 0) {
    std::vector<ThickLeaf> initials, images;
    Graph inclusion_graph;
    bool complete;
    L.find_limit_leaves(initials, images, inclusion_graph, complete, lam_limit_leaves_depth, 0);
    //std::sort(images.begin(), images.end());
    //std::reverse(images.begin(), images.end());
    for (int i=0; i<(int)images.size(); ++i) {
      std::vector<Point2d<float> > polygon_points;
      images[i].polygon_points(polygon_points, 10);
      std::vector<XPoint> draw_points;
      lam_coords_to_pixel(polygon_points, draw_points);
      XSetForeground(display, LP.gc, get_rgb_color(0.5 + (0.5/(double)(i+1)), 0, 0.5 + (0.5/(double)(i+1))) );
      XFillPolygon(display, LP.p, LP.gc, &draw_points[0], draw_points.size(), Nonconvex, CoordModeOrigin);
      XSetForeground(display, LP.gc, get_rgb_color(0,0,0));
      XDrawLines(display, LP.p, LP.gc, &draw_points[0], draw_points.size(), CoordModeOrigin);
    }
    for (int i=0; i<(int)initials.size(); ++i) {
      std::vector<Point2d<float> > polygon_points;
      initials[i].polygon_points(polygon_points, 10);
      std::vector<XPoint> draw_points;
      lam_coords_to_pixel(polygon_points, draw_points);
      XSetForeground(display, LP.gc, get_rgb_color(0, 0.5 + (0.5/(double)(i+1)), 0.5 + (0.5/(double)(i+1))) );
      XFillPolygon(display, LP.p, LP.gc, &draw_points[0], draw_points.size(), Nonconvex, CoordModeOrigin);
      XSetForeground(display, LP.gc, get_rgb_color(0,0,0));
      XDrawLines(display, LP.p, LP.gc, &draw_points[0], draw_points.size(), CoordModeOrigin);
    }
        
  }
  
  /////////////////////////////////////// Main part drawing
  //draw the big black circle
  int ulx, uly, cw;
  lam_coords_to_pixel(-1, 1, ulx, uly);
  lam_size_to_pixel(2, cw);
  XSetForeground(display, LP.gc, BlackPixel(display, screen));
  XDrawArc(display, LP.p, LP.gc, ulx, uly, cw, cw, 0, 360*64);
  
  double epsilon = (2*PI - PI*lam_lambda)*0.99999999;
  int rcol = get_rgb_color(1,0,0);
  int bcol = get_rgb_color(0,0,0);
  int blcol = get_rgb_color(0,0,1);
  int gcol = get_rgb_color(0,1,0);
  
  //draw the domain of f
  XSetForeground(display, LP.gc, blcol);
  XSetLineAttributes(display, LP.gc, 4, LineSolid, CapButt, JoinMiter);
  int a1 = 64*(360.0/TWOPI)*L.f.domain_x;
  int a2 = 64*(360.0/TWOPI)*angle_dist(L.f.domain_x, L.f.domain_y);
  XDrawArc(display, LP.p, LP.gc, ulx, uly, cw, cw, a1, a2);
  
  //draw the domain of g
  XSetForeground(display, LP.gc, gcol);
  XSetLineAttributes(display, LP.gc, 4, LineSolid, CapButt, JoinMiter);
  a1 = 64*(360.0/TWOPI)*L.g.domain_x;
  a2 = 64*(360.0/TWOPI)*angle_dist(L.g.domain_x, L.g.domain_y);
  XDrawArc(display, LP.p, LP.gc, ulx-5, uly-5, cw+10, cw+10, a1, a2);

  
  while (leaf_stack.size() > 0) {
    Leaf ell = leaf_stack.back();
    leaf_stack.pop_back();
    if ((int)ell.word.size() == lam_depth && angle_diff(ell.x, ell.y) > epsilon) {
      //std::cout << "Good leaf: " << ell << "\n";
      XSetForeground(display, LP.gc, rcol);
      XSetLineAttributes(display, LP.gc, 2, LineSolid, CapButt, JoinMiter);
    } else {
      XSetForeground(display, LP.gc, bcol);
      XSetLineAttributes(display, LP.gc, 0, LineSolid, CapButt, JoinMiter);
    }
    double dcenterx, dcentery, dradius, a1, a_extent;
    ell.get_circle_data(dcenterx, dcentery, dradius, a1, a_extent);
    if (dradius < 0) {
      int ix1, iy1, ix2, iy2;
      lam_coords_to_pixel(cos(a1), sin(a1), ix1, iy1);
      lam_coords_to_pixel(cos(a_extent), sin(a_extent), ix2, iy2);
      //std::cout << "Drawing straight leaf at " << a1 << " " << a_extent << "\n";
      //std::cout << "Pixels: " << ix1 << " " << iy1 << " " << ix2 << " " << iy2 << "\n";
      XDrawLine(display, LP.p, LP.gc, ix1, iy1, ix2, iy2);
    } else {
      int icenterx, icentery, iradius;
      lam_coords_to_pixel(dcenterx, dcentery, icenterx, icentery);
      lam_size_to_pixel(dradius, iradius);
      double ia1 = (a1/TWOPI)*360*64;
      double ia_extent = (a_extent/TWOPI)*360*64;
      //std::cout << "Drawing curve leaf at "<< dcenterx << " " << dcentery << " of radius " << dradius << " with angles " << ia1 << " " << ia_extent << "\n";
      XDrawArc(display, LP.p, LP.gc, icenterx-iradius, icentery-iradius, 2*iradius, 2*iradius, ia1, ia_extent);
    }
    if ((int)ell.word.size() >= lam_depth) continue;
    bool can_f_act = L.check_can_act(ell, 0);
    bool can_g_act = L.check_can_act(ell, 1);
    if (can_f_act) leaf_stack.push_back(L.act_on_leaf(ell, 0));
    if (can_g_act) leaf_stack.push_back(L.act_on_leaf(ell, 1));
  }
  XSetLineAttributes(display, LP.gc, 0, LineSolid, CapButt, JoinMiter);
  LP.redraw();
  
  recompute_lam_data();
  
}

void AngleGui::recompute_lam_data() {
  std::stringstream T;
  T.precision(15);
  T.str("");
  T << "Lambda: " << lam_lambda;
  W_lam_lambda_label.update_text(T.str());
  T.str("");
  T << "Theta: " << lam_theta;
  W_lam_theta_label.update_text(T.str());
  Lamination L(lam_lambda, lam_theta);
  LamType ellt;
  int difficulty;
  std::vector<Leaf> viable_leaves;
  if (param_words) {
    L.compute_lam_type_with_words(lam_depth, ellt, difficulty, viable_leaves);
  } else {
    L.compute_lam_type(lam_depth, ellt, difficulty);
  }
  T.str("");
  T << "Lamination type: ";
  T << (ellt == NONE ? "none" : (ellt == PROPER ? "proper" : "cut point"));
  if (ellt == CUT_POINT && param_words) {
    T << " (viable leaves: " << viable_leaves.size() << ")";   
    std::cout << "\n\nViable leaves:\n";
    std::sort(viable_leaves.begin(), viable_leaves.end());
    for (int i=0; i<(int)viable_leaves.size(); ++i) {
      std::cout << i << ": " << viable_leaves[i] << "\n";
    }
  } else {
    T << " (difficulty: " << difficulty << ")";
  }
  W_lam_lam_type.update_text(T.str());
  
  if (lam_limit_leaves_depth > 0) {
    std::vector<ThickLeaf> initials, images;
    Graph inclusion_graph;
    bool complete;
    L.find_limit_leaves(initials, images, inclusion_graph, complete, lam_limit_leaves_depth, 1);
    std::cout << "Leading eigenvalue: " << inclusion_graph.approximate_leading_eigenvalue(inclusion_graph.num_verts) << "\n";
    std::cout << inclusion_graph.mathematica_string() << "\n";
  }
  
}


void AngleGui::reset_param_grid() {
  param_grid = std::vector<std::vector<Point3d<int> > >(param_num_pixel_groups_width);
  for (int i=0; i<(int)param_num_pixel_groups_width; ++i) {
    param_grid[i].resize(param_num_pixel_groups_height);
  }
}



void AngleGui::reset_param_window() {
  param_pixel_group_lambda_size = (param_lambda_u - param_lambda_l)/(double)param_num_pixel_groups_height;
  param_pixel_group_theta_size = (param_theta_u - param_theta_l)/(double)param_num_pixel_groups_width;
  param_pixel_lambda_size = (param_lambda_u - param_lambda_l)/(double)param_pixel_height;
  param_pixel_theta_size = (param_theta_u - param_theta_l)/(double)param_pixel_width;
  draw_param();
}


/*
p.x = 0: lamination doesn't exist
p.x = 1: lamination is proper
p.x = 2: lamination has cut point

p.y: difficulty in computing leaves to the depth, or number of viable 
leaves larger than epsilon at depth, depending on param_words

p.z: maximal rank of free semigroup in inclusion graph (negative if incomplete)

*/
int AngleGui::compute_color_from_grid(const Point3d<int>& p) {
  if (p.x == 0) 
    return get_rgb_color(1,1,1);
  
  if (p.x == 1) {
    if (abs(p.z) > 0) {
      return get_rgb_color(1,1,0); //this means we found limit leaves and no middle leaf limits
    }
    if (param_words) {
      return get_rgb_color(0,0,1);
    } else {
      return get_rgb_color(0,0, 1.0 - 1.0/double(p.y));
    }
  }  
  //otherwise, p.x = 2
  if (p.z==0 && param_inclusion) {
    return get_rgb_color(1,0.5,0); //this means we found middle leaf limits and no limit leaves (error)
  }
  
  if (p.z == 0) {
    return get_rgb_color(1.0 - 1.0/double(p.y/5 + 1) , 0 , 0);
  } else {
    if (p.z > 0) {
      return get_rgb_color(0, 1.0 - 1.0/double(double(p.z+1)/3 + 0.34) , 0);
    } else {
      return get_rgb_color(0, 1.0 - 1.0/double(abs(p.z)) , 1.0 - 1.0/double(abs(p.z)));
    }
  }

}

void AngleGui::draw_param() {
  for (int i=0; i<param_num_pixel_groups_width; ++i) {
    for (int j=0; j<param_num_pixel_groups_height; ++j) {
      double lambda, theta;
      param_pixel_group_to_LT(i, j, lambda, theta);
      if (lambda < 1 || 2 < lambda || theta < 0 || 2*PI < theta) {
        param_grid[i][j].x = param_grid[i][j].y = param_grid[i][j].z = 0;
      } else {
        Lamination L(lambda, theta);
        LamType ellt;
        int difficulty;
        std::vector<Leaf> viable_leaves(0);
        Graph inclusion_graph;
        std::vector<ThickLeaf> images, initials;
        bool inclusion_complete;
        if (param_words) {
          L.compute_lam_type_with_words(param_depth, ellt, difficulty, viable_leaves);
        } else {
          L.compute_lam_type(param_depth, ellt, difficulty);
        }
        param_grid[i][j].x = (ellt == NONE ? 0 : (ellt == PROPER ? 1 : 2));
        param_grid[i][j].y = (param_words ? viable_leaves.size() : difficulty);
        if (param_inclusion && ellt != NONE) {
          L.find_limit_leaves(initials, images, inclusion_graph, inclusion_complete, param_depth, 0);
          param_grid[i][j].z = (inclusion_complete ? 1 : -1)*int(4*(inclusion_graph.approximate_leading_eigenvalue(inclusion_graph.num_verts)));
          //std::cout << i << "," << j << ": " << param_grid[i][j].z << "\n";
          //std::cout << inclusion_graph << "\n";
        } else {     
          param_grid[i][j].z = 0;
        }
      }
      
      int col = compute_color_from_grid(param_grid[i][j]);
      Widget& PP = W_param_plot;
      XSetForeground(display, PP.gc, col);
      XFillRectangle(display, PP.p, PP.gc, i*param_pixel_group_pixels, 
                                           j*param_pixel_group_pixels, 
                                           param_pixel_group_pixels, 
                                           param_pixel_group_pixels);
      XCopyArea(display, PP.p, main_window, PP.gc, i*param_pixel_group_pixels, 
                                                   j*param_pixel_group_pixels, 
                                                   param_pixel_group_pixels, 
                                                   param_pixel_group_pixels,
                                                   PP.ul.x + i*param_pixel_group_pixels,
                                                   PP.ul.y + j*param_pixel_group_pixels);
    }
  }
  
  //if we're supposed to, draw all the points
  if (param_plot_points) {
    std::ifstream f;
    double ell, tee;
    int bcol = get_rgb_color(0,0,0);
    f.open("points_to_plot.txt");
    if (f.is_open()) {
      while (0 != (f >> ell) && 0 != (f>>tee)) {
        draw_param_point(ell, tee, bcol);
      }
    }
  }
  
  draw_highlighted_point();
}

//this updates the highlighted point to agree with whatever the current lambda, theta are
void AngleGui::reset_highlighted_point(double L, double T) {
  //paint over the currently highlighted point
  int ig, jg;
  param_LT_to_pixel_group(lam_lambda, lam_theta, ig, jg);
  int istart = (ig-3 < 0 ? 0 : ig-3);
  int iend = (ig+3 > param_num_pixel_groups_width-1 ? param_num_pixel_groups_width-1 : ig+2);
  int jstart = (jg-3 < 0 ? 0 : jg-3);
  int jend = (jg+3 > param_num_pixel_groups_height-1 ? param_num_pixel_groups_height-1 : jg+2);
  Widget& PP = W_param_plot;
  for (int i=istart; i<iend; ++i) {
    for (int j=jstart; j<jend; ++j) {
      XSetForeground(display, PP.gc, compute_color_from_grid(param_grid[i][j]));
      XFillRectangle(display, PP.p, PP.gc, i*param_pixel_group_pixels, 
                                           j*param_pixel_group_pixels, 
                                           param_pixel_group_pixels, 
                                           param_pixel_group_pixels);
      XCopyArea(display, PP.p, main_window, PP.gc, i*param_pixel_group_pixels, 
                                                   j*param_pixel_group_pixels, 
                                                   param_pixel_group_pixels, 
                                                   param_pixel_group_pixels,
                                                   PP.ul.x + i*param_pixel_group_pixels,
                                                   PP.ul.y + j*param_pixel_group_pixels);
    }
  }
  lam_lambda = L;
  lam_theta = T;
  draw_highlighted_point();
  draw_lam();
}
  

//this just draws the red point on the parameter space
void AngleGui::draw_highlighted_point() {
  int ip,jp;
  Widget& PP = W_param_plot;
  param_LT_to_pixel(lam_lambda, lam_theta, ip, jp);
  XSetForeground(display, PP.gc, get_rgb_color(0,1,1));
  XFillRectangle(display, PP.p, PP.gc, ip-2, jp-2, 5, 5);
  XCopyArea(display, PP.p, main_window, PP.gc, ip-2, jp-2, 5, 5, PP.ul.x + ip-2, PP.ul.y + jp-2);
}

//draw the given point in the given color
void AngleGui::draw_param_point(double L, double T, int col) {
  int ip,jp;
  Widget& PP = W_param_plot;
  param_LT_to_pixel(L, T, ip, jp);
  XSetForeground(display, PP.gc, col);
  XFillRectangle(display, PP.p, PP.gc, ip-2, jp-2, 4, 4);
  XCopyArea(display, PP.p, main_window, PP.gc, ip-2, jp-2, 4, 4, PP.ul.x + ip-2, PP.ul.y + jp-2);
}

/****************************************************************************
 * Widget wrangling
 ****************************************************************************/
 void AngleGui::pack_widget_upper_right(const Widget* w1, Widget* w2) {
  //figure out where it can go
  int desired_x,desired_y;
  if (w1 != NULL) {
    desired_x = w1->ul.x + w1->width;
    desired_y = w1->ul.y;
  } else {
    desired_x = 0;
    desired_y = 0;
  }
  
  
  //std::cout << "Packing widget of size " << w2->width << " " << w2->height << "\n";
  //std::cout << "Desired x: " << desired_x << "\n";
  
  //go through and check the other widgets to see how 
  //far down they obstruct this one
  int greatest_y_obstruction = 0;
  for (int i=0; i<(int)widgets.size(); ++i) {
    if (widgets[i] == w1) continue;
    if (widgets[i]->ul.x == desired_x && 
        widgets[i]->ul.y + widgets[i]->height > greatest_y_obstruction) {
      greatest_y_obstruction = widgets[i]->ul.y + widgets[i]->height;
      //std::cout << "Found widget " << i << " obstructs to height " << greatest_y_obstruction << "\n";
    }
  }
  if (greatest_y_obstruction + w2->height > main_window_height) {
    //std::cout << "Cannot pack widget -- too tall\n";
    return;
  }
  int y = (desired_y > greatest_y_obstruction ? desired_y : greatest_y_obstruction);
  
  //now determine whether we have to shove it over to make room
  int greatest_x_obstruction = desired_x;
  for (int i=0; i<(int)widgets.size(); ++i) {
    if (widgets[i]->ul.y + widgets[i]->height > y && 
        widgets[i]->ul.y < y + w2->height && 
        widgets[i]->ul.x + widgets[i]->width > greatest_x_obstruction) {
      greatest_x_obstruction = widgets[i]->ul.x + widgets[i]->width;
    }
  }
  int x = greatest_x_obstruction;
  
  //find the position
  w2->ul = Point2d<int>(x, y);
  
  //std::cout << "Packed widget at " << w2->ul << "\n";
  
  //record it in the list of widgets
  widgets.push_back(w2);
}


void AngleGui::detach_widget(Widget* w) {
  for (int i=0; i<(int)widgets.size(); ++i) {
    if (widgets[i] == w) {
      widgets.erase(widgets.begin()+i);
      break;
    }
  }
  w->clear();
}

/*****************************************************************************
 * main window creation 
 *****************************************************************************/
void AngleGui::pack_window() {
  //if the window exists, destroy it
  if (main_window_initialized) {
    XDestroyWindow(display, main_window);
  }
  //compute the size of the window
  int display_width = XDisplayWidth(display, screen);
  int display_height = XDisplayHeight(display, screen);
  int x = -1;
  if (2*display_height + 300 > display_width) {
    x = display_width/2 - 150;
  } else {
    x = display_height;
  }
  main_window_height = x + 80;
  main_window_width = 2*x + 300;
  //create the window
  main_window = XCreateSimpleWindow(display, 
                                    RootWindow(display, screen), 20, 20,
                                    main_window_width, main_window_height, 4,
                                    BlackPixel(display, screen), WhitePixel(display, screen));
  XSelectInput(display, main_window, ExposureMask |
                                     PointerMotionMask |
                                     KeyPressMask |
                                     ButtonPressMask |
                                     ButtonReleaseMask );
  XMapWindow(display, main_window);
  //wait until the window is actually mapped
  while (true) {  
    XEvent e;
    XNextEvent(display, &e);
    if (e.type == Expose) break;
  }
  main_window_initialized = true;
  
  //create the wigets
  W_param_plot = WidgetDraw(this, x, x, &AngleGui::S_param_draw);
  W_param_recenter = WidgetButton(this, -1, -1, "Recenter", &AngleGui::S_param_recenter);
  W_param_zoom_in = WidgetButton(this, -1, -1, "In", &AngleGui::S_param_zoom_in);
  W_param_zoom_out = WidgetButton(this, -1, -1, "Out", &AngleGui::S_param_zoom_out);
  W_param_depth = WidgetIntSelector(this, -1, "Depth:", 10, &AngleGui::S_param_depth);
  W_param_mesh = WidgetIntSelector(this, -1, "Mesh:", 4, &AngleGui::S_param_mesh);
  W_param_words = WidgetCheck(this, -1, "Num words colors", param_words, &AngleGui::S_param_words);
  W_param_inclusion = WidgetCheck(this, -1, "Inclusion colors", param_words, &AngleGui::S_param_inclusion);
  W_param_plot_points = WidgetCheck(this, -1, "Plot points", param_plot_points, &AngleGui::S_param_plot_points);
  W_param_mouse_lambda = WidgetText(this, 300, "Mouse lambda: (initializing)");
  W_param_mouse_theta = WidgetText(this, 300, "Mouse theta: (initializing)");
  
  W_lam_plot = WidgetDraw(this, x, x, &AngleGui::S_lam_draw);
  W_lam_depth = WidgetIntSelector(this, -1, "Depth:", lam_depth, &AngleGui::S_lam_depth);
  W_lam_backward_depth = WidgetIntSelector(this, -1, "Backwards:", lam_backward_depth, &AngleGui::S_lam_backward_depth);
  W_lam_limit_leaves_depth = WidgetIntSelector(this, -1, "Limit leaves:", lam_limit_leaves_depth, &AngleGui::S_lam_limit_leaves_depth);
  W_lam_lambda_label = WidgetText(this, 300, "Lambda: (initializing)");
  W_lam_theta_label = WidgetText(this, 300, "Theta: (initializing)");
  W_lam_lam_type = WidgetText(this, 300, "Lamination type: (initializing)");
  
  //record the sizes
  param_pixel_width = W_param_plot.width;
  param_pixel_height = W_param_plot.height;
  param_pixel_group_pixels = 4;
  param_num_pixel_groups_width = param_pixel_width/param_pixel_group_pixels + 1;
  param_num_pixel_groups_height = param_pixel_height/param_pixel_group_pixels + 1;
  param_pixel_group_lambda_size = (param_lambda_u - param_lambda_l)/(double)param_num_pixel_groups_height;
  param_pixel_group_theta_size = (param_theta_u - param_theta_l)/(double)param_num_pixel_groups_width;
  param_pixel_lambda_size = (param_lambda_u - param_lambda_l)/(double)param_pixel_height;
  param_pixel_theta_size = (param_theta_u - param_theta_l)/(double)param_pixel_width;
  reset_param_grid();
  
  lam_pixel_size = W_lam_plot.width;
  
  //pack the widgets
  widgets.resize(0);
  pack_widget_upper_right(NULL, &W_param_plot);
  pack_widget_upper_right(NULL, &W_param_mouse_lambda);
  pack_widget_upper_right(NULL, &W_param_mouse_theta);
  pack_widget_upper_right(&W_param_plot, &W_param_recenter);
  pack_widget_upper_right(&W_param_recenter, &W_param_zoom_in);
  pack_widget_upper_right(&W_param_zoom_in, &W_param_zoom_out);
  pack_widget_upper_right(&W_param_plot, &W_param_depth);
  pack_widget_upper_right(&W_param_plot, &W_param_mesh);
  pack_widget_upper_right(&W_param_plot, &W_param_words);
  pack_widget_upper_right(&W_param_plot, &W_param_inclusion);
  pack_widget_upper_right(&W_param_plot, &W_param_plot_points);
  
  pack_widget_upper_right(&W_param_zoom_out, &W_lam_plot);
  pack_widget_upper_right(&W_lam_plot, &W_lam_depth);
  pack_widget_upper_right(&W_lam_plot, &W_lam_backward_depth);
  pack_widget_upper_right(&W_lam_plot, &W_lam_limit_leaves_depth);
  pack_widget_upper_right(&W_param_words, &W_lam_lambda_label);
  pack_widget_upper_right(&W_param_words, &W_lam_theta_label);
  pack_widget_upper_right(&W_param_words, &W_lam_lam_type);
  
  
  //draw all the widgets
  for (int i=0; i<(int)widgets.size(); ++i) {
    widgets[i]->initial_draw();
  }
  
  //draw the parameter space and lamination
  draw_param();
  draw_lam();
  
}



void AngleGui::main_loop() {
  XEvent e;
  while (true) {
    XNextEvent(display, &e);
    //if it was the keyboard, we deal with it here
    if (e.type == KeyPress) {
      if(XLookupKeysym(&e.xkey, 0) == XK_q){ 
        break;
      }
     
    //if it involves the mouse, we find the appropriate 
    //widget to send it off to
    } else if (e.type == ButtonPress || e.type == MotionNotify) {
      for (int i=0; i<(int)widgets.size(); ++i) {
        if (widgets[i]->contains_pixel( e.xbutton.x, e.xbutton.y) &&
            widgets[i]->event_signal != NULL) {
          //(this->*(widgets[i]->click_signal))(&e);
          widgets[i]->run_event_signal(&e);
          break;
        }
      }
    
    } else if (e.type == Expose) {
      Point2d<int> expose_ul( e.xexpose.x, e.xexpose.y );
      int ewidth = e.xexpose.width;
      int eheight = e.xexpose.height;
      for (int i=0; i<(int)widgets.size(); ++i) {
        if (widgets[i]->intersects_rectangle(expose_ul, ewidth, eheight)) {
          widgets[i]->redraw();
        }
      }
    } 
  }
  
}





void AngleGui::launch() {
  
  param_theta_l = 0.0;
  param_theta_u = TWOPI;
  param_lambda_l = 1.0;
  param_lambda_u = 2.0;
  param_depth = 10;
  param_words = false;
  param_inclusion = false;
  param_plot_points = false;
  
  lam_depth = 10;
  lam_backward_depth = 0;
  lam_limit_leaves_depth = 0;
  lam_lambda = 1.5;
  lam_theta = 0;
  
  display = XOpenDisplay(NULL);
  screen = DefaultScreen(display);
  if (display == NULL) {
    std::cout << "Failed to open display\n";
    return;
  }
  main_window_initialized = false;
  
  pack_window();
  
  main_loop();
  
  XCloseDisplay(display);
}


















