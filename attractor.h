/*
 *  attractor.h
 *  LorenzGL_verHY
 *
 *  Created by Hao Ye on 9/24/10.
 *  Copyright 2010 UCSD. All rights reserved.
 *
 */
#ifndef ATTRACTOR_H
#define ATTRACTOR_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include <math.h>
#include "/usr/X11/include/png.h"
#include "CoreFoundation/CoreFoundation.h"

enum tracer {PROJECT, TRACE, NONE};
enum draw_mode {MANIFOLD, TIME_SERIES, LAGS, 
    RECONSTRUCTION, SHADOW, 
    UNIVARIATE, UNIVARIATE_TS, 
    XMAP, XMAP_TS, GENERIC_RECONSTRUCTION, TAKENS};
enum ode_mode {EULER, RK4};

#define X_LABEL_TEXTURE 0
#define Y_LABEL_TEXTURE 4
#define Z_LABEL_TEXTURE 8
#define M_LABEL_TEXTURE 12
#define TAU_LABEL_TEXTURE 16
#define VIEW_1_LABEL_TEXTURE 18
#define VIEW_2_LABEL_TEXTURE 19
#define VIEW_3_LABEL_TEXTURE 20
#define VIEW_4_LABEL_TEXTURE 23
#define VIEW_5_LABEL_TEXTURE 26
#define VIEW_6_LABEL_TEXTURE 29
#define VIEW_7_LABEL_TEXTURE 32
#define EQUATIONS_LABEL_TEXTURE 38
#define TAKENS_THEOREM_TEXTURE 39


struct texture_2d
{
    GLuint texture_id;
    int width;
    int height;
};

struct texture_label
{
    int index;
    double billboard_matrix[16];
    double z_pos;
    double scale;
    int x_peg;
    int y_peg;
};

bool operator<(const texture_label & a, const texture_label & b);

using namespace std;

class attractor
{
public:
	attractor(const int max_frames);
	~attractor();
	draw_mode VIEW;
	
private:   
    static const double point_color[4];
    static const double neighbor_color[4];
    static const double pred_color[4];
    static const double PI;
    static const double LAG_SAT;
    static const double LAG2_SAT;
    static const double LINE_WIDTH;
    static const double SMALL_LINE_WIDTH;
    static const double POINT_WIDTH;
    
    static const double sigma;
    static const double rho;
    static const double beta;
    static const double dt;
    static const double d;
    static const double draw_fraction;
    static const double init_distance;
    
    static const double xmap_attractor_scale;
    
    // textures
    vector<texture_2d> my_textures;
    vector<texture_label> texture_queue;
    
	// data
	int num_points;
	vector<double> x;
	vector<double> y;
	vector<double> z;
    vector<double> phi_x;
    vector<double> phi_y;
    vector<double> phi_z;
    vector<int> knot;
    int nn_num, nn_skip;
    vector<vector<int> > x_nn_indices;
    vector<vector<double> > x_nn_weights;
    vector<vector<int> > y_nn_indices;
    vector<vector<double> > y_nn_weights;
    vector<vector<int> > z_nn_indices;
    vector<vector<double> > z_nn_weights;
    vector<double> x_xmap_y;
    vector<double> x_xmap_z;
    vector<double> y_xmap_x;
    vector<double> y_xmap_z;
    vector<double> z_xmap_x;
    vector<double> z_xmap_y;
    vector<double> x_forecast;
    vector<double> x_forecast_lag_1;
    vector<double> x_forecast_lag_2;
    vector<double> y_forecast;
    vector<double> y_forecast_lag_1;
    vector<double> y_forecast_lag_2;
    vector<double> z_forecast;
    vector<double> z_forecast_lag_1;
    vector<double> z_forecast_lag_2;
	
	// params
	int lag_dim;
	int pred_dim;
    int x_dim;
    int y_dim;
    int z_dim;
    int x_lag;
    int y_lag;
    int z_lag;
	int tau;
    int tp;
    double curr_x, curr_y, curr_z;
    double x_scale, y_scale, z_scale;
    double texture_scale;
    int window_width;
    int window_height;
    double rot_matrix[16];
	
	// switches
    ode_mode lorenz_sim_mode;
    bool DRAW_CONE;
	int DEBUG;
    bool TSVIEW;
    bool COLOR_METHOD;
    bool MANIFOLD_LABEL;
    bool SPLIT_VIEW;
	tracer x_tracer;
	tracer y_tracer;
	tracer z_tracer;
	int x_start_time;
	int y_start_time;
	int z_start_time;
	
	// rotation point
	double vx, vy, vz;
	double theta;
	double scale;
	
private:
    double f(const double x, const double y, const double z) {return sigma * (y - x);}
    double g(const double x, const double y, const double z) {return rho * x - x * z - y;}
    double h(const double x, const double y, const double z) {return x * y - beta * z;}
	double dist(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2) {return sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2));}
    double dist_to_curr(const double x, const double y, const double z) {return dist(x, y, z, curr_x, curr_y, curr_z);};
    void EULER_sim()
    {
        double xx, yy, zz;
        for(int i = 0; i < num_points-1; i++)
        {
            xx = x[i];
            yy = y[i];
            zz = z[i];
            x[i+1] = xx + f(xx, yy, zz) * dt;
            y[i+1] = yy + g(xx, yy, zz) * dt;
            z[i+1] = zz + h(xx, yy, zz) * dt;
        }
        return;
    }
    void RK4_sim()
    {
        double xx, yy, zz;
        double half_dt = dt/2.0;
        double kx1, kx2, kx3, kx4, ky1, ky2, ky3, ky4, kz1, kz2, kz3, kz4;
        for(int i = 0; i < num_points-1; i++)
        {
            xx = x[i];
            yy = y[i];
            zz = z[i];
            kx1 = f(xx, yy, zz);
            ky1 = g(xx, yy, zz);
            kz1 = h(xx, yy, zz);
            kx2 = f(xx + half_dt*kx1, yy + half_dt*ky1, zz + half_dt*kz1);
            ky2 = g(xx + half_dt*kx1, yy + half_dt*ky1, zz + half_dt*kz1);
            kz2 = h(xx + half_dt*kx1, yy + half_dt*ky1, zz + half_dt*kz1);
            kx3 = f(xx + half_dt*kx2, yy + half_dt*ky2, zz + half_dt*kz2);
            ky3 = g(xx + half_dt*kx2, yy + half_dt*ky2, zz + half_dt*kz2);
            kz3 = h(xx + half_dt*kx2, yy + half_dt*ky2, zz + half_dt*kz2);
            kx4 = f(xx + dt*kx3, yy + dt*ky3, zz + dt*kz3);
            ky4 = g(xx + dt*kx3, yy + dt*ky3, zz + dt*kz3);
            kz4 = h(xx + dt*kx3, yy + dt*ky3, zz + dt*kz3);
            
            x[i] = xx + dt / 6.0 * (kx1 + 2*kx2 + 2*kx3 + kx4);
            y[i] = yy + dt / 6.0 * (ky1 + 2*ky2 + 2*ky3 + ky4);
            z[i] = zz + dt / 6.0 * (kz1 + 2*kz2 + 2*kz3 + kz4);
        }
        return;
    }
    void draw_takens();
	void draw_xmap_ts(const int frame);
    void draw_xmap(const int frame);
    void draw_univariate_ts(const int frame);
	void draw_univariate(const int frame);
	void draw_shadow(const int frame);
	void draw_reconstruction(const int frame);
	void draw_generic_reconstruction(const int frame);
	void draw_lagged_time_series(const int frame);
	void draw_time_series(const int frame);
	void draw_manifold(const int frame);
    vector<double> draw_xmap_manifold(const int frame, const vector<int> & nn_indices, const vector<double> & nn_weights, const int var, bool swap_xy);
    void draw_xmap_generic(const int frame, const int NW_manifold, const int SW_manifold, const int NE_manifold, bool ts_trace);
    void draw_univariate_ts(int frame, vector<double>* ts);
    void draw_xmap_ts(int frame, const int lag, const int skip, const int dim, 
                      const double sat, const double line_width, vector<double>* ts);
    void draw_half_ts(int frame, const int lag, const int skip, const int dim, 
                      const double sat, const double line_width);
	void draw_ts(int frame, const int lag, const int skip, const int dim, 
                 const double sat, const double line_width);
    void draw_embedding(vector<double>::iterator x_i, vector<double>::iterator y_i, vector<double>::iterator z_i, const int frame);
	void draw_tracers(const int frame);
	void draw_axes(bool lag, int lag_dim);
    void draw_lag_axis(int direction, int dim, int lag);
    void draw_axis(const int axis, const double r, const double g, const double b, 
                   const double delta_line_width, const int texture_index, const int lag);
    void draw_labels();
    void enqueue_label(double pos_x, double pos_y, double pos_z, int texture_index, double scale, int x_peg, int y_peg);
    void draw_curve(const double x1, const double y1, const double z1, 
                    const double x2, const double y2, const double z2);
    void col(const double x, const double y, const double z);
    double N(const int i, const int k, const double u);
    void generate_movie();
    void generate_data();
	void transform_data();
    void generate_xmaps();
    void generate_forecasts();
    void load_textures();
    GLuint load_texture(const string filename, int &width, int &height);
    void find_neighbors(const int dim);
    
public:
	void init(bool MOVIE_MODE);
    void reset_rot_matrix();
	void rotate(const double rx, const double ry, const double rz);
	void translate(const double tx, const double ty, const double tz);
	void draw(double & runtime);
	void trace_x(const int frame);
	void trace_y(const int frame);
	void trace_z(const int frame);
	void set_view(const int new_view);
	void set_lagview(const int dim);
	void inc_xview();
	void inc_yview();
	void inc_zview();
	void toggle_predview();
    void toggle_tsview();
    void toggle_color_method();
    void toggle_manifold_label();
    void inc_xtau();
    void inc_ytau();
    void toggle_split_view();
	void change_tau(const int delta);
	void debug_toggle();
	void change_scale(const double new_scale);
    void set_window_size(const int width, const int height);
	
};

#endif