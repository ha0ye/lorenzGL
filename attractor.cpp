/*
 *  attractor.cpp
 *  LorenzGL_verHY
 *
 *  Created by Hao Ye on 9/24/10.
 *  Copyright 2010 UCSD. All rights reserved.
 *
 */

#include "attractor.h"

const double attractor::point_color[4] = {0.0, 0.0, 0.0, 0.8};
const double attractor::neighbor_color[4] = {0.0, 0.4, 0.8, 0.8};
const double attractor::pred_color[4] = {0.5, 0.0, 1.0, 0.8};
const double attractor::PI = 3.1415926535897932384626433832795;
const double attractor::LAG_SAT = 0.7;
const double attractor::LAG2_SAT = 0.4;
const double attractor::LINE_WIDTH = 0.6;
const double attractor::SMALL_LINE_WIDTH = 1.0;
const double attractor::POINT_WIDTH = 9.0;
const double attractor::sigma = 10;
const double attractor::rho = 28;
const double attractor::beta = 8.0/3;
const double attractor::dt = 0.01;
const double attractor::d = 0.85;
const double attractor::draw_fraction = 0.2;
const double attractor::init_distance = 6;
const double attractor::xmap_attractor_scale = 0.60;

bool operator < (const texture_label & a, const texture_label & b)
{
	return a.z_pos < b.z_pos;
}

attractor::attractor(const int max_frames)
{
    // initialize switches
    lorenz_sim_mode = EULER;
    DRAW_CONE = true;
	DEBUG = false;
    TSVIEW = false;
    COLOR_METHOD = false;
    MANIFOLD_LABEL = true;
    SPLIT_VIEW = false;
	x_tracer = NONE;
	y_tracer = NONE;
	z_tracer = NONE;
	scale = 1.0;
    
	// initialize params
    tau = 7;
    tp = 7;
	x_start_time = 0;
	y_start_time = 0;
	z_start_time = 0;
	lag_dim = 1;
	pred_dim = 2;
    x_dim = 1;
    y_dim = 1;
    z_dim = 1;
    x_lag = 0;
    y_lag = 1;
    z_lag = 2;
    nn_num = 4;
    nn_skip = 5;
    
	// initialize vectors
	num_points = max_frames;
	x.resize(num_points);
    y.resize(num_points);
    z.resize(num_points);
    x_nn_indices.resize(num_points);
    x_nn_weights.resize(num_points);
    y_nn_indices.resize(num_points);
    y_nn_weights.resize(num_points);
    z_nn_indices.resize(num_points);
    z_nn_weights.resize(num_points);
    
    x_xmap_y.resize(num_points, 0);
    x_xmap_z.resize(num_points, 0);
    y_xmap_x.resize(num_points, 0);
    y_xmap_z.resize(num_points, 0);
    z_xmap_x.resize(num_points, 0);
    z_xmap_y.resize(num_points, 0);
    
    x_forecast.resize(num_points, 0);
    x_forecast_lag_1.resize(num_points, 0);
    x_forecast_lag_2.resize(num_points, 0);
    y_forecast.resize(num_points, 0);
    y_forecast_lag_1.resize(num_points, 0);
    y_forecast_lag_2.resize(num_points, 0);
    z_forecast.resize(num_points, 0);
    z_forecast_lag_1.resize(num_points, 0);
    z_forecast_lag_2.resize(num_points, 0);
}

attractor::~attractor()
{
}

void attractor::draw_takens()
{
    double depth = -1;
    double low_x, low_y, high_x, high_y;
    double scale;
    
    texture_2d curr_texture = my_textures[TAKENS_THEOREM_TEXTURE];
    
    glBindTexture(GL_TEXTURE_2D, curr_texture.texture_id);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if(0.75*window_width > window_height)
    {
        low_x = -1.5*window_width/window_height;
        high_x = 1.5*window_width/window_height;
        low_y = -1.5;
        high_y = 1.5;
    }
    else
    {
        low_x = -2.0;
        high_x = 2.0;
        low_y = -2.0*window_height/window_width;
        high_y = 2.0*window_height/window_width;
    }
    glOrtho(low_x, high_x, low_y, high_y, -1.5, 1.5);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    if(curr_texture.height > curr_texture.width*0.75)
    {
        scale = 1.5 / curr_texture.height;
    }
    else
    {
        scale = 2.0 / curr_texture.width;
    }
    
    glColor3d(1.0, 1.0, 1.0);
    glEnable(GL_TEXTURE_2D);
    glBegin(GL_QUADS);
    glNormal3f(0.0, 0.0, 1.0);
    glTexCoord2d(0.0, 0.0); glVertex3d(-curr_texture.width*scale, -curr_texture.height*scale, depth);
    glTexCoord2d(0.0, 1.0); glVertex3d(-curr_texture.width*scale, curr_texture.height*scale, depth);
    glTexCoord2d(1.0, 1.0); glVertex3d(curr_texture.width*scale, curr_texture.height*scale, depth);
    glTexCoord2d(1.0, 0.0); glVertex3d(curr_texture.width*scale, -curr_texture.height*scale, depth);
    glEnd();
    glDisable(GL_TEXTURE_2D);
    
    return;
}

void attractor::draw_xmap_ts(const int frame)
{	
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    if(0.75*window_width > window_height)
    {
        glOrtho(-1.5*window_width/window_height, 1.5*window_width/window_height, -1.5, 1.5, -1.5, 1.5);
    }
    else
    {
        glOrtho(-2.0, 2.0, -2.0*window_height/window_width, 2.0*window_height/window_width, -1.5, 1.5);
    }
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMultMatrixd(rot_matrix);
    
    double ts_delta_x = 0.8, delta_x, delta_y = 1.1;
    draw_xmap_generic(frame, pred_dim, lag_dim, -1, true);
    
    vector<double>* prediction_ts;
    switch(lag_dim)
    {
        case 1:
            if(pred_dim == 2)
                prediction_ts = &x_xmap_y;
            else
                prediction_ts = &x_xmap_z;
            break;
        case 2:
            if(pred_dim == 1)
                prediction_ts = &y_xmap_x;
            else
                prediction_ts = &y_xmap_z;
            break;
        case 3:
            if(pred_dim == 2)
                prediction_ts = &z_xmap_y;
            else
                prediction_ts = &z_xmap_x;
            break;
    }
    
    x_scale = 1.0;
    y_scale = xmap_attractor_scale;
    z_scale = xmap_attractor_scale;
    glPushMatrix();
    glTranslated(ts_delta_x*x_scale, (delta_y-d)*y_scale-0.1, 0);
    draw_xmap_ts(frame, 0, 0, pred_dim, 1.0, 1.0, prediction_ts);
    glPopMatrix();
    
	return;
}

void attractor::draw_xmap(const int frame)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    if(0.75*window_width > window_height)
    {
        glOrtho(-1.5*window_width/window_height, 1.5*window_width/window_height, -1.5, 1.5, -1.5, 1.5);
    }
    else
    {
        glOrtho(-2.0, 2.0, -2.0*window_height/window_width, 2.0*window_height/window_width, -1.5, 1.5);
    }
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMultMatrixd(rot_matrix);
    
    draw_xmap_generic(frame, 0, lag_dim, pred_dim, false);
    return;
}

void attractor::draw_univariate_ts(const int frame)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    if(0.75*window_width > window_height)
    {
        glOrtho(-1.5*window_width/window_height, 1.5*window_width/window_height, -1.5, 1.5, -1.5, 1.5);
    }
    else
    {
        glOrtho(-2.0, 2.0, -2.0*window_height/window_width, 2.0*window_height/window_width, -1.5, 1.5);
    }
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMultMatrixd(rot_matrix);
    
    vector<double>* ts;
    double ts_delta_x = 0.8, delta_x = 1.8;
    double project_dist;
    vector<double> point;
    vector<int> nn_indices;
    vector<double>::iterator pred_x_iter, pred_y_iter, pred_z_iter;
    double r = 0.0, g = 0.0, b = 0.0;
    int texture_index;
    
    project_dist = delta_x + ts_delta_x / xmap_attractor_scale;
    
    // ***** get nearest neighbors in reconstructed space *****
    switch(lag_dim)
	{
		case 1:
			ts = &x;
            nn_indices = x_nn_indices[frame];
            pred_y_iter = x_forecast.begin();
            pred_x_iter = x_forecast_lag_1.begin();
            pred_z_iter = x_forecast_lag_2.begin();
            r = 1.0;
            texture_index = X_LABEL_TEXTURE;
            break;
		case 2:
			ts = &y;
            nn_indices = y_nn_indices[frame];
            pred_y_iter = y_forecast.begin();
            pred_x_iter = y_forecast_lag_1.begin();
            pred_z_iter = y_forecast_lag_2.begin();
            g = 1.0;
            texture_index = Y_LABEL_TEXTURE;
			break;
		case 3:
			ts = &z;
            nn_indices = z_nn_indices[frame];
            pred_y_iter = z_forecast.begin();
            pred_x_iter = z_forecast_lag_1.begin();
            pred_z_iter = z_forecast_lag_2.begin();
            b = 1.0;
            texture_index = Z_LABEL_TEXTURE;
			break;
	}
    
    if(nn_indices.size() == 0)
        return;
    
    // ***** ATTRACTOR *****
    x_scale = xmap_attractor_scale;
    y_scale = xmap_attractor_scale;
    z_scale = xmap_attractor_scale;
    glPushMatrix();
    
    glTranslated(-delta_x*x_scale, 0, 0);
    glRotated(theta, 0, 1, 0);
    glTranslated(-d*x_scale, -d*y_scale, -d*z_scale);
    
    // axes
    draw_axis(2, r, g, b, 0, texture_index+1, 0);
    draw_axis(1, r*LAG_SAT, g*LAG_SAT, b*LAG_SAT, 0.5, texture_index+2, 1);
    draw_axis(3, r*LAG2_SAT, g*LAG2_SAT, b*LAG2_SAT, 1, texture_index+3, 2);
    
    glPushMatrix();
    
    glScaled(x_scale, y_scale, z_scale);
    
    if(MANIFOLD_LABEL)
        enqueue_label(d, 2*d, 2*d, M_LABEL_TEXTURE+lag_dim, texture_scale, 1, 0);
    
    // draw attractor
    glLineWidth(scale * LINE_WIDTH);
    draw_embedding(ts->begin()+tau, ts->begin()+2*tau, ts->begin(), frame-2*tau);
    
   	// draw current point
    glPointSize(POINT_WIDTH*scale);
	glColor4dv(point_color);
	glBegin(GL_POINTS);
	glVertex3d(*(ts->begin() + frame-1-tau), *(ts->begin() + frame-1), *(ts->begin() + frame-1-2*tau));
	glEnd();
	
    if(nn_indices.size() == 0)
        return;
    
    // neighbors
    glColor4dv(neighbor_color);
    for(vector<int>::iterator nn_index = nn_indices.begin(); nn_index != nn_indices.end(); nn_index++)
    {
        glBegin(GL_POINTS);
        glVertex3d(*(ts->begin() + *nn_index-1-tau), *(ts->begin() + *nn_index-1), *(ts->begin() + *nn_index-1-2*tau));
        glEnd();
    }
    
    glPopMatrix();
    
    if(DEBUG)
    {
        glPushMatrix();
        glScaled(x_scale, y_scale, z_scale);
        glColor4dv(neighbor_color);
        
        // neighbor trajectories
        glLineWidth(scale * LINE_WIDTH * 1.5);
        for(vector<int>::iterator nn_index = nn_indices.begin(); nn_index != nn_indices.end(); nn_index++)
        {
            glBegin(GL_LINE_STRIP);
            for(int k = 0; k <= tp; k++)
            {
                glVertex3d(*(ts->begin() + *nn_index-1+k-tau), *(ts->begin() + *nn_index-1+k), *(ts->begin() + *nn_index-1-2*tau+k));
            }
            glEnd();
        }
        
        // neighbor forward points
        glColor4d(lag_dim == 1, lag_dim == 2, lag_dim == 3, 0.5);
        for(vector<int>::iterator nn_index = nn_indices.begin(); nn_index != nn_indices.end(); nn_index++)
        {
            glBegin(GL_POINTS);
            glVertex3d(*(ts->begin() + *nn_index-1-tau+tp), *(ts->begin() + *nn_index-1+tp), *(ts->begin() + *nn_index-1-2*tau+tp));
            glEnd();
        }
        
        // draw forecast
        glColor4dv(pred_color);
        glBegin(GL_POINTS);
        glVertex3d(*(pred_x_iter+frame+tp-1), *(pred_y_iter+frame+tp-1), *(pred_z_iter+frame+tp-1));
        glEnd();
        
        glBegin(GL_LINE_STRIP);
        for(int k = 0; k <= tp; k++)
        {
            glVertex3d(*(pred_x_iter+frame+k-1), *(pred_y_iter+frame+k-1), *(pred_z_iter+frame+k-1));
        }
        glEnd();
        
        glPopMatrix();
    }
    
    glPopMatrix();
    
	// ***** TIME SERIES *****
    x_scale = 1.0;
    y_scale = xmap_attractor_scale;
    z_scale = xmap_attractor_scale;
    glPushMatrix();
    glTranslated(ts_delta_x, (-d)*y_scale, 0);
    draw_ts(frame, 0, 2*tau, lag_dim, 1.0, 0);
    if(DEBUG)
    {
        switch(lag_dim)
        {
            case 1:
                draw_univariate_ts(frame, &x_forecast);
                break;
            case 2:
                draw_univariate_ts(frame, &y_forecast);
                break;
            case 3:
                draw_univariate_ts(frame, &z_forecast);
                break;
        }
    }
    //draw_xmap_ts(frame, 0, 0, lag_dim, 1.0, 1.0, &x_forecast);
    
    glPopMatrix();
    
    
    // ***** PROJECTION LINES *****
    if(DEBUG)
    {
        x_scale = xmap_attractor_scale;
        y_scale = xmap_attractor_scale;
        z_scale = xmap_attractor_scale;
        glPushMatrix();
        glScaled(x_scale, y_scale, z_scale);
        glTranslated(-delta_x, 0, 0);
        glRotated(theta, 0, 1, 0);
        glTranslated(-d, -d, -d);
        
        glColor4dv(pred_color);
        glLineWidth(scale * LINE_WIDTH);
        glBegin(GL_LINES);
        glVertex3d(*(pred_x_iter+frame+tp-1), *(pred_y_iter+frame+tp-1), *(pred_z_iter+frame+tp-1));
        glVertex3d(project_dist*cos(theta/180*PI)+d+tp*2.0 / num_points / draw_fraction / x_scale, *(pred_y_iter+frame+tp-1), project_dist*sin(theta/180*PI)+d);
        glEnd();
        
        glPopMatrix();
    }
    
    return;
}


void attractor::draw_univariate(const int frame)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(30.0f,(GLfloat)window_width/(GLfloat)window_height,0.2f,10.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMultMatrixd(rot_matrix);
    glTranslated(0, init_distance, 0);    
    
    vector<double>* ts;
    vector<int> nn_indices;
    vector<double>::iterator pred_x_iter, pred_y_iter, pred_z_iter;
    
    x_scale = 1.0;
    y_scale = 1.0;
    z_scale = 1.0;
	glTranslated(-d, -d, -d);
	
	switch(lag_dim)
	{
		case 1:
			ts = &x;
            nn_indices = x_nn_indices[frame];
            pred_x_iter = x_forecast.begin();
            pred_y_iter = x_forecast_lag_1.begin();
            pred_z_iter = x_forecast_lag_2.begin();
			break;
		case 2:
			ts = &y;
            nn_indices = y_nn_indices[frame];
            pred_x_iter = y_forecast.begin();
            pred_y_iter = y_forecast_lag_1.begin();
            pred_z_iter = y_forecast_lag_2.begin();
			break;
		case 3:
			ts = &z;
            nn_indices = z_nn_indices[frame];
            pred_x_iter = z_forecast.begin();
            pred_y_iter = z_forecast_lag_1.begin();
            pred_z_iter = z_forecast_lag_2.begin();
			break;
	}
	
	// draw attractor
    glLineWidth(scale * LINE_WIDTH);
    draw_embedding(ts->begin()+2*tau, ts->begin()+tau, ts->begin(), frame-2*tau);
    
    if(MANIFOLD_LABEL)
        enqueue_label(d, d, 2*d, M_LABEL_TEXTURE+lag_dim, texture_scale, 1, 0);
    
	// draw current point
    glPointSize(POINT_WIDTH*scale);
	glColor4dv(point_color);
	glBegin(GL_POINTS);
	glVertex3d(*(ts->begin() + frame-1), *(ts->begin() + frame-1-tau), *(ts->begin() + frame-1-2*tau));
	glEnd();
	
	draw_axes(true, lag_dim);
    
    if(nn_indices.size() == 0)
        return;
    
    glColor4dv(neighbor_color);
    for(vector<int>::iterator nn_index = nn_indices.begin(); nn_index != nn_indices.end(); nn_index++)
    {
        glBegin(GL_POINTS);
        glVertex3d(*(ts->begin() + *nn_index-1), *(ts->begin() + *nn_index-1-tau), *(ts->begin() + *nn_index-1-2*tau));
        glEnd();
    }
    
    if(DEBUG)
    {
        // draw trajectories
        glLineWidth(scale * LINE_WIDTH * 1.5);
        for(vector<int>::iterator nn_index = nn_indices.begin(); nn_index != nn_indices.end(); nn_index++)
        {
            glBegin(GL_LINE_STRIP);
            for(int k = 0; k <= tp; k++)
            {
                glVertex3d(*(ts->begin() + *nn_index-1+k), *(ts->begin() + *nn_index-1-tau+k), *(ts->begin() + *nn_index-1-2*tau+k));
            }
            glEnd();
        }
        
        // draw projected neighbors
        glColor4d(lag_dim == 1, lag_dim == 2, lag_dim == 3, 0.5);
        for(vector<int>::iterator nn_index = nn_indices.begin(); nn_index != nn_indices.end(); nn_index++)
        {
            glBegin(GL_POINTS);
            glVertex3d(*(ts->begin() + *nn_index-1+tp), *(ts->begin() + *nn_index-1-tau+tp), *(ts->begin() + *nn_index-1-2*tau+tp));
            glEnd();
        }
        
        // draw forecast
        glColor4dv(pred_color);
        glBegin(GL_POINTS);
        glVertex3d(*(pred_x_iter+frame+tp-1), *(pred_y_iter+frame+tp-1), *(pred_z_iter+frame+tp-1));
        glEnd();
        
        glBegin(GL_LINE_STRIP);
        for(int k = 0; k <= tp; k++)
        {
            glVertex3d(*(pred_x_iter+frame+k-1), *(pred_y_iter+frame+k-1), *(pred_z_iter+frame+k-1));
        }
        glEnd();
    }
    
	return;
}

void attractor::draw_shadow(const int frame)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(30.0f,(GLfloat)window_width/(GLfloat)window_height,0.2f,10.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMultMatrixd(rot_matrix);
    glTranslated(0, init_distance, 0);
    
	vector<double>* ts;
    double sep = 1.4;
    vector<int>nn_indices;
    int nn_index;
    
	glTranslated(-d, -d, -d+.8);
	switch(lag_dim)
	{
		case 1:
			ts = &x;
			break;
		case 2:
			ts = &y;
			break;
		case 3:
			ts = &z;
			break;
	}
	// draw attractor
    glLineWidth(scale * LINE_WIDTH);
    draw_embedding(x.begin(), y.begin(), z.begin(), frame);
    if(MANIFOLD_LABEL)
        enqueue_label(0, 0, d, M_LABEL_TEXTURE, texture_scale, 1, 0);
    
	// draw current point
    glPointSize(POINT_WIDTH*scale);
	glColor4dv(point_color);
	glBegin(GL_POINTS);
	glVertex3d(x[frame-1], y[frame-1], z[frame-1]);
	glEnd();
    if(TSVIEW)
    {
        draw_axes(false, 0);
    }
    
    if(frame-1-2*tau < 0)
        return;
    
    // draw line connecting attractors
    draw_curve(x[frame-1], y[frame-1], z[frame-1], *(ts->begin() + frame-1), *(ts->begin() + frame-1-tau), *(ts->begin() + frame-1-2*tau)-sep);
    
    glPushMatrix();
	glTranslated(0, 0, -sep);
	
	// draw shadow attractor
    glLineWidth(scale * LINE_WIDTH);
    draw_embedding(ts->begin()+2*tau, ts->begin()+tau, ts->begin(), frame-2*tau);
    
    if(TSVIEW)
    {
        draw_axes(true, lag_dim);
    }
    if(MANIFOLD_LABEL)
        enqueue_label(d*0.35, 0, d, M_LABEL_TEXTURE+lag_dim, texture_scale, 1, 0);
	
	// draw current point
    glPointSize(POINT_WIDTH*scale);
	glColor4dv(point_color);
	glBegin(GL_POINTS);
	glVertex3d(*(ts->begin() + frame-1), *(ts->begin() + frame-1-tau), *(ts->begin() + frame-1-2*tau));
	glEnd();
	
    // find index of nearest neighbor
    switch(lag_dim)
	{
		case 1:
            nn_indices = x_nn_indices[frame];
			break;
		case 2:
            nn_indices = y_nn_indices[frame];
			break;
		case 3:
            nn_indices = z_nn_indices[frame];
			break;
	}
    if(nn_indices.size() == 0)
    {
        glPopMatrix();
        return;
    }
    
    nn_index = nn_indices[0];    
    if(DEBUG)
    {
        glColor4dv(neighbor_color);
        glBegin(GL_POINTS);
        glVertex3d(x[nn_index-1], y[nn_index-1], z[nn_index-1]+sep);
        glEnd();
        
        draw_curve(x[nn_index-1], y[nn_index-1], z[nn_index-1]+sep, *(ts->begin() + nn_index-1), *(ts->begin() + nn_index-1-tau), *(ts->begin() + nn_index-1-2*tau));
        
        glColor4dv(neighbor_color);
        glBegin(GL_POINTS);
        glVertex3d(*(ts->begin() + nn_index-1), *(ts->begin() + nn_index-1-tau), *(ts->begin() + nn_index-1-2*tau));
        glEnd();
    }
    
	glPopMatrix();
    
    
    
	return;
}

void attractor::draw_reconstruction(const int frame)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(30.0f,(GLfloat)window_width/(GLfloat)window_height,0.2f,10.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMultMatrixd(rot_matrix);
    glTranslated(0, init_distance, 0);    
    
	vector<double>* ts;
    x_scale = 1.0;
    y_scale = 1.0;
    z_scale = 1.0;
	glTranslated(-d, -d, -d);
	
	switch(lag_dim)
	{
		case 1:
			ts = &x;
			break;
		case 2:
			ts = &y;
			break;
		case 3:
			ts = &z;
			break;
	}
	
	// draw attractor
    glLineWidth(scale * LINE_WIDTH);
    draw_embedding(ts->begin()+2*tau, ts->begin()+tau, ts->begin(), frame-2*tau);
    
    if(MANIFOLD_LABEL)
        enqueue_label(d, d, 2*d, M_LABEL_TEXTURE+lag_dim, texture_scale, 1, 0);
    
	// draw current point
    glPointSize(POINT_WIDTH*scale);
	glColor4dv(point_color);
	glBegin(GL_POINTS);
	glVertex3d(*(ts->begin() + frame-1), *(ts->begin() + frame-1-tau), *(ts->begin() + frame-1-2*tau));
	glEnd();
	
	draw_axes(true, lag_dim);
	return;
}

void attractor::draw_generic_reconstruction(const int frame)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(30.0f,(GLfloat)window_width/(GLfloat)window_height,0.2f,10.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMultMatrixd(rot_matrix);
    glTranslated(0, init_distance, 0);
    
	vector<double>* ts_x;
    vector<double>* ts_y;
    vector<double>* ts_z;
    x_scale = 1.0;
    y_scale = 1.0;
    z_scale = 1.0;
	glTranslated(-d, -d, -d);
	
	switch(x_dim)
	{
		case 1:
			ts_x = &x;
			break;
		case 2:
			ts_x = &y;
			break;
		case 3:
			ts_x = &z;
			break;
	}
    switch(y_dim)
	{
		case 1:
			ts_y = &x;
			break;
		case 2:
			ts_y = &y;
			break;
		case 3:
			ts_y = &z;
			break;
	}
    switch(z_dim)
	{
		case 1:
			ts_z = &x;
			break;
		case 2:
			ts_z = &y;
			break;
		case 3:
			ts_z = &z;
			break;
	}
	
	// draw attractor
    glLineWidth(scale * LINE_WIDTH);
    draw_embedding(ts_x->begin()+x_lag*tau, ts_y->begin()+y_lag*tau, ts_z->begin()+z_lag*tau, frame-2*tau);
    
	// draw current point
    glPointSize(POINT_WIDTH*scale);
	glColor4dv(point_color);
	glBegin(GL_POINTS);
	glVertex3d(*(ts_x->begin() + frame-1-2*tau+x_lag*tau), *(ts_y->begin() + frame-1-2*tau+y_lag*tau), *(ts_z->begin() + frame-1-2*tau+z_lag*tau));
	glEnd();
	
    draw_lag_axis(1, x_dim, x_lag);
    draw_lag_axis(2, y_dim, y_lag);
    draw_lag_axis(3, z_dim, z_lag);
    
	//draw_axes(true, lag_dim);
	return;
}

void attractor::draw_lagged_time_series(const int frame)
{
    vector<double>* ts;
    if(SPLIT_VIEW)
    {
        // draw small attractor
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(30.0f,(GLfloat)window_width/(GLfloat)window_height,0.2f,10.0f);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glMultMatrixd(rot_matrix);
        glTranslated(1, init_distance, 0);
        
        x_scale = 0.7;
        y_scale = 0.7;
        z_scale = 0.7;
        glTranslated(-d*x_scale, -d*y_scale, -d*z_scale);
        
        switch(lag_dim)
        {
            case 1:
                ts = &x;
                break;
            case 2:
                ts = &y;
                break;
            case 3:
                ts = &z;
                break;
        }
        
        // draw attractor
        glLineWidth(scale * LINE_WIDTH);
        glPushMatrix();
        glScaled(x_scale, y_scale, z_scale);
        draw_embedding(ts->begin()+2*tau, ts->begin()+tau, ts->begin(), frame-2*tau);
        glPopMatrix();
        
        // draw current point
        glPointSize(POINT_WIDTH*scale);
        glColor4dv(point_color);
        glBegin(GL_POINTS);
        glVertex3d(*(ts->begin() + frame-1)*x_scale, *(ts->begin() + frame-1-tau)*y_scale, *(ts->begin() + frame-1-2*tau)*z_scale);
        glEnd();
        
        texture_scale *= 1.25;
        if(MANIFOLD_LABEL)
            enqueue_label(d*x_scale, d*y_scale, 2*d*z_scale, M_LABEL_TEXTURE+lag_dim, texture_scale, 1, 0);
        draw_axes(true, lag_dim);
        texture_scale *= 0.8;
        
        // draw labels
        VIEW = RECONSTRUCTION;
        draw_labels();
        VIEW = LAGS;
        texture_queue.clear();
     
        // draw half time series
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(-1.2, 1.2, -1.2, 1.2, -1.5, 1.5);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        
        x_scale = 1.0;
        y_scale = 0.25;
        z_scale = 1.0;
        
        glPushMatrix();
        glTranslated(0.0, 0.4, 0.0);
        draw_half_ts(frame, 0, 16*tau, lag_dim, 1, 0);
        glPopMatrix();
        
        glPushMatrix();
        glTranslated(0.0, -0.3, 0.0);
        draw_half_ts(frame, 8*tau, 16*tau, lag_dim, LAG_SAT, 0.5);
        glPopMatrix();
        
        glPushMatrix();
        glTranslated(0.0, -1.0, 0.0);
        draw_half_ts(frame, 16*tau, 16*tau, lag_dim, LAG2_SAT, 1);
        glPopMatrix();
    }
    else
    {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(-1.2, 1.2, -1.2, 1.2, -1.5, 1.5);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        
        x_scale = 1.0;
        y_scale = 0.25;
        z_scale = 1.0;
        
        glPushMatrix();
        glTranslated(0.0, 0.4, 0.0);
        draw_ts(frame, 0, 16*tau, lag_dim, 1, 0);
        glPopMatrix();
        
        glPushMatrix();
        glTranslated(0.0, -0.3, 0.0);
        draw_ts(frame, 8*tau, 16*tau, lag_dim, LAG_SAT, 0.5);
        glPopMatrix();
        
        glPushMatrix();
        glTranslated(0.0, -1.0, 0.0);
        draw_ts(frame, 16*tau, 16*tau, lag_dim, LAG2_SAT, 1);
        glPopMatrix();
    }
    
	return;
}

void attractor::draw_time_series(const int frame)
{
    if(SPLIT_VIEW)
    {
        // draw small attractor
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(30.0f,(GLfloat)window_width/(GLfloat)window_height,0.2f,10.0f);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glMultMatrixd(rot_matrix);
        glTranslated(1, init_distance, 0);    
        
        x_scale = 0.7;
        y_scale = 0.7;
        z_scale = 0.7;
        glTranslated(-d*x_scale, -d*y_scale, -d*z_scale);
        
        // draw attractor
        glLineWidth(scale * LINE_WIDTH);
        glPushMatrix();
        glScaled(x_scale, y_scale, z_scale);
        draw_embedding(x.begin(), y.begin(), z.begin(), frame);
        glPopMatrix();
        
        // draw current point
        glPointSize(POINT_WIDTH*scale);
        glColor4dv(point_color);
        glBegin(GL_POINTS);
        glVertex3d(x[frame-1]*x_scale, y[frame-1]*y_scale, z[frame-1]*z_scale);
        glEnd();
        
        texture_scale *= 1.25;
        if(MANIFOLD_LABEL)
            enqueue_label(d*x_scale, d*y_scale, 2*d*z_scale, M_LABEL_TEXTURE, texture_scale, 1, 0);
        draw_axes(false, 0);
        texture_scale *= 0.8;
        
        // draw labels
        VIEW = MANIFOLD;
        draw_labels();
        VIEW = TIME_SERIES;
        texture_queue.clear();
        
        // draw half time series
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(-1.2, 1.2, -1.2, 1.2, -1.5, 1.5);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        
        x_scale = 1.0;
        y_scale = 0.25;
        z_scale = 1.0;
        
        glPushMatrix();
        glTranslated(0.0, 0.4, 0.0);
        draw_half_ts(frame, 0, 0, 1, 1, 0);
        glPopMatrix();
        
        glPushMatrix();
        glTranslated(0.0, -0.3, 0.0);
        draw_half_ts(frame, 0, 0, 2, 1, 0);
        glPopMatrix();
        
        glPushMatrix();
        glTranslated(0.0, -1.0, 0.0);
        draw_half_ts(frame, 0, 0, 3, 1, 0);
        glPopMatrix();
        
    }
    else
    {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(-1.2, 1.2, -1.2, 1.2, -1.5, 1.5);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        
        x_scale = 1.0;
        y_scale = 0.25;
        z_scale = 1.0;
        
        glPushMatrix();
        glTranslated(0.0, 0.4, 0.0);
        draw_ts(frame, 0, 0, 1, 1, 0);
        glPopMatrix();
        
        glPushMatrix();
        glTranslated(0.0, -0.3, 0.0);
        draw_ts(frame, 0, 0, 2, 1, 0);
        glPopMatrix();
        
        glPushMatrix();
        glTranslated(0.0, -1.0, 0.0);
        draw_ts(frame, 0, 0, 3, 1, 0);
        glPopMatrix();
    }
    
	return;
}

void attractor::draw_manifold(const int frame)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(30.0f,(GLfloat)window_width/(GLfloat)window_height,0.2f,10.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMultMatrixd(rot_matrix);
    glTranslated(0, init_distance, 0);    
    glRotated(phi_x[frame-1], 1, 0, 0);
    glRotated(phi_y[frame-1], 0, 1, 0);
    glRotated(phi_z[frame-1], 0, 0, 1);
    
    x_scale = 1.0;
    y_scale = 1.0;
    z_scale = 1.0;
	glTranslated(-d, -d, -d);
    
	// draw attractor
    glLineWidth(scale * LINE_WIDTH);
    draw_embedding(x.begin(), y.begin(), z.begin(), frame);
    
    if(MANIFOLD_LABEL)
        enqueue_label(d, d, 2*d, M_LABEL_TEXTURE, texture_scale, 1, 0);
	
	// draw current point
    glPointSize(POINT_WIDTH*scale);
	glColor4dv(point_color);
	glBegin(GL_POINTS);
	glVertex3d(x[frame-1], y[frame-1], z[frame-1]);
	glEnd();
	
	draw_tracers(frame);
	draw_axes(false, 0);
	return;
}

vector<double> attractor::draw_xmap_manifold(const int frame, const vector<int> & nn_indices, 
                                             const vector<double> & nn_weights, const int var, bool swap_xy)
{
    int label_lag, label_other, label_pred;
    int my_frame, delta_frame = 0;
    double total_weight;
    double pred_x = 0.0, pred_y = 0.0, pred_z = 0.0;
    vector<double> point;
    point.resize(3, 0.0);
    vector<double>::iterator my_x, my_y, my_z, temp_iter;
    
    // setup manifold axes
    switch(var)
    {
        case 0: // original manifold
            switch(lag_dim)
        {
            case 1:
                my_x = x.begin();
                label_lag = X_LABEL_TEXTURE;
                if (pred_dim == 2)
                {
                    my_z = z.begin();
                    label_other = Z_LABEL_TEXTURE;
                    my_y = y.begin();
                    label_pred = Y_LABEL_TEXTURE;
                }
                else
                {
                    my_z = y.begin();
                    label_other = Y_LABEL_TEXTURE;
                    my_y = z.begin();
                    label_pred = Z_LABEL_TEXTURE;
                }
                break;
            case 2:
                my_x = y.begin();
                label_lag = Y_LABEL_TEXTURE;
                if (pred_dim == 1)
                {
                    my_z = z.begin();
                    label_other = Z_LABEL_TEXTURE;
                    my_y = x.begin();
                    label_pred = X_LABEL_TEXTURE;
                }
                else
                {
                    my_z = x.begin();
                    label_other = X_LABEL_TEXTURE;
                    my_y = z.begin();
                    label_pred = Z_LABEL_TEXTURE;
                }
                break;
            case 3:
                my_x = z.begin();
                label_lag = Z_LABEL_TEXTURE;
                if (pred_dim == 2)
                {
                    my_z = x.begin();
                    label_other = X_LABEL_TEXTURE;
                    my_y = y.begin();
                    label_pred = Y_LABEL_TEXTURE;
                }
                else
                {
                    my_z = y.begin();
                    label_other = Y_LABEL_TEXTURE;
                    my_y = x.begin();
                    label_pred = X_LABEL_TEXTURE;
                }
                break;
        }
            break;
        case 1: // x reconstruction
            my_x = x.begin()+2*tau;
            my_y = x.begin()+tau;
            my_z = x.begin();
            break;
        case 2: // y reconstruction
            my_x = y.begin()+2*tau;
            my_y = y.begin()+tau;
            my_z = y.begin();
            break;
        case 3: // z reconstruction
            my_x = z.begin()+2*tau;
            my_y = z.begin()+tau;
            my_z = z.begin();
            break;
        default:
            cerr << "ERROR (attractor): out-of-bounds NW_manifold input: " << var << ".\n";
            exit(1);
    }
    
    // swap xy if needed
    if(swap_xy)
    {
        temp_iter = my_x;
        my_x = my_y;
        my_y = temp_iter;
    }
    
    // adjust frame if lags
    if(var != 0)
        delta_frame = 2*tau;
    my_frame = frame - delta_frame;
    
    if(my_frame < 1)
        return point;
    
    // find current point
    point[0] = *(my_x + my_frame-1);
    point[1] = *(my_y + my_frame-1);
    point[2] = *(my_z + my_frame-1);
    
    glPushMatrix();
    glScaled(x_scale, y_scale, z_scale);
    
    // draw attractor
    glLineWidth(scale * SMALL_LINE_WIDTH);
    draw_embedding(my_x, my_y, my_z, my_frame);
    
    // color nearest neighbors
    glPointSize(7.0*scale);
    glColor4dv(neighbor_color);
    glBegin(GL_POINTS);
    total_weight = 0;
    for(int i = 0; i < nn_indices.size(); i++)
    {
        glVertex3d(*(my_x + nn_indices[i]-delta_frame), *(my_y + nn_indices[i]-delta_frame), *(my_z + nn_indices[i]-delta_frame));
        pred_x += *(my_x + nn_indices[i]-delta_frame) * nn_weights[i];
        pred_y += *(my_y + nn_indices[i]-delta_frame) * nn_weights[i];
        pred_z += *(my_z + nn_indices[i]-delta_frame) * nn_weights[i];
        total_weight += nn_weights[i];
    }
    glEnd();
    
    // draw current point
    glPointSize(POINT_WIDTH*scale);
	glColor4dv(point_color);
	glBegin(GL_POINTS);
	glVertex3d(point[0], point[1], point[2]);
	glEnd();
    
    // draw prediction
    glPointSize(POINT_WIDTH*scale);
	glColor4d(lag_dim == 1, lag_dim == 2, lag_dim == 3, 0.5);
	glBegin(GL_POINTS);
    glVertex4d(pred_x, pred_y, pred_z, total_weight);
	glEnd();
    
    glPopMatrix();
    
    // draw axes
    switch(var)
    {
        case 0: // original manifold
            draw_axis(1+swap_xy, lag_dim == 1, lag_dim == 2, lag_dim == 3, 0, label_lag, -1);
            draw_axis(2-swap_xy, pred_dim == 1, pred_dim == 2, pred_dim == 3, 0, label_pred, -1);
            draw_axis(3, lag_dim+pred_dim == 5, lag_dim+pred_dim == 4, lag_dim+pred_dim == 3, 0, label_other, -1);
            break;
        case 1: // x reconstruction
            draw_axis(1+swap_xy, 1, 0, 0, 0, X_LABEL_TEXTURE+1, 0);
            draw_axis(2-swap_xy, LAG_SAT, 0, 0, 0.5, X_LABEL_TEXTURE+2, 1);
            draw_axis(3, LAG2_SAT, 0, 0, 1, X_LABEL_TEXTURE+3, 2);
            break;
        case 2: // y reconstruction
            draw_axis(1+swap_xy, 0, 1, 0, 0, Y_LABEL_TEXTURE+1, 0);
            draw_axis(2-swap_xy, 0, LAG_SAT, 0, 0.5, Y_LABEL_TEXTURE+2, 1);
            draw_axis(3, 0, LAG2_SAT, 0, 1, Y_LABEL_TEXTURE+3, 2);
            break;
        case 3: // z reconstruction
            draw_axis(1+swap_xy, 0, 0, 1, 0, Z_LABEL_TEXTURE+1, 0);
            draw_axis(2-swap_xy, 0, 0, LAG_SAT, 0.5, Z_LABEL_TEXTURE+2, 1);
            draw_axis(3, 0, 0, LAG2_SAT, 1, Z_LABEL_TEXTURE+3, 2);
            break;
        default:
            cerr << "ERROR (attractor): out-of-bounds NW_manifold input: " << var << ".\n";
            exit(1);
    }
    return point;
}

void attractor::draw_xmap_generic(const int frame, const int NW_manifold, const int SW_manifold, const int NE_manifold, bool ts_trace)
{
    double ts_delta_x = 0.8, delta_x, delta_y = 1.1;
    double project_dist;
    vector<double> NW_point, SW_point, NE_point;
    vector<int> nn_indices;
    vector<double> nn_weights;
    
    // setup positions
    if (ts_trace)
        delta_x = 1.8;
    else
        delta_x = 1.4;
    project_dist = delta_x + ts_delta_x / xmap_attractor_scale;
    
    // ***** get nearest neighbors in reconstructed space *****
    switch(lag_dim)
	{
		case 1:
            nn_indices = x_nn_indices[frame];
            nn_weights = x_nn_weights[frame];
			break;
		case 2:
            nn_indices = y_nn_indices[frame];
            nn_weights = y_nn_weights[frame];
			break;
		case 3:
            nn_indices = z_nn_indices[frame];
            nn_weights = z_nn_weights[frame];
			break;
	}
    
    // ***** NW ATTRACTOR *****
    if(NW_manifold >= 0)
    {
        x_scale = xmap_attractor_scale;
        y_scale = xmap_attractor_scale;
        z_scale = xmap_attractor_scale;
        glPushMatrix();
        glTranslated(-delta_x*x_scale, delta_y*y_scale-0.1, 0);
        glRotated(theta, 0, 1, 0);
        glTranslated(-d*x_scale, -d*y_scale, -d*z_scale);
        NW_point = draw_xmap_manifold(frame, nn_indices, nn_weights, NW_manifold, NW_manifold != 0);
        if(MANIFOLD_LABEL)
            enqueue_label(d*x_scale, 2*d*y_scale, 2*d*z_scale, M_LABEL_TEXTURE+NW_manifold, texture_scale, 1, 0);
        glPopMatrix();
    }
    
    // ***** NE TIME SERIES *****
    if(ts_trace)
    {
        x_scale = 1.0;
        y_scale = xmap_attractor_scale;
        z_scale = xmap_attractor_scale;
        glPushMatrix();
        glTranslated(ts_delta_x, (delta_y-d)*y_scale-0.1, 0);
        draw_ts(frame, 0, 0, pred_dim, 1.0, 0);
        glPopMatrix();
    }
    
    // ***** NE ATTRACTOR *****
    if(NE_manifold >= 0)
    {
        x_scale = xmap_attractor_scale;
        y_scale = xmap_attractor_scale;
        z_scale = xmap_attractor_scale;
        glPushMatrix();
        glTranslated(delta_x*x_scale, delta_y*y_scale-0.1, 0);
        glRotated(theta, 0, 1, 0);
        glTranslated(-d*x_scale, -d*y_scale, -d*z_scale);
        NE_point = draw_xmap_manifold(frame, nn_indices, nn_weights, NE_manifold, NE_manifold != 0);
        if(MANIFOLD_LABEL)
            enqueue_label(d*x_scale, 2*d*y_scale, 2*d*z_scale, M_LABEL_TEXTURE+NE_manifold, texture_scale, 1, 0);
        glPopMatrix();
    }
    
    // ***** SW ATTRACTOR *****
    if(SW_manifold >= 0)
    {
        x_scale = xmap_attractor_scale;
        y_scale = xmap_attractor_scale;
        z_scale = xmap_attractor_scale;
        glPushMatrix();
        glTranslated(-delta_x*x_scale, -delta_y*y_scale-0.1, 0);
        glRotated(theta, 0, 1, 0);
        glTranslated(-d*x_scale, -d*y_scale, -d*z_scale);
        SW_point = draw_xmap_manifold(frame, nn_indices, nn_weights, SW_manifold, SW_manifold == 0);
        if(MANIFOLD_LABEL)
            enqueue_label(d*x_scale, 2*d*y_scale, 2*d*z_scale, M_LABEL_TEXTURE+SW_manifold, texture_scale, 1, 0);
        glPopMatrix();
    }
    
	// ***** SE TIME SERIES *****
    if(ts_trace)
    {
        x_scale = 1.0;
        y_scale = xmap_attractor_scale;
        z_scale = xmap_attractor_scale;
        glPushMatrix();
        glTranslated(ts_delta_x, (-delta_y-d)*y_scale-0.1, 0);
        draw_ts(frame, tau, 2*tau, lag_dim, LAG_SAT, 0);
        glPopMatrix();
    }
	
    // ***** PROJECTION LINES *****
    if(DEBUG)
    {
        x_scale = xmap_attractor_scale;
        y_scale = xmap_attractor_scale;
        z_scale = xmap_attractor_scale;
        glPushMatrix();
        glScaled(x_scale, y_scale, z_scale);
        glTranslated(-delta_x, delta_y-0.1/y_scale, 0);
        glRotated(theta, 0, 1, 0);
        glTranslated(-d, -d, -d);
        
        // NW attractor to ts
        if(ts_trace)
        {
            glColor3d(pred_dim == 1, pred_dim == 2, pred_dim == 3);
            glLineWidth(scale * LINE_WIDTH);
            glBegin(GL_LINES);
            glVertex3d(NW_point[0], NW_point[1], NW_point[2]);
            glVertex3d(project_dist*cos(theta/180*PI)+d, NW_point[1], project_dist*sin(theta/180*PI)+d);
            glEnd();
        }
        
        // NW attractor to NE attractor & SW attractor to NE attractor
        if(NE_manifold >= 0)
        {
            draw_curve(NW_point[0], NW_point[1], NW_point[2], 
                       NE_point[0] + 2*delta_x*cos(theta/180*PI), NE_point[1], NE_point[2] + 2*delta_x*sin(theta/180*PI));
            draw_curve(SW_point[0], SW_point[1]-2*delta_y, SW_point[2], 
                       NE_point[0] + 2*delta_x*cos(theta/180*PI), NE_point[1], NE_point[2] + 2*delta_x*sin(theta/180*PI));
        }
        
        // NW attractor to SW attractor
        draw_curve(NW_point[0], NW_point[1], NW_point[2], SW_point[0], SW_point[1]-2*delta_y, SW_point[2]);
        
        // SW attractor to ts
        if(ts_trace)
        {
            glTranslated(0, -2*delta_y, 0);
            glColor3d(SW_manifold == 1, SW_manifold == 2, SW_manifold == 3);
            glLineWidth(scale * LINE_WIDTH);
            glBegin(GL_LINES);
            glVertex3d(SW_point[0], SW_point[1], SW_point[2]);
            glVertex3d(project_dist*cos(theta/180*PI)+d, SW_point[1], project_dist*sin(theta/180*PI)+d);
            glEnd();
        }
        
        glPopMatrix();
    }
    
    return;
}

void attractor::draw_univariate_ts(int frame, vector<double>* ts)
{
    int start_frame = 0;
	int frame_skip = 1;
	double project_t;
	double project_scale;
	vector<double>::iterator ts_i;
    
	project_scale = 2.0 / num_points * frame_skip / draw_fraction;
    project_t = 0.0 - (frame) * project_scale / frame_skip;
	if(project_t < -1.0)
		project_t = -1.0;
    
	if(start_frame < frame - num_points * draw_fraction / 2)
		start_frame = frame - num_points * draw_fraction / 2;
    
    glPushMatrix();
    glScaled(x_scale, y_scale, z_scale);
    
    glLineWidth(scale);
    // draw time series segment
    glColor4dv(pred_color);
	glBegin(GL_LINE_STRIP);
    for(ts_i = ts->begin() + start_frame; project_t < 0.0; ts_i+=frame_skip, project_t += project_scale)
	{
        if(*ts_i != 0.0)
            glVertex2d(project_t, *ts_i);
	}
    glEnd();
    
    glBegin(GL_POINTS);
    glVertex2d(project_t+tp*project_scale, *(ts_i+tp-1));
    glEnd();
    
    // draw rest of time series segment
    glBegin(GL_LINE_STRIP);
	for(vector<double>::iterator ts_j = ts_i; (project_t < 1.0 && ts_j < ts->end()); ts_j+=frame_skip, project_t += project_scale)
	{
		glColor4d(project_t/8+0.75+pred_color[0], project_t/8+0.75+pred_color[1], project_t/8+0.75+pred_color[2], pred_color[3]);
		glVertex2d(project_t, *ts_j);
	}
	glEnd();
    
    glPopMatrix();
    
    return;
}

void attractor::draw_xmap_ts(int frame, const int lag, const int skip, const int dim, 
                             const double sat, const double line_width, vector<double>* ts)
{
    int start_frame = 0;
	int frame_skip = 1;
	int frame_skip_2 = 10;
	double project_t;
	double project_scale;
	double r = lag_dim == 1, g = lag_dim == 2, b = lag_dim == 3;
	vector<double>::iterator ts_i, ts_k1, ts_k2;
    double point_height, lagged_point_height;
	double y_max = d + 0.8;
	double y_min = d - 0.8;
	char label_char;
	
    start_frame = skip - lag;
    
	project_scale = 2.0 / num_points * frame_skip / draw_fraction;
    project_t = 0.0 - (frame-skip) * project_scale / frame_skip;
	if(project_t < -1.0)
		project_t = -1.0;
    
	if(start_frame < frame - lag - num_points * draw_fraction / 2)
		start_frame = frame - lag - num_points * draw_fraction / 2;
    
    glPushMatrix();
    glScaled(x_scale, y_scale, z_scale);
    
    glLineWidth(scale * line_width);
    // draw time series segment
    glColor4d(r, g, b, 0.5);
	glBegin(GL_LINE_STRIP);
    for(ts_i = ts->begin() + start_frame; project_t < 0.0; ts_i+=frame_skip, project_t += project_scale)
	{
        if(*ts_i != 0.0)
            glVertex2d(project_t, *ts_i);
	}
    glEnd();
    
    lagged_point_height = *(ts_i-1+lag);
    point_height = *(ts_i-1);
    
    // draw rest of time series segment
    glBegin(GL_LINE_STRIP);
	for(vector<double>::iterator ts_j = ts_i; (project_t < 1.0 && ts_j < ts->end()); ts_j+=frame_skip, project_t += project_scale)
	{
		glColor3d(project_t/8+0.75*sat+r, project_t/8+0.75*sat+g, project_t/8+0.75*sat+b);
		glVertex2d(project_t, *ts_j);
	}
	glEnd();
    
    glPopMatrix();
    
    return;
}

void attractor::draw_half_ts(int frame, const int lag, const int skip, const int dim, 
                        const double sat, const double line_width)
{
	int start_frame = 0;
	int frame_skip = 2;
	int frame_skip_2 = 10;
	double project_t;
	double project_scale;
	vector<double>* ts;
	double r = 0.0, g = 0.0, b = 0.0;
	vector<double>::iterator ts_i, ts_k1, ts_k2;
    double point_height, lagged_point_height;
	double y_max = 2*d;
	double y_min = 0.0;
    int texture_index;
	
	switch(dim)
	{
		case 1:
			ts = &x;
			r = sat;
            texture_index = X_LABEL_TEXTURE;
			break;
		case 2:
			ts = &y;
			g = sat;
			texture_index = Y_LABEL_TEXTURE;
			break;
		case 3:
			ts = &z;
			b = sat;
			texture_index = Z_LABEL_TEXTURE;
			break;
	}
    if(skip != 0)
    {
        texture_index ++;
        if(lag > 0)
            texture_index++;
        if(lag == 16*tau)
            texture_index ++;
    }
	
    start_frame = skip - lag;
    
	project_scale = 2.0 / num_points * frame_skip / draw_fraction;
    project_t = 0.0 - (frame-skip) * project_scale / frame_skip;
	if(project_t < -1.0)
		project_t = -1.0;
    
	if(start_frame < frame - lag - num_points * draw_fraction / 2)
		start_frame = frame - lag - num_points * draw_fraction / 2;
    
    // draw axis label
    enqueue_label(-1.02*x_scale, d*1.2*y_scale, 0, texture_index, texture_scale, 2, 0);
    
    glPushMatrix();
    glScaled(x_scale, y_scale, z_scale);    
    
	// draw axes
    glLineWidth(scale);
	glColor3d(LAG_SAT, LAG_SAT, LAG_SAT);
	glBegin(GL_LINE_STRIP);
	glVertex2d(-1.0, y_max);
	glVertex2d(-1.0, y_min);
	glVertex2d(0.0, y_min);
	glEnd();
    
    // draw central line
    glLineWidth(2.0*scale);
    glColor4dv(point_color);
	glBegin(GL_LINES);
	glVertex2d(0.0, 2*d);
	glVertex2d(0.0, 0.0);
	glEnd();
    
    glLineWidth(scale * (LINE_WIDTH + line_width));
    // draw time series segment
	glColor3d(r, g, b);
	glBegin(GL_LINE_STRIP);
    for(ts_i = ts->begin() + start_frame; project_t < 0.0; ts_i+=frame_skip, project_t += project_scale)
	{
		glVertex2d(project_t, *ts_i);
	}
    glEnd();
    
    lagged_point_height = *(ts_i-1+lag);
    point_height = *(ts_i-1);
    
    // draw point
    if(lag == 0)
    {
        glPointSize(POINT_WIDTH*scale);
        glColor3d(r, g, b);
        glBegin(GL_POINTS);
        glVertex2d(0, point_height);
        glEnd();
    }
    
    glPopMatrix();
    
	return;
}

void attractor::draw_ts(int frame, const int lag, const int skip, const int dim, 
                        const double sat, const double line_width)
{
	int start_frame = 0;
	int frame_skip = 2;
	int frame_skip_2 = 10;
	double project_t;
	double project_scale;
	vector<double>* ts;
	double r = 0.0, g = 0.0, b = 0.0;
	vector<double>::iterator ts_i, ts_k1, ts_k2;
    double point_height, lagged_point_height;
	double y_max = 2*d;
	double y_min = 0.0;
    int texture_index;
	
	switch(dim)
	{
		case 1:
			ts = &x;
			r = sat;
            texture_index = X_LABEL_TEXTURE;
			break;
		case 2:
			ts = &y;
			g = sat;
			texture_index = Y_LABEL_TEXTURE;
			break;
		case 3:
			ts = &z;
			b = sat;
			texture_index = Z_LABEL_TEXTURE;
			break;
	}
    if(skip != 0)
    {
        texture_index ++;
        if(lag > 0)
            texture_index++;
        if(lag == 16*tau)
            texture_index ++;
    }
	
    start_frame = skip - lag;
    
	project_scale = 2.0 / num_points * frame_skip / draw_fraction;
    project_t = 0.0 - (frame-skip) * project_scale / frame_skip;
	if(project_t < -1.0)
		project_t = -1.0;
    
	if(start_frame < frame - lag - num_points * draw_fraction / 2)
		start_frame = frame - lag - num_points * draw_fraction / 2;
    
    // draw axis label
    enqueue_label(-1.02*x_scale, d*1.2*y_scale, 0, texture_index, texture_scale, 2, 0);
    
    glPushMatrix();
    glScaled(x_scale, y_scale, z_scale);    
    
	// draw axes
    glLineWidth(scale);
	glColor3d(LAG_SAT, LAG_SAT, LAG_SAT);
	glBegin(GL_LINE_STRIP);
	glVertex2d(-1.0, y_max);
	glVertex2d(-1.0, y_min);
	glVertex2d(1.0, y_min);
	glVertex2d(1.0, y_max);
	glEnd();
    
    // draw central line
    glLineWidth(2.0*scale);
    glColor4dv(point_color);
	glBegin(GL_LINES);
	glVertex2d(0.0, 2*d);
	glVertex2d(0.0, 0.0);
	glEnd();
    
    glLineWidth(scale * (LINE_WIDTH + line_width));
    // draw time series segment
	glColor3d(r, g, b);
	glBegin(GL_LINE_STRIP);
    for(ts_i = ts->begin() + start_frame; project_t < 0.0; ts_i+=frame_skip, project_t += project_scale)
	{
		glVertex2d(project_t, *ts_i);
	}
    glEnd();
    
    lagged_point_height = *(ts_i-1+lag);
    point_height = *(ts_i-1);
    
    // draw rest of time series segment
    glBegin(GL_LINE_STRIP);
	for(vector<double>::iterator ts_j = ts_i; (project_t < 1.0 && ts_j < ts->end()); ts_j+=frame_skip, project_t += project_scale)
	{
		glColor3d(project_t/8+0.75*sat+r, project_t/8+0.75*sat+g, project_t/8+0.75*sat+b);
		glVertex2d(project_t, *ts_j);
	}
	glEnd();
    
    // draw point
    if(lag == 0)
    {
        glPointSize(POINT_WIDTH*scale);
        glColor3d(r, g, b);
        glBegin(GL_POINTS);
        glVertex2d(0, point_height);
        glEnd();
    }
    
    if(DEBUG && skip != 0 && lag != 0 && VIEW != XMAP_TS)
    {
        if(lag != 16*tau || DEBUG > 1)
        {
            glLineWidth(2.0*scale);
            glColor3d(0.2, 0.2, 0.2);
            glBegin(GL_LINE_STRIP);
            glVertex2d(0, lagged_point_height);
            glVertex2d(lag*project_scale/frame_skip, lagged_point_height);
            glEnd();
            glPointSize(POINT_WIDTH*scale);
            glColor3d(r/sat, g/sat, b/sat);
            glBegin(GL_POINTS);
            glVertex2d(lag*project_scale/frame_skip, lagged_point_height);
            glEnd();
        }
	}
    if(TSVIEW)
    {
        project_scale = 2.0 / num_points * frame_skip_2;
        project_t = -1.0;
        glTranslated(0.0, -0.6, 0.0);
        glScaled(1.0, 0.4, 1.0);
        glLineWidth(scale);
        glBegin(GL_LINE_STRIP);
        glColor3d(LAG_SAT, LAG_SAT, LAG_SAT);
        vector<double>::iterator ts_k1, ts_k2;
        for(ts_k1 = ts->begin(); ts_k1 < ts->begin()+start_frame; ts_k1+=frame_skip_2, project_t += project_scale)
        {
            glVertex2d(project_t, *ts_k1);
        }
        glColor3d(r, g, b);
        for(ts_k2 = ts_k1; ts_k2 < ts->begin()+frame; ts_k2+=frame_skip_2, project_t += project_scale)
        {
            glVertex2d(project_t, *ts_k2);
        }
        glColor3d(LAG_SAT, LAG_SAT, LAG_SAT);
        for(ts_k1 = ts_k2; ts_k1 < ts->end(); ts_k1+=frame_skip_2, project_t += project_scale)
        {
            glVertex2d(project_t, *ts_k1);
        }
        glEnd();
    }
    
    glPopMatrix();
    
    if(DEBUG && skip != 0 && lag != 0)
    {
        if (lag == 8*tau)
        {
            enqueue_label((double(lag)/320/tau)*x_scale, (lagged_point_height)*y_scale, 0, TAU_LABEL_TEXTURE, texture_scale, 1, 0);
        }
        else if (lag == 16*tau && DEBUG > 1)
        {
            enqueue_label((double(lag)/320/tau)*x_scale, (lagged_point_height)*y_scale, 0, TAU_LABEL_TEXTURE+1, texture_scale, 1, 0);
        }
	}
    
	return;
}

void attractor::draw_embedding(vector<double>::iterator x_i, vector<double>::iterator y_i, vector<double>::iterator z_i, const int frame)
{
    // set current point
    curr_x = *(x_i + frame-1);
    curr_y = *(y_i + frame-1);
    curr_z = *(z_i + frame-1);
    
    glBegin(GL_LINE_STRIP);
	for(int i = 0; i < frame; i++, x_i++, y_i++, z_i++)
	{
		col(*x_i, *y_i, *z_i);
		glVertex3d(*x_i, *y_i, *z_i);
	}
	glEnd();
    
    /*
     glColor3d(1.0, 0.0, 1.0);
     glPointSize(20);
     glBegin(GL_POINTS);
     glVertex3d(d, d, d);
     glEnd();
     */
    
    return;
}

void attractor::draw_tracers(const int frame)
{	
	int frame_skip = 2;
	double project_t;
	double project_scale = 0.0014*frame_skip;
	
	switch(x_tracer)
	{
		case NONE:
			break;
		case TRACE:
			glColor3d(1.0, 0.0, 0.0);
			project_t = -project_scale/frame_skip * (frame-x_start_time-1);
			glBegin(GL_LINE_STRIP);
			for(int i = x_start_time; i < frame; i+=frame_skip, project_t += project_scale)
			{
				glVertex3d(x[i], project_t, 0);
			}
			glEnd();
            glPointSize(POINT_WIDTH*scale);
			glBegin(GL_POINTS);
			glVertex3d(x[frame-1], 0, 0);
			glEnd();
		case PROJECT:
			glColor3d(1.0, 0.0, 0.0);
			glBegin(GL_LINE_STRIP);
			glVertex3d(x[frame-1], y[frame-1], z[frame-1]);
			glVertex3d(x[frame-1], y[frame-1], 0);
			glVertex3d(x[frame-1], 0, 0);
			glEnd();
			break;
	}
	switch(y_tracer)
	{
		case NONE:
			break;
		case TRACE:
			glColor3d(0.0, 1.0, 0.0);
			project_t = -project_scale/frame_skip * (frame-y_start_time-1);
			glBegin(GL_LINE_STRIP);
			for(int i = y_start_time; i < frame; i+=frame_skip, project_t += project_scale)
			{
				glVertex3d(project_t, y[i], 0);
			}
			glEnd();
            glPointSize(POINT_WIDTH*scale);
			glBegin(GL_POINTS);
			glVertex3d(0, y[frame-1], 0);
			glEnd();
		case PROJECT:
			glColor3d(0.0, 1.0, 0.0);
			glBegin(GL_LINE_STRIP);
			glVertex3d(x[frame-1], y[frame-1], z[frame-1]);
			glVertex3d(0, y[frame-1], z[frame-1]);
			glVertex3d(0, y[frame-1], 0);
			glEnd();
			break;
	}
	switch(z_tracer)
	{
		case NONE:
			break;
		case TRACE:
			glColor3d(0.0, 0.0, 1.0);
			project_t = -project_scale/frame_skip * (frame-z_start_time-1);
			glBegin(GL_LINE_STRIP);
			for(int i = z_start_time; i < frame; i+=frame_skip, project_t += project_scale)
			{
				glVertex3d(project_t, 0, z[i]);
			}
			glEnd();
            glPointSize(POINT_WIDTH*scale);
			glBegin(GL_POINTS);
			glVertex3d(0, 0, z[frame-1]);
			glEnd();
		case PROJECT:
			glColor3d(0.0, 0.0, 1.0);
			glBegin(GL_LINE_STRIP);
			glVertex3d(x[frame-1], y[frame-1], z[frame-1]);
			glVertex3d(x[frame-1], 0, z[frame-1]);
			glVertex3d(0, 0, z[frame-1]);
			glEnd();
			break;
	}	
	
	return;
}

void attractor::draw_axes(bool lag, int lag_dim)
{
    double r = 0.0, g = 0.0, b = 0.0;
    int texture_index;
    if(lag)
    {
        switch(lag_dim)
        {
            case 1:
                r = 1.0;
                texture_index = X_LABEL_TEXTURE;
                break;
            case 2:
                g = 1.0;
                texture_index = Y_LABEL_TEXTURE;
                break;
            case 3:
                b = 1.0;
                texture_index = Z_LABEL_TEXTURE;
                break;
        }
        draw_axis(1, r, g, b, 0, texture_index+1, 0);
        draw_axis(2, r*LAG_SAT, g*LAG_SAT, b*LAG_SAT, 0.5, texture_index+2, 1);
        draw_axis(3, r*LAG2_SAT, g*LAG2_SAT, b*LAG2_SAT, 1, texture_index+3, 2);
    }
    else
    {
        draw_axis(1, 1, 0, 0, 0, X_LABEL_TEXTURE, -1);
        draw_axis(2, 0, 1, 0, 0, Y_LABEL_TEXTURE, -1);
        draw_axis(3, 0, 0, 1, 0, Z_LABEL_TEXTURE, -1);
	}
    return;
}

void attractor::draw_lag_axis(int direction, int dim, int lag)
{
    double r = 0.0, g = 0.0, b = 0.0;
    int texture_index;
    double sat;
    sat = 1.0;
    switch(dim)
    {
        case 1:
            r = 1.0;
            texture_index = X_LABEL_TEXTURE;
            break;
        case 2:
            g = 1.0;
            texture_index = Y_LABEL_TEXTURE;
            break;
        case 3:
            b = 1.0;
            texture_index = Z_LABEL_TEXTURE;
            break;
    }
    texture_index = texture_index + lag + 1;
    switch(lag)
    {
        case 1:
            sat = LAG_SAT;
            break;
        case 2:
            sat = LAG2_SAT;
            break;
    }
    draw_axis(direction, r*sat, g*sat, b*sat, 0, texture_index, 0);
}

void attractor::draw_axis(const int axis, const double r, const double g, const double b, 
                          const double delta_line_width, const int texture_index, const int lag)
{
    double x_len = 0, y_len = 0, z_len = 0;
    double k_x = d/30, k_y = d/30, k_z = d/30;
    double axis_scale = 59.0/60.0;
    double label_scale = 1.12;
    
    switch(axis)
    {
        case 1: // x
            x_len = 2*d;
            k_x = 0;
            break;
        case 2: // y
            y_len = 2*d;
            k_y = 0;
            break;
        case 3: // z
            z_len = 2*d;
            k_z = 0;
            break;
        default:
            cerr << "ERROR (attractor): unknown axis input: " << axis << "\n";
            exit(1);
    }
	if(VIEW == XMAP || VIEW == XMAP_TS || VIEW == UNIVARIATE_TS)
        k_x *= 1.5;
    // label
    enqueue_label(x_scale*(x_len*label_scale), y_scale*(y_len*label_scale), z_scale*(z_len*label_scale), texture_index, texture_scale, 1-(axis==2), 1+(axis==2));
    
    glPushMatrix();
    glScaled(x_scale, y_scale, z_scale);
    
	glLineWidth(scale*(LINE_WIDTH+delta_line_width));
    glColor3d(r, g, b);
    
    // draw axis line
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(x_len, y_len, z_len);
	glEnd();
    
    // draw axis end
    if(DRAW_CONE)
    {
        glPushMatrix();
        glTranslated(x_len, y_len, z_len);
        glRotated(90, -y_len, x_len, z_len);
        glutSolidCone(0.018, 0.090, 12, 2);
        glPopMatrix();
    }
    else
    {
        glBegin(GL_LINES);
        glVertex3d(x_len, y_len, z_len);
        glVertex3d(x_len*axis_scale+k_x, y_len*axis_scale, z_len*axis_scale);
        glVertex3d(x_len, y_len, z_len);
        glVertex3d(x_len*axis_scale-k_x, y_len*axis_scale, z_len*axis_scale);
        glVertex3d(x_len, y_len, z_len);
        glVertex3d(x_len*axis_scale, y_len*axis_scale+k_y, z_len*axis_scale);
        glVertex3d(x_len, y_len, z_len);
        glVertex3d(x_len*axis_scale, y_len*axis_scale-k_y, z_len*axis_scale);
        glVertex3d(x_len, y_len, z_len);
        glVertex3d(x_len*axis_scale, y_len*axis_scale, z_len*axis_scale+k_z);
        glVertex3d(x_len, y_len, z_len);
        glVertex3d(x_len*axis_scale, y_len*axis_scale, z_len*axis_scale-k_z);
        glVertex3d(x_len, y_len, z_len);
        glEnd();
    }
    glPopMatrix();
    
    return;
}

void attractor::draw_labels()
{
    texture_2d curr_texture;
    double x_pos;
    double y_pos;
    
    // sort textures by depth
    sort(texture_queue.begin(), texture_queue.end());
    
    // draw textures
    for(vector<texture_label>::iterator iter = texture_queue.begin(); iter != texture_queue.end(); iter++)
    {
        glPushMatrix();
        
        // reset projection matrix
        glLoadIdentity();
        
        switch(VIEW)
        {
            case MANIFOLD:
            case RECONSTRUCTION:
            case GENERIC_RECONSTRUCTION:
            case SHADOW:
            case UNIVARIATE:
            case UNIVARIATE_TS:
            case XMAP:
            case XMAP_TS:
                glMultMatrixd(iter->billboard_matrix);
                break;
            case TIME_SERIES:
            case LAGS:
                glMultMatrixd(iter->billboard_matrix);
                if(window_width > window_height)
                {
                    glScaled(double(window_height)/double(window_width), 1.0, 1.0);
                }
                else
                {
                    glScaled(1.0, double(window_width)/double(window_height), 1.0);
                }
                
                break;
            default:
                glPopMatrix();
                return;
        }
        
        // load texture
        curr_texture = my_textures[iter->index];
        glBindTexture(GL_TEXTURE_2D, curr_texture.texture_id);
        
        // draw
        x_pos = curr_texture.width/2.0*(1-iter->x_peg);
        y_pos = curr_texture.height/2.0*(iter->y_peg-1);
        glColor3d(1.0, 1.0, 1.0);
        glScaled(iter->scale, iter->scale, 1.0);
        glEnable(GL_TEXTURE_2D);
        glBegin(GL_QUADS);
        glNormal3f(0.0, 0.0, -1.0);
        glTexCoord2d(0.0, 0.0); glVertex3d(-curr_texture.width/2.0+x_pos, -curr_texture.height/2.0+y_pos, 0.0);
        glTexCoord2d(0.0, 1.0); glVertex3d(-curr_texture.width/2.0+x_pos, curr_texture.height/2.0+y_pos, 0.0);
        glTexCoord2d(1.0, 1.0); glVertex3d(curr_texture.width/2.0+x_pos, curr_texture.height/2.0+y_pos, 0.0);
        glTexCoord2d(1.0, 0.0); glVertex3d(curr_texture.width/2.0+x_pos, -curr_texture.height/2.0+y_pos, 0.0);
        glEnd();
        glDisable(GL_TEXTURE_2D);
        
        glPopMatrix();
        
    }
    return;
}

void attractor::enqueue_label(double pos_x, double pos_y, double pos_z, int texture_index, double scale, int x_peg, int y_peg)
{
    
    // get camera position and vectors
    double modelview_matrix[16], billboard_matrix[16], curr_matrix[16];
    double camera_pos[3], new_camera_pos[3], pos[3];
    double camera_up[3], look[3];
    double bb_up[3], bb_right[3];
    double temp;
    
    // save in curr_label
    texture_label curr_label;
    
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix);
    glGetDoublev(GL_MODELVIEW_MATRIX, curr_matrix);
    
    // save id
    curr_label.index = texture_index;
    
    // camera position
    switch(VIEW)
    {
        case TIME_SERIES:
        case LAGS:
            camera_pos[0] = pos_x;
            camera_pos[1] = pos_y;
            camera_pos[2] = 1.5;
            break;
        default:
            camera_pos[0] = - modelview_matrix[12];
            camera_pos[1] = - modelview_matrix[13];
            camera_pos[2] = - modelview_matrix[14];
            break;
    }
    
    // camera up vector
    camera_up[0] = modelview_matrix[1];
    camera_up[1] = modelview_matrix[5];
    camera_up[2] = modelview_matrix[9];
    
    // zero the translation in the matrix, so we can use the matrix to transform
	// camera postion to world coordinates using the view matrix
	modelview_matrix[12] = modelview_matrix[13] = modelview_matrix[14] = 0;
    
    // transpose rotation
    temp = modelview_matrix[1];
    modelview_matrix[1] = modelview_matrix[4];
    modelview_matrix[4] = temp;
    temp = modelview_matrix[2];
    modelview_matrix[2] = modelview_matrix[8];
    modelview_matrix[8] = temp;
    temp = modelview_matrix[6];
    modelview_matrix[6] = modelview_matrix[9];
    modelview_matrix[9] = temp;
    
    // correct position of camera in world space
    new_camera_pos[0] = modelview_matrix[0] * camera_pos[0] + 
    modelview_matrix[4] * camera_pos[1] + 
    modelview_matrix[8] * camera_pos[2] + 
    modelview_matrix[12];
    new_camera_pos[1] = modelview_matrix[1] * camera_pos[0] + 
    modelview_matrix[5] * camera_pos[1] + 
    modelview_matrix[9] * camera_pos[2] + 
    modelview_matrix[13];
    new_camera_pos[2] = modelview_matrix[2] * camera_pos[0] + 
    modelview_matrix[6] * camera_pos[1] + 
    modelview_matrix[10] * camera_pos[2] + 
    modelview_matrix[14];
    
    // look vector: pos -> camera_pos
    look[0] = new_camera_pos[0] - pos_x;
    look[1] = new_camera_pos[1] - pos_y;
    look[2] = new_camera_pos[2] - pos_z;
    
    // normalize the look vector
    temp = sqrt(look[0]*look[0] + look[1]*look[1] + look[2]*look[2]);
    if(temp != 0)
    {
        look[0] = look[0] / temp;
        look[1] = look[1] / temp;
        look[2] = look[2] / temp;
    }    
    
    // use camera up vector as billboard up vector
    bb_up[0] = camera_up[0];
    bb_up[1] = camera_up[1];
    bb_up[2] = camera_up[2];
    
    // compute billboard right vector (bb_right = bb_up cross look)
    bb_right[0] = bb_up[1]*look[2] - bb_up[2]*look[1];
    bb_right[1] = bb_up[2]*look[0] - bb_up[0]*look[2];
    bb_right[2] = bb_up[0]*look[1] - bb_up[1]*look[0];
    
    // normalize right vector in case of floating point error
    temp = sqrt(bb_right[0]*bb_right[0] + bb_right[1]*bb_right[1] + bb_right[2]*bb_right[2]);
    if(temp != 0)
    {
        bb_right[0] = bb_right[0] / temp;
        bb_right[1] = bb_right[1] / temp;
        bb_right[2] = bb_right[2] / temp;
    }
    
    // compute correct look vector (look = right cross up)
    look[0] = bb_right[1]*bb_up[2] - bb_right[2]*bb_up[1];
    look[1] = bb_right[2]*bb_up[0] - bb_right[0]*bb_up[2];
    look[2] = bb_right[0]*bb_up[1] - bb_right[1]*bb_up[0];
    
    // make billboard matrix
    billboard_matrix[0] = bb_right[0];
    billboard_matrix[1] = bb_right[1];
    billboard_matrix[2] = bb_right[2];
    billboard_matrix[3] = 0;
    billboard_matrix[4] = bb_up[0];
    billboard_matrix[5] = bb_up[1];
    billboard_matrix[6] = bb_up[2];
    billboard_matrix[7] = 0;
    billboard_matrix[8] = look[0];
    billboard_matrix[9] = look[1];
    billboard_matrix[10] = look[2];
    billboard_matrix[11] = 0;
    billboard_matrix[12] = pos_x;
    billboard_matrix[13] = pos_y;
    billboard_matrix[14] = pos_z;
    billboard_matrix[15] = 1;
    
    // compute combined matrix and depth
    glPushMatrix();
    glLoadIdentity();
    
    // apply both transformations
    glMultMatrixd(curr_matrix);
    glMultMatrixd(billboard_matrix);
    
    // store single transformation matrix
    glGetDoublev(GL_MODELVIEW_MATRIX, curr_label.billboard_matrix);
    
    
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            if(i==j)
                curr_label.billboard_matrix[i*4+j] = 1.0;
            else
                curr_label.billboard_matrix[i*4+j] = 0.0;
    
    
    
    // compute depth of texture
    curr_label.z_pos = curr_label.billboard_matrix[14];
    
    // set other texture drawing properties
    curr_label.scale = scale;
    curr_label.x_peg = x_peg;
    curr_label.y_peg = y_peg;
    
    glPopMatrix();   
    
    texture_queue.push_back(curr_label);
    return;
}

void attractor::draw_curve(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2)
{
    double ctrl_points[4][3];
    ctrl_points[0][0] = x1;
    ctrl_points[0][1] = y1;
    ctrl_points[0][2] = z1;
    
    ctrl_points[3][0] = x2;
    ctrl_points[3][1] = y2;
    ctrl_points[3][2] = z2;
    
    ctrl_points[1][0] = 0.5 * ctrl_points[0][0] + 0.5 * ctrl_points[3][0];
    ctrl_points[1][1] = 0.8 * ctrl_points[0][1] + 0.2 * ctrl_points[3][1];
    ctrl_points[1][2] = 0.5 * ctrl_points[0][2] + 0.5 * ctrl_points[3][2];
    
    ctrl_points[2][0] = 0.2 * ctrl_points[0][0] + 0.8 * ctrl_points[3][0];
    ctrl_points[2][1] = 0.5 * ctrl_points[0][1] + 0.5 * ctrl_points[3][1];
    ctrl_points[2][2] = 0.2 * ctrl_points[0][2] + 0.8 * ctrl_points[3][2];
    
    glMap1d(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 4, &ctrl_points[0][0]);
    glEnable(GL_MAP1_VERTEX_3);
    
    glColor4dv(neighbor_color);
    glLineWidth(scale * LINE_WIDTH);
    glBegin(GL_LINE_STRIP);
    for (int i = 0; i <= 50; i++)
    {
        glEvalCoord1d(double(i)/50.0);
    }
    glEnd();
    glDisable(GL_MAP1_VERTEX_3);
    
    return;
}

void attractor::col(const double x, const double y, const double z)
{
    double sat;
    if(!COLOR_METHOD)
    {
        sat = x - d + (1-d);
        glColor3d(sat, sat, sat);
    }
    else
    {
        sat = exp(-8*d * dist_to_curr(x, y, z));
        glColor4d(0, 0, 0, sat);
    }
    return;
}

double attractor::N(const int i, const int k, const double u)
{
	double num1, denom1, num2, denom2, temp1, temp2;
	// using Cox de Boor algorithm from handout
	if(k == 1)
    {
		if(u >= knot[i] && u < knot[i+1])
        {
			return 1;
        }
		else
        {
			return 0;
        }
    }
    else
    {
        num1 = (u-knot[i])*N(i,k-1,u);
        denom1 = knot[i+k-1]-knot[i];
        if(num1 == 0 && denom1 == 0)
            temp1 = 0;
        else
            temp1 = num1/denom1;
        
        num2 = (knot[i+k]-u)*N(i+1,k-1,u);
        denom2 = knot[i+k]-knot[i+1];
        if(num2 == 0 && denom2 == 0)
            temp2 = 0;
        else
            temp2 = num2/denom2;
        
        return temp1 + temp2;
    }
	return 0;
}

void attractor::generate_movie()
{
    vector<double> p_x;
    vector<double> p_y;
    vector<double> p_z;
    
    // set up control points
    p_x.clear();
    p_y.clear();
    p_z.clear();
    knot.clear();
    
    // control point 0
    p_x.push_back(0);
    p_y.push_back(0);
    p_z.push_back(0);
    
    // control point 1
    p_x.push_back(5000);
    p_y.push_back(2000);
    p_z.push_back(0);
    
    // control point 1
    p_x.push_back(10000);
    p_y.push_back(8000);
    p_z.push_back(-6000);
    
    // control point 2
    p_x.push_back(2000);
    p_y.push_back(4000);
    p_z.push_back(-3000);
    
    // control point 3
    p_x.push_back(0);
    p_y.push_back(0);
    p_z.push_back(0);
    
    knot.push_back(0);
	knot.push_back(0);
	knot.push_back(0);
	knot.push_back(0);
    knot.push_back(5000);
	knot.push_back(num_points);
	knot.push_back(num_points);
	knot.push_back(num_points);
	knot.push_back(num_points);
    
    // computations
    int num_control_points = p_x.size();
    double alpha;
    
    for(int i = 0; i < num_points; i++)
    {
        for(int j = 0; j < num_control_points; j++)
        {
            alpha = N(j, 4, double(i));
            phi_x[i] += p_x[j] * alpha;
            phi_y[i] += p_y[j] * alpha;
            phi_z[i] += p_z[j] * alpha;
        }
    }
    return;
}

void attractor::generate_data()
{
    // generate attractor time series
	x[0] = 20;
	y[0] = 20;
	z[0] = 20;
    
    switch(lorenz_sim_mode)
    {
        case RK4:
            RK4_sim();
            break;
        default:
            EULER_sim();
            break;
    };
	return;
}

void attractor::transform_data()
{
	double xmin, ymin, zmin, xmax, ymax, zmax;
	double tx, ty, tz;
	double scale;
	
	xmin = x[0];
	xmax = x[0];
	ymin = y[0];
	ymax = y[0];
	zmin = z[0];
	zmax = z[0];
	for(int i = 1; i < num_points; i++)
	{
		if(x[i] < xmin)
			xmin = x[i];
		else if(x[i] > xmax)
			xmax = x[i];
		if(y[i] < ymin)
			ymin = y[i];
		else if(y[i] > ymax)
			ymax = y[i];
		if(z[i] < zmin)
			zmin = z[i];
		else if(z[i] > zmax)
			zmax = z[i];
	}
	scale = xmax-xmin;
	if(ymax-ymin > scale)
		scale = ymax-ymin;
	if(zmax-zmin > scale)
		scale = zmax-zmin;
	scale = 1.5 * d / scale;
	
	tx = (xmax+xmin)/2;
	ty = (ymax+ymin)/2;
	tz = (zmax+zmin)/2;
	
	for(int i = 0; i < num_points; i++)
	{
		x[i] = (x[i] - tx) * scale + d;
		y[i] = (y[i] - ty) * scale + d;
		z[i] = (z[i] - tz) * scale + d;
	}
	
	return;
}

void attractor::generate_xmaps()
{
    for(int frame = 0; frame < num_points; frame++)
    {
        x_nn_weights[frame].resize(nn_num, -1);
        y_nn_weights[frame].resize(nn_num, -1);
        z_nn_weights[frame].resize(nn_num, -1);
    }
    
    find_neighbors(1);
    find_neighbors(2);
    find_neighbors(3);
    
    double total_weight;
    vector<double>::iterator weight_iter;
    vector<int>::iterator index_iter;
    double pred_x, pred_y, pred_z;
    
    for(int frame = 2*tau + nn_skip*(nn_num-1); frame < num_points; frame++)
    {
        total_weight = 0;
        pred_y = 0;
        pred_z = 0;
        for(index_iter = x_nn_indices[frame].begin(), weight_iter = x_nn_weights[frame].begin(); 
            (index_iter != x_nn_indices[frame].end()) && (weight_iter != x_nn_weights[frame].end());
            index_iter++, weight_iter++)
        {
            pred_y += y[*index_iter] * (*weight_iter);
            pred_z += z[*index_iter] * (*weight_iter);
            total_weight += (*weight_iter);
        }
        x_xmap_y[frame] = pred_y / total_weight;
        x_xmap_z[frame] = pred_z / total_weight;
        
        total_weight = 0;
        pred_x = 0;
        pred_z = 0;
        for(index_iter = y_nn_indices[frame].begin(), weight_iter = y_nn_weights[frame].begin(); 
            (index_iter != y_nn_indices[frame].end()) && (weight_iter != y_nn_weights[frame].end());
            index_iter++, weight_iter++)
        {
            pred_x += x[*index_iter] * (*weight_iter);
            pred_z += z[*index_iter] * (*weight_iter);
            total_weight += (*weight_iter);
        }
        y_xmap_x[frame] = pred_x / total_weight;
        y_xmap_z[frame] = pred_z / total_weight;
        
        total_weight = 0;
        pred_x = 0;
        pred_y = 0;
        for(index_iter = z_nn_indices[frame].begin(), weight_iter = z_nn_weights[frame].begin(); 
            (index_iter != z_nn_indices[frame].end()) && (weight_iter != z_nn_weights[frame].end());
            index_iter++, weight_iter++)
        {
            pred_x += x[*index_iter] * (*weight_iter);
            pred_y += y[*index_iter] * (*weight_iter);
            total_weight += (*weight_iter);
        }
        z_xmap_x[frame] = pred_x / total_weight;
        z_xmap_y[frame] = pred_y / total_weight;        
    }
    return;
}

void attractor::generate_forecasts()
{    
    double total_weight;
    vector<double>::iterator weight_iter;
    vector<int>::iterator index_iter;
    double pred, pred_lag_1, pred_lag_2;
    
    for(int frame = 2*tau + nn_skip*(nn_num-1); frame < num_points-tp; frame++)
    {
        total_weight = 0;
        pred = 0;
        pred_lag_1 = 0;
        pred_lag_2 = 0;
        for(index_iter = x_nn_indices[frame].begin(), weight_iter = x_nn_weights[frame].begin(); 
            (index_iter != x_nn_indices[frame].end()) && (weight_iter != x_nn_weights[frame].end());
            index_iter++, weight_iter++)
        {
            if(*index_iter < num_points-tp)
            {
                pred += x[(*index_iter)+tp] * (*weight_iter);
                pred_lag_1 += x[(*index_iter)+tp-tau] * (*weight_iter);
                pred_lag_2 += x[(*index_iter)+tp-2*tau] * (*weight_iter);
                total_weight += (*weight_iter);
            }
        }
        x_forecast[frame+tp] = pred / total_weight;
        x_forecast_lag_1[frame+tp] = pred_lag_1 / total_weight;
        x_forecast_lag_2[frame+tp] = pred_lag_2 / total_weight;
        
        total_weight = 0;
        pred = 0;
        pred_lag_1 = 0;
        pred_lag_2 = 0;
        for(index_iter = y_nn_indices[frame].begin(), weight_iter = y_nn_weights[frame].begin(); 
            (index_iter != y_nn_indices[frame].end()) && (weight_iter != y_nn_weights[frame].end());
            index_iter++, weight_iter++)
        {
            pred += y[(*index_iter)+tp] * (*weight_iter);
            pred_lag_1 += y[(*index_iter)+tp-tau] * (*weight_iter);
            pred_lag_2 += y[(*index_iter)+tp-2*tau] * (*weight_iter);
            total_weight += (*weight_iter);
        }
        y_forecast[frame+tp] = pred / total_weight;
        y_forecast_lag_1[frame+tp] = pred_lag_1 / total_weight;
        y_forecast_lag_2[frame+tp] = pred_lag_2 / total_weight;
        
        total_weight = 0;
        pred = 0;
        pred_lag_1 = 0;
        pred_lag_2 = 0;
        for(index_iter = z_nn_indices[frame].begin(), weight_iter = z_nn_weights[frame].begin(); 
            (index_iter != z_nn_indices[frame].end()) && (weight_iter != z_nn_weights[frame].end());
            index_iter++, weight_iter++)
        {
            pred += z[(*index_iter)+tp] * (*weight_iter);
            pred_lag_1 += z[(*index_iter)+tp-tau] * (*weight_iter);
            pred_lag_2 += z[(*index_iter)+tp-2*tau] * (*weight_iter);
            total_weight += (*weight_iter);
        }
        z_forecast[frame+tp] = pred / total_weight;
        z_forecast_lag_1[frame+tp] = pred_lag_1 / total_weight;
        z_forecast_lag_2[frame+tp] = pred_lag_2 / total_weight;
    }
    return;
}

void attractor::load_textures()
{
    texture_2d curr_texture;
    
    my_textures.clear();
    
    curr_texture.texture_id = load_texture("x_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("x_t_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("x_t-tau_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("x_t-2tau_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("y_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("y_t_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("y_t-tau_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("y_t-2tau_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("z_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("z_t_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("z_t-tau_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("z_t-2tau_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("m_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("m_x_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("m_y_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("m_z_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("tau_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("2tau_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_1_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_2_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_3_x_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_3_y_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_3_z_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_4_x_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_4_y_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_4_z_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_5_x_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_5_y_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_5_z_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_6_x_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_6_y_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_6_z_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_7_xy_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_7_xz_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_7_yx_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_7_yz_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_7_zx_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("view_7_zy_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("equations_label.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    curr_texture.texture_id = load_texture("takens_theorem.png", curr_texture.width, curr_texture.height);
    my_textures.push_back(curr_texture);
    
    return;
}

GLuint attractor::load_texture(const string filename, int &width, int &height) 
{
    
    //header for testing if it is a png
    png_byte header[8];
    
    //open file as binary
    FILE *fp = fopen(filename.c_str(), "rb");
    if (!fp) {
        return 0;
    }
    
    //read the header
    fread(header, 1, 8, fp);
    
    //test if png
    int is_png = !png_sig_cmp(header, 0, 8);
    if (!is_png) {
        fclose(fp);
        return 0;
    }
    
    //create png struct
    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL,
                                                 NULL, NULL);
    if (!png_ptr) {
        fclose(fp);
        return 0;
    }
    
    //create png info struct
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        png_destroy_read_struct(&png_ptr, (png_infopp) NULL, (png_infopp) NULL);
        fclose(fp);
        return 0;
    }
    
    //create png info struct
    png_infop end_info = png_create_info_struct(png_ptr);
    if (!end_info) {
        png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp) NULL);
        fclose(fp);
        return 0;
    }
    
    //png error stuff, not sure libpng man suggests this.
    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
        fclose(fp);
        return 0;
    }
    
    //init png reading
    png_init_io(png_ptr, fp);
    
    //let libpng know you already read the first 8 bytes
    png_set_sig_bytes(png_ptr, 8);
    
    // read all the info up to the image data
    png_read_info(png_ptr, info_ptr);
    
    //variables to pass to get info
    int bit_depth, color_type;
    png_uint_32 twidth, theight;
    
    // get info about png
    png_get_IHDR(png_ptr, info_ptr, &twidth, &theight, &bit_depth, &color_type,
                 NULL, NULL, NULL);
    
    //update width and height based on png info
    width = twidth;
    height = theight;
    
    // Update the png info struct.
    png_read_update_info(png_ptr, info_ptr);
    
    // Row size in bytes.
    int rowbytes = png_get_rowbytes(png_ptr, info_ptr);
    
    // Allocate the image_data as a big block, to be given to opengl
    png_byte *image_data = new png_byte[rowbytes * height];
    if (!image_data) {
        //clean up memory and close stuff
        png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
        fclose(fp);
        return 0;
    }
    
    //row_pointers is for pointing to image_data for reading the png with libpng
    png_bytep *row_pointers = new png_bytep[height];
    if (!row_pointers) {
        //clean up memory and close stuff
        png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
        delete[] image_data;
        fclose(fp);
        return 0;
    }
    // set the individual row_pointers to point at the correct offsets of image_data
    for (int i = 0; i < height; ++i)
        row_pointers[height - 1 - i] = image_data + i * rowbytes;
    
    //read the png into image_data through row_pointers
    png_read_image(png_ptr, row_pointers);
    
    //Now generate the OpenGL texture object
    GLuint texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, (GLvoid*) image_data);
    glGenerateMipmap(GL_TEXTURE_2D);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    
    //clean up memory and close stuff
    png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
    delete[] image_data;
    delete[] row_pointers;
    fclose(fp);
    
    return texture;
}

void attractor::find_neighbors(const int dim)
{
    vector<double>::iterator x_i, y_i, z_i;
    vector<int> nn_indices;
    vector<double> nn_distances;
	double temp_distance;
	int temp_rank, temp_index;
    double curr_x, curr_y, curr_z;
    int index;
    
    // select right time series
    switch(dim)
    {
        case 1:
            x_i = x.begin();
            y_i = x.begin() + tau;
            z_i = x.begin() + 2*tau;
            break;
        case 2:
            x_i = y.begin();
            y_i = y.begin() + tau;
            z_i = y.begin() + 2*tau;
            break;
        case 3:
            x_i = z.begin();
            y_i = z.begin() + tau;
            z_i = z.begin() + 2*tau;
            break;
        default:
            cerr << "ERROR (attractor): invalid dimension given to find_neighbors, dim = " << dim << ".\n";
            exit(1);
    }
    
    for(int frame = 2*tau, my_frame = 0; frame < num_points; frame++, my_frame++)
    {
        if(my_frame >= nn_skip*nn_num)
        {
            nn_indices.clear();
            nn_distances.clear();
            
            curr_x = *(x_i + my_frame-1);
            curr_y = *(y_i + my_frame-1);
            curr_z = *(z_i + my_frame-1);
            
            // initialize neighbors
            index = my_frame - nn_skip;
            for(int i = 0; i < nn_num; i++, index -= nn_skip)
            {
                nn_indices.push_back(index + 2*tau);
                nn_distances.push_back(dist(*(x_i + index), *(y_i + index), *(z_i + index), curr_x, curr_y, curr_z));
            }
            
            // sort distances
            for(int i = 0; i < nn_num; i++)
            {
                for(int j = nn_num-1; j > i; j--)
                {
                    if(nn_distances[j] < nn_distances[j-1])
                    {
                        temp_distance = nn_distances[j];
                        nn_distances[j] = nn_distances[j-1];
                        nn_distances[j-1] = temp_distance;
                        temp_index = nn_indices[j];
                        nn_indices[j] = nn_indices[j-1];
                        nn_indices[j-1] = temp_index;
                    }
                } 
            }
            
            // search for nn_num nearest neighbors
            for(int i = index; i >= 0; i-=nn_skip)
            {
                temp_distance = dist(*(x_i + i), *(y_i + i), *(z_i + i), curr_x, curr_y, curr_z);
                temp_rank = nn_num;
                for(vector<double>::reverse_iterator j = nn_distances.rbegin(); j < nn_distances.rend(); j++, temp_rank--)
                {
                    if(temp_distance > *j)
                    {
                        break;
                    }	
                }
                if(temp_rank < nn_num)
                {
                    nn_distances.insert(nn_distances.begin()+temp_rank, temp_distance);
                    nn_indices.insert(nn_indices.begin()+temp_rank, i + 2*tau);
                    nn_distances.pop_back();
                    nn_indices.pop_back();				
                }
            }
            
            // compute simplex weights
            switch(dim)
            {
                case 1:
                    for(int i = 0; i < nn_num; i++)
                    {
                        x_nn_weights[frame][i] = exp(-nn_distances[i] / nn_distances[0]);
                        if(x_nn_weights[frame][i] < 0.00001)
                            x_nn_weights[frame][i] = 0.00001;
                    }
                    x_nn_indices[frame] = nn_indices;
                    break;
                case 2:
                    for(int i = 0; i < nn_num; i++)
                    {
                        y_nn_weights[frame][i] = exp(-nn_distances[i] / nn_distances[0]);
                        if(y_nn_weights[frame][i] < 0.00001)
                            y_nn_weights[frame][i] = 0.00001;
                    }
                    y_nn_indices[frame] = nn_indices;
                    break;
                case 3:
                    for(int i = 0; i < nn_num; i++)
                    {
                        z_nn_weights[frame][i] = exp(-nn_distances[i] / nn_distances[0]);
                        if(z_nn_weights[frame][i] < 0.00001)
                            z_nn_weights[frame][i] = 0.00001;
                    }
                    z_nn_indices[frame] = nn_indices;
                    break;
            }
        }
    }
    
    return;
}

void attractor::init(bool MOVIE_MODE)
{
    // modify path to use local resource directory
    CFBundleRef mainBundle = CFBundleGetMainBundle();
    CFURLRef resourcesURL = CFBundleCopyResourcesDirectoryURL(mainBundle);
    char path[300];
    if (!CFURLGetFileSystemRepresentation(resourcesURL, TRUE, (UInt8 *)path, PATH_MAX))
    {
        // error!
    }
    CFRelease(resourcesURL);
    chdir(path);
    
    cerr << "generating movie path...";
    phi_x.resize(num_points, 0);
    phi_y.resize(num_points, 0);
    phi_z.resize(num_points, 0);
    if(MOVIE_MODE)
    {
        generate_movie();
    }
    cerr << "done!\n";    
    
    cerr << "generating attractor...";
	generate_data();
	transform_data();
    cerr << "done!\n";
    
    cerr << "generating cross maps...";
	generate_xmaps();
    cerr << "done!\n";
    
    cerr << "generating forecasts...";
    generate_forecasts();
    cerr << "done!\n";
    
    cerr << "loading textures...";
    load_textures();
    cerr << "done!\n";
    
	return;
}

void attractor::reset_rot_matrix()
{
    switch(VIEW)
	{
		case MANIFOLD:
		case RECONSTRUCTION:
        case GENERIC_RECONSTRUCTION:
		case SHADOW:
        case UNIVARIATE:
        case TIME_SERIES:
        case LAGS:
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			gluPerspective(30.0f,(GLfloat)window_width/(GLfloat)window_height,0.2f,10.0f);
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			glRotated(-90, 1, 0, 0);
            glGetDoublev(GL_MODELVIEW_MATRIX, rot_matrix);
			break;
        case XMAP:
        case XMAP_TS:
        case UNIVARIATE_TS:
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			
			if(0.75*window_width > window_height)
			{
				glOrtho(-1.5*window_width/window_height, 1.5*window_width/window_height, -1.5, 1.5, -1.5, 1.5);
			}
			else
			{
				glOrtho(-2.0, 2.0, -2.0*window_height/window_width, 2.0*window_height/window_width, -1.5, 1.5);
			}
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
            glGetDoublev(GL_MODELVIEW_MATRIX, rot_matrix);
			break;
        default:
            break;
	}
    
    return;
}


void attractor::rotate(const double rx, const double ry, const double rz)
{	
	if(VIEW == XMAP || VIEW == XMAP_TS || VIEW == UNIVARIATE_TS)
	{
		theta += ry;
	}
	else if(VIEW == MANIFOLD || VIEW == RECONSTRUCTION || VIEW == SHADOW || 
            VIEW == UNIVARIATE || VIEW == TIME_SERIES || VIEW == LAGS || VIEW == GENERIC_RECONSTRUCTION)
	{
		glGetDoublev(GL_MODELVIEW_MATRIX, rot_matrix);
		glLoadIdentity();
		glTranslated(-vx, -vy, -vz);
		glRotated(rx, 1, 0, 0);
		glTranslated(vx, vy, vz);
		glMultMatrixd(rot_matrix);
		
		glGetDoublev(GL_MODELVIEW_MATRIX, rot_matrix);
		glLoadIdentity();
		glTranslated(-vx, -vy, -vz);
		glRotated(ry, 0, 1, 0);
		glTranslated(vx, vy, vz);
		glMultMatrixd(rot_matrix);
        
        glGetDoublev(GL_MODELVIEW_MATRIX, rot_matrix);
	}
	
	return;
}

void attractor::translate(const double tx, const double ty, const double tz)
{
	double temp_matrix[16];
	
	glGetDoublev(GL_MODELVIEW_MATRIX, temp_matrix);
    glLoadIdentity();
    glTranslated(tx, ty, tz);
    glMultMatrixd(temp_matrix);
	
	return;
}

void attractor::draw(double & runtime)
{
    double depth = -1;
    double x_scale = 0.002;
    double y_scale = 0.002;
    double x_pos = 0.0;
    double y_pos = 1.4;
    double scale;
    texture_2d curr_texture;
    
	int frame;    
	if(runtime > num_points)
	{
		frame = num_points;
	}
	else
	{
		frame = int(runtime);
	}
    texture_queue.clear();
    
    // draw
	switch(VIEW)
	{
        case TAKENS:
            draw_takens();
            break;
		case MANIFOLD:
			draw_manifold(frame);
            curr_texture = my_textures[VIEW_1_LABEL_TEXTURE];
			break;
		case TIME_SERIES:
			draw_time_series(frame);
            curr_texture = my_textures[VIEW_2_LABEL_TEXTURE];
			break;
		case LAGS:
			draw_lagged_time_series(frame);
            curr_texture = my_textures[VIEW_3_LABEL_TEXTURE+lag_dim-1];
			break;
		case RECONSTRUCTION:
			draw_reconstruction(frame);
            curr_texture = my_textures[VIEW_4_LABEL_TEXTURE+lag_dim-1];
			break;
        case GENERIC_RECONSTRUCTION:
            draw_generic_reconstruction(frame);
            curr_texture = my_textures[VIEW_4_LABEL_TEXTURE+lag_dim-1];
            break;
		case SHADOW:
			draw_shadow(frame);
            curr_texture = my_textures[VIEW_5_LABEL_TEXTURE+lag_dim-1];
            x_pos = -1.5;
            y_pos = 0.0;
			break;
		case UNIVARIATE:
			draw_univariate(frame);
            curr_texture = my_textures[VIEW_6_LABEL_TEXTURE+lag_dim-1];
			break;
        case UNIVARIATE_TS:
            draw_univariate_ts(frame);
            curr_texture = my_textures[VIEW_6_LABEL_TEXTURE+lag_dim-1];
            break;
        case XMAP:
            draw_xmap(frame);
            curr_texture = my_textures[VIEW_7_LABEL_TEXTURE+lag_dim*2-2+(pred_dim > 6-lag_dim-pred_dim)];
            break;
		case XMAP_TS:
			draw_xmap_ts(frame);
            curr_texture = my_textures[VIEW_7_LABEL_TEXTURE+lag_dim*2-2+(pred_dim > 6-lag_dim-pred_dim)];
			break;
	}
	
    // draw texture labels
    draw_labels();
    
    // set up ortho view for fixed location textures
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if(0.75*window_width > window_height)
    {
        glOrtho(-1.5*window_width/window_height, 1.5*window_width/window_height, -1.5, 1.5, -1.5, 1.5);
    }
    else
    {
        glOrtho(-2.0, 2.0, -2.0*window_height/window_width, 2.0*window_height/window_width, -1.5, 1.5);
    }
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    if(false) // draw slide titles
    {
        glBindTexture(GL_TEXTURE_2D, curr_texture.texture_id);

        glColor3d(1.0, 1.0, 1.0);
        glEnable(GL_TEXTURE_2D);
        glBegin(GL_QUADS);
        glNormal3f(0.0, 0.0, 1.0);
        glTexCoord2d(0.0, 0.0); glVertex3d(-curr_texture.width/2.0*x_scale+x_pos, -curr_texture.height/2.0*y_scale+y_pos, depth);
        glTexCoord2d(0.0, 1.0); glVertex3d(-curr_texture.width/2.0*x_scale+x_pos, curr_texture.height/2.0*y_scale+y_pos, depth);
        glTexCoord2d(1.0, 1.0); glVertex3d(curr_texture.width/2.0*x_scale+x_pos, curr_texture.height/2.0*y_scale+y_pos, depth);
        glTexCoord2d(1.0, 0.0); glVertex3d(curr_texture.width/2.0*x_scale+x_pos, -curr_texture.height/2.0*y_scale+y_pos, depth);
        glEnd();
        glDisable(GL_TEXTURE_2D);
    }
    if (VIEW == MANIFOLD && DEBUG)
    {
        x_pos = 1.3;
        y_pos = 1.0;
        
        curr_texture = my_textures[EQUATIONS_LABEL_TEXTURE];
        glBindTexture(GL_TEXTURE_2D, curr_texture.texture_id);
        
        glColor3d(1.0, 1.0, 1.0);
        glEnable(GL_TEXTURE_2D);
        glBegin(GL_QUADS);
        glNormal3f(0.0, 0.0, 1.0);
        glTexCoord2d(0.0, 0.0); glVertex3d(-curr_texture.width/2.0*x_scale+x_pos, -curr_texture.height/2.0*y_scale+y_pos, depth);
        glTexCoord2d(0.0, 1.0); glVertex3d(-curr_texture.width/2.0*x_scale+x_pos, curr_texture.height/2.0*y_scale+y_pos, depth);
        glTexCoord2d(1.0, 1.0); glVertex3d(curr_texture.width/2.0*x_scale+x_pos, curr_texture.height/2.0*y_scale+y_pos, depth);
        glTexCoord2d(1.0, 0.0); glVertex3d(curr_texture.width/2.0*x_scale+x_pos, -curr_texture.height/2.0*y_scale+y_pos, depth);
        glEnd();
        glDisable(GL_TEXTURE_2D);
        
        /*
        glLineWidth(1.0);
        glBegin(GL_LINES);
        glColor3d(0.8, 0.8, 0.8);
        glVertex2d(-1.5, -1.45);
        glVertex2d(1.5, -1.45);
        glColor4d(0.0, 0.0, 0.0, 1.0);
        glVertex2d(-1.5, -1.4);
        glVertex2d(-1.5, -1.5);
        glVertex2d(1.5, -1.4);
        glVertex2d(1.5, -1.5);
        glVertex2d(-1.5, -1.45);
        glVertex2d(3.0*double(frame)/double(num_points) - 1.5, -1.45);
        glEnd();
        */
    }
    
	while(runtime > num_points)
	{
		runtime -= num_points;
		x_start_time = 0;
		y_start_time = 0;
		z_start_time = 0;
	}
	
	return;
}

void attractor::trace_x(const int frame)
{
	switch(x_tracer)
	{
		case NONE:
			x_tracer = PROJECT;
			break;
		case PROJECT:
			x_tracer = TRACE;
			x_start_time = frame;
			break;
		case TRACE:
			x_tracer = NONE;
			break;
	}
	return;
}

void attractor::trace_y(const int frame)
{
	switch(y_tracer)
	{
		case NONE:
			y_tracer = PROJECT;
			break;
		case PROJECT:
			y_tracer = TRACE;
			y_start_time = frame;
			break;
		case TRACE:
			y_tracer = NONE;
			break;
	}
	return;
}

void attractor::trace_z(const int frame)
{
	switch(z_tracer)
	{
		case NONE:
			z_tracer = PROJECT;
			break;
		case PROJECT:
			z_tracer = TRACE;
			z_start_time = frame;
			break;
		case TRACE:
			z_tracer = NONE;
			break;
	}
	return;
}

void attractor::set_view(const int new_view)
{
    double persp_scale = 0.002;
    double ortho_scale = 0.0016;
    
	switch(new_view)
	{
        case -1:
            vx = 0;
			vy = -.2;
			vz = 6;	
			VIEW = TAKENS;
			break;
        case 0:
            texture_scale = persp_scale;
			vx = 0;
			vy = 0;
			vz = init_distance;
			VIEW = GENERIC_RECONSTRUCTION;
            break;
		case 1:
            texture_scale = persp_scale;
			vx = 0;
            vy = 0;
			vz = init_distance;
			VIEW = MANIFOLD;
			break;
		case 2:
            texture_scale = ortho_scale;
            vx = -1;
            vy = 0;
            vz = init_distance;
			VIEW = TIME_SERIES;
			break;
		case 3:
            texture_scale = ortho_scale;
			vx = -1;
			vy = 0;
			vz = init_distance;	
			VIEW = LAGS;
			break;
		case 4:
            texture_scale = persp_scale;
			vx = 0;
			vy = 0;
			vz = init_distance;
			VIEW = RECONSTRUCTION;
			break;
		case 5:
            texture_scale = persp_scale;
			vx = 0;
			vy = 0;
			vz = init_distance;	
			VIEW = SHADOW;
			break;
		case 6:
            texture_scale = persp_scale;
			vx = 0;
			vy = 0;
			vz = init_distance;	
			VIEW = UNIVARIATE;
			break;
        case 7:
            texture_scale = ortho_scale;
			vx = 0;
			vy = 0;
			vz = 0;	
			theta = 0;
			VIEW = UNIVARIATE_TS;
			break;
        case 8:
            texture_scale = ortho_scale;
			vx = 0;
			vy = 0;
			vz = 0;	
			theta = 0;
			VIEW = XMAP;
			break;
		case 9:
            texture_scale = ortho_scale;
			vx = 0;
			vy = 0;
			vz = 0;
			theta = 0;
			VIEW = XMAP_TS;
			break;
		default:
			cerr << "ERROR: UNKNOWN VIEWING MODE: " << new_view << "\n";
			cerr << "USING manifold VIEW INSTEAD.\n";
			break;
	}
    reset_rot_matrix();
    switch(new_view)
	{
        case -1:
            break;
        case 0:
			break;
		case 1:
            //rotate(10, 0, 0);
            //rotate(0, -5, 0);
			break;
		case 2:
            rotate(0, -20, 0);
			break;
		case 3:
            rotate(0, -20, 0);
			break;
		case 4:
			break;
		case 5:
            rotate(-11, -40, 0);
            translate(-0.25, -0.25, -0.25);
			break;
		case 6:
			break;
        case 7:
			break;
        case 8:
			break;
		case 9:
			break;
		default:
			cerr << "ERROR: UNKNOWN VIEWING MODE: " << new_view << "\n";
			cerr << "USING manifold VIEW INSTEAD.\n";
			break;
	}

	return;
}

void attractor::set_lagview(const int dim)
{
	lag_dim = dim;
	if(lag_dim == pred_dim)
	{
		pred_dim = lag_dim%3+1;
	}
	return;
}

void attractor::inc_xview()
{
	x_dim = x_dim % 3 + 1;
	return;
}

void attractor::inc_yview()
{
	y_dim = y_dim % 3 + 1;
	return;
}

void attractor::inc_zview()
{
	z_dim = z_dim % 3 + 1;
	return;
}

void attractor::toggle_predview()
{
    pred_dim = 6 - lag_dim - pred_dim;
    return;
}

void attractor::toggle_tsview()
{
    TSVIEW = !TSVIEW;
    return;
}

void attractor::toggle_color_method()
{
    COLOR_METHOD = !COLOR_METHOD;
    return;
}

void attractor::toggle_manifold_label()
{
    MANIFOLD_LABEL = !MANIFOLD_LABEL;
    return;
}

void attractor::inc_xtau()
{
    x_lag = (++x_lag) % 3;
    return;
}

void attractor::inc_ytau()
{
    y_lag = (++y_lag) % 3;
    return;
}

void attractor::toggle_split_view()
{
    if(VIEW == GENERIC_RECONSTRUCTION)
        z_lag = (++z_lag) % 3;
    else
        SPLIT_VIEW = !SPLIT_VIEW;
    return;
}

void attractor::change_tau(const int delta)
{
	tau += delta;
	if(tau < 0)
		tau = 0;
	return;
}

void attractor::debug_toggle()
{
    if(VIEW == LAGS)
    {
        DEBUG = (DEBUG + 1) % 3;
    }
    else
    {
        DEBUG = !DEBUG;
	}
	return;
}

void attractor::change_scale(const double new_scale)
{
	scale = new_scale;
	return;
}

void attractor::set_window_size(const int width, const int height)
{
    window_width = width;
    window_height = height;
    return;
}