/*
 *  main.cpp
 *  LorenzGL_verHY
 *
 *  Created by Hao Ye on 9/24/10.
 *  Copyright 2010 UCSD. All rights reserved.
 *
 */

#include <cstdlib>
#include <iostream>
#include <GLUT/glut.h>
#include "attractor.h"

using namespace std;

enum button { LEFT, MIDDLE, RIGHT, OFF };
bool PAUSED, PREV_PAUSED;
bool TAKENS_VIEW;
draw_mode PREV_VIEW;
bool FULLSCREEN;
bool MOVIE;

attractor* a;

button my_button; // enum var for which button is pressed
int px, py; // previous mouse position
double rotx, roty; // rotation
double tx, ty, tz; // translation
int my_width, my_height;
double speed;
double p_time;

const double ROTATION_SCALE = 0.5;
const double TRANSLATION_SCALE = 0.01;

const int max_frames = 10000;
double runtime;

void display();
void reshape(int width, int height);
void mouse(int b, int s, int x, int y);
void motion(int x, int y);
void keyboard(unsigned char k, int x, int y);
void idle();
void initGL();

void reset_window_title(int param)
{
    glutSetWindowTitle("Lorenz Simulator ver. HY");
    return;
}

int main(int argc, char* argv[])
{
    if(argc > 1 && strcmp(argv[1], "-m") == 0)
        MOVIE = true;
    
    glutInit(&argc, argv);
    
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    
	glutInitWindowSize(800, 600);
    glutCreateWindow("");
    reset_window_title(0);
    glutPositionWindow(80,80);
    FULLSCREEN = false;
	
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);
    
	PAUSED = false;
    PREV_PAUSED = PAUSED;
    TAKENS_VIEW = false;
	speed = 16;
	
	initGL();
	a = new attractor(max_frames);
	a->init(MOVIE);
	a->set_view(1);
	runtime = 0;
    
    glutMainLoop();
    return 0;
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glPushMatrix();
	a->draw(runtime);
	glPopMatrix();
	
    glutSwapBuffers();
}

void reshape(int width, int height)
{
	if (height == 0)
	{
		height = 1;
	}
    glViewport(0, 0, width, height);
	my_width = width;
	my_height = height;
	a->change_scale(min(width/800.0, height/600.0));
	a->set_window_size(width, height);
	if(a->VIEW == XMAP || a->VIEW == XMAP_TS || a->VIEW == UNIVARIATE_TS)
        a->reset_rot_matrix();
	
	return;
}

void mouse(int b, int s, int x, int y)
{
	switch(b)
	{
		case GLUT_LEFT_BUTTON:
			if(s == GLUT_DOWN)
				my_button = LEFT;
			break;
		case GLUT_MIDDLE_BUTTON:
			if(s == GLUT_DOWN)
				my_button = MIDDLE;
			break;
		case GLUT_RIGHT_BUTTON:
			if(s == GLUT_DOWN)
				my_button = RIGHT;
			break;
		default:
			my_button = OFF;
			break;
	}
	px = x;
	py = y;
	return;
}


void motion(int x, int y)
{
	if (my_button == LEFT)
	{
		rotx = ROTATION_SCALE * double(x - px);
		roty = ROTATION_SCALE * double(y - py);
	}
	else if(my_button == RIGHT)
	{
		tz = TRANSLATION_SCALE * (x - px);
		tx = TRANSLATION_SCALE * (y - py);
	}
	else if(my_button == MIDDLE)
	{
		ty = TRANSLATION_SCALE * (y - py);
	}
	else
	{
	}
	px = x;
	py = y;
	return;
}

void keyboard(unsigned char k, int x, int y)
{
    if(MOVIE)
    {
        switch(k)
        {
            case 27:
            case 'q':
            case 'Q':
                exit (0);
                break;
            case 'f':
            case 'F':
                FULLSCREEN = !FULLSCREEN;
                if (FULLSCREEN)
                {
                    glutFullScreen();
                }
                else
                {
                    glutReshapeWindow(800, 600);
                    glutPositionWindow(80,80);
                    glutTimerFunc(0, reset_window_title, 0);
                }
                break;
            case 'r':
                runtime = 0;
                break;
            case '`':
                a->debug_toggle();
                break;
            case '[':
                speed *= 1.25;
                break;
            case ']':
                speed *= 0.8;
                break;
        }
    }
    else
    {
        switch(k)
        {
            case 27:
            case 'q':
            case 'Q':
                exit (0);
                break;
            case 'f':
            case 'F':
                FULLSCREEN = !FULLSCREEN;
                if (FULLSCREEN)
                {
                    glutFullScreen();
                }
                else
                {
                    glutReshapeWindow(800, 600);
                    glutPositionWindow(80,80);
                    glutTimerFunc(0, reset_window_title, 0);
                }
                break;
            case 'r':
                runtime = 0;
                break;
            case '`':
                a->debug_toggle();
                break;
            case ' ':
                PAUSED = !PAUSED;
                PREV_PAUSED = PAUSED;
                break;
            case 'x':
            case 'X':
                switch(a->VIEW)
			{
				case MANIFOLD:
					a->trace_x(int(runtime));
					break;
                case GENERIC_RECONSTRUCTION:
                    a->inc_xview();
                    break;
                default:
					a->set_lagview(1);
					break;
			}
                break;
            case 'y':
            case 'Y':
                switch(a->VIEW)
            {
                case MANIFOLD:
                    a->trace_y(int(runtime));
                    break;
                case GENERIC_RECONSTRUCTION:
                    a->inc_yview();
                    break;
                default:
                    a->set_lagview(2);
                    break;
            }
                break;
            case 'z':
            case 'Z':
                switch(a->VIEW)
            {
                case MANIFOLD:
                    a->trace_z(int(runtime));
                    break;
                case GENERIC_RECONSTRUCTION:
                    a->inc_zview();
                    break;
                default:
                    a->set_lagview(3);
                    break;
            }
                break;
            case 'm':
            case 'M':
                a->toggle_manifold_label();
                break;
            case 'n':
            case 'N':
                a->toggle_predview();
                break;
            case 'k':
            case 'K':
                a->toggle_color_method();
                break;
            case 'l':
            case 'L':
                a->toggle_tsview();
                break;
            case '[':
                speed *= 1.25;
                break;
            case ']':
                speed *= 0.8;
                break;
            case ',':
                a->inc_xtau();
                break;
            case '.':
                a->inc_ytau();
                break;
            case '/':
                a->toggle_split_view();
                break;
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
            case '0':
                PAUSED = PREV_PAUSED;
                a->set_view(int(k-'0'));
                break;
            case 't':
            case 'T':
                TAKENS_VIEW = !TAKENS_VIEW;
                if(TAKENS_VIEW)
                {
                    PREV_PAUSED = PAUSED;
                    PAUSED = true;
                    PREV_VIEW = a->VIEW;
                    a->set_view(-1);
                }
                else
                {
                    PAUSED = PREV_PAUSED;
                    a->set_view(PREV_VIEW +1);
                }
                
                break;
        }
    }
    return;
}

void idle()
{
	double c_time = glutGet(GLUT_ELAPSED_TIME);
	
	// increment frame appropriately
	if(!PAUSED)
	{
		runtime += (c_time - p_time) / speed;
	}
	p_time = c_time;
	
	// handle motion
    if(!MOVIE)
    {
        a->rotate(roty, rotx, 0);
        a->translate(tz, -tx, ty);
    }
    else
    {
        //a->set_angle(phi, theta, nu);
    }
	tx = 0;
	ty = 0;
	tz = 0;
    rotx = 0;
    roty = 0;
	
	// draw
    glutPostRedisplay();
	return;
}

void initGL()
{
	glShadeModel(GL_SMOOTH);
    
    glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glEnable(GL_MULTISAMPLE);
    
    glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glPolygonMode(GL_FRONT, GL_FILL);
	//glEnable(GL_DEPTH_TEST);
  	
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
	tx = 0;
	ty = 0;
	tz = 0;
	rotx = 0;
	roty = 0;
	p_time = 0;
	return;
}
