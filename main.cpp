#ifdef WIN32
#include <windows.h>
#endif

#if defined (__APPLE__) || defined(MACOSX)
#include <OpenGL/gl.h>
//#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#else //linux
#include <GL/gl.h>
#include <GL/glut.h>
#endif

//other includes
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>

#include <map>
int grid_width;
int grid_height;
int pixel_size;
int win_width;
int win_height;
std::string input;
int numOfPoly;
int numOfPoints;
//double x, y, z;
int e1, e2;
bool key1 = false;
bool key2 = false;
bool key3 = false;
int fresh = 1;
std::vector<std::vector<double>> cList;


void check();
void init();
void lineConnectionXY();
void lineConnectionYZ();
void lineConnectionXZ();
void reshape(int width, int height);
void reshape(int width, int height);
void key(unsigned char ch, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void idle();
void displayMain();
void displayXY();
void displayXZ();
void displayYZ();
void scaling(std::vector<double> factor, int id);
void rotation(double angle, int id, std::vector<double> p0, std::vector<double> p);
void translation(std::vector<double> movement, int id);


class Edge {
public:
    void setEdge(int s, int e) {
        startP = s;
        endP = e;
    }
    int getStartP() {return startP;}
    int getEndP() {return endP;}
private:
    int startP;
    int endP;
    std::vector<std::vector<double>> between;
};


class Points {
public:
    Points(int num) {numOfPoints = num;}
    std::vector<double>& getPoint(int index);
    void setNumOfEdges(int num) {numOfEdges = num;}
    void setCoordinates(std::vector<double> pix);
    int& getNumOfPoints() {return numOfPoints;}
    Edge& getEdges(int key) {
        return edges[key];
    }
    int getEdgesNum() {
        return numOfEdges;
    }
    void addEdge(Edge e) {
        edges[edges.size()] = e;
    }
private:
    std::vector<std::vector<double>> pixels;
    std::map<int, Edge> edges;
    int numOfPoints;
    int numOfEdges;
};

class Poly {
public:
    void setPoly(int num) {numOfPolys = num;}
    int getNumOfPolys() {return numOfPolys;}
    void setList(Points id);
    Points& getPoly(int index);
private:
    int numOfPolys;
    std::vector<Points> polyList;
};

void Poly::setList(Points id) {
    polyList.push_back(id);
}

std::vector<double>& Points::getPoint(int index) {
    return pixels[index];
}

void Points::setCoordinates(std::vector<double> pix) {
    pixels.push_back(pix);
}

Points& Poly::getPoly(int index) {
    return polyList[index];
}

Poly P1;
int main(int argc, char** argv) {
    //the number of pixels in the grid
    grid_width = 700;
    grid_height = 700;

    //the size of pixels sets the inital window height and width
    //don't make the pixels too large or the screen size will be larger than
    //your display size
    pixel_size = 1;

    /*Window information*/
    win_height = grid_height*pixel_size;
    win_width = grid_width*pixel_size;
    std::cout << "Please type the name of the file you want to read. :)" << std::endl;
    std::cin >> input;
    std::ifstream infile;
    infile.open(input);
    if (!infile) {
        std::cout << "unable to open file test_scene" << std::endl;
        exit(1);
    }
    infile >> numOfPoly;
    std::vector<double> coordin;
    int edgesNum;
    P1.setPoly(numOfPoly);
    for (int i = 0; i < numOfPoly; i++) {
        infile >> numOfPoints;
        Points P(numOfPoints);
        double a, b, c;
        for (int j = 0; j < numOfPoints; j++) {
            infile >> a >> b >> c;
            coordin.push_back(round(a * 700));
            coordin.push_back(round(b * 700));
            coordin.push_back(round(c * 700));
            P.setCoordinates(coordin);
            coordin.clear();
        }
        infile >> edgesNum;
        P.setNumOfEdges(edgesNum);
        std::map<int, Edge> storeEdge;
        for (int k = 0; k < P.getEdgesNum(); k++) {
            infile >> e1 >> e2;
            Edge readEdge;
            readEdge.setEdge(e1, e2);
            P.addEdge(readEdge);
        }
        P1.setList(P);
    }




    infile.close();
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    /*initialize variables, allocate memory, create buffers, etc. */
    //create window of size (win_width x win_height
    glutInitWindowSize(win_width, win_height);
    //windown title is "glut demo"
    glutCreateWindow("Project 2");
    glutDisplayFunc(displayMain);
    glutReshapeFunc(reshape); //update GL on window size change
    glutMouseFunc(mouse);     //mouse button events
    glutMotionFunc(motion);   //mouse movement events
    glutKeyboardFunc(key);    //Keyboard events
    glutIdleFunc(idle);

    //initialize opengl variables
    init();
    //start glut event loop
    glutMainLoop();
    return 0;
}

std::vector<double> centroid(Points poly) {
    double cx = 0;
    double cy = 0;
    double cz = 0;
    for (int i = 0; i < poly.getNumOfPoints(); i++) {
        cx += poly.getPoint(i)[0] * 1.0 / poly.getNumOfPoints();
        cy += poly.getPoint(i)[1] * 1.0 / poly.getNumOfPoints();
        cz += poly.getPoint(i)[2] * 1.0 / poly.getNumOfPoints();
    }
    std::vector<double> center;
    center.push_back(cx);
    center.push_back(cy);
    center.push_back(cz);
    return center;
}

void scaling(std::vector<double> factor, int id) {
    std::vector<double> center = centroid(P1.getPoly(id - 1));
    for (int i = 0; i < P1.getPoly(id - 1).getNumOfPoints(); i++) {
        P1.getPoly(id - 1).getPoint(i)[0] = factor[0] * P1.getPoly(id - 1).getPoint(i)[0] + center[0]
                    - center[0] * factor[0];
        P1.getPoly(id - 1).getPoint(i)[1] = factor[1] * P1.getPoly(id - 1).getPoint(i)[1] + center[1]
                    - center[1] * factor[1];
        P1.getPoly(id - 1).getPoint(i)[2] = factor[2] * P1.getPoly(id - 1).getPoint(i)[2] + center[2]
                    - center[2] * factor[2];
    }
}

void idle()
{
    //displayDDA();

    glutPostRedisplay();
}

void displayMain() {
    glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
    //clears the opengl Modelview transformation matrix
    glLoadIdentity();
    //checks for opengl errors
    displayXY();
    displayXZ();
    displayYZ();
    glutSwapBuffers();
    check();
}

void rotation(double angle, int id, std::vector<double> p0, std::vector<double> p) {
    double radian = angle * M_PI / 180.0;
    std::vector<double> axis;
    double del_x = p[0] - p0[0];
    double del_y = p[1] - p0[1];
    double del_z = p[2] - p0[2];
    double uvX = del_x / sqrt(del_x*del_x + del_y*del_y + del_z*del_z);
    double uvY = del_y / sqrt(del_x*del_x + del_y*del_y + del_z*del_z);
    double uvZ = del_z / sqrt(del_x*del_x + del_y*del_y + del_z*del_z);
    axis.push_back(uvX);
    axis.push_back(uvY);
    axis.push_back(uvZ);
    double ux = axis[0];
    double uy = axis[1];
    double uz = axis[2];
    for (int i = 0; i < P1.getPoly(id - 1).getNumOfPoints(); i++) {
        double x = P1.getPoly(id - 1).getPoint(i)[0];
        double y = P1.getPoly(id - 1).getPoint(i)[1];
        double z = P1.getPoly(id - 1).getPoint(i)[2];
        double tempX = (cos(radian) + ux*ux*(1 - cos(radian)))*x + (ux*uy*(1 - cos(radian)) - uz*sin(radian))*y
                + (ux*uz*(1 - cos(radian)) + uy*sin(radian))*z;
        double tempY = (uy*ux*(1 - cos(radian)) + uz*sin(radian))*x + (cos(radian) + uy*uy*(1 - cos(radian)))*y
                + (uy*uz*(1 - cos(radian)) - ux*sin(radian))*z;
        double tempZ = (uz*ux*(1 - cos(radian)) - uy*sin(radian))*x + (uz*uy*(1 - cos(radian)) + ux*sin(radian))*y
                + (cos(radian) + uz*uz*(1 - cos(radian)))*z;
        P1.getPoly(id - 1).getPoint(i)[0] = tempX;
        P1.getPoly(id - 1).getPoint(i)[1] = tempY;
        P1.getPoly(id - 1).getPoint(i)[2] = tempZ;
    }
}

void translation(std::vector<double> movement, int id) {
    double mx = movement[0];
    double my = movement[1];
    double mz = movement[2];
    for (int i = 0; i < P1.getPoly(id - 1).getNumOfPoints(); i++) {
        double &x = P1.getPoly(id - 1).getPoint(i)[0];
        double &y = P1.getPoly(id - 1).getPoint(i)[1];
        double &z = P1.getPoly(id - 1).getPoint(i)[2];
        x += mx;
        y += my;
        z += mz;
    }
}

void key(unsigned char ch, int x, int y)
{
    if (ch == 's') {
        key1 = true;
    } else if (ch == 'r') {
        key2 = true;
    } else if (ch == 't') {
        key3 = true;
    } else if (ch == 'e') {
        std::ofstream writefile;
        writefile.open(input);
        writefile << P1.getNumOfPolys() << std::endl << std::endl;
        for (int i = 0; i < P1.getNumOfPolys(); i++) {
            writefile << P1.getPoly(i).getNumOfPoints() << std::endl;
            for (int j = 0; j < P1.getPoly(i).getNumOfPoints(); j++) {
                writefile << P1.getPoly(i).getPoint(j)[0] << " " << P1.getPoly(i).getPoint(j)[1]
                    << " " << P1.getPoly(i).getPoint(j)[2] << std::endl;
            }
            writefile << P1.getPoly(i).getEdgesNum() << std::endl;
            for (int k = 0; k < P1.getPoly(i).getEdgesNum(); k++) {
                writefile << P1.getPoly(i).getEdges(k).getStartP() << " " << P1.getPoly(i).getEdges(k).getEndP() << std::endl;
            }
            writefile << std::endl << std::endl;
        }
        exit(1);
    }
    //redraw the scene after keyboard input
    glutPostRedisplay();
}

void motion(int x, int y)
{
    //redraw the scene after mouse movement
    glutPostRedisplay();
}

void mouse(int button, int state, int x, int y)
{
    //print the pixel location, and the grid location
    printf ("MOUSE AT PIXEL: %d %d, GRID: %d %d\n",x,y,(int)(x/pixel_size),(int)((win_height-y)/pixel_size));
//    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN && mouseInput) {
//        mouseInput1 = true;
//        mousefunc(state, x, y);
//        mouseInput = false;
//    }

    //redraw the scene after mouse click
    glutPostRedisplay();
}

void reshape(int width, int height)
{
    /*set up projection matrix to define the view port*/
    //update the ne window width and height
    win_width = width;
    win_height = height;

    //creates a rendering area across the window
    glViewport(0,0,width,height);
    // up an orthogonal projection matrix so that
    // the pixel space is mapped to the grid space
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0,grid_width,0,grid_height,-10,10);

    //clear the modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    //set pixel size based on width, if the aspect ratio
    //changes this hack won't work as well
    pixel_size = width/(float)grid_width;

    //set pixel size relative to the grid cell size
    glPointSize(pixel_size);
    //check for opengl errors
    check();
}

void displayXY() {
//    glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
//    //clears the opengl Modelview transformation matrix
//    glLoadIdentity();
    lineConnectionXY();
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex2f(0, grid_height);
    glVertex2f(grid_width, grid_height);
    glEnd();

    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex2f(grid_width, grid_height);
    glVertex2f(grid_width, 0);
    glEnd();
    for (int i = 0; i < P1.getNumOfPolys(); i++) {
        cList.push_back(centroid(P1.getPoly(i)));
    }

    std::vector<double> factor1;
    std::vector<double> unitV;
    if (fresh % 10 == 0) {
        if (key1) {
            double factor;
            int id;
            std::cout << "Please type scaling Factor and id of the polygon in \"factor id\" form" << std::endl;
            std::cin >> factor >> id;
            factor1.push_back(factor);
            factor1.push_back(factor);
            factor1.push_back(factor);
            scaling(factor1, id);
            key1 = false;
        } else if (key2) {
            double degree, x0, y0, z0, x, y, z;
            int id;
            std::cout << "Please type the angle, two points in the axis, and the id of the polygon in \"degree x0 y0 z0 x y z id\" form" << std::endl;
            std::cin >> degree >> x0 >> y0 >> z0 >> x >> y >> z >> id;
            std::vector<double> p0;
            std::vector<double> p;
            p0.push_back(x0);
            p0.push_back(y0);
            p0.push_back(z0);
            p.push_back(x);
            p.push_back(y);
            p.push_back(z);
            rotation(degree, id, p0, p);
            key2 = false;
        } else if (key3) {
            double x, y, z;
            int id;
            std::vector<double> movement;
            std::cout << "Please type the displacement vector and id of the polygon in the form \"x y z id\"." << std::endl;
            std::cin >> x >> y >> z >> id;
            movement.push_back(x);
            movement.push_back(y);
            movement.push_back(z);
            translation(movement, id);
            key3 = false;
        }
    }
    fresh++;
//    glutSwapBuffers();
    //checks for opengl errors
//    check();
}

void displayYZ() {
//    glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
//    //clears the opengl Modelview transformation matrix
//    glLoadIdentity();
    lineConnectionYZ();
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex2f(0, 0);
    glVertex2f(grid_width, 0);
    glEnd();

    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex2f(grid_width, grid_height);
    glVertex2f(grid_width, 0);
    glEnd();
    for (int i = 0; i < P1.getNumOfPolys(); i++) {
        cList.push_back(centroid(P1.getPoly(i)));
    }

    std::vector<double> factor1;
//    std::vector<double> unitV;
    if (fresh % 10 == 0) {
        if (key1) {
            double factor;
            int id;
            std::cout << "Please type scaling Factor and id of the polygon in \"x y z id\" form" << std::endl;
            std::cin >> factor >> id;
            factor1.push_back(factor);
            factor1.push_back(factor);
            factor1.push_back(factor);
            scaling(factor1, id);
            key1 = false;
        } else if (key2) {
            double degree, x0, y0, z0, x, y, z;
            int id;
            std::cout << "Please type the angle, two points in the axis, and the id of the polygon in \"degree x0 y0 z0 x y z id\" form" << std::endl;
            std::cin >> degree >> x0 >> y0 >> z0 >> x >> y >> z >> id;
            std::vector<double> p0;
            std::vector<double> p;
            p0.push_back(x0);
            p0.push_back(y0);
            p0.push_back(z0);
            p.push_back(x);
            p.push_back(y);
            p.push_back(z);

            rotation(degree, id, p0, p);
            key2 = false;
        } else if (key3) {
            double x, y, z;
            int id;
            std::vector<double> movement;
            std::cout << "Please type the displacement vector and id of the polygon in the form \"x y z id\"." << std::endl;
            std::cin >> x >> y >> z >> id;
            movement.push_back(x);
            movement.push_back(y);
            movement.push_back(z);
            translation(movement, id);
            key3 = false;
        }
    }
    fresh++;
//    glutSwapBuffers();
//    //checks for opengl errors
//    check();
}

void displayXZ() {
//    glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
//    //clears the opengl Modelview transformation matrix
//    glLoadIdentity();
    lineConnectionXZ();
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex2f(0, 0);
    glVertex2f(grid_width, 0);
    glEnd();
    for (int i = 0; i < P1.getNumOfPolys(); i++) {
        cList.push_back(centroid(P1.getPoly(i)));
    }

    std::vector<double> factor1;
//    std::vector<double> unitV;
    if (fresh % 10 == 0) {
        if (key1) {
            double factor;
            int id;
            std::cout << "Please type scaling Factor and id of the polygon in \"x y z id\" form" << std::endl;
            std::cin >> factor >> id;
            factor1.push_back(factor);
            factor1.push_back(factor);
            factor1.push_back(factor);
            scaling(factor1, id);
            key1 = false;
        } else if (key2) {
            double degree, x0, y0, z0, x, y, z;
            int id;
            std::cout << "Please type the angle, two points in the axis, and the id of the polygon in \"degree x0 y0 z0 x y z id\" form" << std::endl;
            std::cin >> degree >> x0 >> y0 >> z0 >> x >> y >> z >> id;
            std::vector<double> p0;
            std::vector<double> p;
            p0.push_back(x0);
            p0.push_back(y0);
            p0.push_back(z0);
            p.push_back(x);
            p.push_back(y);
            p.push_back(z);

            rotation(degree, id, p, p0);
            key2 = false;
        } else if (key3) {
            double x, y, z;
            int id;
            std::vector<double> movement;
            std::cout << "Please type the displacement vector and id of the polygon in the form \"x y z id\"." << std::endl;
            std::cin >> x >> y >> z >> id;
            movement.push_back(x);
            movement.push_back(y);
            movement.push_back(z);
            translation(movement, id);
            key3 = false;
        }
    }
    fresh++;
//    glutSwapBuffers();
//    //checks for opengl errors
//    check();
}

void init()
{
    //set clear color (Default background to white)
    glClearColor(1.0,1.0,1.0,1.0);
    //checks for OpenGL errors
    check();
}

void check()
{
    GLenum err = glGetError();
    if(err != GL_NO_ERROR)
    {
        printf("GLERROR: There was an error %s\n",gluErrorString(err) );
        exit(1);
    }
}

void lineConnectionXY() {
    glViewport(0, 0, win_width / 2, win_height / 2);
    for (int i = 0; i < P1.getNumOfPolys(); i++) {
        for (int j = 0; j < P1.getPoly(i).getEdgesNum(); j++) {
            glBegin(GL_LINES);
            glColor3f(0.7, 0.4, 0.2);
            glVertex2f(P1.getPoly(i).getPoint(P1.getPoly(i).getEdges(j).getStartP() - 1)[0],
                    P1.getPoly(i).getPoint(P1.getPoly(i).getEdges(j).getStartP() - 1)[1]);
            glVertex2f(P1.getPoly(i).getPoint(P1.getPoly(i).getEdges(j).getEndP() - 1)[0],
                       P1.getPoly(i).getPoint(P1.getPoly(i).getEdges(j).getEndP() - 1)[1]);
            glEnd();
        }
    }
}

void lineConnectionYZ() {
    glViewport(0, win_height / 2, win_width / 2, win_height / 2);
    for (int i = 0; i < P1.getNumOfPolys(); i++) {
        for (int j = 0; j < P1.getPoly(i).getEdgesNum(); j++) {
            glBegin(GL_LINES);
            glColor3f(0.7, 0.4, 0.2);
            glVertex2f(P1.getPoly(i).getPoint(P1.getPoly(i).getEdges(j).getStartP() - 1)[1],
                       P1.getPoly(i).getPoint(P1.getPoly(i).getEdges(j).getStartP() - 1)[2]);
            glVertex2f(P1.getPoly(i).getPoint(P1.getPoly(i).getEdges(j).getEndP() - 1)[1],
                       P1.getPoly(i).getPoint(P1.getPoly(i).getEdges(j).getEndP() - 1)[2]);
            glEnd();
        }
    }
}

void lineConnectionXZ() {
    glViewport(win_width / 2, win_height / 2, win_width / 2, win_height / 2);
    for (int i = 0; i < P1.getNumOfPolys(); i++) {
        for (int j = 0; j < P1.getPoly(i).getEdgesNum(); j++) {
            glBegin(GL_LINES);
            glColor3f(0.7, 0.4, 0.2);
            glVertex2f(P1.getPoly(i).getPoint(P1.getPoly(i).getEdges(j).getStartP() - 1)[0],
                       P1.getPoly(i).getPoint(P1.getPoly(i).getEdges(j).getStartP() - 1)[2]);
            glVertex2f(P1.getPoly(i).getPoint(P1.getPoly(i).getEdges(j).getEndP() - 1)[0],
                       P1.getPoly(i).getPoint(P1.getPoly(i).getEdges(j).getEndP() - 1)[2]);
            glEnd();
        }
    }
}
