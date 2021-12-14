/*  */
#define GLEW_STATIC
#define STB_IMAGE_IMPLEMENTATION

#include <iostream>
#include <random>
#include <cmath>
#include <windows.h>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

using namespace std;

void squareShapeDataSet(double** data_vertices, int size) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> disx(-0.9, 0.9);
    for (int i = 0; i < size; i++) {
        data_vertices[i] = (double*)malloc(2 * sizeof(double));
        data_vertices[i][0] = disx(gen);
        data_vertices[i][1] = disx(gen);
    }
}

void circleShapeDataSet(double** data_vertices, int size) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> disx(-0.9, 0.9);
    for (int i = 0; i < size; i++) {
        data_vertices[i] = (double*)malloc(2 * sizeof(double));
        double x = disx(gen);
        double y = disx(gen);
        while (x * x + y * y > 0.81) {
            x = disx(gen);
            y = disx(gen);
        }
        data_vertices[i][0] = x;
        data_vertices[i][1] = y;
    }
}

void ringsShapeDataSet(double** data_vertices, int size, int rings) {
    if (rings <= 0)rings = (int)sqrt(size);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> disx(-0.9, 0.9);
    double interval = 0.9 / rings;
    int point = 0;
    for (int i = 1; i <= rings; i++) {
        int num = i * size / ((1 + rings) * rings / 2);
        if (num < 3)num=3;
        if (i == rings)num = size - point;
        double r = i * 0.9 / rings;
        double angle = 6.28 / num;
        for (int j = 0; j < num; j++) {
            data_vertices[point] = (double*)malloc(2 * sizeof(double));
            data_vertices[point][0] = r * cos(j * angle);
            data_vertices[point][1] = r * sin(j * angle);
            point++;
        }
    }
}

void generatePointSet(double** data_vertices, int size, int shape, int rings = 0) {
    if (shape == 0)squareShapeDataSet(data_vertices, size);
    if (shape == 1)circleShapeDataSet(data_vertices, size);
    if (shape == 2)ringsShapeDataSet(data_vertices, size, rings);
}

//because P1, P2 and P3 are in xy plane, the cross product of vector(P1 P2) and vector(P1 P3) will on z-axis, thus the function only return the z value.
double crossProduct(double* point1, double* point2, double* point3) {
    return (point2[0] - point1[0]) * (point3[1] - point1[1]) - (point2[1] - point1[1]) * (point3[0] - point1[0]);
}

//test if the segment P1-P2 and P3-P4 are intersect
bool isIntersect(double* point1, double* point2, double* point3, double* point4) {
    if (crossProduct(point1, point2, point3) * crossProduct(point1, point2, point4) < 0 &&
        crossProduct(point3, point4, point1) * crossProduct(point3, point4, point2) < 0) {
        return true;
    }
    return false;
}

//return the x and y position of the intersect point of line (P1 P2) and line (P3 P4)
double* twoLineIntersect(double* point1, double* point2, double* point3, double* point4) {
    double* res = (double*)malloc(2 * sizeof(double));
    if (point1[0] == point2[0])point1[0] += 0.0000001;
    if (point3[0] == point4[0])point3[0] += 0.0000001;
    double a1 = (point1[1] - point2[1]) / (point1[0] - point2[0]);
    double b1 = (point2[1] * point1[0] - point1[1] * point2[0]) / (point1[0] - point2[0]);
    double a2 = (point3[1] - point4[1]) / (point3[0] - point4[0]);
    double b2 = (point4[1] * point3[0] - point3[1] * point4[0]) / (point3[0] - point4[0]);
    res[0] = (b2 - b1) / (a1 - a2);
    res[1] = (a1 * b2 - a2 * b1) / (a1 - a2);
    return res;
}

double cosine(double* point1, double* point2, double* point3, double* point4) {
    double x1 = point2[0] - point1[0];
    double x2 = point4[0] - point3[0];
    double y1 = point2[1] - point1[1];
    double y2 = point4[1] - point3[1];
    return (x1 * x2 + y1 * y2) / (sqrt(x1 * x1 + y1 * y1) * sqrt(x2 * x2 + y2 * y2));
}

/*
* This function will reorder the input point set, the result is the spiral line.
* Parameter:
* data_vertices: the point set data, the shape is size*2, data_vertices[i][0] is the x coordinate of i th point, data_vertices[i][1] is the y coordinate.
* size: means how many points in the point set.
* h: the index from 0 to h is the convex hull of the point set.
* y: the index for the inner starshape, please see document for detail.
*/
void convexLayer(double** data_vertices, int size, int& y, int& h) {
    int min = 0;
    for (int i = 1; i < size; i++) {
        if (data_vertices[i][0] < data_vertices[min][0] || (data_vertices[i][0] == data_vertices[min][0] && data_vertices[i][1] < data_vertices[min][1])) {
            min = i;
        }
    }
    double* t = data_vertices[min];
    data_vertices[min] = data_vertices[0];
    data_vertices[0] = t;
    for (int i = 1; i < size; i++) {
        int temp = i;
        for (int j = i + 1; j < size; j++) {
            if (crossProduct(data_vertices[i - 1], data_vertices[temp], data_vertices[j]) > 0) {
                temp = j;
            }
        }
        double* t = data_vertices[temp];
        data_vertices[temp] = data_vertices[i];
        data_vertices[i] = t;
    }
    y = size - 3;
    while (crossProduct(data_vertices[size - 2], data_vertices[size - 1], data_vertices[y]) < 0)y--;
    h = 1;
    for (int i = 1; i < size; i++) {
        if (crossProduct(data_vertices[0], data_vertices[h], data_vertices[i]) < 0)h = i;
    }
}

/*
* This function will do the path triangulation of the input spiral line, and save the result in the 2D array triangles by each points' index.
* Parameter:
* data_vertices, size, y, h: data get form the function convexLayer
*/
void pathTriangulation(double **data_vertices, int size, int &y, int &h, int* triangles) {
    int outer = 0;
    int inner = h;
    int tri = 0;
    while (outer < y || inner < size - 1) {
        if (inner < size - 1 &&
            cosine(data_vertices[inner], data_vertices[outer], data_vertices[inner], data_vertices[inner + 1]) > cosine(data_vertices[inner], data_vertices[outer], data_vertices[outer], data_vertices[outer + 1]) &&
            crossProduct(data_vertices[inner], data_vertices[outer], data_vertices[inner + 1]) < 0 ) {
            triangles[3 * tri] = outer;
            triangles[3 * tri + 1] = inner;
            triangles[3 * tri + 2] = inner + 1;
            inner++;
        }
        else {
            triangles[3 * tri] = outer;
            triangles[3 * tri + 1] = inner;
            triangles[3 * tri + 2] = outer + 1;
            outer++;
        }
        tri++;
    }
    for (int i = y; i < size - 2; i++) {
        triangles[3 * tri] = i;
        triangles[3 * tri + 1] = size - 1;
        triangles[3 * tri + 2] = i + 1;
        tri++;
    }
    cout <<"tri: "<< tri << " " << 2*size-h-3;
}

void generateHexagons(int* triangles, int triangles_num, vector<int>& hexagons) {
    for (int i = triangles_num - 4; i >= 0; i -= 4) {
        int outer0 = triangles[3 * i + 0];
        int inner0 = triangles[3 * i + 1];
        int outer1 = 0, inner1 = 0;
        if (triangles[3 * i + 11] == triangles[3 * i + 10] + 1) {
            outer1 = triangles[3 * i + 9];
            inner1 = triangles[3 * i + 11];
        }
        else {
            outer1 = triangles[3 * i + 11];
            inner1 = triangles[3 * i + 10];
        }
        for (int j = inner0; j <= inner1; j++) {
            hexagons.push_back(j);
        }
        for (int k = outer1; k >= outer0; k--) {
            hexagons.push_back(k);
        }
    }
}

void quadStarshape(int* hex, int& offset, vector<int>& quads, int steiner) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            quads.push_back(hex[(j + 2 * i + offset) % 6]);
        }
        quads.push_back(steiner);
    }
}

void quadSteiner(int* hex, int& offset, vector<int>&quads, int steiner, int start) {
    for (int i = start; i < start + 3; i++) {
        quads.push_back(hex[(i + offset) % 6]);
    }
    quads.push_back(steiner);
}

void quadCCCCCC(int* hex, double** hexagon, int& offset, vector<int>& quads, vector<double>& data_vertices) {
    cout << "quad CCCCCC" << endl;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 4; j++) {
            quads.push_back(hex[(j + 3 * i + offset) % 6]);
        }
    }
}

void quadRCCCCC(int* hex, double** hexagon, int& offset, vector<int>& quads, vector<double>& data_vertices) {
    //verify if d is in wedge(a)
    if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 5) % 6]) > 0 &&
        crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 3) % 6]) > 0)
    {
        //this is the RCCCCC case 1
        cout << "quad RCCCCC case 1" << endl;
        quadCCCCCC(hex, hexagon, offset, quads, data_vertices);
    }
    //if c-e intersect with a-b, c must not see e
    else if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 4) % 6])>0) {
        //this is the RCCCCC case 2, c sees e
        cout << "quad RCCCCC case 2" << endl;
        double* intersect1 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 4) % 6]);
        if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 4) % 6]) > 0) {
            intersect1 = hexagon[(offset + 4) % 6];
        }
        double* intersect2 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 4) % 6]);
        if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 5) % 6]) > 0) {
            intersect2 = hexagon[(offset + 2) % 6];
        }
        int s = data_vertices.size()/2;
        double steiner[] = { (intersect1[0] + intersect2[0] + hexagon[(offset + 0) % 6][0]) / 3 ,(intersect1[1] + intersect2[1] + hexagon[(offset + 0) % 6][1]) / 3 };
        data_vertices.push_back(steiner[0]);
        data_vertices.push_back(steiner[1]);
        quadStarshape(hex, offset, quads, s);
    }
    else {
        //the c could not see e case
        if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 5) % 6]) <= 0) {
            //if d lies L(a,b) RCCCCC case 3
            cout << "quad RCCCCC case 3" << endl;
            double* intersect1 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6]);
            double* intersect2 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6]);
            int s = data_vertices.size() / 2;
            double steiner[] = { (intersect1[0] + intersect2[0] + hexagon[(offset + 0) % 6][0]) / 3 ,(intersect1[1] + intersect2[1] + hexagon[(offset + 0) % 6][1]) / 3 };
            data_vertices.push_back(steiner[0]);
            data_vertices.push_back(steiner[1]);
            quadSteiner(hex, offset, quads, s, 4);
            offset = (offset + 5) % 6;
            hexagon[(offset + 0) % 6][0] = steiner[0];
            hexagon[(offset + 0) % 6][1] = steiner[1];
            hex[(offset + 0) % 6] = s;
            quadRCCCCC(hex, hexagon, offset, quads, data_vertices);
        }
        else {
            //if d lies R(a,b) RCCCCC case 4
            cout << "quad RCCCCC case 4" << endl;
            double* intersect1 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 3) % 6]);
            double* intersect2 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 3) % 6]);
            int s = data_vertices.size() / 2;
            double steiner[] = { (intersect1[0] + intersect2[0] + hexagon[(offset + 0) % 6][0]) / 3,(intersect1[1] + intersect2[1] + hexagon[(offset + 0) % 6][1]) / 3 };
            data_vertices.push_back(steiner[0]);
            data_vertices.push_back(steiner[1]);
            quadSteiner(hex, offset, quads, s, 0);
            offset = (offset + 1) % 6;
            hexagon[(offset + 0) % 6][0] = steiner[0];
            hexagon[(offset + 0) % 6][1] = steiner[1];
            hex[(offset + 0) % 6] = s;
            quadRCCCCC(hex, hexagon, offset, quads, data_vertices);
        }
    }
}

void quadRCRCCC(int* hex, double** hexagon, int& offset, vector<int>& quads, vector<double>& data_vertices) {
    if (!isIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 3) % 6]) &&
        !isIntersect(hexagon[(offset + 2) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6])) {
        //a and c can see e, the starshaped case RCRCCC case 1
        cout << "quad RCRCCC case 1" << endl;
        double* triangleTop = hexagon[(offset + 4) % 6];
        double* triangleLeft = hexagon[(offset + 0) % 6];
        double* triangleRight = hexagon[(offset + 2) % 6];
        if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], hexagon[(offset + 2) % 6]) > 0) {
            triangleRight = twoLineIntersect(hexagon[(offset + 5) % 6], hexagon[(offset + 0) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 4) % 6]);
        }
        if (crossProduct(hexagon[(offset + 2) % 6], hexagon[(offset + 0) % 6], hexagon[(offset + 3) % 6]) > 0) {
            triangleLeft= twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 3) % 6]);
        }
        if (crossProduct(hexagon[(offset + 2) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 4) % 6]) > 0) {
            triangleTop= twoLineIntersect(hexagon[(offset + 1) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 0) % 6], hexagon[(offset + 4) % 6]);
            triangleRight = twoLineIntersect(hexagon[(offset + 1) % 6], hexagon[(offset + 2) % 6], triangleLeft, triangleRight);
        }
        if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 1) % 6]) > 0) {
            triangleTop = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 4) % 6]);
            triangleLeft = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], triangleLeft, triangleRight);
        }
        int s = data_vertices.size() / 2;
        double steiner[] = { (triangleLeft[0] + triangleRight[0] + triangleTop[0]) / 3,(triangleLeft[1] + triangleRight[1] + triangleTop[1]) / 3 };
        data_vertices.push_back(steiner[0]);
        data_vertices.push_back(steiner[1]);
        quadStarshape(hex, offset, quads, s);
    }
    else {
        double steiner[] = { 0.0,0.0 };
        if (isIntersect(hexagon[(offset + 2) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6])) {
            //RCRCCC case 2
            cout << "quad RCRCCC case 2" << endl;
            double* midPoint1 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6]);
            double* midPoint2 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6]);          
            if (crossProduct(hexagon[(offset + 2) % 6], hexagon[(offset + 1) % 6], midPoint2) > 0) {
                midPoint2 = twoLineIntersect(hexagon[(offset + 1) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6]);
            }
            steiner[0] = (midPoint1[0] + midPoint2[0] + hexagon[(offset + 0) % 6][0]) / 3;
            steiner[1] = (midPoint1[1] + midPoint2[1] + hexagon[(offset + 0) % 6][1]) / 3;
        }
        else {
            //RCRCCC case 3
            cout << "quad RCRCCC case 3" << endl;
            double* midPoint1 = twoLineIntersect(hexagon[(offset + 2) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 5) % 6]);
            double* midPoint2 = twoLineIntersect(hexagon[(offset + 2) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 5) % 6]);
            if (crossProduct(hexagon[(offset + 0) % 6], midPoint2, hexagon[(offset + 1) % 6]) > 0) {
                midPoint2= twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 5) % 6]);
            }
            steiner[0] = (midPoint1[0] + midPoint2[0] + hexagon[(offset + 2) % 6][0]) / 3;
            steiner[1] = (midPoint1[1] + midPoint2[1] + hexagon[(offset + 2) % 6][1]) / 3;
        }
        int s = data_vertices.size() / 2;
        quadSteiner(hex, offset, quads, s, 0);
        data_vertices.push_back(steiner[0]);
        data_vertices.push_back(steiner[1]);
        offset = (offset + 1) % 6;
        hexagon[(offset + 0) % 6][0] = steiner[0];
        hexagon[(offset + 0) % 6][1] = steiner[1];
        hex[(offset + 0) % 6] = s;
        quadRCCCCC(hex, hexagon, offset, quads, data_vertices);
    }
}

void quadRRCCCC(int* hex, double** hexagon, int& offset, vector<int>& quads, vector<double>& data_vertices) {
    cout << "quad RRCCCC" << endl;
    double* intersect0 = hexagon[(offset + 0) % 6];
    if (crossProduct(hexagon[(offset + 1) % 6], hexagon[(offset + 0) % 6], hexagon[(offset + 3) % 6]) > 0) {
        intersect0 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 3) % 6]);
    }
    double* intersect1 = hexagon[(offset + 4) % 6];
    if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 1) % 6])>0) {
        intersect1 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6]);
    }
    double* intersect2 = hexagon[(offset + 3) % 6];
    if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 5) % 6]) > 0) {
        intersect2 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 3) % 6]);
    }else intersect2 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6]);
    int s = data_vertices.size() / 2;
    double steiner[] = { (intersect0[0] + intersect1[0] + intersect2[0]) / 3,(intersect0[1] + intersect1[1] + intersect2[1]) / 3 };
    data_vertices.push_back(steiner[0]);
    data_vertices.push_back(steiner[1]);
    quadSteiner(hex, offset, quads, s, 4);
    offset = (offset + 5) % 6;
    hexagon[(offset + 0) % 6][0] = steiner[0];
    hexagon[(offset + 0) % 6][1] = steiner[1];
    hex[(offset + 0) % 6] = s;
    quadRCRCCC(hex, hexagon, offset, quads, data_vertices);
}

void quadRCCRCC(int* hex, double** hexagon, int& offset, vector<int>& quads, vector<double>& data_vertices) {
    if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 5) % 6]) > 0 &&
        crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 3) % 6]) > 0) {
        //RCCRCC case 1
        cout << "quad RCCRCC case 1" << endl;
        quadCCCCCC(hex, hexagon, offset, quads, data_vertices);
    }
    else if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 3) % 6]) < 0) {
        //RCCRCC case 2
        cout << "quad RCCRCC case 2" << endl;
        double steiner[] = { 0.0,0.0 };
        if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 5) % 6]) > 0) {
            double * intersect1= twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6]);
            steiner[0] = (hexagon[(offset + 3) % 6][0] + hexagon[(offset + 5) % 6][0] + intersect1[0]) / 3;
            steiner[1] = (hexagon[(offset + 3) % 6][1] + hexagon[(offset + 5) % 6][1] + intersect1[1]) / 3;
        }
        else {
            double* intersect1 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 5) % 6]);
            double* intersect2 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6]);
            steiner[0] = (hexagon[(offset + 3) % 6][0] + intersect1[0] + intersect2[0]) / 3;
            steiner[1] = (hexagon[(offset + 3) % 6][1] + intersect1[1] + intersect2[1]) / 3;
        }
        int s = data_vertices.size() / 2;
        quadSteiner(hex, offset, quads, s, 3);
        data_vertices.push_back(steiner[0]);
        data_vertices.push_back(steiner[1]);
        offset = (offset + 4) % 6;
        hexagon[(offset + 0) % 6][0] = steiner[0];
        hexagon[(offset + 0) % 6][1] = steiner[1];
        hex[(offset + 0) % 6] = s;
        quadRCRCCC(hex, hexagon, offset, quads, data_vertices);
    }
    else {
        //RCCRCC case 3
        cout << "quad RCCRCC case 3" << endl;
        double steiner[] = { 0.0,0.0 };
        if (crossProduct(hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 1) % 6]) > 0) {
            double* intersect1 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6]);
            steiner[0] = (hexagon[(offset + 0) % 6][0] + hexagon[(offset + 4) % 6][0] + intersect1[0]) / 3;
            steiner[1] = (hexagon[(offset + 0) % 6][1] + hexagon[(offset + 4) % 6][1] + intersect1[1]) / 3;
        }
        else {
            double* intersect1 = twoLineIntersect(hexagon[(offset + 1) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6]);
            double* intersect2 = twoLineIntersect(hexagon[(offset + 1) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 0) % 6], hexagon[(offset + 4) % 6]);
            steiner[0] = (hexagon[(offset + 3) % 6][0] + intersect1[0] + intersect2[0]) / 3;
            steiner[1] = (hexagon[(offset + 3) % 6][1] + intersect1[1] + intersect2[1]) / 3;
        }
        int s = data_vertices.size() / 2;
        data_vertices.push_back(steiner[0]);
        data_vertices.push_back(steiner[1]);
        quadSteiner(hex, offset, quads, s, 4);
        offset = (offset + 3) % 6;
        hexagon[(offset + 2) % 6][0] = steiner[0];
        hexagon[(offset + 2) % 6][1] = steiner[1];
        hex[(offset + 2) % 6] = s;
        quadRCRCCC(hex, hexagon, offset, quads, data_vertices);
    }
}

void quadRCRCRC(int* hex, double** hexagon, int& offset, vector<int>& quads, vector<double>& data_vertices) {
    bool isStarShape = true;
    for (int i = 0; i <= 4; i += 2) {
        if (isIntersect(hexagon[(offset + 0 + i) % 6], hexagon[(offset + 4 + i) % 6], hexagon[(offset + 2 + i) % 6], hexagon[(offset + 3 + i) % 6])) {
            isStarShape = false;
            offset = (offset + i) % 6;
        }
    }
    if (isStarShape) {
        //RCRCRC case 1
        cout << "quad RCRCRC case 1" << endl;
        double* triangle[] = { hexagon[(offset + 0) % 6],hexagon[(offset + 2) % 6],hexagon[(offset + 4) % 6] };
        for (int i = 0; i < 3; i++) {
            if (crossProduct(hexagon[(offset + 0 + 2 * i) % 6], hexagon[(offset + 5 + 2 * i) % 6], hexagon[(offset + 2 + 2 * i) % 6]) > 0) {
                triangle[(i + 1) % 3] = twoLineIntersect(hexagon[(offset + 0 + 2 * i) % 6], hexagon[(offset + 5 + 2 * i) % 6], triangle[(i + 1) % 3], triangle[(i + 2) % 3]);
                triangle[(i + 0) % 3] = twoLineIntersect(hexagon[(offset + 0 + 2 * i) % 6], hexagon[(offset + 5 + 2 * i) % 6], triangle[(i + 0) % 3], triangle[(i + 2) % 3]);
            }
            if (crossProduct(hexagon[(offset + 0 + 2 * i) % 6], hexagon[(offset + 4 + 2 * i) % 6], hexagon[(offset + 1 + 2 * i) % 6]) > 0) {
                triangle[(i + 2) % 3] = twoLineIntersect(hexagon[(offset + 0 + 2 * i) % 6], hexagon[(offset + 1 + 2 * i) % 6], triangle[(i + 1) % 3], triangle[(i + 2) % 3]);
                triangle[(i + 0) % 3] = twoLineIntersect(hexagon[(offset + 0 + 2 * i) % 6], hexagon[(offset + 1 + 2 * i) % 6], triangle[(i + 1) % 3], triangle[(i + 0) % 3]);
            }
        }
        int s = data_vertices.size() / 2;
        double steiner[] = { (triangle[0][0] + triangle[1][0] + triangle[2][0]) / 3, (triangle[0][1] + triangle[1][1] + triangle[2][1]) / 3 };
        data_vertices.push_back(steiner[0]);
        data_vertices.push_back(steiner[1]);
        quadStarshape(hex, offset, quads, s);
    }
    else {
        //RCRCRC case 2
        cout << "quad RCRCRC case 2" << endl;
        double* intersect1 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 5) % 6]);
        double* intersect2 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 3) % 6]);
        double steiner[] = { (intersect1[0] + intersect2[0] + hexagon[(offset + 2) % 6][0]) / 3,(intersect1[1] + intersect2[1] + hexagon[(offset + 2) % 6][1]) / 3 };
        int s = data_vertices.size() / 2;
        data_vertices.push_back(steiner[0]);
        data_vertices.push_back(steiner[1]);
        quadSteiner(hex, offset, quads, s, 0);
        offset = (offset + 1) % 6;
        hexagon[(offset + 0) % 6][0] = steiner[0];
        hexagon[(offset + 0) % 6][1] = steiner[1];
        hex[(offset + 0) % 6] = s;
        quadRCCRCC(hex, hexagon, offset, quads, data_vertices);
    }
}

void quadRCRRCC(int* hex, double** hexagon, int& offset, vector<int>& quads, vector<double>& data_vertices);

void quadRRCRCC(int* hex, double** hexagon, int& offset, vector<int>& quads, vector<double>& data_vertices) {
    if (crossProduct(hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 0) % 6]) > 0) {
        //RRCRCC case 1
        cout << "quad RRCRCC case 1" << endl;
        double* intersect1 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 3) % 6]);
        double* intersect2 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 1) % 6]);
        double* intersect3 = hexagon[(offset + 3) % 6];
        if (crossProduct(hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 1) % 6]) > 0) {
            intersect3 = twoLineIntersect(hexagon[(offset + 1) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 3) % 6]);
        }
        if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], intersect3) > 0) {
            intersect3 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], intersect1, intersect3);
        }
        double steiner[] = { (intersect1[0] + intersect2[0] + intersect3[0]) / 3,(intersect1[1] + intersect2[1] + intersect3[1]) / 3 };
        int s = data_vertices.size() / 2;
        data_vertices.push_back(steiner[0]);
        data_vertices.push_back(steiner[1]);
        quadSteiner(hex, offset, quads, s, 4);
        offset = (offset + 5) % 6;
        hexagon[(offset + 0) % 6][0] = steiner[0];
        hexagon[(offset + 0) % 6][1] = steiner[1];
        hex[(offset + 0) % 6] = s;
        quadRCRCRC(hex, hexagon, offset, quads, data_vertices);
    }
    else {
        //RRCRCC case 2
        cout << "quad RRCRCC case 2" << endl;
        double* intersect1 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6]);
        double steiner[] = { (hexagon[(offset + 3) % 6][0] + hexagon[(offset + 5) % 6][0] + intersect1[0]) / 3,
            (hexagon[(offset + 3) % 6][1] + hexagon[(offset + 5) % 6][1] + intersect1[1]) / 3 };
        int s = data_vertices.size() / 2;
        data_vertices.push_back(steiner[0]);
        data_vertices.push_back(steiner[1]);
        quadSteiner(hex, offset, quads, s, 3);
        offset = (offset + 4) % 6;
        hexagon[(offset + 0) % 6][0] = steiner[0];
        hexagon[(offset + 0) % 6][1] = steiner[1];
        hex[(offset + 0) % 6] = s;
        quadRCRRCC(hex, hexagon, offset, quads, data_vertices);
    }
}

void quadRCRRCC(int* hex, double** hexagon, int& offset, vector<int>& quads, vector<double>& data_vertices) {
    if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 5) % 6]) > 0) {
        //RRCRCC case 1
        cout << "quad RCRRCC case 1" << endl;
        double* intersect1 = twoLineIntersect(hexagon[(offset + 1) % 6], hexagon[(offset + 0) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 5) % 6]);
        double* intersect2 = twoLineIntersect(hexagon[(offset + 1) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 5) % 6]);
        double* intersect3 = hexagon[(offset + 0) % 6];
        if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 5) % 6]) > 0) {
            intersect3 = twoLineIntersect(hexagon[(offset + 1) % 6], hexagon[(offset + 0) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 5) % 6]);
        }
        if (crossProduct(hexagon[(offset + 3) % 6], intersect3, hexagon[(offset + 4) % 6]) > 0) {
            intersect3 = twoLineIntersect(hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6], intersect1, intersect3);
        }
        double steiner[] = { (intersect1[0] + intersect2[0] + intersect3[0]) / 3,(intersect1[1] + intersect2[1] + intersect3[1]) / 3 };
        int s = data_vertices.size() / 2;
        data_vertices.push_back(steiner[0]);
        data_vertices.push_back(steiner[1]);
        quadSteiner(hex, offset, quads, s, 3);
        offset = (offset + 4) % 6;
        hexagon[(offset + 0) % 6][0] = steiner[0];
        hexagon[(offset + 0) % 6][1] = steiner[1];
        hex[(offset + 0) % 6] = s;
        quadRCRCRC(hex, hexagon, offset, quads, data_vertices);
    }
    else {
        //RRCRCC case 2
        cout << "quad RCRRCC case 1" << endl;
        double* intersect1 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6]);
        double steiner[] = { (hexagon[(offset + 0) % 6][0] + hexagon[(offset + 4) % 6][0] + intersect1[0]) / 3,
                            (hexagon[(offset + 0) % 6][1] + hexagon[(offset + 4) % 6][1] + intersect1[1]) / 3 };
        int s = data_vertices.size() / 2;
        data_vertices.push_back(steiner[0]);
        data_vertices.push_back(steiner[1]);
        quadSteiner(hex, offset, quads, s, 4);
        offset = (offset + 2) % 6;
        hexagon[(offset + 3) % 6][0] = steiner[0];
        hexagon[(offset + 3) % 6][1] = steiner[1];
        hex[(offset + 3) % 6] = s;
        quadRRCRCC(hex, hexagon, offset, quads, data_vertices);
    }   
}

void quadRRRCCC(int* hex, double** hexagon, int& offset, vector<int>& quads, vector<double>data_vertices) {
    cout << "quad RRRCCC" << endl;
    double* intersect1 = twoLineIntersect(hexagon[(offset + 1) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6]);
    double steiner[] = { (hexagon[(offset + 0) % 6][0] + hexagon[(offset + 4) % 6][0] + intersect1[0]) / 3,
                        (hexagon[(offset + 0) % 6][1] + hexagon[(offset + 4) % 6][1] + intersect1[1]) / 3 };
    int s = data_vertices.size() / 2;
    data_vertices.push_back(steiner[0]);
    data_vertices.push_back(steiner[1]);
    quadSteiner(hex, offset, quads, s, 4);
    offset = (offset + 5) % 6;
    hexagon[(offset + 0) % 6][0] = steiner[0];
    hexagon[(offset + 0) % 6][1] = steiner[1];
    hex[(offset + 0) % 6] = s;
    quadRCRRCC(hex, hexagon, offset, quads, data_vertices);
}

/*
* This function will figure out the input hexagon is which kind of shape and call the corresponding quad function.
* Parameter:
* hexagon: the input 6*2 array, should be ordered in counter-clockwised
* quads: to receive result. each 8 double represent a quadrangle
*/
void quadHexagon(int* hex, vector<int>& quads, vector<double>& data_vertices) {
    double** hexagon = new double* [6];
    for (int i = 0; i < 6; i++) {
        hexagon[i] = new double[2];
        hexagon[i][0] = data_vertices[2 * hex[i]];
        hexagon[i][1] = data_vertices[2 * hex[i] + 1];
    }
    int angle[] = {0,0,0,0,0,0};
    for (int i = 0; i < 6; i++) {
        if (crossProduct(hexagon[i], hexagon[(i + 5) % 6], hexagon[(i + 1) % 6]) > 0) {
            angle[i] = 1;
        }
    }
    int offset = 0;
    int score = 64;
    for (int i = 0; i < 6; i++) {
        int temp = 0;
        int helper = 1;
        for (int j = 0; j < 6; j++) {
            temp += helper * angle[(i + j) % 6];
            helper *= 2;
        }
        if (temp < score) {
            score = temp;
            offset = i;
        }
    }
    cout << "===== score : " << score << endl;
    if (score == 0)  quadCCCCCC(hex, hexagon, offset, quads, data_vertices);
    if (score == 1)  quadRCCCCC(hex, hexagon, offset, quads, data_vertices);
    if (score == 5)  quadRCRCCC(hex, hexagon, offset, quads, data_vertices);
    if (score == 3)  quadRRCCCC(hex, hexagon, offset, quads, data_vertices);
    if (score == 9)  quadRCCRCC(hex, hexagon, offset, quads, data_vertices);
    if (score == 21) quadRCRCRC(hex, hexagon, offset, quads, data_vertices);
    if (score == 11) quadRRCRCC(hex, hexagon, offset, quads, data_vertices);
    if (score == 13) quadRCRRCC(hex, hexagon, offset, quads, data_vertices);
    if (score == 7)  quadRRRCCC(hex, hexagon, offset, quads, data_vertices);
}

/*
* This function will quadrangulate the input quadrangle with an inner point
*/
void quadQuad(int* q, int inner, vector<int>& quads, vector<double>& data_vertices) {
    double** quad = new double* [4];
    for (int i = 0; i < 4; i++) {
        quad[i] = new double[2];
        quad[i][0] = data_vertices[2 * q[i]];
        quad[i][1] = data_vertices[2 * q[i] + 1];
    }
    int offset = -1;
    for (int i = 0; i < 4; i++) {
        if (crossProduct(quad[i], quad[(i + 3) % 4], quad[(i + 1) % 4]) > 0) {
            offset = i;
        }
    }
    int* hex = new int[6];
    double innerPoint[] = { data_vertices[2 * inner],data_vertices[2 * inner + 1] };
    if (offset == -1) {
        double* temp = quad[3];
        if (crossProduct(innerPoint, quad[1], quad[3]) > 0) {
            temp = twoLineIntersect(quad[1], innerPoint, quad[0], quad[3]);
        }
        
        double steiner[] = { (temp[0] + quad[0][0] + innerPoint[0]) / 3,(temp[1] + quad[1][0] + innerPoint[1]) / 3 };
        for (int i = 0; i < 4; i++) {
            hex[i] = q[(i+1) % 4];
        }
        int s = data_vertices.size() / 2;
        data_vertices.push_back(steiner[0]);
        data_vertices.push_back(steiner[1]);
        hex[4] = s;
        hex[5] = inner;
        quads.push_back(q[0]);
        quads.push_back(q[1]);
        quads.push_back(inner);
        quads.push_back(s);
    }
    else {
        if (crossProduct(quad[(offset + 0) % 4], quad[(offset + 2) % 4], innerPoint) < 0) {
            double* intersect1 = twoLineIntersect(quad[(offset + 1) % 4], innerPoint, quad[(offset + 2) % 4], quad[(offset + 3) % 4]);
            double* intersect2 = innerPoint;
            if (crossProduct(quad[(offset + 0) % 4], quad[(offset + 3) % 4], innerPoint) > 0)
                intersect2 = twoLineIntersect(quad[(offset + 1) % 4], innerPoint, quad[(offset + 0) % 4], quad[(offset + 3) % 4]);
            double steiner[] = { (intersect1[0] + intersect2[0] + quad[(offset + 0) % 4][0]) / 3,(intersect1[1] + intersect2[1] + quad[(offset + 0) % 4][1]) / 3 };
            int s = data_vertices.size() / 2;
            data_vertices.push_back(steiner[0]);
            data_vertices.push_back(steiner[1]);
            quads.push_back(q[(offset + 0) % 4]);
            quads.push_back(q[(offset + 1) % 4]);
            quads.push_back(inner);
            quads.push_back(s);
            for (int i = 0; i < 4; i++) {
                hex[i] = q[(i + 1 + offset) % 4];
            }
            hex[4] = s;
            hex[5] = inner;
        }
        else {
            cout << offset << endl;
            double* intersect1 = twoLineIntersect(quad[(offset + 3) % 4], innerPoint, quad[(offset + 1) % 4], quad[(offset + 2) % 4]);
            double* intersect2 = innerPoint;
            if (crossProduct(quad[(offset + 0) % 4], quad[(offset + 1) % 4], innerPoint) < 0)
                intersect2 = twoLineIntersect(quad[(offset + 3) % 4], innerPoint, quad[(offset + 0) % 4], quad[(offset + 1) % 4]);
            double steiner[] = { (intersect1[0] + intersect2[0] + quad[(offset + 0) % 4][0]) / 3,(intersect1[1] + intersect2[1] + quad[(offset + 0) % 4][1]) / 3 };
            int s = data_vertices.size() / 2;
            data_vertices.push_back(steiner[0]);
            data_vertices.push_back(steiner[1]);
            quads.push_back(q[(offset + 0) % 4]);
            quads.push_back(q[(offset + 3) % 4]);
            quads.push_back(inner);
            quads.push_back(s);
            for (int i = 0; i <  4; i++) {
                hex[i] = q[(i + offset) % 4];
            }
            hex[4] = inner;
            hex[5] = s;
        }
    }
    quadHexagon(hex, quads, data_vertices);
}

void statistic(vector<int>& quads, int quad_num, double* vertices_array) {
    double max_max_angle = 0.0;
    double min_min_angle = 180.0;
    double min_max_angle = 180.0;
    double max_min_angle = 0.0;
    double avg_max_angle = 0.0;
    double avg_min_angle = 0.0;
    double squ_max = 0.0;
    double squ_min = 0.0;
    double squ = 0.0;
    double pi = 3.1415926536;
    double rad_to_deg = 180 / pi;
    for (int i = 0; i < quad_num; i++) {
        double max_angle = 0.0;
        double min_angle = 180.0;
        double** p = new double* [4];
        for (int j = 0; j < 4; j++) {
            p[j] = new double[2];
            p[j][0] = vertices_array[2 * quads[4 * i + j]];
            p[j][1] = vertices_array[2 * quads[4 * i + j] + 1];
        }
        for (int j = 0; j < 4; j ++) {
            double angle = acos(cosine(p[j], p[(j + 3) % 4], p[j], p[(j + 1) % 4]));
            max_angle = max(max_angle, angle);
            min_angle = min(min_angle, angle);
            squ += (angle * angle);
        }
        max_max_angle = max(max_max_angle, max_angle);
        min_min_angle = min(min_min_angle, min_angle);
        max_min_angle = max(max_min_angle, min_angle);
        min_max_angle = min(min_max_angle, max_angle);
        avg_max_angle += max_angle;
        avg_min_angle += min_angle;
        squ_max += max_angle * max_angle;
        squ_min += min_angle * min_angle;
    }
    cout << "=======  Statistic of quadrangles  =======" << endl;
    cout << "the maximum of maximum angle of each quadrangle is: " << max_max_angle * rad_to_deg << endl;
    cout << "the minimum of minimum angle of each quadrangle is: " << min_min_angle * rad_to_deg << endl;
    cout << "the minimum of maximum angle of each quadrangle is: " << min_max_angle * rad_to_deg << endl;
    cout << "the maximum of minimum angle of each quadrangle is: " << max_min_angle * rad_to_deg << endl;
    cout << "the average of maximum angle of each quadrangle is: " << avg_max_angle * rad_to_deg / quad_num << endl;
    cout << "the average of minimum angle of each quadrangle is: " << avg_min_angle * rad_to_deg / quad_num << endl;
    cout << "the standard deviation of maximum angle of each quadrangle is: " << sqrt((squ_max - avg_max_angle * avg_max_angle / quad_num) / quad_num) * rad_to_deg << endl;
    cout << "the standard deviation of minimum angle of each quadrangle is: " << sqrt((squ_min - avg_min_angle * avg_min_angle / quad_num) / quad_num) * rad_to_deg << endl;
    cout << "the standard deviation of all the angle of each quadrangle is: " << sqrt(squ / (4 * quad_num) -pi * pi / 4) * rad_to_deg << endl;
}

void outputMesh(const char* fileName, vector<double>& vertices, int n_points, vector<vector<int> >& meshs) {
    FILE* file = fopen(fileName, "w");
    fprintf(file, "OFF\n");
    fprintf(file, "%d %d %d\n", n_points, meshs.size(), 0);
    for (int i = 0; i < n_points; i++) {
        fprintf(file, "%f %f %f\n", vertices[2 * i], vertices[2 * i + 1], 0.0);
    }
    for (int i = 0; i < meshs.size(); i ++) {
        vector<int>mesh = meshs[i];
        fprintf(file, "%d ", mesh.size());
        for (int j = 0; j < mesh.size(); j++) {
            fprintf(file, "%d ", mesh[j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

int main() {
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(800, 800, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }
    glViewport(0, 0, 800, 800);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
;
    int size = 100;
    double** vertices = new double* [size];
    generatePointSet(vertices, size, 2,9);
    
    int y = 0, h = 0;
    convexLayer(vertices, size, y, h);
    int* spiral_line = new int[size + 1];
    spiral_line[0] = h;
    for (int i = 1; i <= size; i++)spiral_line[i] = i - 1;

    int triangles_num = (2 * size - 3 - h);
    int* triangles = new int[3 * triangles_num];
    for (int i = 0; i < 3 * triangles_num; i++)triangles[i] = 0;
    pathTriangulation(vertices, size, y, h, triangles);
    vector<vector<int> >triangle_mesh;
    for (int i = 0; i < triangles_num; i++) {
        vector<int>triangle;
        for (int j = 0; j < 3; j++)triangle.push_back(triangles[3 * i + j]);
        triangle_mesh.push_back(triangle);
    }
    vector<double> data_vertices;
    for (int i = 0; i < size; i++) {
        data_vertices.push_back(vertices[i][0]);
        data_vertices.push_back(vertices[i][1]);
    }

    //combine each 4 consceutive triangles to hexagons
    vector<int>hexagons;
    generateHexagons(triangles, triangles_num, hexagons);
    int hexagon_num = triangles_num / 4;
    vector<vector<int> >hexagon_mesh;

    //first classify hexagon or quadrangle with inner point, then quadrangulate it
    vector<int>quads;
    for (int i = 0; i < hexagon_num; i++) {
        if (hexagons[6 * i] == hexagons[6 * i + 2]) {
            int* quad = new int[4];
            vector<int>mesh;
            for (int j = 0; j < 4; j++) {
                quad[j] = hexagons[6 * i + 2 + j];
                mesh.push_back(hexagons[6 * i + 2 + j]);
            }
            int inner = hexagons[6 * i + 1];
            hexagon_mesh.push_back(mesh);
            quadQuad(quad, inner, quads, data_vertices);           
        }
        else {
            int* hex = new int[6];
            vector<int>mesh;
            for (int j = 0; j < 6; j++) {
                hex[j] = hexagons[6 * i + j];
                mesh.push_back(hexagons[6 * i + j]);
            }
            hexagon_mesh.push_back(mesh);
            quadHexagon(hex, quads, data_vertices);
        }
    }

    //When after quadrangulating each hexagons, the remain triangles are 2 or 3, and 2 of them form a reflex quadrangle.
    int temp = triangles_num % 4;
    if (temp > 1) {
        temp--;
        double** quad = new double* [4];
        int* q = new int[4];
        if (triangles[3 * temp] == triangles[3 * temp - 1]) {
            quad[0] = vertices[triangles[3 * temp - 3]]; q[0] = triangles[3 * temp - 3];
            quad[1] = vertices[triangles[3 * temp - 2]]; q[1] = triangles[3 * temp - 2];
            quad[2] = vertices[triangles[3 * temp + 2]]; q[2] = triangles[3 * temp + 2];
            quad[3] = vertices[triangles[3 * temp]]; q[3] = triangles[3 * temp];
        }
        else {
            quad[0] = vertices[triangles[3 * temp - 2]]; q[0] = triangles[3 * temp - 2];
            quad[1] = vertices[triangles[3 * temp - 1]]; q[1] = triangles[3 * temp - 1];
            quad[2] = vertices[triangles[3 * temp + 2]]; q[2] = triangles[3 * temp + 2];
            quad[3] = vertices[triangles[3 * temp]]; q[3] = triangles[3 * temp];
        }
        if (!isIntersect(quad[0], quad[2], quad[1], quad[3])) {
            double innerpoint[2] = { (quad[0][0] + quad[1][0] + quad[2][0]) / 3,(quad[0][1] + quad[1][1] + quad[2][1]) / 3 };
            int inner = data_vertices.size() / 2;
            data_vertices.push_back(innerpoint[0]);
            data_vertices.push_back(innerpoint[1]);
            quadQuad(q, inner, quads, data_vertices);
        }
        else {
            vector<int>mesh;
            for (int i = 0; i < 4; i++) {
                mesh.push_back(q[i]);
                quads.push_back(q[i]);
            }
            hexagon_mesh.push_back(mesh);
        }
    }
    vector<vector<int> >quad_meah;
    for (int i = 0; i < quads.size() / 4; i++) {
        vector<int>mesh;
        for (int j = 0; j < 4; j++) {
            mesh.push_back(quads[4 * i + j]);
        }
        quad_meah.push_back(mesh);
    }
    if (triangles_num % 2 == 1) {
        vector<int>triangle;
        for (int i = 0; i < 3; i++) {
            triangle.push_back(triangles[i]);
        }
        hexagon_mesh.push_back(triangle);
        quad_meah.push_back(triangle);
    }
    outputMesh("triangle_mesh.off", data_vertices, size, triangle_mesh);
    outputMesh("hexagon_mesh.off", data_vertices,size, hexagon_mesh);
    outputMesh("quadrangle_mesh.off", data_vertices, data_vertices.size() / 2, quad_meah);
    int quad_num = quads.size() / 4;
    double* vertiecs_array = new double[data_vertices.size()];
    for (int i = 0; i < data_vertices.size(); i++)vertiecs_array[i] = data_vertices[i];
    statistic(quads, quad_num, vertiecs_array);

    //code below this line are GLFW part.
    const char* vertexShaderSource = "#version 330 core\n"
        "layout (location = 0) in vec2 aPos;\n"
        "void main()\n"
        "{\n"
        "   gl_Position = vec4(aPos, 0.0, 1.0);\n"
        "}\0";

    unsigned int vertexShader;
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);

    const char* fragmentShaderSource = "#version 330 core\n"
        "out vec4 FragColor;\n"
        "uniform vec3 color;\n"
        "void main()\n"
        "{\n"
        "   FragColor = vec4(color,1.0);\n"
        "}\n\0";

    unsigned int fragmentShader;
    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);

    unsigned int shaderProgram;
    shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    unsigned int VAO;
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    unsigned int VBO;
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);

    unsigned int EBO;
    glGenBuffers(1, &EBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);

    glVertexAttribPointer(0, 2, GL_DOUBLE, GL_FALSE, 2 * sizeof(double), (void*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    int triangle_start = 0, triangle_end = 0;
    bool triangleFinish = false;
    int hexagon_start = 0, hexagon_end = -1;
    int quad_end = 0;
    while (!glfwWindowShouldClose(window))
    {
        Sleep(40);
        processInput(window);
        glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glUseProgram(shaderProgram);
        glBindVertexArray(VAO);
        glUniform3f(glGetUniformLocation(shaderProgram, "color"), 1.0f, 1.0f, 1.0f);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, data_vertices.size() * sizeof(double), vertiecs_array, GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * triangles_num * sizeof(int), triangles , GL_STATIC_DRAW);
        glDrawElements(GL_TRIANGLES, 3 * (triangle_end - triangle_start), GL_UNSIGNED_INT, 0);
        if (triangle_end < triangles_num && !triangleFinish)triangle_end++;

        if (hexagon_end == hexagon_num && quad_end < quad_num) quad_end++;
        glUniform3f(glGetUniformLocation(shaderProgram, "color"), 1.0f, 1.0f, 1.0f);
        for (int i = hexagon_start; i < quad_end; i++) {
            int q[] = { quads[4 * i],quads[4 * i + 1],quads[4 * i + 2],quads[4 * i + 3],quads[4 * i] };
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, 5 * sizeof(int), q, GL_STATIC_DRAW);
            glDrawElements(GL_LINE_STRIP, 5, GL_UNSIGNED_INT, 0);
        }

        if (triangle_end == triangles_num)triangleFinish = true;

        if (triangleFinish && hexagon_end < hexagon_num) {
            triangle_end -= 4;
            if (triangle_end < 0 && triangles_num % 2 == 1)triangle_end = 1;
            hexagon_end++;
        }
        glUniform3f(glGetUniformLocation(shaderProgram, "color"), 1.0f, 1.0f, 1.0f);
        for (int i = hexagon_start; i < hexagon_end; i++) {
            int hexagon[] = { hexagons[6 * i],hexagons[6 * i + 1],hexagons[6 * i + 2],hexagons[6 * i + 3],hexagons[6 * i + 4],hexagons[6 * i + 5],hexagons[6 * i] };
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, 7 * sizeof(int), hexagon, GL_STATIC_DRAW);
            glDrawElements(GL_LINE_STRIP, 7, GL_UNSIGNED_INT, 0);
        }

        glUniform3f(glGetUniformLocation(shaderProgram, "color"), 1.0f, 1.0f, 1.0f);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, (size + 1) * sizeof(int), spiral_line, GL_STATIC_DRAW);
        glDrawElements(GL_LINE_STRIP, size + 1, GL_UNSIGNED_INT, 0);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glfwTerminate();

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteProgram(shaderProgram);
    return 0;
}