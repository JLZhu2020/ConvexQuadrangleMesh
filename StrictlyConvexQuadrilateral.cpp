/*  */
#define GLEW_STATIC
#define STB_IMAGE_IMPLEMENTATION

#include <iostream>
#include <random>
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

class Edge;
class Vertex;
class Face;

class Vertex {
public:
    double x_coordinate;
    double y_coordinate;
    Edge* pre;
    Edge* next;
    Vertex() {}
    Vertex(double x_coordinate, double y_coordinate) {
        this->x_coordinate = x_coordinate;
        this->y_coordinate = y_coordinate;
    }
    Vertex(Edge *pre, Edge *next, double x_coordinate, double y_coordinate) {
        this->x_coordinate = x_coordinate;
        this->y_coordinate = y_coordinate;
        this->pre = pre;
        this->next = next;
    }
    Vertex copy(Edge *pre, Edge *next) {
        return Vertex(pre, next, this->x_coordinate, this->y_coordinate);
    }
    Vertex copy() {
        return Vertex(this->x_coordinate, this->y_coordinate);
    }
    void setPre(Edge *pre) {
        this->pre = pre;
    }
    void setNext(Edge *next) {
        this->next = next;
    }
};

class Edge{
public:
    Face* face;
    Edge* twin;
    Edge* next;
    Edge* pre;
    double* coordinate;
    Edge() {
        this->pre = nullptr;
        this->next = nullptr;
        this->coordinate = new double[2];
        this->coordinate[0] = 0;
        this->coordinate[1] = 0;
    }
    Edge(Edge& pre, Edge& next, double* coord){
        this->pre = &pre;
        this->next = &next;
        this->coordinate = new double[2];
        this->coordinate[0] = coord[0];
        this->coordinate[1] = coord[1];
    }
    void setTwin(Edge& twin) {
        this->twin = &twin;
        twin.twin = this;
    }
    void setFace(Face& face) {
        this->face = &face;
    }
    void setPre(Edge& pre) {
        this->pre = &pre;
        pre.next = this;
    }
    void setNext(Edge& next) {
        this->next = &next;
        next.pre = this;
    }
};

//class Face {
//public:
//    Edge* edge;
//    int size;
//    Face(Edge& edge, int size) {
//        this->edge = &edge;
//        this->size = size;
//    }
//
//    Face& addSteiner(int begin, int end, double x_coordinate, double y_coordinate) {
//        Vertex steinerQ(x_coordinate, y_coordinate);
//        Vertex steinerR(x_coordinate, y_coordinate);
//        Vertex beginR = this->list[begin];
//        Vertex beginQ = beginR.copy();
//        beginQ.next = beginR.next;
//        Vertex endR = this->list[end];
//        Vertex endQ = endR.copy();
//        endQ.pre = endR.pre;
//        Edge edge1R(steinerR, endR);
//        Edge edge2R(beginR, steinerR);
//        Edge edge1Q(endQ, steinerQ);
//        Edge edge2Q(steinerQ, beginQ);
//        edge1Q.setTwin(&edge1R);
//        edge2Q.setTwin(&edge2R);
//        edge1R.setFace(this);
//        edge2R.setFace(this);
//        vector<Vertex> quadVertice = { beginQ,this->list[(begin + 1) % size],endQ,steinerQ };
//        Face quad(quadVertice, 4);
//        this->list[(begin + 1) % this->size] = steinerR;
//        return quad;
//    }
//
//    void split3(int begin, double x_coordinate, double y_coordinate, vector<Face> quads) {
//        if (this->size != 6)return;
//        for (int i = 0; i <= 4; i+=2) {
//            Vertex v0 = this->list[(begin + i + 0) % 6];
//            Vertex v1 = this->list[(begin + i + 1) % 6];
//            Vertex v2 = this->list[(begin + i + 1) % 6].copy();
//            v2.setPre(this->list[(begin + i + 1) % 6].pre);
//            Vertex steiner(x_coordinate, y_coordinate);
//            Edge e1(v2, steiner);
//            Edge e2(steiner, v0);
//            vector<Vertex> quadVertice;
//            quadVertice.push_back(v0);
//            quadVertice.push_back(v1);
//            quadVertice.push_back(v2);
//            quadVertice.push_back(steiner);
//            Face quad(quadVertice, 4);
//            quads.push_back(quad);
//        }
//        for (int i = 0; i < 3; i++) {
//            Edge edge1 = *(quads[i].list[2].next);
//            Edge edge2 = *(quads[(i + 1) % 3].list[0].pre);
//            edge1.setTwin(&edge2);
//        }
//    }
//
//    bool combine(Face* another) {
//        int begin = -1;
//        int end = -1;
//        for (int i=0; i < this->size; i++) {
//            if (this->list[i].next->twin->face == another) {
//                if (begin < 0)begin = i;
//                end = max(i, end);
//            }
//        }
//        if (begin < 0)return false;
//        end++;
//        if (end == this->size && begin==0) {
//            end = 0;
//            begin = size - 1;
//            while (this->list[end].next->twin->face == another) {
//                end++;
//            }
//            while (this->list[begin].pre->twin->face == another) {
//                begin--;
//            }
//        }
//        vector<Vertex> vertices;
//        Vertex pin = this->list[end];
//        Vertex stop = this->list[begin];
//        while (&pin != &stop) {
//            vertices.push_back(pin);
//            pin = *(pin.next->to);
//        }
//        pin = *(stop.next->twin->to);
//        stop = *(this->list[end].pre->twin->from);
//        while (&pin != &stop) {
//            vertices.push_back(pin);
//            pin = *(pin.next->to);
//        }
//        this->list = vertices;
//        this->size = vertices.size();
//    }
//};

/*
* This function will reorder the input point set, the result is the spiral line.
* Parameter:
* data_vertices: the point set data, the shape is size*2, data_vertices[i][0] is the x coordinate of i th point, data_vertices[i][1] is the y coordinate.
* size: means how many points in the point set.
* h: the index from 0 to h is the convex hull of the point set.
* y: the index for the inner starshape, please see document for detail.
*/
void polygonChain(double** data_vertices, int size, int& y, int& h) {
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
* data_vertices, size, y, h: data get form the function polygonChain
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

void quadStarshape(double** hexagon, int& offset, vector<double>&quads, double* steiner) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            quads.push_back(hexagon[(j + 2 * i + offset) % 6][0]);
            quads.push_back(hexagon[(j + 2 * i + offset) % 6][1]);
        }
        quads.push_back(steiner[0]);
        quads.push_back(steiner[1]);
        quads.push_back(hexagon[(2 * i + offset) % 6][0]);
        quads.push_back(hexagon[(2 * i + offset) % 6][1]);
    }
}

void quadSteiner(double** hexagon, int& offset, vector<double>&quads, double* steiner, int start) {
    for (int i = start; i < start + 3; i++) {
        quads.push_back(hexagon[(i + offset) % 6][0]);
        quads.push_back(hexagon[(i + offset) % 6][1]);
    }
    quads.push_back(steiner[0]);
    quads.push_back(steiner[1]);
    quads.push_back(hexagon[(offset + start) % 6][0]);
    quads.push_back(hexagon[(offset + start) % 6][1]);
}

void quadCCCCCC(double** hexagon, int &offset, vector<double>&quads) {
    cout << "quad CCCCCC" << endl;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 4; j++) {
            quads.push_back(hexagon[(j + 3 * i + offset) % 6][0]);
            quads.push_back(hexagon[(j + 3 * i + offset) % 6][1]);
        }
        quads.push_back(hexagon[(3 * i + offset) % 6][0]);
        quads.push_back(hexagon[(3 * i + offset) % 6][1]);
    }
}

void quadRCCCCC(double** hexagon, int &offset, vector<double>&quads) {
    //verify if d is in wedge(a)
    if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 5) % 6]) > 0 &&
        crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 3) % 6]) > 0)
    {
        //this is the RCCCCC case 1
        cout << "quad RCCCCC case 1" << endl;
        quadCCCCCC(hexagon, offset, quads);
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
        double steiner[] = { (intersect1[0] + intersect2[0] + hexagon[(offset + 0) % 6][0]) / 3,(intersect1[1] + intersect2[1] + hexagon[(offset + 0) % 6][1]) / 3 };
        quadStarshape(hexagon, offset, quads, steiner);
    }
    else {
        //the c could not see e case
        if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 5) % 6]) <= 0) {
            //if d lies L(a,b) RCCCCC case 3
            cout << "quad RCCCCC case 3" << endl;
            double* intersect1 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6]);
            double* intersect2 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6]);
            double steiner[] = { (intersect1[0] + intersect2[0] + hexagon[(offset + 0) % 6][0]) / 3,(intersect1[1] + intersect2[1] + hexagon[(offset + 0) % 6][1]) / 3 };
            quadSteiner(hexagon, offset, quads, steiner, 4);
            offset = (offset + 5) % 6;
            hexagon[(offset + 0) % 6][0] = steiner[0];
            hexagon[(offset + 0) % 6][1] = steiner[1];
            quadRCCCCC(hexagon, offset, quads);
        }
        else {
            //if d lies R(a,b) RCCCCC case 4
            cout << "quad RCCCCC case 4" << endl;
            double* intersect1 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 3) % 6]);
            double* intersect2 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 3) % 6]);
            double steiner[] = { (intersect1[0] + intersect2[0] + hexagon[(offset + 0) % 6][0]) / 3,(intersect1[1] + intersect2[1] + hexagon[(offset + 0) % 6][1]) / 3 };
            quadSteiner(hexagon, offset, quads, steiner, 0);
            offset = (offset + 1) % 6;
            hexagon[(offset + 0) % 6][0] = steiner[0];
            hexagon[(offset + 0) % 6][1] = steiner[1];
            quadRCCCCC(hexagon, offset, quads);
        }
    }
}

void quadRCRCCC(double** hexagon, int& offset, vector<double>&quads) {
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
        double steiner[] = { (triangleLeft[0] + triangleRight[0] + triangleTop[0]) / 3,(triangleLeft[1] + triangleRight[1] + triangleTop[1]) / 3 };
        quadStarshape(hexagon, offset, quads, steiner);
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
        quadSteiner(hexagon, offset, quads, steiner, 0);
        offset = (offset + 1) % 6;
        hexagon[(offset + 0) % 6][0] = steiner[0];
        hexagon[(offset + 0) % 6][1] = steiner[1];
        quadRCCCCC(hexagon, offset, quads);
    }
}

void quadRRCCCC(double** hexagon, int& offset, vector<double>&quads) {
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
    double steiner[] = { (intersect0[0] + intersect1[0] + intersect2[0]) / 3,(intersect0[1] + intersect1[1] + intersect2[1]) / 3 };
    quadSteiner(hexagon, offset, quads, steiner, 4);
    offset = (offset + 5) % 6;
    hexagon[(offset + 0) % 6][0] = steiner[0];
    hexagon[(offset + 0) % 6][1] = steiner[1];
    quadRCRCCC(hexagon, offset, quads);
}

void quadRCCRCC(double** hexagon, int& offset, vector<double>&quads) {
    if (crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 5) % 6]) > 0 &&
        crossProduct(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 3) % 6]) > 0) {
        //RCCRCC case 1
        cout << "quad RCCRCC case 1" << endl;
        quadCCCCCC(hexagon, offset, quads);
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
        quadSteiner(hexagon, offset, quads, steiner, 3);
        offset = (offset + 4) % 6;
        hexagon[(offset + 0) % 6][0] = steiner[0];
        hexagon[(offset + 0) % 6][1] = steiner[1];
        quadRCRCCC(hexagon, offset, quads);
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
        quadSteiner(hexagon, offset, quads, steiner, 4);
        offset = (offset + 3) % 6;
        hexagon[(offset + 2) % 6][0] = steiner[0];
        hexagon[(offset + 2) % 6][1] = steiner[1];
        quadRCRCCC(hexagon, offset, quads);    
    }
}

void quadRCRCRC(double** hexagon, int& offset, vector<double>&quads) {
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
        double steiner[] = { (triangle[0][0] + triangle[1][0] + triangle[2][0]) / 3, (triangle[0][1] + triangle[1][1] + triangle[2][1]) / 3 };
        for (int i = 0; i < 3; i++)cout << triangle[i][0] << " " << triangle[i][1] << endl;
        quadStarshape(hexagon, offset, quads, steiner);
    }
    else {
        //RCRCRC case 2
        cout << "quad RCRCRC case 2" << endl;
        double* intersect1 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 5) % 6]);
        double* intersect2 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 1) % 6], hexagon[(offset + 2) % 6], hexagon[(offset + 3) % 6]);
        double steiner[] = { (intersect1[0] + intersect2[0] + hexagon[(offset + 2) % 6][0]) / 3,(intersect1[1] + intersect2[1] + hexagon[(offset + 2) % 6][1]) / 3 };
        quadSteiner(hexagon, offset, quads, steiner, 0);
        offset = (offset + 1) % 6;
        hexagon[(offset + 0) % 6][0] = steiner[0];
        hexagon[(offset + 0) % 6][1] = steiner[1];
        quadRCCRCC(hexagon, offset, quads);
    }
}

void quadRCRRCC(double** hexagon, int& offset, vector<double>&quads);

void quadRRCRCC(double** hexagon, int& offset, vector<double>&quads) {
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
        quadSteiner(hexagon, offset, quads, steiner, 4);
        offset = (offset + 5) % 6;
        hexagon[(offset + 0) % 6][0] = steiner[0];
        hexagon[(offset + 0) % 6][1] = steiner[1];
        quadRCRCRC(hexagon, offset, quads);
    }
    else {
        //RRCRCC case 2
        cout << "quad RRCRCC case 2" << endl;
        double* intersect1 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6]);
        double steiner[] = { (hexagon[(offset + 3) % 6][0] + hexagon[(offset + 5) % 6][0] + intersect1[0]) / 3,
            (hexagon[(offset + 3) % 6][1] + hexagon[(offset + 5) % 6][1] + intersect1[1]) / 3 };
        quadSteiner(hexagon, offset, quads, steiner, 3);
        offset = (offset + 4) % 6;
        hexagon[(offset + 0) % 6][0] = steiner[0];
        hexagon[(offset + 0) % 6][1] = steiner[1];
        quadRCRRCC(hexagon, offset, quads);
    }
}

void quadRCRRCC(double** hexagon, int& offset, vector<double>&quads) {
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
        quadSteiner(hexagon, offset, quads, steiner, 3);
        offset = (offset + 4) % 6;
        hexagon[(offset + 0) % 6][0] = steiner[0];
        hexagon[(offset + 0) % 6][1] = steiner[1];
        quadRCRCRC(hexagon, offset, quads);
    }
    else {
        //RRCRCC case 2
        cout << "quad RCRRCC case 1" << endl;
        double* intersect1 = twoLineIntersect(hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6], hexagon[(offset + 3) % 6], hexagon[(offset + 4) % 6]);
        double steiner[] = { (hexagon[(offset + 0) % 6][0] + hexagon[(offset + 4) % 6][0] + intersect1[0]) / 3,
                            (hexagon[(offset + 0) % 6][1] + hexagon[(offset + 4) % 6][1] + intersect1[1]) / 3 };
        quadSteiner(hexagon, offset, quads, steiner, 4);
        offset = (offset + 2) % 6;
        hexagon[(offset + 3) % 6][0] = steiner[0];
        hexagon[(offset + 3) % 6][1] = steiner[1];
        quadRRCRCC(hexagon, offset, quads);
    }   
}

void quadRRRCCC(double** hexagon, int& offset, vector<double>&quads) {
    cout << "quad RRRCCC" << endl;
    double* intersect1 = twoLineIntersect(hexagon[(offset + 1) % 6], hexagon[(offset + 4) % 6], hexagon[(offset + 0) % 6], hexagon[(offset + 5) % 6]);
    double steiner[] = { (hexagon[(offset + 0) % 6][0] + hexagon[(offset + 4) % 6][0] + intersect1[0]) / 3,
                        (hexagon[(offset + 0) % 6][1] + hexagon[(offset + 4) % 6][1] + intersect1[1]) / 3 };
    quadSteiner(hexagon, offset, quads, steiner, 4);
    offset = (offset + 5) % 6;
    hexagon[(offset + 0) % 6][0] = steiner[0];
    hexagon[(offset + 0) % 6][1] = steiner[1];
    quadRCRRCC(hexagon, offset, quads);
}

/*
* This function will figure out the input hexagon is which kind of shape and call the corresponding quad function.
* Parameter:
* hexagon: the input 6*2 array, should be ordered in counter-clockwised
* quads: to receive result. each 8 double represent a quadrangle
*/
void quadHexagon(double **hexagon, vector<double>&quads) {
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
    if (score == 0 ) quadCCCCCC(hexagon, offset, quads);
    if (score == 1 ) quadRCCCCC(hexagon, offset, quads);
    if (score == 5 ) quadRCRCCC(hexagon, offset, quads);
    if (score == 3 ) quadRRCCCC(hexagon, offset, quads);
    if (score == 9 ) quadRCCRCC(hexagon, offset, quads);
    if (score == 21) quadRCRCRC(hexagon, offset, quads);
    if (score == 11) quadRRCRCC(hexagon, offset, quads);
    if (score == 13) quadRCRRCC(hexagon, offset, quads);
    if (score == 7 ) quadRRRCCC(hexagon, offset, quads);
}

/*
* This function will quadrangulate the input quadrangle with an inner point
*/
void quadQuad(double** quad, double* innerPoint, vector<double>& quads) {
    int offset = -1;
    for (int i = 0; i < 4; i++) {
        if (crossProduct(quad[i], quad[(i + 3) % 4], quad[(i + 1) % 4]) > 0) {
            offset = i;
        }
    }
    double** hexagon = new double* [6];
    if (offset == -1) {
        double* temp = quad[3];
        if (crossProduct(innerPoint, quad[1], quad[3]) > 0) {
            temp = twoLineIntersect(quad[1], innerPoint, quad[0], quad[3]);
        }
        double steiner[] = { (temp[0] + quad[0][0] + innerPoint[0]) / 3,(temp[1] + quad[1][0] + innerPoint[1]) / 3 };
        for (int i = 0; i < 4; i++) {
            hexagon[i] = quad[(i+1) % 4];
        }
        hexagon[4] = steiner;
        hexagon[5] = innerPoint;
        quads.push_back(quad[0][0]); quads.push_back(quad[0][1]);
        quads.push_back(quad[1][0]); quads.push_back(quad[1][1]);
        quads.push_back(innerPoint[0]); quads.push_back(innerPoint[1]);
        quads.push_back(steiner[0]); quads.push_back(steiner[1]);
        quads.push_back(quad[0][0]); quads.push_back(quad[0][1]);
    }
    else {
        if (crossProduct(quad[(offset + 0) % 4], quad[(offset + 2) % 4], innerPoint) < 0) {
            double* intersect1 = twoLineIntersect(quad[(offset + 1) % 4], innerPoint, quad[(offset + 2) % 4], quad[(offset + 3) % 4]);
            double* intersect2 = innerPoint;
            if (crossProduct(quad[(offset + 0) % 4], quad[(offset + 3) % 4], innerPoint) > 0)
                intersect2 = twoLineIntersect(quad[(offset + 1) % 4], innerPoint, quad[(offset + 0) % 4], quad[(offset + 3) % 4]);
            double steiner[] = { (intersect1[0] + intersect2[0] + quad[(offset + 0) % 4][0]) / 3,(intersect1[1] + intersect2[1] + quad[(offset + 0) % 4][1]) / 3 };
            quads.push_back(quad[(offset + 0) % 4][0]); quads.push_back(quad[(offset + 0) % 4][1]);
            quads.push_back(quad[(offset + 1) % 4][0]); quads.push_back(quad[(offset + 1) % 4][1]);
            quads.push_back(innerPoint[0]); quads.push_back(innerPoint[1]);
            quads.push_back(steiner[0]); quads.push_back(steiner[1]);
            quads.push_back(quad[(offset + 0) % 4][0]); quads.push_back(quad[(offset + 0) % 4][1]);
            for (int i = 0; i < 4; i++) {
                hexagon[i] = quad[(i + 1 + offset) % 4];
            }
            hexagon[4] = steiner;
            hexagon[5] = innerPoint;
        }
        else {
            cout << offset << endl;
            double* intersect1 = twoLineIntersect(quad[(offset + 3) % 4], innerPoint, quad[(offset + 1) % 4], quad[(offset + 2) % 4]);
            double* intersect2 = innerPoint;
            if (crossProduct(quad[(offset + 0) % 4], quad[(offset + 1) % 4], innerPoint) < 0)
                intersect2 = twoLineIntersect(quad[(offset + 3) % 4], innerPoint, quad[(offset + 0) % 4], quad[(offset + 1) % 4]);
            double steiner[] = { (intersect1[0] + intersect2[0] + quad[(offset + 0) % 4][0]) / 3,(intersect1[1] + intersect2[1] + quad[(offset + 0) % 4][1]) / 3 };
            quads.push_back(quad[(offset + 0) % 4][0]); quads.push_back(quad[(offset + 0) % 4][1]);
            quads.push_back(quad[(offset + 3) % 4][0]); quads.push_back(quad[(offset + 3) % 4][1]);
            quads.push_back(innerPoint[0]); quads.push_back(innerPoint[1]);
            quads.push_back(steiner[0]); quads.push_back(steiner[1]);
            quads.push_back(quad[(offset + 0) % 4][0]); quads.push_back(quad[(offset + 0) % 4][1]);
            for (int i = 0; i <  4; i++) {
                hexagon[i] = quad[(i + offset) % 4];
            }
            hexagon[4] = innerPoint;
            hexagon[5] = steiner;
        }
    }
    quadHexagon(hexagon, quads);
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
 //   double data_vertices[] = { -1,0,  0,-1,  1,0,  0,1,  -1,0,  -1,0,  -1,0 };
    int size = 500;
    double** vertices = new double* [size];
    generatePointSet(vertices, size, 1,25);
    int y = 0, h = 0;
    polygonChain(vertices, size, y, h);
    int triangles_num = (2 * size - 3 - h);
    int* triangles = new int[3 * triangles_num];
    for (int i = 0; i < 3 * triangles_num; i++)triangles[i] = 0;
    pathTriangulation(vertices, size, y, h, triangles);
    double* data_vertices = new double[size * 2];
    for (int i = 0; i < size; i++) {
        data_vertices[i * 2] = vertices[i][0];
        data_vertices[i * 2 + 1] = vertices[i][1];
    }

    int* spiral_line = new int[size + 1];
    spiral_line[0] = h;
    for (int i = 1; i <= size; i++)spiral_line[i] = i - 1;
    vector<int>hexagons;
    for (int i = triangles_num -4; i >=0; i -= 4) {
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
    int hexagon_num = triangles_num / 4;
    double* data_hexagon = new double[hexagon_num * 12];
    for (int i = 0; i < hexagon_num * 6 ; i++) {
        data_hexagon[2 * i] = data_vertices[2 * hexagons[i]];
        data_hexagon[2 * i + 1] = data_vertices[2 * hexagons[i] + 1];
    }
    int hexagon[] = { 0,1,2,3,4,5,0 };
    vector<double>quads;
    for (int i = 0; i < hexagon_num; i++) {
        if (hexagons[6 * i] == hexagons[6 * i + 2]) {
            double** quad = new double* [4];
            for (int j = 0; j < 4; j++) {
                quad[j] = new double[2];
                quad[j][0] = vertices[hexagons[6 * i + 2 + j]][0];
                quad[j][1] = vertices[hexagons[6 * i + 2 + j]][1];
            }
            double inner[2] = { vertices[hexagons[6 * i + 1]][0],vertices[hexagons[6 * i + 1]][1] };
            quadQuad(quad, inner, quads);
            
        }
        else {
            double** hex = new double* [6];
            for (int j = 0; j < 6; j++) {
                int k = hexagons[6 * i + j];
                hex[j] = new double[2];
                hex[j][0] = data_vertices[2 * k];
                hex[j][1] = data_vertices[2 * k + 1];
            }
            quadHexagon(hex, quads);
        }
    }
    int temp = triangles_num % 4;
    if (temp > 1) {
        temp--;
        double** quad = new double* [4];
        if (triangles[3 * temp] == triangles[3 * temp - 1]) {
            quad[0] = vertices[triangles[3 * temp - 3]];
            quad[1] = vertices[triangles[3 * temp - 2]];
            quad[2] = vertices[triangles[3 * temp + 2]];
            quad[3] = vertices[triangles[3 * temp]];
            for (int i = 0; i < 4; i++)cout << quad[i][0] << " " << quad[i][1] << endl;
        }
        else {
            quad[0] = vertices[triangles[3 * temp - 2]];
            quad[1] = vertices[triangles[3 * temp - 1]];
            quad[2] = vertices[triangles[3 * temp + 2]];
            quad[3] = vertices[triangles[3 * temp]];
        }
        if (!isIntersect(quad[0], quad[2], quad[1], quad[3])) {
            double innerpoint[2] = { (quad[0][0] + quad[1][0] + quad[2][0]) / 3,(quad[0][1] + quad[1][1] + quad[2][1]) / 3 };
            quadQuad(quad, innerpoint, quads);
        }
    }
    double* data_quad = new double[quads.size()];
    for (int i = 0; i < quads.size(); i++) {
        data_quad[i] = quads[i];
    }
    int quad_num = quads.size() / 10;
    int quad[] = { 0,1,2,3,0 };

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
        Sleep(80);
        processInput(window);
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glUseProgram(shaderProgram);
        glBindVertexArray(VAO);
        glUniform3f(glGetUniformLocation(shaderProgram, "color"), 0.0f, 0.0f, 1.0f);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, 2 * size * sizeof(double), data_vertices, GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * triangles_num * sizeof(int), triangles , GL_STATIC_DRAW);
        glDrawElements(GL_TRIANGLES, 3 * (triangle_end - triangle_start), GL_UNSIGNED_INT, 0);
        if (triangle_end < triangles_num && !triangleFinish)triangle_end++;

        if (hexagon_end == hexagon_num && quad_end < quad_num) quad_end++;
        glUniform3f(glGetUniformLocation(shaderProgram, "color"), 1.6f, -0.3f, -0.3f);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 5 * sizeof(int), quad, GL_STATIC_DRAW);
        for (int i = hexagon_start; i < quad_end; i++) {
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferData(GL_ARRAY_BUFFER, 10 * sizeof(double), data_quad + 10 * i, GL_STATIC_DRAW);
            glDrawElements(GL_LINE_STRIP, 5, GL_UNSIGNED_INT, 0);
        }

        if (triangle_end == triangles_num)triangleFinish = true;

        if (triangleFinish && hexagon_end < hexagon_num) {
            triangle_end -= 4;
            if (triangle_end < 0 && triangles_num % 2 == 1)triangle_end = 1;
            hexagon_end++;
        }
        glUniform3f(glGetUniformLocation(shaderProgram, "color"), 0.0f, 1.0f, 0.0f);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 7 * sizeof(int), hexagon, GL_STATIC_DRAW);
        for (int i = hexagon_start; i < hexagon_end; i++) {
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferData(GL_ARRAY_BUFFER, 12 * sizeof(double), data_hexagon + 12 * i, GL_STATIC_DRAW);
            glDrawElements(GL_LINE_STRIP, 7, GL_UNSIGNED_INT, 0);
        }

        glUniform3f(glGetUniformLocation(shaderProgram, "color"), 1.0f, 1.0f, 1.0f);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, 2 * size * sizeof(double), data_vertices, GL_STATIC_DRAW);
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


