#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <iostream>
#include <chrono>
#include <sstream>
#include <stdio.h>
#include <functional>
#include "lbfgs/lbfgs.h"
#include "lbfgs/arithmetic_ansi.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//Vector class from project 1 
class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    Vector normalize() {
        double n = norm();
        return Vector(data[0]/n, data[1]/n, data[2]/n);
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

//Lecture 6 slide 15
Vector intersect(Vector &A, Vector &B, std::vector<Vector> &edge){
    Vector v1 = edge[0], v2 = edge[1];
    Vector N = Vector(v2[1] - v1[1], v1[0] - v2[0]).normalize();
    double t = dot(v1 - A, N) / dot(B - A, N);
    if (t >= 0 && t <= 1){
         return A + t * (B - A);
    }
    return Vector();
}

//Lecture 6 slide 61
bool inside(Vector &P, std::vector<Vector> &edge) {
    Vector v1 = edge[0], v2 = edge[1];
    Vector N = Vector(v2[1] - v1[1], v1[0] - v2[0]).normalize();
    if (dot(P - v1, N) <= 0){
        return true;
    }
    return false;
}


//Polygon class
class Polygon {
public:
    Polygon() {};
    explicit Polygon(const std::vector<Vector> &vertices) {
        this->vertices = vertices;
    };
    void add_vertex(Vector& a) {
        this->vertices.push_back(a);
    };
    double get_area();
    Vector get_centroid();
	std::vector<Vector> vertices = {};
};

//Lecture 7 slide 9 
double Polygon::get_area(){
    if(vertices.size() < 3){ 
        return 0;
    }
    double A = 0;
    for (int i = 0; i < vertices.size() - 1; i++){
        A = A + vertices[i][0]*vertices[i + 1][1] - vertices[i][1]*vertices[i + 1][0];
    }
    A += vertices[vertices.size() - 1][0] * vertices[0][1] - vertices[vertices.size() - 1][1]*vertices[0][0];
    return std::abs(A)/2;
};

Vector Polygon::get_centroid() {
    double Cx = 0, Cy = 0, A=0;
    
    for (int i = 0; i < vertices.size() - 1; i++){
        Cx += (vertices[i][0] + vertices[i+1][0]) * (vertices[i][0] * vertices[i+1][1] - vertices[i+1][0] * vertices[i][1]);
        Cy += (vertices[i][1] + vertices[i+1][1]) * (vertices[i][0] * vertices[i+1][1] - vertices[i+1][0] * vertices[i][1]);
    }

    Cx += (vertices[vertices.size()-1][0] + vertices[0][0]) * (vertices[vertices.size()-1][0] * vertices[0][1] - vertices[0][0] * vertices[vertices.size()-1][1]);
    Cy += (vertices[vertices.size()-1][1] + vertices[0][1]) * (vertices[vertices.size()-1][0] * vertices[0][1] - vertices[0][0] * vertices[vertices.size()-1][1]);
    for (int i = 0; i < vertices.size() - 1; i++){
        A = A + vertices[i][0]*vertices[i + 1][1] - vertices[i][1]*vertices[i + 1][0];
    }
    A += vertices[vertices.size() - 1][0] * vertices[0][1] - vertices[vertices.size() - 1][1]*vertices[0][0];
    A=A/2;
    return Vector(-Cx / (6* A), -Cy / (6 * A), 0.0);
}

//Lecture notes page 88
Polygon clip(Polygon subjectPolygon, Polygon clipPolygon) {
    Polygon outPolygon;
    for(int i = 0; i < clipPolygon.vertices.size(); i++){
        Vector v1 = clipPolygon.vertices[i];
        Vector v2 = clipPolygon.vertices[(i > 0) ? (i - 1) : (clipPolygon.vertices.size() - 1)];
        std::vector<Vector> edge = {v1, v2};
        outPolygon = Polygon();
        for(int j = 0; j < subjectPolygon.vertices.size(); j++){
            Vector curVertex = subjectPolygon.vertices[j];
            Vector prevVertex = subjectPolygon.vertices[(j > 0) ? (j - 1) : (subjectPolygon.vertices.size() - 1)];
            Vector intersection = intersect(prevVertex, curVertex, edge);
            if(inside(curVertex, edge)){
                if(!inside(prevVertex, edge)){
                    outPolygon.add_vertex(intersection);
                }
                outPolygon.add_vertex(curVertex);
            } 
            else if(inside(prevVertex, edge)){
                outPolygon.add_vertex(intersection);
            }
        }
        subjectPolygon = outPolygon;
    }
    return outPolygon;
}


// SVG
//from https://pastebin.com/bEYVtqYy
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
    FILE* f = fopen(filename.c_str(), "w+"); 
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (int i=0; i<polygons.size(); i++) {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \""); 
        for (int j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}

//from https://pastebin.com/jVcNAE5Q
void save_frame(const std::vector<Polygon> &cells, std::string filename, int frameid,int N){
    int W = 1000, H = 1000;
    std::vector<unsigned char> image(W*H * 3, 255);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < cells.size(); i++) {

        double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
        for (int j = 0; j < cells[i].vertices.size(); j++) {
            bminx = std::min(bminx, cells[i].vertices[j][0]);
            bminy = std::min(bminy, cells[i].vertices[j][1]);
            bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
            bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
        }
        bminx = std::min(W-1., std::max(0., W * bminx));
        bminy = std::min(H-1., std::max(0., H * bminy));
        bmaxx = std::max(W-1., std::max(0., W * bmaxx));
        bmaxy = std::max(H-1., std::max(0., H * bmaxy));

        for (int y = bminy; y < bmaxy; y++) {
            for (int x = bminx; x < bmaxx; x++) {
                int prevSign = 0;
                bool isInside = true;
                double mindistEdge = 1E9;
                for (int j = 0; j < cells[i].vertices.size(); j++) {
                    double x0 = cells[i].vertices[j][0] * W;
                    double y0 = cells[i].vertices[j][1] * H;
                    double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
                    double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
                    double det = (x - x0)*(y1-y0) - (y - y0)*(x1-x0);
                    int sign = (det > 0) - (det < 0);
                    if (prevSign == 0) prevSign = sign; else
                        if (sign == 0) sign = prevSign; else
                        if (sign != prevSign) {
                            isInside = false;
                            break;
                        }
                    prevSign = sign;
                    double edgeLen = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
                    double distEdge = std::abs(det)/ edgeLen;
                    double dotp = (x - x0)*(x1 - x0) + (y - y0)*(y1 - y0);
                    if (dotp<0 || dotp>edgeLen*edgeLen) distEdge = 1E9;
                    mindistEdge = std::min(mindistEdge, distEdge);
                }
                if (isInside) {
                    // if (i < N) { // the N first particles may represent fluid, displayed in blue
                    //     image[((H - y - 1)*W + x) * 3] = 0;
                    //     image[((H - y - 1)*W + x) * 3 + 1] = 0;
                    //     image[((H - y - 1)*W + x) * 3 + 2] = 255;
                    // }
                    if (mindistEdge <= 2) {
                        image[((H - y - 1)*W + x) * 3] = 0;
                        image[((H - y - 1)*W + x) * 3 + 1] = 0;
                        image[((H - y - 1)*W + x) * 3 + 2] = 0;
                    }
                }
            }
        }
    }
    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}



// Voronoi diagrams
Vector voronoi_intersect(Vector &A, Vector &B, std::vector<Vector> &edge){
    Vector v1 = edge[0], v2 = edge[1];
    Vector M = (v1 + v2) / 2 ;
    double t = dot(M - A, v1 - v2)/dot(B - A, v1 - v2);
    if(t >= 0 && t <= 1){
        return A + t*(B - A);
    }
    return Vector();
}

bool voronoi_inside(Vector &P, std::vector<Vector> &edge){
    Vector v1 = edge[0], v2 = edge[1];
    Vector M = (v1 + v2) / 2;
    if (dot(P - M, v2 - v1) < 0){
        return true;
    } 
    return false;
}

std::vector<Polygon> voronoi_diagram(std::vector<Vector> &points, Polygon &frame){
    std::vector<Polygon> ans(points.size());
    Polygon outPolygon;
    for(int i = 0; i < points.size(); i++){
        Vector P1 = points[i];
        Polygon cell = frame;
        for(int j = 0; j < points.size(); j++){
            if(i != j){
                Vector P2 = points[j];
                std::vector<Vector> edge = {P1, P2};
                outPolygon = Polygon();
                for(int k = 0; k < cell.vertices.size(); k++){
                    Vector curVertex = cell.vertices[k];
                    Vector prevVertex = cell.vertices[(k > 0)? (k - 1) : (cell.vertices.size() - 1)];
                    Vector intersection = voronoi_intersect(prevVertex, curVertex, edge);

                    if(voronoi_inside(curVertex, edge)){
                        if(!voronoi_inside(prevVertex, edge)){
                            outPolygon.add_vertex(intersection);
                        }
                        outPolygon.add_vertex(curVertex);
                    } 

                    else if(voronoi_inside(prevVertex, edge)){
                        outPolygon.add_vertex(intersection);
                    }
                }
                cell = outPolygon;
            }
            ans[i] = cell;
        }
    }
    return ans;
}

// Power diagrams

Vector power_intersect(Vector &A, Vector &B, std::vector<Vector> &edge, std::vector<double> &weight){
    Vector v1 = edge[0], v2 = edge[1];
    Vector M = (v1 + v2) / 2 + (weight[0] - weight[1]) / (2*(v1 - v2).norm2()) * (v2 - v1);
    double t = dot(M - A, v1 - v2)/dot(B - A, v1 - v2);
    if (t >= 0 && t <= 1){
        return A + t*(B - A);
    }
    return Vector();
}

bool power_inside(Vector &P, std::vector<Vector> &edge, std::vector<double> &weight){
    Vector v1 = edge[0], v2 = edge[1];
    Vector M = (v1 + v2)/2 + (weight[0] - weight[1])/(2*(v1 - v2).norm2())*(v2 - v1);
    if (dot(P - M, v2 - v1) < 0){
        return true;
    } 
    return false;
}

std::vector<Polygon> power_diagram(std::vector<Vector> &points, Polygon &frame, std::vector<double> &weights){
    std::vector<Polygon> ans(points.size());
    for(int i = 0; i < points.size(); i++) {
        Vector P1 = points[i];
        Polygon cell = frame;
        for(int j = 0; j < points.size(); j++) {
            if(i != j){
                Vector P2 = points[j];
                std::vector<Vector> edge = {P1, P2};
                std::vector<double> weight = {weights[i], weights[j]};
                Polygon outPolygon;
                for(int k = 0; k < cell.vertices.size(); k++) {
                    Vector curVertex = cell.vertices[k];
                    Vector prevVertex = cell.vertices[(k > 0)? (k - 1) : (cell.vertices.size() - 1)];
                    Vector intersection = power_intersect(prevVertex, curVertex, edge, weight);

                    if(power_inside(curVertex, edge, weight)){
                        if(!power_inside(prevVertex, edge, weight)){
                            outPolygon.add_vertex(intersection);
                        }
                        outPolygon.add_vertex(curVertex);
                    } 
                    
                    else if(power_inside(prevVertex, edge, weight)){
                        outPolygon.add_vertex(intersection);
                    }
                }
                cell = outPolygon;
            }
        }
        ans[i] = cell;
    }
    return ans;
}

Polygon Circle(Vector C, double R){
    Polygon circle=Polygon();
    for(double theta=0; theta<2*M_PI; theta+=2*M_PI/100){
        double x=C[0] + R*cos(theta);
        double y=C[1] + R*sin(theta);
        x=std::min(1.,std::max(x,0.));
        y=std::min(1.,std::max(y,0.));
        Vector vertex=Vector(x,y);
        circle.add_vertex(vertex);
        //std::cout<<x<<" "<<y<<'\n';
    }
    return circle;
}

std::vector<Polygon> fluid_diagram(std::vector<Vector> &points, Polygon &frame, std::vector<double> &weights, double wair){
    std::vector<Polygon> ans(points.size());
    for(int i = 0; i < points.size(); i++) {
        Vector P1 = points[i];
        //std::cout<<"weight: "<<weights[i]<<" wair:"<<wair<<'\n'; 
        double R = sqrt(std::max(0.,weights[i] - wair));
        //std::cout<<"radius: "<<R<<'\n';
        Polygon cell = Circle(P1, R);
        for(int j = 0; j < points.size(); j++) {
            if(i != j){
                Vector P2 = points[j];
                std::vector<Vector> edge = {P1, P2};
                std::vector<double> weight = {weights[i], weights[j]};
                Polygon outPolygon;
                for(int k = 0; k < cell.vertices.size(); k++) {
                    Vector curVertex = cell.vertices[k];
                    Vector prevVertex = cell.vertices[(k > 0)? (k - 1) : (cell.vertices.size() - 1)];
                    Vector intersection = power_intersect(prevVertex, curVertex, edge, weight);

                    if(power_inside(curVertex, edge, weight)){
                        if(!power_inside(prevVertex, edge, weight)){
                            outPolygon.add_vertex(intersection);
                        }
                        outPolygon.add_vertex(curVertex);
                    } 
                    
                    else if(power_inside(prevVertex, edge, weight)){
                        outPolygon.add_vertex(intersection);
                    }
                }
                cell = outPolygon;
            }
        }
        ans[i] = cell;
    }
    return ans;
}
// LBFGS optimization based on provided sample: https://github.com/chokkan/liblbfgs/blob/master/sample/sample.cpp
class objective_function {
protected:
    lbfgsfloatval_t *m_x;
    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        ) 
    {
        return reinterpret_cast<objective_function*>(instance)->evaluate(x, g, n, step);
    }

    //Lecture 8 slide 34
    lbfgsfloatval_t evaluate(
            const lbfgsfloatval_t *x,
            lbfgsfloatval_t *g,
            const int n,
            const lbfgsfloatval_t step
        ) {
        lbfgsfloatval_t fx = 0;
        std::vector<double> weights(x, x + n);
        polygons = power_diagram(points, frame, weights);
        if(iterations == 10 || iterations == 100 || iterations == 1000){
            std::string filename = "imgs/iter_" + std::to_string(iterations) + ".svg";
            save_svg(polygons, filename);
        }
        for(int i = 0; i < n; i++){
            std::vector<Vector> vertices = polygons[i].vertices;
            size_t n = polygons[i].vertices.size();
            Vector point = points[i];
            double area = polygons[i].get_area();
            double tmp = 0;
            if(n > 0){
                Vector c1 = vertices[0];
                for(int i = 0; i < n - 2; i++){
                    Vector c2 = vertices[i + 1];
                    Vector c3 = vertices[i + 2];
                    tmp += (Polygon({c1,c2,c3}).get_area()/6) * (dot(c1 - point, c1 - point) + dot(c1 - point, c2 - point) + dot(c1 - point, c3 - point) + dot(c2 - point, c2 - point) +dot(c2 - point, c3 - point) +dot(c3 - point, c3 - point));
                }
            }
            fx += tmp - x[i]*area + lambdas[i]*x[i];
            g[i] = area - lambdas[i];
            //g[i] = 
        }
        iterations += 1;
        return -fx;
    }
    // lbfgsfloatval_t evaluate(
    //         const lbfgsfloatval_t *x,
    //         lbfgsfloatval_t *g,
    //         const int n,
    //         const lbfgsfloatval_t step
    //     ) {
    //     lbfgsfloatval_t fx = 0;
    //     std::vector<double> weights(x, x + (n-1));
    //     std::cout<<"weights size: "<<weights.size()<<'\n';
    //     for(int i=0; i<weights.size(); i++){
    //         std::cout<<weights[i]<<" ";
    //     }
    //     std::cout<<'\n';
    //     for(int i=0; i<points.size(); i++){
    //         std::cout<<points[i][0]<<" "<<points[i][1];
    //     }
    //     std::cout<<"\n wair: "<<x[n-1]<<'\n';
    //     std::cout<<'\n';
    //     polygons = fluid_diagram(points, frame, weights, x[n-1]);
    //     save_svg(polygons, "imgs/test.svg");
    //     //abort();
    //     std::cout<<"polygons: "<<polygons.size()<<'\n';
    //     for(int i=0; i<polygons.size(); i++){
    //         //std::cout<<polygons[i].vertices[0][0]<<" ";
    //     }
    //     std::cout<<'\n';
    //     if(iterations == 10 || iterations == 100 || iterations == 1000){
    //         //std::string filename = "imgs/optimized_" + std::to_string(iterations) + ".svg";
    //         //save_svg(polygons, filename);
    //     }
    //     for(int i = 0; i <n-1; i++){
    //         std::vector<Vector> vertices = polygons[i].vertices;
    //         size_t n = polygons[i].vertices.size();
    //         //std::cout<<"n: "<<n<<'\n';
    //         Vector point = points[i];
    //         double area = polygons[i].get_area();
    //         //std::cout<<"area: "<<area<<'\n';
    //         double tmp = 0;
    //         if(n > 0){
    //             Vector c1 = vertices[0];
    //             for(int i = 0; i < n - 2; i++){
    //                 Vector c2 = vertices[i + 1];
    //                 Vector c3 = vertices[i + 2];
    //                 tmp += (Polygon({c1,c2,c3}).get_area()/6) * (dot(c1 - point, c1 - point) + dot(c1 - point, c2 - point) + dot(c1 - point, c3 - point) + dot(c2 - point, c2 - point) +dot(c2 - point, c3 - point) + dot(c3 - point, c3 - point));
    //             }
    //         }
    //         fx += tmp - x[i]*area + lambdas[i]*x[i];
    //         //std::cout<<"fx: "<<fx<<'\n';
    //         g[i] = area - lambdas[i];
    //         std::cout<<"gradient: "<<g[i]<<'\n';
    //         //g[i] = 
    //     }
    //     double dair, eair=0.;
    //     for(int i=0; i<n-1; i++){
    //         //std::cout<<i<<'\n';
    //         eair += polygons[i].get_area();
    //     }
    //     eair=1-eair;
    //     std::cout<<"eair: "<<eair<<'\n';

    //     dair = 1 - lambdas[0]*(n-1);
    //     //std::cout<<n<<" "<<lambdas[0]<<'\n';
    //     //std::cout<<dair<<'\n';
    //     fx += x[n-1]*(dair-eair);
    //     //std::cout<<"fx: "<<fx<<'\n';
    //     g[n-1]=eair-dair;
    //     std::cout<<"gradient: "<<g[n-1]<<'\n';
    //     iterations += 1;
    //     return -fx;
    // }
    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls) 
    {
        return reinterpret_cast<objective_function*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        ) 
    {
        printf("Iteration %d:\n", k);
        printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
        printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
        printf("\n");
        return 0;
    }

    std::vector<double> lambdas;
    int iterations;
    std::vector<Vector> points;
    Polygon frame;
    std::vector<Polygon> polygons;
public:
    objective_function(std::vector<double> &lambdas, std::vector<Vector> &points,Polygon &frame) 
    {
        this->m_x = NULL;
        this->lambdas = lambdas;
        this->points = points;
        this->frame = frame;
        this->iterations = 0;
    }
    virtual ~objective_function(){
        if (m_x != NULL) {
            lbfgs_free(m_x);
            m_x = NULL;
        }
    }   
    std::vector<Polygon> run(int N)
    {
        lbfgsfloatval_t fx;
        lbfgsfloatval_t *m_x = lbfgs_malloc(N);

        if (m_x == NULL) {
            printf("ERROR: Failed to allocate a memory block for variables.\n");
            exit(1);
        }

        /* Initialize the variables. */
        for (int i = 0; i < N-1; i ++) {
            m_x[i] = 1;
        }
        m_x[N-1]=0;

        /*
            Start the L-BFGS optimization; this will invoke the callback functions
            evaluate() and progress() when necessary.
            */
        int ret = lbfgs(N, m_x, &fx, _evaluate, _progress, this, NULL);

        /* Report the result. */
        printf("L-BFGS optimization terminated with status code = %d\n", ret);
        printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, m_x[0], m_x[1]);

        return this->polygons;
    }
};


// Fluids

void simulation_step(
        std::vector<Vector> &positions,
        std::vector<Vector> &velocities,
        std::vector<double> &masses, int nbframes,
        Polygon &bounds, std::string filename,
        int N
    ) {
    double eps = 1e-5;
    double dt = 1e-5;
    Vector gravity{0, -0.1, 0};
    for (int i = 0; i < nbframes; i++) {        
        objective_function solver(
            masses, positions, bounds
        );
        std::vector<Polygon> fluid_polygons = solver.run(N+1);
        //std::cout<<fluid_polygons[0].vertices.size()<<'\n';
        for (size_t k = 0; k < N; k++) {
            std::cout<<fluid_polygons[k].vertices[0][0]<<'\n';
            Vector F_spring = (fluid_polygons[k].get_centroid() - positions[k])/eps;
            Vector F = F_spring + masses[k]*gravity;
            velocities[k] = velocities[k] + (dt/masses[k])*F;
            positions[k] = positions[k] + dt*velocities[k];
            if(positions[k][0]>1){
                positions[k][0]=1;
                velocities[k][0]=0;
            }
            if(positions[k][1]>1){
                positions[k][1]=1;
                velocities[k][1]=0;
            }
            if(positions[k][0]<0){
                positions[k][0]=0;
                velocities[k][0]=0;
            }
            if(positions[k][1]<0){
                positions[k][1]=0;
                velocities[k][1]=0;
            }
        }
        std::cout<<positions.size()<<'\n';
        for(int j=0; j<positions.size();j++){
            std::cout<<j<<" "<<positions[j][0]<<" "<<positions[j][1]<<'\n';
        }
        
    save_frame(fluid_polygons, filename, i, N);
    }
}

int main() {

    Polygon subjectPolygon = Circle(Vector(0.3, 0.3), 0.2);
    Polygon clipPolygon({
        Vector(0.3, 0.3), Vector(0.3, 0.7),
        Vector(0.7, 0.7), Vector(0.7, 0.3)
    });
    save_svg({subjectPolygon, clipPolygon}, "pics/before_clipping.svg");
    save_svg({clip(subjectPolygon, clipPolygon)}, "pics/after_clipping.svg");

    int n = 1000;
    Polygon bounds({
        Vector(0, 0), Vector(0, 1),
        Vector(1, 1), Vector(1, 0)
    });
    std::vector<Vector> points(n);
    for (int i = 0; i < n; i++) {
        double x = (double) rand() / RAND_MAX;
        double y = (double) rand() / RAND_MAX;
        points[i] = Vector(x, y, 0);
    }
    save_svg(voronoi_diagram(points, bounds), "pics/voronoi_diagram.svg");

    std::vector<double> weights(n);
    for (int i = 0; i < n; i++) {
        if (points[i][0] < 0.1 || points[i][0] > 0.9 || points[i][1] < 0.1 || points[i][1] > 0.9 ) {
            weights[i] = 0.99;
        } 
        else{
            weights[i] = 1;
        }
    }
    save_svg(power_diagram(points, bounds, weights), "pics/power_diagram.svg");

    std::vector<double> lambdas(n);
    double T;
    Vector C(0.5,0.5,0);
    for (int i = 0; i < n; i++) {
        lambdas[i] = std::exp(-(points[i] - C).norm2()/0.02);
        T += lambdas[i];
    }
    for (int i = 0; i < n; i++){
        lambdas[i] /= T;
    } 
    objective_function solver(lambdas, points, bounds);
    save_svg(solver.run(n), "pics/optimized_final.svg");

    // double desired = 0.02;
    // int N = 10; // fluid particles
    // std::vector<Vector> particles(N);
    // for (int i = 0; i < N; i++) {
    //     //std::cout<<i<<"\n";
    //     double x = (double) rand() / RAND_MAX;
    //     double y = (double) rand() / RAND_MAX;
    //     std::cout<<x<<" "<<y<<'\n';
    //     particles[i] = Vector(x, y, 0);
    // }
    
    // std::vector<double> lambdas(N, desired/N);
    // std::vector<Vector> velocities(N, Vector(0.,0.,0.));
    // simulation_step(
    //     particles, velocities, lambdas, 100, bounds, "imgs/frame", N
    // );
}
