#include "trace.H"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <getopt.h>
#include <math.h>
#ifdef __APPLE__
#define MAX std::numeric_limits<double>::max()
#else
#include <values.h>
#define MAX 999999999
#endif

// return the determinant of the matrix with columns a, b, c.
double det(const SlVector3 &a, const SlVector3 &b, const SlVector3 &c) {
    return a[0] * (b[1] * c[2] - c[1] * b[2]) +
        b[0] * (c[1] * a[2] - a[1] * c[2]) +
            c[0] * (a[1] * b[2] - b[1] * a[2]);
}

float max_f(double x, double y) {

    if (x >= y) {

        return x;

    }

    else {

        return y;

    }

}

inline double sqr(double x) {return x*x;} 

double multi(const SlVector3 &a, const SlVector3 &b){
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

bool Triangle::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {
    SlVector3 ba = a-b;
    SlVector3 ca = a-c;
    SlVector3 ea = a-r.e;
    double detA = det(ba, ca, r.d);
    double t = det(ba, ca, ea)/detA;
    if (t < t0 || t > t1) return false;
    double beta = det(ea, ca, r.d)/detA;
    if (beta < 0 || beta > 1) return false;
    double gamma = det(ba, ea, r.d)/detA;
    if (gamma < 0.0 || gamma > 1.0-beta) return false;
    hr.t = t;
    hr.p = r.e + t * r.d;
    hr.n = cross(ba,ca);
    normalize(hr.n);
    hr.alpha = 1.0 - beta - gamma;
    hr.beta = beta;
    hr.gamma = gamma;
    return true;
    

}

bool TrianglePatch::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {
    //std::cout<<"triangle"<<std::endl;
    bool temp = Triangle::intersect(r,t0,t1,hr);
    if (temp) {
        hr.n = hr.alpha * n1 + hr.beta * n2 + hr.gamma * n3;
        normalize(hr.n);
        return true;
    }
    return false;
}


bool Sphere::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {
    //std::cout<<"sphere"<<std::endl;
    double a1, a2, a3;
    double A = multi(r.d, r.d);
    double B = multi(r.d, (r.e-c));
    double C = multi((r.e-c),(r.e-c)) - sqr(rad);
    double delta = B*B-A*C;
    if(delta>=0){
        a1 = (-B+sqrt(delta))/A;
        a2 = (-B-sqrt(delta))/A;
        a3 = a1;
        if (a1 < 0 || (a2 > 0 && a2 < a1)){
            a3 = a2;
        }   
        if(a3 < t0 || a3 > t1){
            return false;
        }
        hr.t = a3;
        hr.p = r.e + a3 * r.d;
        hr.n = (hr.p - c) / rad;
        return true;
    }
    
    return false;
}




Tracer::Tracer(const std::string &fname) {
    std::ifstream in(fname.c_str(), std::ios_base::in);
    std::string line;
    char ch;
    Fill fill;
    bool coloredlights = false;
    while (in) {
        getline(in, line);
        switch (line[0]) {
            case 'b': {
                std::stringstream ss(line);
                ss>>ch>>bcolor[0]>>bcolor[1]>>bcolor[2];
                break;
            }

            case 'v': {
                getline(in, line);
                std::string junk;
                std::stringstream fromss(line);
                fromss>>junk>>eye[0]>>eye[1]>>eye[2];

                getline(in, line);
                std::stringstream atss(line);
                atss>>junk>>at[0]>>at[1]>>at[2];

                getline(in, line);
                std::stringstream upss(line);
                upss>>junk>>up[0]>>up[1]>>up[2];

                getline(in, line);
                std::stringstream angless(line);
                angless>>junk>>angle;

                getline(in, line);
                std::stringstream hitherss(line);
                hitherss>>junk>>hither;

                getline(in, line);
                std::stringstream resolutionss(line);
                resolutionss>>junk>>res[0]>>res[1];
                break;
            }

            case 'p': {
                bool patch = false;
                std::stringstream ssn(line);
                unsigned int nverts;
                if (line[1] == 'p') {
                    patch = true;
                    ssn>>ch;
                }
                ssn>>ch>>nverts;
                std::vector<SlVector3> vertices;
                std::vector<SlVector3> normals;
                for (unsigned int i=0; i<nverts; i++) {
                    getline(in, line);
                    std::stringstream ss(line);
                    SlVector3 v,n;
                    if (patch) ss>>v[0]>>v[1]>>v[2]>>n[0]>>n[1]>>n[2];
                    else ss>>v[0]>>v[1]>>v[2];
                    vertices.push_back(v);
                    normals.push_back(n);
                }
                bool makeTriangles = false;
                if (vertices.size() == 3) {
                    if (patch) {
                        surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2], 
                        normals [0], normals [1], normals [2]), fill));
                    } else {
                        surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]), fill));
                    }
                } else if (vertices.size() == 4) {
                    SlVector3 n0 = cross(vertices[1] - vertices[0], vertices[2] - vertices[0]);
                    SlVector3 n1 = cross(vertices[2] - vertices[1], vertices[3] - vertices[1]);
                    SlVector3 n2 = cross(vertices[3] - vertices[2], vertices[0] - vertices[2]);
                    SlVector3 n3 = cross(vertices[0] - vertices[3], vertices[1] - vertices[3]);
                    if (dot(n0,n1) > 0 && dot(n0,n2) > 0 && dot(n0,n3) > 0) {
                        makeTriangles = true;
                        if (patch) {
                            surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2], 
                            normals[0], normals[1], normals[2]), fill));
                            surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[2], vertices[3], 
                            normals[0], normals[2], normals[3]), fill));
                        } else {
                            surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]), fill));
                            surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[2], vertices[3]), fill));
                        }
                    }
                    if (!makeTriangles) {
                        std::cerr << "I didn't make triangles.  Poly not flat or more than quad.\n";
                    }
                }
                break;
            }

            case 's' : {
                std::stringstream ss(line);
                SlVector3 c;
                double r;
                ss>>ch>>c[0]>>c[1]>>c[2]>>r;
                surfaces.push_back(std::pair<Surface *, Fill>(new Sphere(c,r), fill));
                break;
            }
	  
            case 'f' : {
                std::stringstream ss(line);
                ss>>ch>>fill.color[0]>>fill.color[1]>>fill.color[2]>>fill.kd>>fill.ks>>fill.shine>>fill.t>>fill.ior;
                break;
            }

            case 'l' : {
                std::stringstream ss(line);
                Light l;
                ss>>ch>>l.p[0]>>l.p[1]>>l.p[2];
                if (!ss.eof()) {
                    ss>>l.c[0]>>l.c[1]>>l.c[2];
                    coloredlights = true;
                }
                lights.push_back(l);
                break;
            }

            default:
            break;
        }
    }
    if (!coloredlights) for (unsigned int i=0; i<lights.size(); i++) lights[i].c = 1.0/sqrt(lights.size());
    im = new SlVector3[res[0]*res[1]];
    shadowbias = 1e-6;
    samples = 1;
    aperture = 0.0;
}

Tracer::~Tracer() {
    if (im) delete [] im;
    for (unsigned int i=0; i<surfaces.size(); i++) delete surfaces[i].first;
}


SlVector3 Tracer::shade(const HitRecord &hr) const {
    //std::cout<<"shade"<<std::endl;
    if (color) return hr.f.color;

    SlVector3 color(0.0);
    SlVector3 pos = hr.p;
    SlVector3 LightPos(0.0);
    SlVector3 normal = hr.n;
    SlVector3 diffuse(0.0);
    SlVector3 reflection1(0.0);
    SlVector3 reflection2(0.0);
    SlVector3 viewer = hr.v;
    SlVector3 specular(0.0);
    SlVector3 surfaceTolight = (0.0);
    SlVector3 Totalcolor(0.0);
    SlVector3 Normalcolor = (0.1, 0.1, 0.1);
    bool reflected = false;
    HitRecord dummy;
    double t1 = MAX;
    HitRecord rhr;
    normalize(normal);
    for (unsigned int i=0; i<lights.size(); i++) {
        const Light &light = lights[i];
        bool shadow = false;

        LightPos = light.p;
        SlVector3 V1 = LightPos-pos;
        SlVector3 V2 = pos-LightPos;
        normalize(V1);
        normalize(V2);
        for(unsigned int j=0; j < surfaces.size(); j++){
            Ray surfaceToray = Ray(pos, V1);
            normalize(surfaceToray.d);
            if (surfaces[j].first->intersect(surfaceToray, 0.01, MAX, dummy)){
                shadow = true;
            }
        }

        reflection1 = 2*normal*multi(normal, V1) - V1;
        normalize(reflection1);

        if (!shadow) {
            diffuse = hr.f.kd * max_f(0,multi(V1,normal));
            specular = hr.f.ks * pow(max_f(0,multi(reflection1, viewer)),hr.f.shine);
            Totalcolor = hr.f.color*light.c;
            color = (diffuse + specular + Normalcolor) * Totalcolor + color;
        }
        else{
            color += Normalcolor * Totalcolor;
        }
    }
    reflection2 = 2*normal*multi(normal, viewer) - viewer;
    normalize(reflection2);
    Ray reflectionRay = Ray(hr.p, reflection2);
    reflectionRay.depth = hr.raydepth + 1;
    normalize(reflectionRay.d);
    for(unsigned int k=0; k < surfaces.size(); k++){
        if (surfaces[k].first->intersect(reflectionRay, hither, MAX, dummy)){
            reflected = true;
        }
    }
    if(reflected){
        SlVector3 reflectColor = trace(reflectionRay, hither, MAX);
        color = color + reflectColor * hr.f.ks;
    }

    return color;
}

SlVector3 Tracer::trace(const Ray &r, double t0, double t1) const {
    //std::cout<<"trace"<<std::endl;
    //std::cout<<r.depth<<std::endl;
    HitRecord hr;
    SlVector3 color(bcolor);
  
    bool hit = false;
   for (unsigned int k=0; k<surfaces.size(); k++) {
        if (surfaces[k].first->intersect(r, t0, t1, hr)) { 
            t1 = hr.t;
            hr.f = surfaces[k].second;
            hr.raydepth = r.depth;
            hr.v = r.e - hr.p;
            normalize(hr.v);
            normalize(hr.n);
            hit = true;
        }
    }

    if (hit && hr.raydepth < 3) {
        color = shade(hr);
    }
    return color;
}

void Tracer::traceImage() {
    //std::cout<<"image"<<std::endl;
    // set up coordinate system
    SlVector3 w = eye - at;
    w /= mag(w);
    SlVector3 u = cross(up,w);
    normalize(u);
    SlVector3 v = cross(w,u);
    normalize(v);

    double d = mag(eye - at);
    double h = tan((M_PI/180.0) * (angle/2.0)) * d;
    double l = -h;
    double r = h;
    double b = h;
    double t = -h;

    SlVector3 *pixel = im;

    for (unsigned int j=0; j<res[1]; j++) {
        for (unsigned int i=0; i<res[0]; i++, pixel++) {

            SlVector3 result(0.0,0.0,0.0);

            for (int k = 0; k < samples; k++) {

                double rx = 1.1 * rand() / RAND_MAX;
                double ry = 1.1 * rand() / RAND_MAX;

                double x = l + (r-l)*(i+rx)/res[0];
                double y = b + (t-b)*(j+ry)/res[1];
                SlVector3 dir = -d * w + x * u + y * v;
	
                Ray r(eye, dir);
                normalize(r.d);
                //std::cout<<"ijk:"<<i<<" "<<j<<" "<<k<<std::endl;
                result += trace(r, hither, MAX);

            }
            (*pixel) = result / samples;
        }
    }
}

void Tracer::writeImage(const std::string &fname) {
#ifdef __APPLE__
    std::ofstream out(fname, std::ios::out | std::ios::binary);
#else
    std::ofstream out(fname.c_str(), std::ios_base::binary);
#endif
    out<<"P6"<<"\n"<<res[0]<<" "<<res[1]<<"\n"<<255<<"\n";
    SlVector3 *pixel = im;
    char val;
    for (unsigned int i=0; i<res[0]*res[1]; i++, pixel++) {
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[0])) * 255.0);
        out.write (&val, sizeof(unsigned char));
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[1])) * 255.0);
        out.write (&val, sizeof(unsigned char));
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[2])) * 255.0);
        out.write (&val, sizeof(unsigned char));
    }
    out.close();
}


int main(int argc, char *argv[]) {
    int c;
    double aperture = 0.0;
    int samples = 1;
    int maxraydepth = 5;
    bool color = false;
    while ((c = getopt(argc, argv, "a:s:d:c")) != -1) {
        switch(c) {
            case 'a':
            aperture = atof(optarg);
            break;
            case 's':
            samples = atoi(optarg);
            break;
            case 'c':
            color = true;
            break;
            case 'd':
            maxraydepth = atoi(optarg);
            break;
            default:
            abort();
        }
    }

    if (argc-optind != 2) {
        std::cout<<"usage: trace [opts] input.nff output.ppm"<<std::endl;
        for (unsigned int i=0; i<argc; i++) std::cout<<argv[i]<<std::endl;
        exit(0);
    }

    Tracer tracer(argv[1]);
    tracer.aperture = aperture;
    tracer.samples = samples;
    tracer.color = color;
    tracer.maxraydepth = maxraydepth;
    tracer.traceImage();
    tracer.writeImage(argv[2]);
};
