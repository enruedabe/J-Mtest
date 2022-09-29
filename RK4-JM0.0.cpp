#include <iostream>
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <vector>

using alg=double (double, double, double, double, double);

double dif1(double x, double y, double z, double w, double t);
double dif2(double x, double y, double z, double w, double t);
double dif3(double x, double y, double z, double w, double t);
double dif4(double x, double y, double z, double w, double t);
double A(double t);
void RK(double x0, double y0, double z0, double w0, double t0, int N);

const double h=0.001;
const double omega=1.01;
const double phi=0.695;
const double f=0.183;
const double g=0.1965;
const double k=0.07;
const double gamma=0.0009;
const double omega_a=1.0;
const double omega_s=1.1;





int main(){
    int N=1000000;
    double x0=4.0;
    double y0=0.8;
    double z0=2.0;
    double w0=0.4;
    double t0=0;

    RK(x0, y0, z0, w0, t0, N);
    return 0;
}

double dif1(double x, double y, double z, double w, double t){
    return phi*(1+f*std::cos(A(t))+omega_a*z)+g*w-0.5*k*x;
}
double dif2(double x, double y, double z, double w, double t){
    return omega_s*w-g*z+2*g*z*y*y+2*g*z*w*w-0.5*gamma*y;
}
double dif3(double x, double y, double z, double w, double t){
    return -omega_a*x-g*y-0.5*k*z;
}
double dif4(double x, double y, double z, double w, double t){
    return -omega_s*w+g*x-2*g*x*y*y-2*g*x*w*w-0.5*gamma*w;
}
double A(double t){
    return omega*t;
}

void RK(double x0, double y0, double z0, double w0, double t0, int N){
    double t=t0;
    double x=x0;
    double y=y0;
    double z=z0;
    double w=w0;
    std::vector<double>data_x(N);
    std::vector<double>data_y(N);
    std::vector<double>data_z(N);
    std::vector<double>data_w(N);
    double k1, k2, k3, k4, l1, l2, l3, l4, m1, m2, m3, m4, n1, n2, n3, n4;
        for(int i=0;i<N;i++){
        k1=dif1(x,y,z,w,t);
        l1=dif2(x,y,z,w,t);
        m1=dif3(x,y,z,w,t);
        n1=dif4(x,y,z,w,t);

        k2=dif1(x+0.5*h*k1, y+0.5*h*l1, z+0.5*h*m1, w+0.5*h*n1, t+0.5*h);
        l2=dif2(x+0.5*h*k1, y+0.5*h*l1, z+0.5*h*m1, w+0.5*h*n1, t+0.5*h);
        m2=dif3(x+0.5*h*k1, y+0.5*h*l1, z+0.5*h*m1, w+0.5*h*n1, t+0.5*h);
        n2=dif4(x+0.5*h*k1, y+0.5*h*l1, z+0.5*h*m1, w+0.5*h*n1, t+0.5*h);


        k3=dif1(x+0.5*h*k2, y+0.5*h*l2, z+0.5*h*m2, w+0.5*h*n2, t+0.5*h);
        l3=dif2(x+0.5*h*k2, y+0.5*h*l2, z+0.5*h*m2, w+0.5*h*n2, t+0.5*h);
        m3=dif3(x+0.5*h*k2, y+0.5*h*l2, z+0.5*h*m2, w+0.5*h*n2, t+0.5*h);
        n3=dif4(x+0.5*h*k2, y+0.5*h*l2, z+0.5*h*m2, w+0.5*h*n2, t+0.5*h);


        k4=dif1(x+h*k3, y+h*l3, z+h*m3, w+h*n3, t+h);
        l4=dif2(x+h*k3, y+h*l3, z+h*m3, w+h*n3, t+h);
        m4=dif3(x+h*k3, y+h*l3, z+h*m3, w+h*n3, t+h);
        n4=dif4(x+h*k3, y+h*l3, z+h*m3, w+h*n3, t+h);

        data_x[i]=x;
        data_y[i]=y;
        data_z[i]=z;
        data_w[i]=w;

        t=t+h;

        x=x+h*(k1+2*k2+2*k3+k4)/6;
        y=y+h*(l1+2*l2+2*l3+l4)/6;
        z=z+h*(m1+2*m2+2*m3+m4)/6;
        w=w+h*(n1+2*n2+2*n3+n4)/6;
    }
    std::ofstream outfile;

    outfile.open("data_x.txt");
    for(int i=0; i<N; i++){
        outfile<<data_x[i]<<std::endl;
    }
    outfile.close();

    outfile.open("data_y.txt");
    for(int i=0; i<N; i++){
        outfile<<data_y[i]<<std::endl;
    }
    outfile.close();

    outfile.open("data_z.txt");
    for(int i=0; i<N; i++){
        outfile<<data_z[i]<<std::endl;
    }
    outfile.close();

    outfile.open("data_w.txt");
    for(int i=0; i<N; i++){
        outfile<<data_w[i]<<std::endl;
    }
    outfile.close();

}