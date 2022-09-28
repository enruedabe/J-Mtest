#include <cmath>


double h(double& a_r,double& a_i,double& sigma_i, double& t, const double phi, const double f, const double Omega, const double omega_a,const double g, const double k); //dar/dt
double i(double& sigma_r,double& a_i,double& sigma_i, double& t,const double omega_sigma, const double g,const double gamma); //d(sigma_r)/dt
double j(double& a_r,double& a_i,double& sigma_i, double& t, const double omega_a, const double g, const double k); // d(ai)/dt
double k(double& a_r,double& a_i,double& sigma_i, double& t, const double omega_sigma, const double g, const double gamma); //d(sigma_i)/dt
/*d(a_r)/dt = (\phi)(1+fcos(lambda)) + (omega_a)*a_i + g (sigma_i) - k/2(a_r)
		  d(sigma_r)/dt = (omega_sigma)(sigma_i) - g(a_i)- 2g(a_i)(sigma_r)^2 + 2g(a_i)(sigma_i)^2 - (gamma)/2 (sigma_r)
		  d(a_i)/dt = -(omega_a)(a_r) - g(sigma_r) - k/2 (a_i)
		  d(sigma_i)/dt = -(omega_sigma)(sigma_i) + g(a_r) - 2g(a_r)(sigma_r)^2 - 2g(a_r)(sigma_i)^2 - (gamma)/2 (sigma_i)
		  d(lambda)/dt = Omega => lambda=Omega*t
*/
void RK4(std::vector<double>& t, std::vector<double>& a_r, std::vector<double>& a_i, std::vector<double>& sigma_r,std::vector<double>& sigma_i,const double& h, const int& n); 



int main(int argc, char const *argv[])
{
	//Consts
	double f=0;
	double g=0;
	double k=0;
	double phi=0;
	double gamma=0;
	double Omega=0;
	double omega_a=0;
	double omega_sigma=0;

		
	int n= 100; //steps
	double h=M_PI/5000; //step size
	std::vector<double> t;
	std::vector<double> a_r;
	std::vector<double> a_i;
	std::vector<double> sigma_r;
	std::vector<double> sigma_i;

	//ini cond
	t[0]=0.0;
	a_r[0]=4.0;
	a_i[0]=4.0;
	sigma_r[0]=0.2;
	sigma_i[0]=0.2;
	




	return 0;
}


double h(double& a_r,double& a_i,double& sigma_i, double& t, const double phi, const double f, const double Omega, const double omega_a,const double g, const double k){

	return phi*(1.0+f*std::cos(Omega*t)) + omega_a*a_i + g*sigma_i - (k/2.0)*(a_r);
}//ar

double i(double& sigma_r,double& a_i,double& sigma_i, double& t,const double omega_sigma, const double g,const double gamma){

	return omega_sigma*sigma_i - g*a_i + 2.0*g*a_i*std::pow(sigma_r,2.0) + (2.0)*g*a_i*std::pow(sigma_i,2.0) - (gamma/2.0)*(sigma_r);
}//sigmar

double j(double& a_r,double& a_i,double& sigma_i, double& t, const double omega_a, const double g, const double k){

	return (-1.0)*omega_a*a_r - g*sigma_r - (k/2.0)*a_i;
}//ai

double k(double& a_r,double& a_i,double& sigma_i, double& t, const double omega_sigma, const double g, const double gamma){

	return (-1.0)*omega_sigma*sigma_i + g*a_r - (2.0)*g*a_r*std::pow(sigma_r,2.0) - (2.0)*g*a_r*std::pow(sigma_i,2.0) - (gamma/2.0)*(sigma_i);
}//sigmai

void RK4(std::vector<double>& t, std::vector<double>& a_r, std::vector<double>& a_i, std::vector<double>& sigma_r,std::vector<double>& sigma_i,const double& h, const int& n){

	double H1, H2, H3, H4, I1, I2, I3, I4, J1, J2, J3, J4, K1, K2, K3, K4=0;
	for (int i = 0; i < n; ++i)
	{
		H1 = h(a_r[i],a_i[i],sigma_i[i],t[i],phi,f,Omega,omega_a,g,k);
		I1 = i(sigma_r[i],a_i[i],sigma_i[i],t[i],omega_sigma,g,gamma);
		J1 = j(a_r[i],a_i[i],sigma_i[i],t[i],omega_a,g,k);
		K1 = k(a_r[i],a_i[i],sigma_i[i],t[i],omega_sigma,g,gamma);

		H2 = h(a_r[i]+H1*(h/2.0),a_i[i]+J1*(h/2.0),sigma_i[i]+K1*(h/2.0),t[i]+(h/2.0),phi,f,Omega,omega_a,g,k);
		I2 = i(sigma_r[i]+I1*(h/2.0),a_i[i]+J1*(h/2.0),sigma_i[i]+K1*(h/2.0),t[i]+(h/2.0),omega_sigma,g,gamma);
		J2 = j(a_r[i]+H1*(h/2.0),a_i[i]+J1*(h/2.0),sigma_i[i]+K1*(h/2.0),t[i]+(h/2.0),omega_a,g,k);
		K2 = k(a_r[i]+H1*(h/2.0),a_i[i]+J1*(h/2.0),sigma_i[i]+K1*(h/2.0),t[i]+(h/2.0),omega_sigma,g,gamma);

		
	}
}

