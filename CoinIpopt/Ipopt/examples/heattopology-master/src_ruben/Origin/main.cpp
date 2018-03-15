#include <dolfin.h>
#include "MixedPoissonDual.h"

using namespace dolfin;

double V=0.4;		//Volume bound
double k_plastic=0.1;	//Thermal conductivity of plastic in Watt per meter Kelvin
double k_metal=80;	//Thermal conductivity of metal in Watt per meter Kelvin
double size_x=0.01;	//Plate size in m
double size_y=0.01;	//Plate size in m
double size_z=0.001;	//Plate size in m
double T_out=300;	//Output temperature in Kelvin
double Q_in=2.5;	//Input heat in Watt
int p=5;		//Penalization factor
int n=250;		//Mesh size

//Source term
class Source:public Expression{
	void eval(Array<double>& values,const Array<double>& x)const{
		values[0]=Q_in/(size_x*size_y);
	}
};

//Boundary source for Neumann boundary condition
class BoundarySource:public Expression{
	void eval(Array<double>& values,const Array<double>& x)const{
		values[0]=0;
	}
};

//Subdomain for Dirichlet boundary condition
class DirichletBoundary:public SubDomain{
	bool inside(const Array<double>& x,bool on_boundary)const{
		return x[0]>size_x-DOLFIN_EPS&&x[1]<=0.7*size_y+DOLFIN_EPS&&x[1]>=0.3*size_y-DOLFIN_EPS;
	}
};

int main(){
	// Create mesh
	UnitSquareMesh mesh(n,n);
	
	// Construct function space
	MixedPoissonDual::FunctionSpace W(mesh);
	MixedPoissonDual::BilinearForm a(W,W);
	MixedPoissonDual::LinearForm L(W);
	
	// Create sources and assign to L
	Source f;
	BoundarySource g;
	L.f=f;
	L.g=g;
	
	// Define boundary condition
	SubSpace W1(W,1);
	DirichletBoundary boundary;
	DirichletBC bc(W1,T_out,boundary);
	
	// Compute solution
	Function w(W);
	solve(a==L,w,bc);
	
	// Extract sub functions (function views)
	Function& sigma=w[0];
	Function& u=w[1];
	
	// Plot solutions
	plot(u);
	plot(sigma);
	interactive();
	
	return 0;
}
