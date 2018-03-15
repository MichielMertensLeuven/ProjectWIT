#include <iostream>
#include <IpIpoptApplication.hpp>
#include "generate_matrices.hpp"
#include "optimization.hpp"

int main(int argc, char** argv){
	const int n=100; // Amount of grid points in each direction
	
	Ipopt::SmartPtr<Ipopt::TNLP> mynlp=new HeatNLP<n>();
	SmartPtr<IpoptApplication> app=IpoptApplicationFactory();
	app->Options()->SetNumericValue("tol",1e-12);
	app->Options()->SetStringValue("mu_strategy","adaptive");
	app->Options()->SetStringValue("output_file","ipopt.out");
	app->Options()->SetStringValue("hessian_approximation","limited-memory");
	Ipopt::ApplicationReturnStatus status=app->Initialize();
	if(status!=Solve_Succeeded){
		printf("\n\n*** Error during initialization!\n");
		return (int) status;
	}
	status=app->OptimizeTNLP(mynlp);
	
	if(status==Ipopt::Solve_Succeeded){
		printf("\n\n*** The problem solved!\n");
	}else{
		printf("\n\n*** The problem FAILED!\n");
	}
	return 0;
}
