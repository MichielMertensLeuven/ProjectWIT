#include <iostream>
#include <IpIpoptApplication.hpp>
#include "generate_matrices.hpp"
#include "optimization_functions.hpp"

#define EIGEN_MALLOC_ALREADY_ALIGNED = 0

// How to run:
// g++ --std=c++14 -I /localhost/packages/math_eng_soft/CoinIpopt/build/include/coin -I Eigen/ -L/localhost/packages/math_eng_soft/CoinIpopt/build/lib -o test test_heat.cpp print.cpp -lipopt -lm -ldl -g3 && ./test && xdg-open plate_heat.ppm && xdg-open plate_material.ppm && notify-send "done"
// export LD_LIBRARY_PATH=/localhost/packages/math_eng_soft/CoinIpopt/build/lib:$LD_LIBRARY_PATH
// ./test



int main(int argc, char** argv){
	const int n=50; // Amount of grid points in each direction
	
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
