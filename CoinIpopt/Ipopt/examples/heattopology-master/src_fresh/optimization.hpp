#ifndef optimization_hpp
#define optimization_hpp

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iostream>
#include <random>
#include "print.hpp"
#include "generate_matrices.hpp"
#include "IpIpoptApplication.hpp"

using namespace Ipopt;

template <int N>
class HeatNLP : public TNLP{
	Eigen::SparseMatrix<double> S;
	Eigen::VectorXd u;
	Eigen::MatrixXd k;
	double size=0.01; // Size of the plate in x and y direction in meters
	int p=4; // Penalization factor
	double k_plastic=0.1; // Thermal conductivity of plastic in Watt per Kelvin per meter per meter
	double k_steel=80; // Thermal conductivity of metal in Watt per Kelvin per meter per meter
	double T=300; // External temperature in Kelvin
	double Q=2.5; // Input heat in Watt
	
public:
	
	// constructor
	HeatNLP(){};
	
	//destructor
	~HeatNLP(){};
	
	// returns the size of the problem
	bool get_nlp_info(Index& nn, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style){
		nn=N*N; // The problem N^2 variables
		m=1; // one inequality constraint
		nnz_jac_g=N*N; // jacobian is dense and contains N^2 nonzeros
		index_style=TNLP::C_STYLE;// use the C style indexing (0-based)
		return true;
	}
	
	// returns the variable bounds
	bool get_bounds_info(Index nn, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u){
		for(Index i=0;i<nn;++i){
			x_l[i]=0.0; // the variables have lower bounds of 0
			x_u[i]=1.0; // the variables have upper bounds of 1
		}
		
		g_l[0]=0.39; // the first constraint g1 has a lower bound of 0
		g_u[0]=0.4; // the first constraint g1 has an upper bound of 0.4
		
		return true;
	}
	
	// returns the initial point for the problem
	bool get_starting_point(Index nn, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda){
		assert(init_x==true);
		assert(init_z==false);
		assert(init_lambda==false);
		
		double lower_bound=0.39;
		double upper_bound=0.4;
		std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
		std::default_random_engine re;
		
		for(int i=0;i<nn;i++){ // initialize to the given starting point
			x[i]=unif(re);
		}
		Eigen::MatrixXd A(N,N);
		for(int i=0;i<N;++i){
			A.row(i)=Eigen::VectorXd::Map(&x[i*N],N);
		}
		double material;
		for(int i=0;i<nn;i++){
			material+=x[i];
		}
		material=material/(nn);
		
		wit::print_material(A);
		return true;
	}
	
	// returns the value of the objective function
	bool eval_f(Index nn, const Number* x, bool new_x, Number& obj_value){
		assert(nn==N*N);
		//if(!new_x)std::cout<<"X did not change"<<std::endl;
		Eigen::MatrixXd A(N,N);
		for(int i=0;i<N;++i){
			A.row(i)=Eigen::VectorXd::Map(&x[i*N],N);
		}
		k=generate_k(N,A,p,k_plastic,k_steel);
		S=generate_s(N,size,k);
		u=generate_u(N,size,p,k_plastic,k_steel,T,Q,A);
		double w=0.5; // weight parameter 0<=w<=1 to optimise objective function towards minimal average or variance. Also to be changed in generate_matrices.hpp
		double J1=u.dot(Eigen::VectorXd::Ones(nn))/nn;
		double J21=u.dot(u)/nn;
		double J22=(u.dot(Eigen::VectorXd::Ones(nn))*Eigen::VectorXd::Ones(nn).dot(u))/nn/nn;
		double obj_value_double=w*J1+(1-w)*(J21-J22);
		obj_value=(Number)obj_value_double;
		double material;
		for(int i=0;i<nn;i++){
			material+=x[i];
		}
		material=material/(nn);
		std::cout<<"obj:	"<<obj_value<<"	; J1:	"<<J1<<"	; J21-J22:	"<<J21-J22<<"	; Avg material:	"<<material<<"	; Max temp:	"<<u.maxCoeff()<<std::endl;
		return true;
	}
	
	// return the gradient of the objective function grad_{x} f(x)
	bool eval_grad_f(Index nn, const Number* x, bool new_x, Number* grad_f){
		assert(nn==N*N);
		Eigen::MatrixXd A(N,N);
		for(int i=0;i<N*N;i+=N){
			A.row((int)i/N)=Eigen::VectorXd::Map(&x[i],N);
		}
		k=generate_k(N,A,p,k_plastic,k_steel);
		if(S.nonZeros()==0){
			S=generate_s(N,size,k);
			u=generate_u(N,size,p,k_plastic,k_steel,T,Q,A);
		}
		grad_f=generate_dydp(N,size,k,A,p,k_plastic,k_steel,u,generate_s(N,size,k)).data();
		return true;
	}
	
	// return the value of the constraints: g(x)
	bool eval_g(Index nn, const Number* x, bool new_x, Index m, Number* g){
		assert(nn==N*N);
		assert(m==1);
		
		for(int i=0;i<nn;i++){
			g[0]+=x[i];
		}
		g[0]=g[0]/(nn);
		return true;
	}
	
	// return the structure or values of the jacobian of g
	bool eval_jac_g(Index nn, const Number* x, bool new_x, Index m, Index n_elem_jac, Index* iRow, Index *jCol, Number* values){
		if(values==NULL){
			for(int i=0;i<nn;i++){
				iRow[i]=0;
				jCol[i]=i;
			}
		}else{
			for(int i=0;i<nn;i++){
				values[i]=1/(nn);
			}
		}
		assert(nn==n_elem_jac);
		return true;
	}
	
	//return the structure or values of the hessian of the lagrangian:
	// sigma_f * hes(f(x)) + sum_i(lambda_i*g_i(x))
	bool eval_h(Index nn, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values){
		return false;
	}
	
	void finalize_solution(SolverReturn status, Index nn, const Number* x, const Number* z_L, const Number* z_U, Index m, const Number* g, const Number* lambda, Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq){
		if(status!=SUCCESS){
			std::cout<<std::endl<<"ended with non-succes status: "<<status<<std::endl;
		}
		Eigen::MatrixXd A(N,N);
		for(int i=0;i<N;++i){
			A.row(i)=Eigen::VectorXd::Map(&x[i*N],N);
		}
		k=generate_k(N,A,p,k_plastic,k_steel);
		u=generate_u(N,size,p,k_plastic,k_steel,T,Q,A);
		//u=generate_dydp(N,size,k,A,p,k_plastic,k_steel,generate_u(N,size,p,k_plastic,k_steel,T,Q,A),generate_s(N,size,k));
		Eigen::MatrixXd U(N,N); // Reform solution to matrix, for easy viewing and understanding
		for(int i=0;i<N;++i){
			for(int j=0;j<N;++j){
				U.coeffRef(i,j)=u(i+N*j);
			}
		}
		wit::print_heat(U);
		wit::print_material(A);
		double material;
		for(int i=0;i<nn;i++){
			material+=x[i];
		}
		material=material/(nn);
		std::cout<<"Average temperature in Kelvin: "<<(Eigen::VectorXd::Ones(nn).transpose()*u)*1/nn<<std::endl;
		std::cout<<"Temperature variance : "<<(u.dot(u))*1/nn-(u.dot(Eigen::VectorXd::Ones(nn))*Eigen::VectorXd::Ones(nn).dot(u))*1/nn/nn<<std::endl;
		std::cout<<"Maximum temperature in Kelvin: "<<u.maxCoeff()<<std::endl;
		std::cout<<"Percentage of metal surface: "<<material<<std::endl;
	}
};

#endif
