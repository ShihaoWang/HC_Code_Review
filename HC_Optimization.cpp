#include "HC_Header.h"
#include <string>
#include <math.h>
#include <iostream>
#include <time.h>
#include <cmath>
#include <stdlib.h>     /* srand, rand */
#include <dlib/matrix.h>
#include <dlib/optimization.h>

using namespace std;
using namespace dlib;

extern RobotState Robot_State_In_Solver;

NormalDistribution Impact_Gauss(Impact_Avg, Impact_Std);
NormalDistribution Mu_A_Gauss(Mu_A_Avg, Mu_A_Std);
NormalDistribution Mu_E_Gauss(Mu_D_Avg, Mu_D_Std);

// extern NormalDistribution Impact_Gauss;
// extern NormalDistribution Mu_A_Gauss;
// extern NormalDistribution Mu_E_Gauss;

extern Mu_Database Mu_A_Lib;
extern Mu_Database Mu_D_Lib;

// The cumulative density function and probability density function of a normal distribution

std::vector<double> Hand_Contact_Optimization(RobotState &Initial_State_i, Optimal_Control_Class &Optimal_Control_Seq)
{
	// This the main function to conduct the optimization to get the optimal control sequence

	// The first step is to conduct an integration based on the current initial condition
	Collision_Struct  Collision_Condition_i;

	Robot_State_In_Solver = Initial_State_i;  //optimal_initial_condition_in_solver will be transmitted to the solver

	// Here the initial control sequence has already been initialized
	Init_Opt_Ctrl_Generator(Initial_State_i, Optimal_Control_Seq, Collision_Condition_i);

	// This is the solver to find the pre-impact optimal control sequence
	find_min_box_constrained(bfgs_search_strategy(), objective_delta_stop_strategy(1e-9, 5).be_verbose(), myObjfn, derivative(myObjfn), Optimal_Control_Seq, -0.1, 0.1);

	std::vector<double> Optimal_Control_Seq_Vec;
	for (int i = 0; i < 2*(grids-1); i++)
	{
		Optimal_Control_Seq_Vec.push_back(Optimal_Control_Seq(i));
	}
	return Optimal_Control_Seq_Vec;
}

void Init_Opt_Ctrl_Generator(RobotState optimal_initial_condition, Optimal_Control_Class &Init_Ctrl_Sequence, Collision_Struct &Collision_Condition_i)
{
	// This function is used to generate the initial control sequence used for the later optimization process

	double kP_alpha = 250; double kD_alpha = 20;
	double kP_beta = 250;  double kD_beta = 20;

	double ddt = duration_time/(grids*1.0 - 1.0);

	// For the vertical wall and cubic parabola case
	// These two values are chosen heuristically
	double alpha_ref = 5 * dlib::pi  / 6;
	double beta_ref = 5 * dlib::pi  / 8;

	//// For the flat contact case
	//double alpha_ref = dlib::pi  / 2.0;
	//double beta_ref = 3.6*dlib::pi  /5.0;

	RobotState Pre_Impact_State_i = optimal_initial_condition;

	for (int i = 0; i < grids-1; i++)
	{
		double Theta =    Pre_Impact_State_i.theta;
		double Alpha =    Pre_Impact_State_i.alpha;
		double Beta =     Pre_Impact_State_i.beta;
		double Thetadot = Pre_Impact_State_i.thetadot;
		double Alphadot = Pre_Impact_State_i.alphadot;
		double Betadot =  Pre_Impact_State_i.betadot;

		double R_00 = ((cos(Beta*2.0)*4.543259751217974E36 + cos(Alpha)*1.305781530283095E38 - cos(Alpha + Beta*2.0)*8.112963841460668E36 - 1.399373244558452E39)*(cos(Alpha + Theta)*(-2.1672252E-1) + cos(Alpha + Beta + Theta)*2.4525E-2 + cos(Theta)*7.864554375E-1 - (Alphadot*Alphadot)*sin(Alpha + Beta)*3.125E-4 - (Betadot*Betadot)*sin(Alpha + Beta)*3.125E-4 + (Alphadot*Alphadot)*sin(Alpha)*2.7615E-3 + (Betadot*Betadot)*sin(Beta)*1.75E-4 - Alphadot*Betadot*sin(Alpha + Beta)*6.25E-4 - Alphadot*Thetadot*sin(Alpha + Beta)*6.25E-4 - Betadot*Thetadot*sin(Alpha + Beta)*6.25E-4 + Alphadot*Betadot*sin(Beta)*3.5E-4 + Alphadot*Thetadot*sin(Alpha)*5.523E-3 + Betadot*Thetadot*sin(Beta)*3.5E-4)*8.070450532247929E22) / (cos(Alpha*2.0)*1.36466322867013E58 + cos(Beta*2.0)*1.810096145662263E57 + cos(Alpha*2.0 + Beta*2.0)*3.087528119013518E58 - 1.068113514259467E60) - ((cos(Alpha + Beta + Theta)*(-2.4525E-2) - (Thetadot*Thetadot)*sin(Alpha + Beta)*3.125E-4 + (Alphadot*Alphadot)*sin(Beta)*1.75E-4 + (Thetadot*Thetadot)*sin(Beta)*1.75E-4 + Alphadot*Thetadot*sin(Beta)*3.5E-4)*(cos(Alpha + Beta)*2.322293885282703E41 + cos(Alpha)*1.24267240715314E40 + cos(Beta)*4.590862064203024E40 - cos(Alpha + Beta)*cos(Alpha)*2.294164443975028E40 - cos(Alpha + Beta)*cos(Beta)*1.453843120389752E39 + pow(cos(Alpha + Beta), 2.0)*2.596148429267414E39 - cos(Alpha)*cos(Beta)*1.284732088626016E40 - 4.440565984720187E40)*9.007199254740992E20) / (cos(Alpha*2.0)*1.36466322867013E58 + cos(Beta*2.0)*1.810096145662263E57 + cos(Alpha*2.0 + Beta*2.0)*3.087528119013518E58 - 1.068113514259467E60) - ((cos(Alpha)*1.24267240715314E41 + pow(cos(Beta), 2.0)*4.070760737091305E39 - cos(Alpha + Beta)*cos(Beta)*1.453843120389752E40 + pow(cos(Alpha + Beta), 2.0)*1.298074214633707E40 - 8.509828931667414E41)*(cos(Alpha + Theta)*(-2.1672252E-1) + cos(Alpha + Beta + Theta)*2.4525E-2 + (Thetadot*Thetadot)*sin(Alpha + Beta)*3.125E-4 + (Betadot*Betadot)*sin(Beta)*1.75E-4 - (Thetadot*Thetadot)*sin(Alpha)*2.7615E-3 + Alphadot*Betadot*sin(Beta)*3.5E-4 + Betadot*Thetadot*sin(Beta)*3.5E-4)*1.801439850948198E20) / (cos(Alpha*2.0)*1.36466322867013E58 + cos(Beta*2.0)*1.810096145662263E57 + cos(Alpha*2.0 + Beta*2.0)*3.087528119013518E58 - 1.068113514259467E60);
		double R_10 = ((cos(Alpha + Theta)*(-2.1672252E-1) + cos(Alpha + Beta + Theta)*2.4525E-2 + (Thetadot*Thetadot)*sin(Alpha + Beta)*3.125E-4 + (Betadot*Betadot)*sin(Beta)*1.75E-4 - (Thetadot*Thetadot)*sin(Alpha)*2.7615E-3 + Alphadot*Betadot*sin(Beta)*3.5E-4 + Betadot*Thetadot*sin(Beta)*3.5E-4)*(cos(Alpha + Beta)*2.322293885282703E41 + cos(Alpha)*1.24267240715314E40 + cos(Beta)*4.590862064203024E40 - cos(Alpha + Beta)*cos(Alpha)*2.294164443975028E40 - cos(Alpha + Beta)*cos(Beta)*1.453843120389752E39 + pow(cos(Alpha + Beta), 2.0)*2.596148429267414E39 - cos(Alpha)*cos(Beta)*1.284732088626016E40 - 4.440565984720187E40)*9.007199254740992E20) / (cos(Alpha*2.0)*1.36466322867013E58 + cos(Beta*2.0)*1.810096145662263E57 + cos(Alpha*2.0 + Beta*2.0)*3.087528119013518E58 - 1.068113514259467E60) - ((cos(Alpha + Beta)*1.611415088366918E24 + cos(Alpha)*8.622771990548647E22 - cos(Alpha + Beta)*cos(Beta)*1.008806316530991E22 - cos(Alpha)*cos(Beta)*8.914619657921062E22)*(cos(Alpha + Theta)*(-2.1672252E-1) + cos(Alpha + Beta + Theta)*2.4525E-2 + cos(Theta)*7.864554375E-1 - (Alphadot*Alphadot)*sin(Alpha + Beta)*3.125E-4 - (Betadot*Betadot)*sin(Alpha + Beta)*3.125E-4 + (Alphadot*Alphadot)*sin(Alpha)*2.7615E-3 + (Betadot*Betadot)*sin(Beta)*1.75E-4 - Alphadot*Betadot*sin(Alpha + Beta)*6.25E-4 - Alphadot*Thetadot*sin(Alpha + Beta)*6.25E-4 - Betadot*Thetadot*sin(Alpha + Beta)*6.25E-4 + Alphadot*Betadot*sin(Beta)*3.5E-4 + Alphadot*Thetadot*sin(Alpha)*5.523E-3 + Betadot*Thetadot*sin(Beta)*3.5E-4)*1.298074214633707E38) / (cos(Alpha*2.0)*1.36466322867013E58 + cos(Beta*2.0)*1.810096145662263E57 + cos(Alpha*2.0 + Beta*2.0)*3.087528119013518E58 - 1.068113514259467E60) + ((cos(Beta)*2.295431032101512E42 + pow(cos(Alpha), 2.0)*5.068268089629632E42 - cos(Alpha + Beta)*cos(Alpha)*1.147082221987514E42 + pow(cos(Alpha + Beta), 2.0)*6.490371073168535E40 - 1.844402881666282E44)*(cos(Alpha + Beta + Theta)*(-2.4525E-2) - (Thetadot*Thetadot)*sin(Alpha + Beta)*3.125E-4 + (Alphadot*Alphadot)*sin(Beta)*1.75E-4 + (Thetadot*Thetadot)*sin(Beta)*1.75E-4 + Alphadot*Thetadot*sin(Beta)*3.5E-4)*3.602879701896397E19) / (cos(Alpha*2.0)*1.36466322867013E58 + cos(Beta*2.0)*1.810096145662263E57 + cos(Alpha*2.0 + Beta*2.0)*3.087528119013518E58 - 1.068113514259467E60);


		double H_00 = (cos(Beta*2.0)*3.666615307735769E59 + cos(Alpha*2.0 + Beta*2.0)*1.169201309864722E60 + cos(Alpha)*2.107649049214543E61 - cos(Alpha + Beta*2.0)*1.309505467048489E60 - 1.517635867819378E62) / (cos(Alpha*2.0)*1.36466322867013E58 + cos(Beta*2.0)*1.810096145662263E57 + cos(Alpha*2.0 + Beta*2.0)*3.087528119013518E58 - 1.068113514259467E60);
		double H_01 = (cos(Alpha*2.0 + Beta*2.0)*-1.169201309864722E60 - cos(Alpha + Beta)*2.033877185724722E62 - cos(Alpha)*1.053824524607272E61 - cos(Beta)*3.101881122829559E61 + cos(Alpha - Beta)*5.785918955607044E60 + cos(Alpha + Beta*2.0)*6.547527335242445E59 + cos(Alpha*2.0 + Beta)*1.033199813501258E61 + 3.882786131833514E61) / (cos(Alpha*2.0)*1.36466322867013E58 + cos(Beta*2.0)*1.810096145662263E57 + cos(Alpha*2.0 + Beta*2.0)*3.087528119013518E58 - 1.068113514259467E60);
		double H_10 = (cos(Alpha*2.0 + Beta*2.0)*-1.169201309864722E60 - cos(Alpha + Beta)*2.033877185724722E62 - cos(Alpha)*1.053824524607272E61 - cos(Beta)*3.101881122829559E61 + cos(Alpha - Beta)*5.785918955607044E60 + cos(Alpha + Beta*2.0)*6.547527335242445E59 + cos(Alpha*2.0 + Beta)*1.033199813501258E61 + 3.882786131833514E61) / (cos(Alpha*2.0)*1.36466322867013E58 + cos(Beta*2.0)*1.810096145662263E57 + cos(Alpha*2.0 + Beta*2.0)*3.087528119013518E58 - 1.068113514259467E60);
		double H_11 = (cos(Alpha*2.0)*9.130180111947915E61 + cos(Alpha*2.0 + Beta*2.0)*1.169201309864722E60 + cos(Beta)*6.203762245659118E61 - cos(Alpha*2.0 + Beta)*2.066399627002516E61 - 6.552690702045327E63) / (cos(Alpha*2.0)*1.36466322867013E58 + cos(Beta*2.0)*1.810096145662263E57 + cos(Alpha*2.0 + Beta*2.0)*3.087528119013518E58 - 1.068113514259467E60);

		double alphaddot_ref = -kP_alpha * (Alpha - alpha_ref) - kD_alpha * Alphadot;
		double betaddot_ref = -kP_beta * (Beta - beta_ref) - kD_beta * Betadot;

		double D_00 = alphaddot_ref;
		double D_10 = betaddot_ref;

		double Inv_Coef = H_00 * H_11 - H_01 * H_10;
		// Now it is time to do the inverse matrix multiplication

		//u_vec = H\(D - R);
		//u_alpha = u_vec(1);
		//u_beta = u_vec(2);

		// These two are the first element in each control sequence
		double u_alpha_i = (H_11 * (D_00 - R_00) - H_01 * (D_10 - R_10)) / Inv_Coef;
		double u_beta_i = (-H_10 * (D_00 - R_00) + H_00 * (D_10 - R_10)) / Inv_Coef;

		Init_Ctrl_Sequence(2 * i, 0) = u_alpha_i;
		Init_Ctrl_Sequence(2 * i + 1, 0) = u_beta_i;

		Pre_Impact_State_i = Dynamics_Integrator_Onestep(Pre_Impact_State_i, ddt, u_alpha_i, u_beta_i, 1, Collision_Condition_i);
	}
	return;
}

double myObjfn(const Optimal_Control_Class &Optimal_Control_Seq)
{
	// This function calculates the value of the objective function

	RobotState Pre_Impact_State_i = Robot_State_In_Solver;
	RobotState Post_Impact_State;

	Collision_Struct  Collision_Detection_Condition;

	int integration_index_i = 0;                // Number of Iterations, should be an integer between 0 to grids -1
	double ddt = duration_time / (grids * 1.0 - 1.0), u_alpha_i, u_beta_i, obj_per_index;

	// ddt is the time step for the each step of Euler integration
	// u_alpha_i, u_beta_i: the control at each sequence of time
	// obj_per_index: the performance index of the given control sequence

	double Fai, rDx, rDy, Fric_A, Fric_D, impulse_mag, Impact_Prob, Fric_A_Prob, Fric_D_Prob, Post_Impact_Beta, Post_Impact_Betadot;

	while (Collision_Detection_Condition.CollisionOrNot == false)  // The contact must be established to make sure the control sequence can be optimized
	{
		if (integration_index_i < grids-1)
		{
			u_alpha_i = Optimal_Control_Seq(2 * integration_index_i);
			u_beta_i = Optimal_Control_Seq(2 * integration_index_i + 1);
		}
		else
		{
			// In this case, the control alpha/beta will be the last control element
			u_alpha_i = Optimal_Control_Seq(2 * (grids - 2));
			u_beta_i = Optimal_Control_Seq(2 * (grids - 2) + 1);
		}
		Pre_Impact_State_i = Dynamics_Integrator_Onestep(Pre_Impact_State_i, ddt, u_alpha_i, u_beta_i, 1, Collision_Detection_Condition);

		integration_index_i++;
	}

	impulse_mag = Impact_MapNMag_fn(Pre_Impact_State_i, Post_Impact_State);

	Fai = Collision_Detection_Condition.fai;
	rDx = Collision_Detection_Condition.rD_x;
	rDy = Collision_Detection_Condition.rD_y;
	Post_Impact_Beta = Post_Impact_State.beta;
	Post_Impact_Betadot = Post_Impact_State.betadot;

	Fric_A = FiveD_Intep(Fai, rDx, rDy, Post_Impact_Beta, Post_Impact_Betadot, Mu_A_Lib);
	Fric_D = FiveD_Intep(Fai, rDx, rDy, Post_Impact_Beta, Post_Impact_Betadot, Mu_D_Lib);

	Impact_Prob = Impact_Gauss.cdf(impulse_mag);
	Fric_A_Prob = Mu_A_Gauss.cdf(Fric_A);
	Fric_D_Prob = Mu_E_Gauss.cdf(Fric_D);

	obj_per_index = 1 - (1 - Impact_Prob) * (1 - Fric_A_Prob) * (1 - Fric_D_Prob);

	return obj_per_index;
}
