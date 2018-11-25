#pragma once
#include <vector>
#include <fstream>
#include <dlib/matrix.h>
#include <dlib/optimization.h>

// Some data structures

struct RobotState
{
	double theta, alpha, beta;
	double thetadot, alphadot, betadot;
};

struct Stateddot
{
	double thetaddot, alphaddot, betaddot;
};

struct Force_Class
{
	double FA_x, FA_y, FD_x, FD_y;
};

struct Control_Class
{
	double u_alpha;
	double u_beta;
};

struct Mu_A_Mu_D_Class
{
	double mu_A, mu_D;  // K is the derivative gain for the post-impact system,
											// mu_A/mu_D: upper bound of the required friction coefficient
};

struct Impact_Moment
{
	double Impact_Moment_x, Impact_Moment_y;
};

struct Collision_Struct
{
	double rD_x, rD_y, fai;				// rD_x, rD_y saves the position of the hand contact point while fai saves the normal angle of the hand contact surface
																// where the vertical wall is said to have a 0 fai angle
	bool CollisionOrNot;   // The default collision_detection_flag is set to be false
};

struct Point_Pos
{
	// Just a structure to save the information of a point: x, y coordinate
	double x, y;
};

class NormalDistribution
{
public:
	NormalDistribution(double _mu, double _sigma) : mu(_mu), sigma(_sigma) {}
	inline double pdf(double x);
	double cdf(double x);
private:
	double mu;
	double sigma;
};

// Some pre-defined constrants
const int N_Nor_Angl = 46;          // This is the number of points that the normal direction could take
const int N_Hz_Dist = 16;           // This is the number of points in the horizontal direction
const int N_Vt_Dist = 28;           // This is the number of points in the vertical direction
const int N_Beta = 41;							// This is the number of points for the discretization of beta angles
const int N_Betadot = 41;           // This is the number of points for the discretization of betadot

const int N_Database_Len = N_Nor_Angl * N_Hz_Dist *N_Vt_Dist * N_Beta * N_Betadot;

typedef dlib::matrix<float, N_Database_Len, 1> Mu_Database;					// This is actually the pre-computed database used to accelerate the optimiztion process

const float Impact_Avg = 0.35,	Impact_Std = 0.1;
const float Mu_A_Avg = 0.5,			Mu_A_Std = 0.1;
const float Mu_D_Avg = 0.4,			Mu_D_Std = 0.1;

const double duration_time = 0.25;																	// The total time of duration
const int grids = 11;																								// The total number of grids within that given duration

typedef dlib::matrix<double, 2 * (grids-1), 1> Optimal_Control_Class;					// This type is used to define the optimal control class which is a piecewise linear function

using namespace std;

#define Eps_Value std::numeric_limits<double>::epsilon()

// Dynamics functions
Stateddot Pre_Impact_EOM(RobotState Pre_Impact_State, double u_alpha, double u_beta);
double Impact_MapNMag_fn(RobotState Pre_Impact_State_i, RobotState &Post_Impact_State);
std::vector<double> End_Link_Pos_Fn(RobotState robot_state_i);

// Database functions
double FiveD_Intep(double Fai, double Wall_x, double Wall_y, double Beta, double Betadot, Mu_Database &Mu_Lib);
double Three_Dim_Intep(int Fai_index, int Hz_index, double rEy, double Beta, double Betadot, Mu_Database &Mu_Lib);
int Index_Selection_fn(double value_i, double value_low, double value_high, double value_unit, int value_number);
int Index_Five2One(int i, int j, int k, int l, int m);

// Integration functions
RobotState Dynamics_Integrator_Onestep(RobotState initial_state, double dt, double u_alpha, double u_beta, int EulerOrRK4, Collision_Struct &collision_condition);
RobotState Euler_Step(RobotState Pre_Impact_State, double ddt, double u_alpha, double u_beta);
RobotState RK4_Step(RobotState state, double ddt, double u_alpha, double u_beta);
RobotState State_Update_with_Acc(RobotState state, Stateddot stateddot, double ddt);
Collision_Struct Collision_Condition_Update(RobotState &pre_impact_state_i, Collision_Struct &Collision_Condition_i);
bool Intersection_Cal(Point_Pos p1, Point_Pos p2, Point_Pos p3, Point_Pos p4, Point_Pos &Inters_Point);

// Loader function
bool Robotstate_Reader(string txt_name,  RobotState &Initial_State_i);
bool Envi_Map_Reader(vector<double> &map_vector_x, vector<double> &map_vector_y);
bool Database_Reader(Mu_Database &Mu_A_Lib, Mu_Database &Mu_D_Lib);

// Optimization function
std::vector<double> Hand_Contact_Optimization(RobotState &Initial_State_i, Optimal_Control_Class &Optimal_Control_Seq);
void Init_Opt_Ctrl_Generator(RobotState optimal_initial_condition, Optimal_Control_Class &Init_Ctrl_Sequence, Collision_Struct &Collision_Condition_i);
double myObjfn(const Optimal_Control_Class &Optimal_Control_Seq);
