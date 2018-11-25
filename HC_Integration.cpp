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

extern std::vector<double> Map_Vector_x, Map_Vector_y;

RobotState Dynamics_Integrator_Onestep(RobotState initial_state, double dt, double u_alpha, double u_beta, int EulerOrRK4, Collision_Struct &collision_condition)
{
	// This function is used to integrate the system dynamics forwards
	// Inputs:
	//									RobotState initial_state			is the initial state for this integration
	//									double dt											is the total integration time duration
	//									double u_alpha, u_beta_i			the control lasts within this time duration
	//									collision_structure_i					the structure used to detect the collision at the robot end effector
	//									EulerOrRK4										the choice of the integrator type: Euler or RK4?
	//									collision_condition						a user-defined structure used to save the information of the collision


	// The default idea is to divide the dt into a certain number of grids and conduct the integration step by step
	double ddt = dt/(grids * 1.0 - 1.0);
	RobotState State_i = initial_state;

	for (int i = 0; i < grids; i++)
	{
		// Integration
		if (EulerOrRK4 == 1)
		{
			 // In this case, the integration is conduct with the Euler rule
			 State_i = Euler_Step(State_i, ddt, u_alpha, u_beta);
		}
		else
		{
			// In this case, the integration is conducted with the RK4 rule
			State_i = RK4_Step(State_i, ddt, u_alpha, u_beta);
		}
	}

	// Then it is to update the collision struct
	collision_condition = Collision_Condition_Update(State_i, collision_condition);

	return State_i;
}

RobotState Euler_Step(RobotState Pre_Impact_State, double ddt, double u_alpha, double u_beta)
{
	RobotState State_i;

	Stateddot Pre_Impact_Stateddot = Pre_Impact_EOM(Pre_Impact_State, u_alpha, u_beta);

	State_i.theta = 		Pre_Impact_State.theta + 			ddt * Pre_Impact_State.thetadot;
	State_i.alpha = 		Pre_Impact_State.alpha + 			ddt * Pre_Impact_State.alphadot;
	State_i.beta = 			Pre_Impact_State.beta + 			ddt * Pre_Impact_State.betadot;
	State_i.thetadot = 	Pre_Impact_State.thetadot + 	ddt * Pre_Impact_Stateddot.thetaddot;
	State_i.alphadot = 	Pre_Impact_State.alphadot + 	ddt * Pre_Impact_Stateddot.alphaddot;
	State_i.betadot = 	Pre_Impact_State.betadot + 		ddt * Pre_Impact_Stateddot.betaddot;

	return State_i;
};

RobotState RK4_Step(RobotState state, double ddt, double u_alpha, double u_beta)
{
	RobotState State_i;

	Stateddot Pre_Stateddot_1 = Pre_Impact_EOM(state, u_alpha, u_beta);
	Stateddot Pre_Stateddot_2 = Pre_Impact_EOM(State_Update_with_Acc(state, Pre_Stateddot_1, ddt / 2), u_alpha, u_beta);
	Stateddot Pre_Stateddot_3 = Pre_Impact_EOM(State_Update_with_Acc(state, Pre_Stateddot_2, ddt / 2), u_alpha, u_beta);
	Stateddot Pre_Stateddot_4 = Pre_Impact_EOM(State_Update_with_Acc(state, Pre_Stateddot_3, ddt), u_alpha, u_beta);

	RobotState state1 = State_Update_with_Acc(state, Pre_Stateddot_1, ddt / 2);
	RobotState state2 = State_Update_with_Acc(state, Pre_Stateddot_2, ddt / 2);
	RobotState state3 = State_Update_with_Acc(state, Pre_Stateddot_3, ddt);

	State_i.theta = 		state.theta + ddt / 6 * 		(state.thetadot + 2 * state1.thetadot 											+ 2 * state2.thetadot 					+ state3.thetadot);
	State_i.alpha = 		state.alpha + ddt / 6 * 		(state.alphadot + 2 * state1.alphadot 											+ 2 * state2.alphadot 					+ state3.alphadot);
	State_i.beta = 			state.beta + ddt / 6 * 			(state.betadot + 2 * state1.betadot 												+ 2 * state2.betadot 						+ state3.betadot);
	State_i.thetadot = 	state.thetadot + ddt / 6 * 	(Pre_Stateddot_1.thetaddot + 2 * Pre_Stateddot_2.thetaddot 	+ 2 * Pre_Stateddot_3.thetaddot + Pre_Stateddot_4.thetaddot);
	State_i.alphadot = 	state.alphadot + ddt / 6 * 	(Pre_Stateddot_1.alphaddot + 2 * Pre_Stateddot_2.alphaddot 	+ 2 * Pre_Stateddot_3.alphaddot + Pre_Stateddot_4.alphaddot);
	State_i.betadot = 	state.betadot + ddt / 6 * 	(Pre_Stateddot_1.betaddot + 2 * Pre_Stateddot_2.betaddot 		+ 2 * Pre_Stateddot_3.betaddot 	+ Pre_Stateddot_4.betaddot);

	return State_i;
};

RobotState State_Update_with_Acc(RobotState state, Stateddot stateddot, double ddt)
{
	// This function updates the robotstate using the acceleration and the time duration
	RobotState state_add_i;

	state_add_i.theta = 		state.theta + 		state.thetadot*ddt;
	state_add_i.alpha = 		state.alpha + 		state.alphadot*ddt;
	state_add_i.beta = 			state.beta + 			state.betadot*ddt;
	state_add_i.thetadot = 	state.thetadot + 	stateddot.thetaddot*ddt;
	state_add_i.alphadot = 	state.alphadot + 	stateddot.alphaddot*ddt;
	state_add_i.betadot = 	state.betadot + 	stateddot.betaddot*ddt;

	return state_add_i;
};

Collision_Struct Collision_Condition_Update(RobotState &pre_impact_state_i, Collision_Struct &Collision_Condition_i)
{
	// This function is used to detect the collision given the current state and the Map_Vector
	// If the collision has been detected, then the collision condition will be updated


	std::vector<double> End_Link_Pos = End_Link_Pos_Fn(pre_impact_state_i);

	double rC_x = End_Link_Pos[0];						double rC_y = End_Link_Pos[1];
	double rD_x = End_Link_Pos[2];						double rD_y = End_Link_Pos[3];

	// The comparison of the collision is through the detection of the intersection between the end link of the Map_Vector_x/y

	Point_Pos rC, rD;													rC.x = rC_x;  rC.y = rC_y;				rD.x = rD_x;  rD.y = rD_y;
	Point_Pos Inters_Point, Edge_Left, Edge_Right;

	for (int i = 0; i < Map_Vector_x.size()/2; i++)
	{
		Edge_Left.x = Map_Vector_x[2*i];
		Edge_Left.y = Map_Vector_y[2*i];

		Edge_Right.x = Map_Vector_x[2*i+1];
		Edge_Right.y = Map_Vector_y[2*i+1];

		bool CollisionOrNot = Intersection_Cal(rC, rD, Edge_Left, Edge_Right, Inters_Point);
		if(Collision_Condition_i.CollisionOrNot == true)
		{
			// In this case, the given robot state has a collision condition
			break;
		}
		else
		{
			// In this case, the given robot state does not have a collision condition
			if(CollisionOrNot == true)
			{
				// However, a new collision has been found
				Collision_Condition_i.rD_x = Inters_Point.x;
				Collision_Condition_i.rD_y = Inters_Point.y;
				Collision_Condition_i.fai = dlib::pi /2.0 - atan2(Edge_Right.y - Edge_Left.y, Edge_Right.x - Edge_Left.x);
				Collision_Condition_i.CollisionOrNot = true;
				break;
			}
			else
			{
				continue;
			}
		}
		return Collision_Condition_i;
	}
}

bool Intersection_Cal(Point_Pos p1, Point_Pos p2, Point_Pos p3, Point_Pos p4, Point_Pos &Inters_Point)
{
	// Store the values for fast access and easy
	// equations-to-code conversion
	float x1 = p1.x, x2 = p2.x, x3 = p3.x, x4 = p4.x;
	float y1 = p1.y, y2 = p2.y, y3 = p3.y, y4 = p4.y;

	float d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);			// The only case for this d to be 0 is that line:p1_p2 is parallel to link:p3_p4

	if (d <= Eps_Value)
	{
		return false;
	}
	else
	{
		// Get the x and y
		float pre = (x1*y2 - y1*x2), post = (x3*y4 - y3*x4);
		float x = (pre * (x3 - x4) - (x1 - x2) * post) / d;
		float y = (pre * (y3 - y4) - (y1 - y2) * post) / d;

		// Check if the x and y coordinates are within both lines
		if (x < min(x1, x2) || x > max(x1, x2) ||
		x < min(x3, x4) || x > max(x3, x4)) return false;
		if (y < min(y1, y2) || y > max(y1, y2) ||
		y < min(y3, y4) || y > max(y3, y4)) return false;

		// Return the point of intersection
		Inters_Point.x = x;
		Inters_Point.y = y;
		return true;
	}
}
