#include <string>
#include <math.h>
#include <iostream>
#include <time.h>
#include <cmath>
#include <stdlib.h>     /* srand, rand */
#include <dlib/matrix.h>
#include <dlib/optimization.h>
#include "HC_Header.h"

using namespace std;
using namespace dlib;

std::vector<double> Map_Vector_x;
std::vector<double> Map_Vector_y;

Mu_Database Mu_A_Lib, Mu_D_Lib;   									// This is the offline database

Optimal_Control_Class Optimal_Control_Seq;					// This is the optimal control sequence

RobotState Robot_State_In_Solver;



int main()
{
	// This is the main function to run this optimization with hand contact

	// 1. Read-in the initial disturbed robot state
	RobotState Initial_State_i;
	string Init_State_File = "Initial_State";
	bool Init_State_Res = Robotstate_Reader(Init_State_File, Initial_State_i);
	if (Init_State_Res == 0)
	{
		std::cout << "Initial Robot state reading failure!" << endl;
		return 0;
	}

	// 2. Read-in the shape of the environment shape
	bool Init_Map_Res = Envi_Map_Reader(Map_Vector_x, Map_Vector_y);
	if (Init_Map_Res == 0)
	{
		std::cout << "Initial Map reading failure!" << endl;
		return 0;
	}

	// 3. Read-in the pre-computed database of the friction coefficients for both hand contact and foot contact
	bool Init_Database_Res = Database_Reader(Mu_A_Lib, Mu_D_Lib);
	if (Init_Database_Res == 0)
	{
		std::cout << "Initial Database reading failure!" << endl;
		return 0;
	}

	// Here the initial read-in functions have already been completed!

	clock_t start, end;
	start = clock();

	std::vector<double> Optimal_Control_Sequence = Hand_Contact_Optimization(Initial_State_i, Optimal_Control_Seq);

	end = clock();
	double running_time = double(end - start);
	cout << Optimal_Control_Seq << endl;
	cout << running_time << "ms" << endl;

	ofstream Optimal_Control_file;
	Optimal_Control_file.open("Optimal_Control_Seq.txt");
	Optimal_Control_file << Optimal_Control_Seq;
	Optimal_Control_file.close();

	return 0;
}
