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

// Initial reading function
bool Robotstate_Reader(string txt_name,  RobotState &Initial_State_i)
{
	// The default name is "Initial_State.txt"
	// This function is used to read in the inital disturbed robot state into the pre-impact state from the file with "name"
	// If the output is 0, this means that the reading is failure. If the return value is 1, this menas that the reading is successful!

	ifstream InitFiles;
	string txt_name_file = txt_name + ".txt";
	InitFiles.open(txt_name_file);			// Here the Init_State should be a 6 by 1 vector saving the initial disturbed robot state
	if (!InitFiles)
	{
		std::cout << "Open File Error\n";
		std::cout<< txt_name_file << "does not exist!\n";
		return 0;
	}
	RobotState Initial_State;
	double theta, alpha, beta, thetadot, alphadot, betadot;
	InitFiles >> theta >> alpha >> beta >> thetadot >> alphadot >> betadot;
	InitFiles.close();

	Initial_State_i.theta = dlib::pi *theta;
	Initial_State_i.alpha = dlib::pi *alpha;
	Initial_State_i.beta = dlib::pi *beta;
	Initial_State_i.thetadot = thetadot;
	Initial_State_i.alphadot = alphadot;
	Initial_State_i.betadot = betadot;
	return 1;
}

bool Envi_Map_Reader(vector<double> &map_vector_x, vector<double> &map_vector_y)
{
	// The default name is "Map_File.txt"
	// This function is used to read in the map profile needed for the hand contact optimization
	// Since the hand contact is happened inside the falling plane. As a result, the coordinate of the vertex is of 2 dimension.

	// The preferred way to define this map is to list the x, y coordinates of the vertex tuples: x1, y1, x1', y1',...
	// Since the number of segments is not known in advance, we would have to read all of them in and then figure out the number of segments inside the environment system.

	ifstream MapFiles;
	string map_name_file = "Map_File.txt";
	MapFiles.open(map_name_file);
	if (!MapFiles)
	{
		std::cout << "Open File Error\n";
		std::cout << map_name_file <<"does not exist!\n";
		return 0;
	}
	// Otherwise, the map information will be read into a global vector
	// extracting words form the file

	double Map_i = 0.0;

	std::vector<double> map_vector;

	while (MapFiles >> Map_i)
	{
		map_vector.push_back(Map_i);
	}

	// Here the operation is to assign the x, y coordinate of map_vector into map_vector_x, map_vector_y
	for (int i = 0; i < map_vector.size()/2; i++)
	{
		map_vector_x.push_back(map_vector[2 * i]);
		map_vector_y.push_back(map_vector[2 * i + 1]);
	}

	return 1;
}

bool Database_Reader(Mu_Database &Mu_A_Lib, Mu_Database &Mu_D_Lib)
{
	// This function is used to read-in the pre-computed database for the hand contact and foot contact optimization
	FILE* f_A = fopen("Mu_A_Lib.bin", "rb");
	int File_A_Res, File_D_Res;
	File_A_Res = fread(&Mu_A_Lib(0), sizeof(float), N_Database_Len, f_A);
	fclose(f_A);

	FILE* f_D = fopen("Mu_D_Lib.bin", "rb");
	File_D_Res = fread(&Mu_D_Lib(0), sizeof(float), N_Database_Len, f_D);
	fclose(f_D);

	int File_Result = 1;

	if (File_A_Res == 0)
	{
		std::cout << "Mu_A_Lib.bin Open File Error\n";
		File_Result = 0;
	}

	if (File_D_Res == 0)
	{
		std::cout << "Mu_D_Lib.bin Open File Error\n";
		File_Result = 0;
	}
	return File_Result;
}
