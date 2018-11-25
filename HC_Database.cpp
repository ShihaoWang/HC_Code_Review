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

double FiveD_Intep(double Fai, double Wall_x, double Wall_y, double Beta, double Betadot, Mu_Database &Mu_Lib)
{
	/*
	Given to this function are five variables determined by the environment

	%     Fai:   the angle of the normal vector at the contact point

	%  Wall_x : the horizontal position of the contact point with respect to
	%           the feet of the robot
	%  Wall_y : the vertical position of the contact point with respect to the
	%           feet of the robot

	%    Beta : the post impact angle
	% Betadot : the post impaat angular velocity

	%  Accoridng to the intepolation pattern,
	% 3 dimensional interpolation == == = > 7 times
	% 4 dimensional interpolation == == = > 15 times
	% 5 dimensional interpolation == == = > 31 times
	*/

	//Five dimensional array

	double Fai_low = 0.0;
	double Hz_Dist_low = 0.1;

	double Fai_hgh = dlib::pi / 2.0;
	double Hz_Dist_hgh = 0.25;

	// The units for all the five variables
	double Fai_unit = (Fai_hgh - Fai_low) / (N_Nor_Angl - 1.0);
	double Hz_Dist_unit = (Hz_Dist_hgh - Hz_Dist_low) / (N_Hz_Dist - 1.0);

	int Fai_index = Index_Selection_fn(Fai, Fai_low, Fai_hgh, Fai_unit, N_Nor_Angl);
	int Hz_index = Index_Selection_fn(Wall_x, Hz_Dist_low, Hz_Dist_hgh, Hz_Dist_unit, N_Hz_Dist);

	double Cost_A = Three_Dim_Intep(Fai_index - 1, Hz_index - 1, Wall_y, Beta, Betadot, Mu_Lib);
	double Cost_B = Three_Dim_Intep(Fai_index, Hz_index - 1, Wall_y, Beta, Betadot, Mu_Lib);
	double Cost_C = Three_Dim_Intep(Fai_index, Hz_index, Wall_y, Beta, Betadot, Mu_Lib);
	double Cost_D = Three_Dim_Intep(Fai_index - 1, Hz_index, Wall_y, Beta, Betadot, Mu_Lib);

	// This stage is a pure bi - linear intepolation

	double Fai_in_lib = Fai_low + (Fai_index - 1.0) * Fai_unit;
	double Wall_x_in_lib = Hz_Dist_low + (Hz_index - 1.0) *Hz_Dist_unit;

	double Cost_MAB = (Cost_B - Cost_A) / Hz_Dist_unit * (Wall_x - (Wall_x_in_lib - Hz_Dist_unit)) + Cost_A;
	double Cost_MCD = (Cost_C - Cost_D) / Hz_Dist_unit * (Wall_x - (Wall_x_in_lib - Hz_Dist_unit)) + Cost_D;

	double Cost_MABCD = (Cost_MCD - Cost_MAB) / Fai_unit * (Fai - (Fai_in_lib - Fai_unit)) + Cost_MAB;

	//Reserved for test

	double Cost_MAD = (Cost_D - Cost_A) / Fai_unit * (Fai - (Fai_in_lib - Fai_unit)) + Cost_A;
	double Cost_MBC = (Cost_C - Cost_B) / Fai_unit * (Fai - (Fai_in_lib - Fai_unit)) + Cost_B;
	double Cost_MADCB = (Cost_MBC - Cost_MAD) / Hz_Dist_unit * (Wall_x - (Wall_x_in_lib - Hz_Dist_unit)) + Cost_MAD;

	double Inteplt_Val = Cost_MABCD;

	return Inteplt_Val;
}

double Three_Dim_Intep(int Fai_index, int Hz_index, double rEy, double Beta, double Betadot, Mu_Database &Mu_Lib)
{
	// This function is used to do the interpolation in a three dimensional space.
	// The upper bound and lower bound of all the four variables

	double Vt_Dist_low = 0.0;
	double Beta_low = dlib::pi / 4.0;
	double Betadot_low = -6.0;

	double Vt_Dist_hgh = 0.27;
	double Beta_hgh = dlib::pi / 1.0;
	double Betadot_hgh = 6.0;

	// The units for all the four variables
	double Vt_Dist_unit = (Vt_Dist_hgh - Vt_Dist_low) / (N_Vt_Dist - 1.0);
	double Beta_unit = (Beta_hgh - Beta_low) / (N_Beta - 1.0);
	double Betadot_unit = (Betadot_hgh - Betadot_low) / (N_Betadot - 1.0);

	int Vt_index = Index_Selection_fn(rEy, Vt_Dist_low, Vt_Dist_hgh, Vt_Dist_unit, N_Vt_Dist);
	int Beta_index = Index_Selection_fn(Beta, Beta_low, Beta_hgh, Beta_unit, N_Beta);
	int Betadot_index = Index_Selection_fn(Betadot, Betadot_low, Betadot_hgh, Betadot_unit, N_Betadot);

	int Vt_index_1 = Vt_index - 1;
	int Beta_index_1 = Beta_index - 1;
	int Betadot_index_1 = Betadot_index - 1;

	// First layer of the four points
	int Point_A_Index = Index_Five2One(Fai_index, Hz_index, Vt_index_1, Beta_index_1, Betadot_index_1);
	int Point_B_Index = Index_Five2One(Fai_index, Hz_index, Vt_index, Beta_index_1, Betadot_index_1);
	int Point_C_Index = Index_Five2One(Fai_index, Hz_index, Vt_index, Beta_index, Betadot_index_1);
	int Point_D_Index = Index_Five2One(Fai_index, Hz_index, Vt_index_1, Beta_index, Betadot_index_1);

	//Second layer of the four points
	int Point_E_Index = Index_Five2One(Fai_index, Hz_index, Vt_index_1, Beta_index, Betadot_index);
	int Point_F_Index = Index_Five2One(Fai_index, Hz_index, Vt_index_1, Beta_index_1, Betadot_index);
	int Point_G_Index = Index_Five2One(Fai_index, Hz_index, Vt_index, Beta_index_1, Betadot_index);
	int Point_H_Index = Index_Five2One(Fai_index, Hz_index, Vt_index, Beta_index, Betadot_index);

	//float Cost_A = Mu_Lib(Point_A_Index);
	float Cost_A = Mu_Lib(Point_A_Index);
	float Cost_B = Mu_Lib(Point_B_Index);
	float Cost_C = Mu_Lib(Point_C_Index);
	float Cost_D = Mu_Lib(Point_D_Index);

	float Cost_E = Mu_Lib(Point_E_Index);
	float Cost_F = Mu_Lib(Point_F_Index);
	float Cost_G = Mu_Lib(Point_G_Index);
	float Cost_H = Mu_Lib(Point_H_Index);

	double Betadot_in_lib = Betadot_low + (Betadot_index - 1.0) * Betadot_unit;
	double Beta_in_lib = Beta_low + (Beta_index - 1.0) * Beta_unit;
	double rEy_in_lib = Vt_Dist_low + (Vt_index - 1.0) * Vt_Dist_unit;

	//There are totally 7 interpolations to make
	double Cost_MAF = (Cost_F - Cost_A) / (Betadot_unit)*(Betadot - (Betadot_in_lib - Betadot_unit)) + Cost_A;
	double Cost_MBG = (Cost_G - Cost_B) / (Betadot_unit)*(Betadot - (Betadot_in_lib - Betadot_unit)) + Cost_B;
	double Cost_MCH = (Cost_H - Cost_C) / (Betadot_unit)*(Betadot - (Betadot_in_lib - Betadot_unit)) + Cost_C;
	double Cost_MDE = (Cost_E - Cost_D) / (Betadot_unit)*(Betadot - (Betadot_in_lib - Betadot_unit)) + Cost_D;

	double Cost_MAF_MDE = (Cost_MDE - Cost_MAF) / Beta_unit * (Beta - (Beta_in_lib - Beta_unit)) + Cost_MAF;
	double Cost_MBG_MCH = (Cost_MCH - Cost_MBG) / Beta_unit * (Beta - (Beta_in_lib - Beta_unit)) + Cost_MBG;

	double Cost_MAF_MDE_MBG_MCH = (Cost_MBG_MCH - Cost_MAF_MDE) / Vt_Dist_unit*(rEy - (rEy_in_lib - Vt_Dist_unit)) + Cost_MAF_MDE;

	double Three_D_val = Cost_MAF_MDE_MBG_MCH;

	return Three_D_val;
}

int Index_Selection_fn(double value_i, double value_low, double value_high, double value_unit, int value_number)
{
	int Index_i;

	if (value_i <= value_low + 1e-5)
	{
		Index_i = 2;
		return Index_i;
	}

	if (value_i >= value_high)
	{
		Index_i = value_number;
		return Index_i;
	}

	if ((value_i > value_low) && (value_i < value_high))
	{
		Index_i = ceil((value_i - value_low) / value_unit);
		if (Index_i == 1)
			Index_i = 2;
		return Index_i;
	}
}

int Index_Five2One(int i, int j, int k, int l, int m)
{
	// This function is used to generate the index of the (i,j,k,l,m) element in the given database

	int Index = (i - 1) * (N_Hz_Dist * N_Vt_Dist * N_Beta * N_Betadot) + (j - 1) * (N_Vt_Dist * N_Beta *  N_Betadot) + (k - 1) * (N_Beta * N_Betadot) + (l - 1) * N_Betadot + m - 1;

	return Index;
}
