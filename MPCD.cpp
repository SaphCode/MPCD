// MPCD.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <Eigen/Dense>
#include "Particle.h"
#include "MPCD.h"

using namespace Eigen;
using namespace MPCD;

int main()
{
	Vector2d pos0(Grid::cell_dim / 10, Grid::cell_dim / 10);
	Vector2d vel0(0, 0);
	Particle p(pos0, vel0);

	Vector2d shift1(Grid::max_shift, Grid::max_shift);
	Vector2i index_shift1(0, 0);
	Vector2d shift2(Grid::max_shift, -Grid::max_shift);
	Vector2i index_shift2(1, -1);

	Vector2d shift_zero(0, 0);
	Vector2i zero_index = p.shift(shift_zero);

	Vector2i zeros(0, 0);

	//ASSERT_TRUE(areVectorsEqual(zero_index, zeros)) << "Particle with position inside first cell should have (0,0) index (first cell should be (0,0))";

	Vector2i shift1_index = p.shift(shift1);
	Vector2i shift1_index_calc = zero_index + index_shift1;
	Vector2i r_shift1_index = p.shift(-shift1);
	Vector2i r_shift1_index_calc = zero_index - index_shift1;

	//EXPECT_EQ(r_shift1_index_calc(0), r_shift1_index(0));
	//EXPECT_EQ(r_shift1_index_calc(1), r_shift1_index(1));

	//ASSERT_TRUE(areVectorsEqual(shift1_index_calc, shift1_index)) << "Shifting once should not change index in this case.";
	//ASSERT_TRUE(areVectorsEqual(r_shift1_index_calc, r_shift1_index)) << "Shifting by + and - same vector should get us back to the beginning index.";
	//ASSERT_TRUE(areVectorsEqual(pos0, p.getPosition())) << "Position does not change in a shift. This would take unnecessary storage and computation time.";

	Vector2i shift2_index = p.shift(shift2);
	Vector2i shift2_index_calc = zero_index + index_shift2;
	Vector2i r_shift2_index = p.shift(-shift2);
	Vector2i r_shift2_index_calc = zero_index - index_shift2;

	//ASSERT_TRUE(areVectorsEqual(r_shift2_index_calc, r_shift2_index)) << "Shifting once should change index in this case.";;
	//ASSERT_TRUE(areVectorsEqual(shift2_index_calc, shift2_index)) << "Shifting by + and - same vector should get us back to the beginning index.";;
	//ASSERT_TRUE(areVectorsEqual(pos0, p.getPosition())) << "Position does not change in a shift. This would take unnecessary storage and computation time.";

	//ASSERT_TRUE(areVectorsEqual(p.getVelocity(), vel0)) << "Velocity does not change in a shift.";
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
