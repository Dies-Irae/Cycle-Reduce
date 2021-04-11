/*
 * main.cpp
 * Cycle Reduce Alogorithm.
 * Alogorithm 2 from "Random Redundant Soft-In Soft-Out Decoding of Linear Block Codes"
 * Thomas R. Halford and Keith M. Chugg

 * Copyright (c) 2020 Dies-Irae
 * All rights reserved.
 * https://github.com/Dies-Irae

 * Implemented by Dies-Irae in Chongqing University.
 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * with the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is furnished
 * to do so, subject to the following conditions:
 
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimers.
 * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimers in the documentation
 * and/or other materials provided with the distribution.
 * Neither the names of Thomas R. Halford, the Communication Sciences Institute,
 * the University of Southern California nor the names of its contributors may
 * be used to endorse or promote products derived from this Software without
 * specific prior written permission. 

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
 * AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.
 */
#include <iostream>
#include <fstream>
#include <string>
#include "Short_Cycle_Matrix.h"
#include "Short_Cycle_Counter.h"
#include <vector>
#include <string>
#include <sstream>
using namespace std;
vector<vector<int>> readHMatrix(string filename)
{
	fstream fin(filename);
	string line;
	vector<vector<int>> result;
	while (getline(fin, line))
	{
		int tmp;
		vector<int> vline;
		istringstream iss(line);
		while (iss >> tmp)
			vline.push_back(tmp);
		result.push_back(vline);
	}
	return result;
}
vector<int> binary_vector_sum(vector<int> v1, vector<int> v2)
{
	auto iter2 = v2.begin();
	for (auto iter = v1.begin(); iter != v1.end(); iter++)
	{
		*iter = *iter ^ *iter2;
		iter2++;
	}
	return v1;
}

vector<vector<int>> cycle_reduce(Short_Cycle_Matrix E, vector<vector<int>> H)
{


	Short_Cycle_Counter E_counter(E);
	E_counter.count();
	size_t r1_star = 1, r2_star = 1, g_star = E_counter.girth();
	size_t Ng_star = E_counter.Ng();
	size_t Ng2_star = E_counter.Ng2();
	do {
		// if r1* != r2*, replace row
		if (r1_star != r2_star)
		{
			H[r2_star] = binary_vector_sum(H[r1_star], H[r2_star]);
		}
		r1_star = -1, r2_star = -1;
		for (int r1 = 0; r1 < H.size(); r1++)
		{
			for (int r2 = 0; r2 < H.size(); r2++)
			{
				if (r1 != r2)
				{
					// replace row
					vector<int> temp(H[r2]);
					int r2_s_tmp = r2;
					H[r2] = binary_vector_sum(H[r1], H[r2]);
					// re calculate cycles
					E.read_incidence_from_vector(H[0].size(), H.size(), H);// convert vector to matrix for counting
					E_counter.initialize(E);// counter should be re-intialized after modifing matrix
					E_counter.count();

					size_t g = E_counter.girth();
					size_t Ng = E_counter.Ng();
					size_t Ng2 = E_counter.Ng2();

					if (g < g_star)
					{
						g_star = g, r1_star = r1, r2_star = r2, Ng2_star = Ng;
						Ng2_star = Ng2;
					}
					else if (Ng < Ng_star)
					{
						r1_star = r1, r2_star = r2;
						Ng_star = Ng, Ng2_star = Ng2;
					}
					else if (Ng = Ng_star)
					{
						if (Ng2 < Ng2_star)
						{
							r1_star = r1, r2_star = r2;
							Ng2_star = Ng2;
						}
					}
					H[r2_s_tmp] = temp;
				}
			}
		}
	} while (r1_star != -1 or r2_star != -1);
	return H;
}

int main()
{
	Short_Cycle_Matrix E;
	//E.read_alist_file("BCH63_45.txt");
	E.read_incidence_matrix_file(73, 36, "QR73_37.txt");
	Short_Cycle_Counter E_counter(E);
	E_counter.count();
	int g = E_counter.girth();

	cout << endl
		<< "Cycle Count:" << endl
		<< "girth = " << g << endl
		<< "N_" << g << " = " << E_counter.Ng() << endl
		<< "N_" << g + 2 << " = " << E_counter.Ng2() << endl
		<< "N_" << g + 4 << " = " << E_counter.Ng4() << endl;
	double mg = 0.0, sg = 0.0, mg2 = 0.0, sg2 = 0.0, mg4 = 0.0, sg4 = 0.0;
	E_counter.cycle_dist(&mg, &sg, &mg2, &sg2, &mg4, &sg4);
	cout << endl
		<< "Cycle Distribution:" << endl
		<< "mu_" << g << " = " << mg << ", sigma_" << g << " = " << sg << endl
		<< "mu_" << g + 2 << " = " << mg2 << ", sigma_" << g + 2 << " = " << sg2 << endl
		<< "mu_" << g + 4 << " = " << mg4 << ", sigma_" << g + 4 << " = " << sg4 << endl
		<< endl;
	vector<vector<int>> H;
	H = readHMatrix("QR73_37.txt");
	E.read_incidence_from_vector(H[0].size(), H.size(), H);
	auto H2 = cycle_reduce(E, H);
	E.read_incidence_from_vector(H2[0].size(), H2.size(), H2);
	E_counter.initialize(E);
	E_counter.count();
	g = E_counter.girth();

	cout << endl
		<< "Cycle Count:" << endl
		<< "girth = " << g << endl
		<< "N_" << g << " = " << E_counter.Ng() << endl
		<< "N_" << g + 2 << " = " << E_counter.Ng2() << endl
		<< "N_" << g + 4 << " = " << E_counter.Ng4() << endl;
	E_counter.cycle_dist(&mg, &sg, &mg2, &sg2, &mg4, &sg4);
	cout << endl
		<< "Cycle Distribution:" << endl
		<< "mu_" << g << " = " << mg << ", sigma_" << g << " = " << sg << endl
		<< "mu_" << g + 2 << " = " << mg2 << ", sigma_" << g + 2 << " = " << sg2 << endl
		<< "mu_" << g + 4 << " = " << mg4 << ", sigma_" << g + 4 << " = " << sg4 << endl
		<< endl;
	for (auto i = H2.begin(); i != H2.end(); i++)
	{
		for (auto j = i->begin(); j != i->end(); j++)
		{
			cout << *j << " ";
		}
		cout << "\n";
	}
	return 0;

}