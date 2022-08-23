#include "Solver.h"

void TravelTimeNode::add_value(const double& value) {
	if (tt_wave == nullptr) {
		tt_wave = new double[1];
		tt_wave[0] = value;
		count = 1;
	}
	else {
		double* temp = new double[count + 1];
		for (int i = 0; i < count; ++i) {
			temp[i] = tt_wave[i];
		}

		temp[count] = value;
		++count;
		delete[] tt_wave;
		tt_wave = temp;
		temp = nullptr;
	}
}

void TravelTimeNode::clear() {
	if (tt_wave != nullptr) {
		delete[] tt_wave;
	}
	count = 0;
}



TravelTimeFld::TravelTimeFld() {
	size = Index();
	mesh = nullptr;
}

TravelTimeFld::~TravelTimeFld() {
	long int mesh_pos = 0;
	for (int i = 0; i < size.i; ++i) {
		for (int j = 0; j < size.j; ++j) {
			for (int k = 0; k < size.k; ++k) {
				mesh[mesh_pos].clear();
				++mesh_pos;
			}
		}
	}
	delete[] mesh;
}

int TravelTimeFld::set_size(int& I, int& J, int& K) {
	if ((I <= 0) || (J <= 0) || (K <= 0))
		return -1;

	size.i = I;
	size.j = J;
	size.k = K;
	mesh = new TravelTimeNode[I * J * K];
	return 0;
}

int TravelTimeFld::add_traveltime(int& i, int& j, int& k, double& value) {
	mesh[(i * size.j + j) * size.k + k].add_value(value);
	return 0;
}

int TravelTimeFld::add_traveltime(double* grid, bool* frozen_cells) {
	return 0;
}