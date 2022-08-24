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
	mesh_ = nullptr;
}

TravelTimeFld::~TravelTimeFld() {
	long int mesh_pos = 0;
	for (int i = 0; i < size_.i; ++i) {
		for (int j = 0; j < size_.j; ++j) {
			for (int k = 0; k < size_.k; ++k) {
				mesh_[mesh_pos].clear();
				++mesh_pos;
			}
		}
	}
	delete[] mesh_;
}

int TravelTimeFld::set_size(const Index &size) {
	if ((size.i <= 0) || (size.j <= 0) || (size.k <= 0))
		return -1;

	size_ = size;
	mesh_ = new TravelTimeNode[size_.i * size.j * size_.k];
	return 0;
}

int TravelTimeFld::add_traveltime(int& i, int& j, int& k, double& value) {
	mesh_[(i * size_.j + j) * size_.k + k].add_value(value);
	return 0;
}

int TravelTimeFld::add_traveltime(double* grid, bool* frozen_cells) {
	return 0;
}