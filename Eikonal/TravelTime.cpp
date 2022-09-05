#include "Solver.h"

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
	mesh_ = new std::list<double>[size_.i * size.j * size_.k];
	return 0;
}

int TravelTimeFld::add_traveltime(double* grid, bool* frozen_cells) {
	for (long int i = 0; i < size_.i * size_.j * size_.k; ++i) {
		if (!frozen_cells[i]) {
			auto it = std::find(mesh_[i].begin(), mesh_[i].end(), grid[i]);
			if (it == mesh_[i].end()) {
				mesh_[i].push_back(grid[i]);
			}
		}
	}
	return 0;
}

void TravelTimeFld::write() {
	std::ofstream out("Size.bin", std::ios::binary | std::ios::out);
	for (int i = 0; i < size_.i * size_.j * size_.k; i++) {
		int s = mesh_[i].size();
		out.write((char*)&s, sizeof(s));
	}
	out.close();

	std::ofstream tt("TravelTime_Data.bin", std::ios::binary | std::ios::out);
	for (int i = 0; i < size_.i * size_.j * size_.k; i++) {
		for (auto it = mesh_[i].begin(); it != mesh_[i].end(); it++) {
			double val = *it;
			tt.write((char*)&val, sizeof(double));
		}
	}
	tt.close();
}
