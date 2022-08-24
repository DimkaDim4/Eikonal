#pragma once
#include <iostream>
#include "ContainerList.h"
#include "Containers.h"
#include "EnviromentModel.h"

struct TravelTimeNode {
	double* tt_wave = nullptr;
	int count = 0;

	void add_value(const double& value);
	void clear();
};


class TravelTimeFld {
public:
	TravelTimeFld();
	~TravelTimeFld();

	int set_size(const Index &size);
	int add_traveltime(int &i, int &j, int &k, double &value);
	int add_traveltime(double* grid, bool* frozen_cells);

private:
	Index size_;
	TravelTimeNode* mesh_;
};


struct Task {
	std::list<int> calculation_layers;
	std::list<int> no_calculation_layers;
	WaveFront init_front;
	void add_front(ContactBoundary*) {};
};


class Eikonal {
public:
	Eikonal();

	void SetModel(EnviromentModel* env_model);
	void SetSourse(const double& x, const double& y, const double& z);

	void Calculate(const double& time_limit);

	void WriteToFile();

private:
	std::list<Task> tasks_;
	TravelTimeFld ttf_;
	EnviromentModel *env_model_;

	void AnalysCurrentSolution();
};

void FastSweep3D(double* grid, bool* frozen_cells, EnviromentModel* model);
