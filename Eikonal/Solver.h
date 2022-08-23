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

	int set_size(int &I, int &J, int &K);
	int add_traveltime(int &i, int &j, int &k, double &value);
	int add_traveltime(double* grid, bool* frozen_cells);

private:
	Index size;
	TravelTimeNode* mesh;
};


struct Task {
	std::list<int> calculation_layers_;
	WaveFront init_front;
};


class TaskManager {
public:
	TaskManager();

	int add_task(std::list<int>& calculation_layers, WaveFront& init_front);
	bool empty();

private:
	std::list<Task> tasks_;
};


class Eikonal {
public:
	Eikonal();


private:
	TaskManager tasks_;
	TravelTimeFld ttf_;
	EnviromentModel env_model_;
};

void FastSweep3D(double* grid, bool* frozen_cells, EnviromentModel* model);
