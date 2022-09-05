#pragma once
#include <iostream>
#include "ContainerList.h"
#include "Containers.h"
#include "EnviromentModel.h"
#include <fstream>
#include <string>

class TravelTimeFld {
public:
	TravelTimeFld();
	~TravelTimeFld();

	int set_size(const Index &size);
	int add_traveltime(double* grid, bool* frozen_cells);
	void write();

private:
	Index size_;
	std::list<double>* mesh_;
};


struct Task {
	std::list<int> calculation_layers;
	std::list<int> no_calculation_layers;

	int size_wave_front;
	Index* wave_front_coords;
	double* travel_time_values;

	// ���� � �������� ��������� ������ ������� ���������� �������, �� ��� ������� �� ������ �����
	// �.�. ���� wave_front_coords ���������� �������� ��������� -> ��� ���������� ���� ������, ������� ������ �� �����
	// wave_front_coords == nullptr
	// ���� ��������� ������ �������� ���� �������, 
	bool is_contact_bound;

	std::list<int> travel_layers; // ������������������ �����, ����� ������� ������ ����� �� �����
	std::list<std::string> filenames; // ������� ������ - �������� �����, � ������� ����� ��������� ���� ������ � ����
};


class Eikonal {
public:
	Eikonal();

	void SetModel(EnviromentModel* env_model);
	void AddSourse(const double& x, const double& y, const double& z);

	void Calculate(const double& time_limit);

private:
	std::list<Task> tasks_;
	EnviromentModel *env_model_;

	std::list<Index> source_;
};

void FastSweep3D(double* grid, bool* frozen_cells, EnviromentModel* model, int type_wave);

std::list<std::list<int>> gen(const std::list<int>& list);
