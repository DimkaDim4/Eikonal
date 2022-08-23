#pragma once
#include "Containers.h"
#include "ContainerList.h"
#include <iostream>
#include <list>
#include <map>
#include <vector>

const int count_neighbors = 26;
const int neighboring_nodes[count_neighbors][3] = // i, j, k
{ { -1, -1, -1 }, { -1, -1, 0 }, { -1, -1, 1 },
  { -1,  0, -1 }, { -1,  0, 0 }, { -1,  0, 1 },
  { -1,  1, -1 }, { -1,  1, 0 }, { -1,  1, 1 },

  {  0, -1, -1 }, {  0, -1, 0 }, {  0, -1, 1 },
  {  0,  0, -1 },                {  0,  0, 1 },
  {  0,  1, -1 }, {  0,  1, 0 }, {  0,  1, 1 },

  {  1, -1, -1 }, {  1, -1, 0 }, {  1, -1, 1 },
  {  1,  0, -1 }, {  1,  0, 0 }, {  1,  0, 1 },
  {  1,  1, -1 }, {  1,  1, 0 }, {  1,  1, 1 }, };

// храниться скоростная модель
// хранятся контактные границы

class EnviromentModel {
public:
	EnviromentModel();
	~EnviromentModel();

	int InitSize(const int& I, const int& J, const int& K);
	int InitEnviromentFromFile(const char* File_Vp, const char* File_Vs, const char* File_Rho);
	int InitEnviromentFromArray(double* Vp, double* Vs, double* Rho);

	std::list<ContactBoundary>* Bounds();
	int* Layers();

	double* Vp();
	double* Vs();
	double* Rho();
	double h();
	Index Size();

	void Clear();

private:
	Index size_;

	double* vp_;
	double* vs_;
	double* rho_;

	std::list<ContactBoundary> bounds_;

	int* layers_;
	int count_layers_;

	std::map<std::pair<int, int>, bool> connect_layers_;
	std::map<std::pair<int, int>, ContactBoundary*> bound_on_connect_layers_;

	void FindLayers();
	void FindLayer(int i, int j, int k, int layer);
	void ConnectLayers();

	void FindBounds();
	bool SplitBound(ContactBoundary* bound);
	void ConnectLayersWithBounds();
};