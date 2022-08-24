#include "Solver.h"

void FastSweep3D(double* grid, bool* frozen_cells, EnviromentModel* model) {
	const int n_sweep = 8;
	const int I = model->Size().i;
	const int J = model->Size().i;
	const int K = model->Size().i;
	const double h = model->h();

	double tmp;
	double aa[3];
	double f;
	double d_curr, d_new;
	long int grid_pos;

	// sweep directions { start, end, step }
	const int dir_x[n_sweep][3] = { {0, I - 1, 1}, {0, I - 1, 1},  {0, I - 1, 1},  {0, I - 1, 1},  {I - 1, 0, -1}, {I - 1, 0, -1}, {I - 1, 0, -1}, {I - 1, 0, -1} };
	const int dir_y[n_sweep][3] = { {0, J - 1, 1}, {0, J - 1, 1},  {J - 1, 0, -1}, {J - 1, 0, -1}, {0, J - 1, 1},  {0, J - 1, 1},  {J - 1, 0, -1}, {J - 1, 0, -1} };
	const int dir_z[n_sweep][3] = { {0, K - 1, 1}, {K - 1, 0, -1}, {0, K - 1, 1},  {K - 1, 0, -1}, {0, K - 1, 1},  {K - 1, 0, -1}, {0, K - 1, 1},  {K - 1, 0, -1} };

	for (int s = 0; s < n_sweep; s++) {
		for (int i = dir_x[s][0]; dir_x[s][2] * i <= dir_x[s][1]; i += dir_x[s][2]) {
			for (int j = dir_y[s][0]; dir_y[s][2] * j <= dir_y[s][1]; j += dir_y[s][2]) {
				for (int k = dir_z[s][0]; dir_z[s][2] * k <= dir_z[s][1]; k += dir_z[s][2]) {
					grid_pos = (i * J + j) * K + k;

					if (!frozen_cells[grid_pos]) {
						if ((i == 0) || (i == I - 1)) {
							if (i == 0) {
								aa[0] = grid[grid_pos] < grid[grid_pos + J * K] ? grid[grid_pos] : grid[grid_pos + J * K];
							}
							if (i == I - 1) {
								aa[0] = grid[grid_pos - J * K] < grid[grid_pos] ? grid[grid_pos - J * K] : grid[grid_pos];
							}
						} else {
							aa[0] = grid[grid_pos + J * K] < grid[grid_pos - J * K] ? grid[grid_pos + J * K] : grid[grid_pos - J * K];
						}

						if ((j == 0) || (j == J - 1)) {
							if (j == 0) {
								aa[1] = grid[grid_pos] < grid[grid_pos + K] ? grid[grid_pos] : grid[grid_pos + K];
							}
							if (j == J - 1) {
								aa[1] = grid[grid_pos - K] < grid[grid_pos] ? grid[grid_pos - K] : grid[grid_pos];
							}
						} else {
							aa[1] = grid[grid_pos + K] < grid[grid_pos - K] ? grid[grid_pos + K] : grid[grid_pos - K];
						}

						if ((k == 0) || (k == K - 1)) {
							if (k == 0) {
								aa[2] = grid[grid_pos] < grid[grid_pos + 1] ? grid[grid_pos] : grid[grid_pos + 1];
							}
							if (k == K - 1) {
								aa[2] = grid[grid_pos - 1] < grid[grid_pos] ? grid[grid_pos - 1] : grid[grid_pos];
							}
						} else {
							aa[2] = grid[grid_pos + 1] < grid[grid_pos - 1] ? grid[grid_pos + 1] : grid[grid_pos - 1];
						}

						if (aa[0] > aa[1]) { tmp = aa[0]; aa[0] = aa[1]; aa[1] = tmp; }
						if (aa[1] > aa[2]) { tmp = aa[1]; aa[1] = aa[2]; aa[2] = tmp; }
						if (aa[0] > aa[1]) { tmp = aa[0]; aa[0] = aa[1]; aa[1] = tmp; }

						f = model->Vp()[grid_pos];
						d_curr = aa[0] + h * f;
						if (d_curr < aa[1]) {
							d_new = d_curr;
						} else {
							double a = 2.0;
							double b = -2.0 * (aa[0] + aa[1]);
							double c = aa[0] * aa[0] + aa[1] * aa[1] - h * h * f * f;
							double D = sqrt(b * b - 4.0 * a * c);

							d_curr = ((-b + D) > (-b - D) ? (-b + D) : (-b - D)) / (2.0 * a);

							if (d_curr < aa[2]) {
								d_new = d_curr;
							} else {
								a = 3.0;
								b = -2.0 * (aa[0] + aa[1] + aa[2]);
								c = aa[0] * aa[0] + aa[1] * aa[1] + aa[2] * aa[2] - h * h * f * f;
								D = sqrt(b * b - 4.0 * a * c);
								d_new = ((-b + D) > (-b - D) ? (-b + D) : (-b - D)) / (2.0 * a);
							}

						}
						grid[grid_pos] = grid[grid_pos] < d_new ? grid[grid_pos] : d_new;
					}
				}
			}
		}
	}
}


Eikonal::Eikonal() {
	this->env_model_ = nullptr;
}

void Eikonal::SetModel(EnviromentModel* env_model) {
	env_model_ = env_model;
	this->ttf_.set_size(env_model_->Size());
}

void Eikonal::SetSourse(const double& x, const double& y, const double& z) {
	Task new_task;
	new_task.calculation_layers.push_back(1);
	tasks_.push_back(new_task);
}

void Eikonal::Calculate(const double& time_limit) {
	Task* current_task = nullptr;
	bool* frozen_cells = nullptr;
	double* grid = nullptr;
	while (!tasks_.empty()) {
		current_task = &tasks_.front();

		// *** //
		//FastSweep3D(grid, frozen_cells, env_model_);
		// *** //

		// analys solution task
		std::list<std::pair<int, int>> pairs;
		for (auto it_1 = current_task->calculation_layers.begin(); it_1 != current_task->calculation_layers.end(); it_1++) {
			for (auto it_2 = env_model_->ConnectLayer()->begin(); it_2 != env_model_->ConnectLayer()->end(); it_2++) {
				if (it_2->first.first == *it_1) {
					pairs.push_back(it_2->first);
				}
			}
		}

		for (auto it_1 = pairs.begin(); it_1 != pairs.end();) {
			for (auto it_2 = current_task->calculation_layers.begin(); it_2 != current_task->calculation_layers.end(); it_2++) {
				if (it_1->second == *it_2) {
					it_1 = pairs.erase(it_1);
					break;
				}
				it_1++;
			}
		}

		std::map<int, std::list<int>> layer_neigbor;
		std::list<std::list<int>> group_layer;
		for (auto it_1 = pairs.begin(); it_1 != pairs.end(); it_1++) {
			for (auto it_2 = pairs.begin(); it_2 != pairs.end(); it_2++) {
				if (env_model_->ConnectLayer()->count(std::make_pair(it_1->second, it_2->second))) {
					layer_neigbor[it_1->second].push_back(it_2->second);
				}
			}
		}

		for (auto it = layer_neigbor.begin(); it != layer_neigbor.end(); it++) {
			it->second.sort();
			it->second.unique();
		}

		if (!layer_neigbor.empty()) {
			bool is_end = false;
			group_layer.resize(1);
			group_layer.begin()->push_back(*layer_neigbor.at(0).begin());
			//while (!is_end) {
			for (auto it_1 = group_layer.begin(); it_1 != group_layer.end(); it_1++) {
				for (auto it_2 = it_1->begin(); it_2 != it_1->end(); it_2++) {
					for (auto it_3 = layer_neigbor[*it_2].begin(); it_3 != layer_neigbor[*it_2].end(); it_3++) {
						auto it = std::find(it_1->begin(), it_1->end(), *it_3);
						if (it == it_1->end()) {
							it_1->push_back(*it_3);
						}
					}
				}
				for (auto it_2 = layer_neigbor.begin(); it_2 != layer_neigbor.end(); it_2++) {
					auto it_3 = std::find(it_1->begin(), it_1->end(), it_2->first);
					if (it_3 == it_1->end()) {
						std::list<int> new_group;
						new_group.push_back(it_2->first);
						group_layer.push_back(new_group);
					}
				}
			}
			//}
		}

		// получили слой - список его соседей, слои и соседи являются соседями current_task->calculation_layers
		// нужно сгрупировать соседей, одна группа может не касаться другой!
		// по аналогии с разделением контактных границ
		// после этого необходимо расчитать всевозможные варианты расчетов
		// т.е. если в группе 3 слоя, нужно посчитать по отдельности, затем 1 с 2, 2 с 3, 3 с 1, и все три сразу,
		// при этом нужно учитывать, что в группе все слои не контактируют друг с другом сразу!


		std::map<int, Task> tasks;
		for (auto it_1 = pairs.begin(); it_1 != pairs.end(); it_1++) {
			tasks[it_1->second].add_front(env_model_->BoundOnConnectLayers()->at(*it_1));
		}
		for (auto it = tasks.begin(); it != tasks.end(); it++) {
			it->second.no_calculation_layers.insert(it->second.no_calculation_layers.end(), current_task->calculation_layers.begin(), current_task->calculation_layers.end());
			it->second.calculation_layers.push_back(it->first);
		}
		
		


		// нашли слои, с которыми контактирует наша область
		
		// зная, будет ли простое отражение или нет, добавляем новые задачи для расчета в найденных слоях

		// *** //
	}
	// - решаем уравнение в текущем слое
	// - смотрим, с какими областями контактирует рассматриваемый слой
	// - по очереди рассчитываем, как из текущего слоя распространится фронт в соседнем слое
	// - - если скорость в расмматриваемом слое больше, то просто рассчитываем в этой среде
}

void Eikonal::WriteToFile()
{

}
