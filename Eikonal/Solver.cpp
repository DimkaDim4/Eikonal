#include "Solver.h"

void FastSweep3D(double* grid, bool* frozen_cells, EnviromentModel* model, int type_wave) {
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

						//if (type_wave == 0) {
							f = 1. / model->Vp()[grid_pos];
						//} else {
							//f = 1. / model->Vs()[grid_pos];
						//}

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

std::list<std::list<int>> gen(const std::list<int>& list) {
	std::vector<std::vector<bool>> gen_comb;
	std::vector<bool> element;
	element.resize(list.size());
	gen_comb.resize(pow(2, list.size()) - 1);
	for (int i = 0; i < list.size(); i++) {
		element[i] = false;
	}

	for (int i = 0; i < pow(2, list.size()) - 1; i++) {
		int j = list.size() - 1;
		while (element[j]) {
			element[j] = false;
			j--;
		}
		element[j] = true;
		gen_comb[i] = element;
	}

	std::list<std::list<int>> result;
	result.resize(pow(2, list.size()) - 1);

	int i = 0;
	for (auto it_1 = result.begin(); it_1 != result.end(); it_1++) {
		int j = 0;
		for (auto it = list.begin(); it != list.end(); it++) {
			if (gen_comb[i][j]) {
				it_1->push_back(*it);
			}
			j++;
		}
		i++;
	}

	return result;
}

Eikonal::Eikonal() {
	this->env_model_ = nullptr;
}

void Eikonal::SetModel(EnviromentModel* env_model) {
	env_model_ = env_model;
}

void Eikonal::AddSourse(const double& x, const double& y, const double& z) {
	Index coord;
	Index size = env_model_->Size();
	double h = env_model_->h();
	for (int i = 0; i < size.i; i++) {
		if (x <= i * h) {
			coord.i = i;
			break;
		}
	}
	for (int i = 0; i < size.j; i++) {
		if (y <= i * h) {
			coord.j = i;
			break;
		}
	}
	for (int i = 0; i < size.k; i++) {
		if (z <= i * h) {
			coord.k = i;
			break;
		}
	}

	Task new_task;
	new_task.calculation_layers.push_back(env_model_->Layers()[(coord.i * size.j + coord.j) * size.k + coord.k]);
	new_task.travel_time_values = new double[1];
	new_task.travel_time_values[0] = 0.;
	new_task.wave_front_coords = new Index[1];
	new_task.wave_front_coords[0] = coord;
	new_task.size_wave_front = 1;
	tasks_.push_back(new_task);
}

// наложить сетку источников
// именновать типы волн - отраженная, прямая, преломленная
// пользователь должен иметь возможность выбирать, какие волны мы хотим расчитывать


void Eikonal::Calculate(const double& time_limit) {
	int num_file = 0;
	Task* current_task = nullptr;
	Index size = env_model_->Size();
	Index* coords = nullptr;
	double* values = nullptr;

	bool* frozen_cells = new bool[size.i * size.j * size.k];
	double* grid = new double[size.i * size.j * size.k];
	const int* layers = env_model_->Layers();

	std::ofstream result_calculation;
	result_calculation.open("fronts.txt");

	int num_task = 0;
	while (!tasks_.empty()) {
		num_task++;
		current_task = &tasks_.front();
		for (long int i = 0; i < size.i * size.j * size.k; i++) {
			frozen_cells[i] = true;
		}

		for (auto it = current_task->calculation_layers.begin(); it != current_task->calculation_layers.end(); it++) {
			for (long int i = 0; i < size.i * size.j * size.k; i++) {
				if (layers[i] == *it) {
					frozen_cells[i] = false;
				}
			}
		}

		for (long int i = 0; i < size.i * size.j * size.k; i++) {
			grid[i] = time_limit * 100.;
		}
		
		Index coord;
		long int pos;
		for (long int i = 0; i < current_task->size_wave_front; ++i) {
			coord = current_task->wave_front_coords[i];
			pos = (coord.i * size.j + coord.j) * size.k + coord.k;
			grid[pos] = current_task->travel_time_values[i];
			frozen_cells[pos] = true;
		}

		FastSweep3D(grid, frozen_cells, env_model_, 0);
		for (auto it = current_task->calculation_layers.begin(); it != current_task->calculation_layers.end(); it++) {
			current_task->travel_layers.push_back(*it);
			std::string filename = "data" + std::to_string(num_file + 1) + ".bin";
			std::ofstream out(filename, std::ios::binary | std::ios::out);
			for (long int i = 0; i < size.i * size.j * size.k; i++) {
				if (layers[i] == *it) {
					out.write((char*)&grid[i], sizeof(double));
				}
			}
			out.close();
			num_file++;
			current_task->filenames.push_back(filename);
		}

		auto it_filenames = current_task->filenames.begin();
		auto it_travel_layers = current_task->travel_layers.begin();

		result_calculation << "namefile:\t";
		for (int i = 0; i < current_task->filenames.size(); i++) {
			result_calculation << *it_filenames << "\t";
			it_filenames++;
		}
		result_calculation << "\n";

		result_calculation << "layers:\t";
		for (int i = 0; i < current_task->travel_layers.size(); i++) {
			result_calculation << *it_travel_layers << "\t";
			it_travel_layers++;
		}
		result_calculation << "\n\n";

		double max_value_in_calc_layers = 0.;
		for (auto it = current_task->calculation_layers.begin(); it != current_task->calculation_layers.end(); it++) {
			for (long int i = 0; i < size.i * size.j * size.k; i++) {
				if (layers[i] == *it) {
					max_value_in_calc_layers = max_value_in_calc_layers > grid[i] ? max_value_in_calc_layers : grid[i];
				}
			}
		}

		std::cout << "task:\t" << num_task << "\n";
		std::cout << "\tcalculation layers:\n";
		for (auto it = current_task->calculation_layers.begin(); it != current_task->calculation_layers.end(); it++) {
			std::cout << "\t\t" << *it << "\n";
		}

		if (current_task->no_calculation_layers.size() != 0) {
			std::cout << "\tno calculation layers:\n";
			for (auto it = current_task->no_calculation_layers.begin(); it != current_task->no_calculation_layers.end(); it++) {
				std::cout << "\t\t" << *it << "\n";
			}
		}
		std::cout << "\tmax value of calc layers: " << max_value_in_calc_layers << '\n';

		// ищем следующий слой, в котором мы будем расчитывать фронт
		const std::pair<int, int>* pair = nullptr;
		for (auto it = env_model_->ConnectLayer()->begin(); it != env_model_->ConnectLayer()->end(); it++) {
			auto it_1 = std::find(current_task->calculation_layers.begin(), current_task->calculation_layers.end(), it->first.first);
			auto it_2 = std::find(current_task->calculation_layers.begin(), current_task->calculation_layers.end(), it->first.second);
			auto it_3 = std::find(current_task->no_calculation_layers.begin(), current_task->no_calculation_layers.end(), it->first.second);
			if ((it_1 != current_task->calculation_layers.end()) &&
				(it_2 == current_task->calculation_layers.end()) &&
				(it_3 == current_task->no_calculation_layers.end())) {
				pair = &it->first;
				break;
			}
		}

		if (pair) {
			Task new_task_1;
			new_task_1.travel_layers = current_task->travel_layers;
			new_task_1.no_calculation_layers = current_task->no_calculation_layers;

			bool is_refracted = env_model_->ConnectLayer()->find(*pair)->second;
			if (is_refracted) {
				for (auto it = current_task->calculation_layers.begin(); it != current_task->calculation_layers.end(); it++) {
					if (*it != pair->first) {
						new_task_1.no_calculation_layers.push_back(*it);
					}
				}

				new_task_1.calculation_layers.push_back(pair->first);
				new_task_1.calculation_layers.push_back(pair->second);

				new_task_1.filenames = current_task->filenames;
				new_task_1.filenames.pop_back();
				new_task_1.travel_layers.pop_back();

				new_task_1.size_wave_front = current_task->size_wave_front;
				new_task_1.travel_time_values = current_task->travel_time_values;
				new_task_1.wave_front_coords = current_task->wave_front_coords;
			} else {
				new_task_1.calculation_layers.push_back(pair->second);
				new_task_1.no_calculation_layers.push_back(pair->first);

				new_task_1.filenames = current_task->filenames;

				ContactBoundary* front = env_model_->BoundOnConnectLayers()->at(*pair);
				new_task_1.size_wave_front = front->Count();
				new_task_1.wave_front_coords = front->Coords();
				new_task_1.travel_time_values = new double[new_task_1.size_wave_front];
				for (long int i = 0; i < new_task_1.size_wave_front; ++i) {
					coord = new_task_1.wave_front_coords[i];
					new_task_1.travel_time_values[i] = grid[(coord.i * size.j + coord.j) * size.k + coord.k];
				}
			}

			tasks_.push_back(new_task_1);

			// отраженная волна
			Task new_task_2;
			new_task_2.calculation_layers.push_back(pair->first);
			new_task_2.no_calculation_layers.push_back(pair->second);
			new_task_2.filenames = current_task->filenames;
			new_task_2.travel_layers = current_task->travel_layers;

			ContactBoundary* front = env_model_->BoundOnConnectLayers()->at(*pair);
			new_task_2.size_wave_front = front->Count();
			new_task_2.wave_front_coords = front->Coords();
			new_task_2.travel_time_values = new double[new_task_2.size_wave_front];
			double min = time_limit * 100.;
			for (int i = 0; i < new_task_2.size_wave_front; i++) {
				coord = new_task_2.wave_front_coords[i];
				double value = grid[(coord.i * size.j + coord.j) * size.k + coord.k];
				new_task_2.travel_time_values[i] = value;
				min = min < value ? min : value;
			}

			if (min < time_limit) {
				tasks_.push_back(new_task_2);
			}

			//FastSweep3D(grid, frozen_cells, env_model_, 1);
			//ttf_.add_traveltime(grid, frozen_cells);

			/*min = 0.;
			new_task_2.travel_time_values.clear();
			for (int i = 0; i < front->Count(); i++) {
				Index coord = coords[i];
				double value = grid[(coord.i * size.j + coord.j) * size.k + coord.k];
				new_task_2.travel_time_values.push_back(value);
				min = min > value ? min : value;
			}

			if (min < time_limit) {
				tasks_.push_back(new_task_2);
			}*/

			//front = nullptr;
			//coords = nullptr;
		}
		tasks_.pop_front();
		//system("CLS");

		// расчитать отражение в текущих слоях
		// расчитать прямую волну в текущими и новом слоях

		// analys solution task
		// нахоим слои, контактирующие с current_task->calculation_layers
		// игнорируем слои из no_calculation_layers;
		/*std::list<std::pair<int, int>> pairs; 
		for (auto it_1 = current_task->calculation_layers.begin(); it_1 != current_task->calculation_layers.end(); it_1++) {
			for (auto it_2 = env_model_->ConnectLayer()->begin(); it_2 != env_model_->ConnectLayer()->end(); it_2++) {
				if (it_2->first.first == *it_1) {
					auto it = std::find(current_task->no_calculation_layers.begin(), current_task->no_calculation_layers.end(), it_2->first.second);
					if (it == current_task->no_calculation_layers.end()) {
						pairs.push_back(it_2->first);
					}
				}
			}
		}*/

		/*for (auto it_1 = pairs.begin(); it_1 != pairs.end();) {
			bool is_deleted = false;
			for (auto it_2 = current_task->calculation_layers.begin(); it_2 != current_task->calculation_layers.end(); it_2++) {
				if (it_1->second == *it_2) {
					it_1 = pairs.erase(it_1);
					is_deleted = true;
					break;
				}
			}
			if (!is_deleted) it_1++;
		}*/

		//////////////////////
		// слой, контактирующий с calculation_layers - список его соседец, так же контактирующих с calculation_layers
		/*std::map<int, std::list<int>> layer_neigbor;
		std::list<std::list<int>> group_layer;
		for (auto it_1 = pairs.begin(); it_1 != pairs.end(); it_1++) {
			for (auto it_2 = pairs.begin(); it_2 != pairs.end(); it_2++) {
				if (env_model_->ConnectLayer()->count(std::make_pair(it_1->second, it_2->second))) {
					layer_neigbor[it_1->second].push_back(it_2->second);
				}
			}
		}*/

		/*for (auto it = layer_neigbor.begin(); it != layer_neigbor.end(); it++) {
			it->second.sort();
			it->second.unique();
		}*/

		// получили слой - список его соседей, слои и соседи являются соседями current_task->calculation_layers
		// нужно сгрупировать соседей, одна группа может не касаться другой!
		// по аналогии с разделением контактных границ

		//if (!layer_neigbor.empty()) {
		//	//bool is_end = false;
		//	group_layer.resize(1);
		//	group_layer.begin()->push_back(*layer_neigbor.begin()->second.begin());
		//	for (auto it_1 = group_layer.begin(); it_1 != group_layer.end(); it_1++) {
		//		for (auto it_2 = it_1->begin(); it_2 != it_1->end(); it_2++) {
		//			for (auto it_3 = layer_neigbor[*it_2].begin(); it_3 != layer_neigbor[*it_2].end(); it_3++) {
		//				auto it = std::find(it_1->begin(), it_1->end(), *it_3);
		//				if (it == it_1->end()) {
		//					it_1->push_back(*it_3);
		//				}
		//			}
		//		}
		//		for (auto it_2 = layer_neigbor.begin(); it_2 != layer_neigbor.end(); it_2++) {
		//			auto it_3 = std::find(it_1->begin(), it_1->end(), it_2->first);
		//			if (it_3 == it_1->end()) {
		//				std::list<int> new_group;
		//				new_group.push_back(it_2->first);
		//				group_layer.push_back(new_group);
		//			}
		//		}
		//	}
		//}

		
		// после этого необходимо расчитать всевозможные варианты расчетов
		// т.е. если в группе 3 слоя, нужно посчитать по отдельности, затем 1 с 2, 2 с 3, 3 с 1, и все три сразу,
		// при этом нужно учитывать, что в группе все слои не контактируют друг с другом сразу!
		
		// для каждой группы необходимо сгенерировать сочетания без повторений из n по m для каждой группы
		// n = 1..size, m = size

		/*std::list<std::list<std::list<int>>> gen_comb;
		for (auto it = group_layer.begin(); it != group_layer.end(); it++) {
			gen_comb.push_back(gen(*it));
		}*/

		// сгенерировали комбинации
		// теперь нужно убрать комбинации слоев, которые между собой не контактируют

		/*for (auto it_1 = gen_comb.begin(); it_1 != gen_comb.end(); it_1++) {
			for (auto it_2 = it_1->begin(); it_2 != it_1->end();) {
				if (it_2->size() == 1) { it_2++; continue; }
				std::list<int> comb_connect;
				comb_connect.push_back(it_2->front());

				for (auto it_3 = comb_connect.begin(); it_3 != comb_connect.end(); it_3++) {
					for (auto it_4 = layer_neigbor[*it_3].begin(); it_4 != layer_neigbor[*it_3].end(); it_4++) {
						auto it_5 = std::find(it_2->begin(), it_2->end(), *it_4);
						auto it_6 = std::find(comb_connect.begin(), comb_connect.end(), *it_4);
						if ((it_5 != it_2->end()) && (it_6 == comb_connect.end())) {
							comb_connect.push_back(*it_4);
						}
					}
				}
				if (comb_connect.size() != it_2->size()) {
					it_2 = it_1->erase(it_2);
				} else {
					it_2++;
				}
			}
		}*/

		// для сгенерированных областей нуобходимо создать задачи
		// ининциализировать начальный фронт
		// слои, в которых будут вестись расчеты
		// слои, в которых мы уже все посчитали
		// добавить инит фронт, 
		// обойти лист pairs, найти те пары, у которых второй равен одному и списка calculation_layers

		/*for (auto it_1 = gen_comb.begin(); it_1 != gen_comb.end(); it_1++) {
			for (auto it_2 = it_1->begin(); it_2 != it_1->end(); it_2++) {
				Task new_task;
				new_task.calculation_layers = *it_2;
				new_task.no_calculation_layers = current_task->no_calculation_layers;
				new_task.no_calculation_layers.insert(new_task.no_calculation_layers.end(), current_task->calculation_layers.begin(), current_task->calculation_layers.end());

				for (auto it = layer_neigbor.begin(); it != layer_neigbor.end(); it++) {
					auto it_ = std::find(it_2->begin(), it_2->end(), it->first);
					if (it_ == it_2->end()) {
						new_task.no_calculation_layers.push_back(it->first);
					}
				}

				for (auto it_3 = pairs.begin(); it_3 != pairs.end(); it_3++) {
					for (auto it_4 = it_2->begin(); it_4 != it_2->end(); it_4++) {
						if (it_3->second == *it_4) {
							new_task.add_front(env_model_->BoundOnConnectLayers()->at(*it_3));
						}
					}
				}
				tasks_.push_back(new_task);
			}
		}*/

		//tasks_.pop_front();
		// нашли слои, с которыми контактирует наша область
		
		// зная, будет ли простое отражение или нет, добавляем новые задачи для расчета в найденных слоях

		// *** //
	}
	// - решаем уравнение в текущем слое
	// - смотрим, с какими областями контактирует рассматриваемый слой
	// - по очереди рассчитываем, как из текущего слоя распространится фронт в соседнем слое
	// - - если скорость в расмматриваемом слое больше, то просто рассчитываем в этой среде

	current_task = nullptr;
	coords = nullptr;
	values = nullptr;

	delete[] frozen_cells;
	delete[] grid;
	layers = nullptr;

	result_calculation.close();

	//ttf_.write();
}
