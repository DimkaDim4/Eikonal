#include "EnviromentModel.h"

EnviromentModel::EnviromentModel() {
	size_ = Index();
	vp_ = nullptr;
	vs_ = nullptr;
	rho_ = nullptr;
	layers_ = nullptr;
	count_layers_ = 0;
}

EnviromentModel::~EnviromentModel() {
	Clear();
}

void EnviromentModel::FindLayers() {
	if (layers_ != nullptr) {
		delete[] layers_;
	}

	layers_ = new (std::nothrow) int[size_.i * size_.j * size_.k];
	if (layers_ == nullptr) {
		std::cout << "Failed! EnviromentModel::layers_ = new (std::nothrow) int[size_.i * size_.j * size_.k]";
		return;
	}

	for (long int i = 0; i < size_.i * size_.j * size_.k; ++i) {
		layers_[i] = 0;
	}
	count_layers_ = 0;

	bool is_end = false;
	long int grid_pos;

	for (int i = 0; i < size_.i; ++i) {
		for (int j = 0; j < size_.j; ++j) {
			for (int k = 0; k < size_.k; ++k) {
				grid_pos = (i * size_.j + j) * size_.k + k;
				if (layers_[grid_pos] == 0) {
					++count_layers_;
					layers_[grid_pos] = count_layers_;
					FindLayer(i, j, k, count_layers_);
				}
			}
		}
	}
	ConnectLayers();
}

void EnviromentModel::FindLayer(int i, int j, int k, int layer) {
	Index coord(i, j, k);
	Index coord_neighbor;
	long int grid_pos;
	long int grid_pos_neighbor;

	std::list<Index> nodes_to_view;
	nodes_to_view.push_back(coord);

	while (!nodes_to_view.empty()) {
		coord = nodes_to_view.front();
		grid_pos = (coord.i * size_.j + coord.j) * size_.k + coord.k;
		for (int s = 0; s < count_neighbors; ++s) {
			if ((coord.i == 0) && (neighboring_nodes[s][0] == -1)) { continue; }
			if ((coord.j == 0) && (neighboring_nodes[s][1] == -1)) { continue; }
			if ((coord.k == 0) && (neighboring_nodes[s][2] == -1)) { continue; }
			if ((coord.i == size_.i - 1) && (neighboring_nodes[s][0] == 1)) { continue; }
			if ((coord.j == size_.j - 1) && (neighboring_nodes[s][1] == 1)) { continue; }
			if ((coord.k == size_.k - 1) && (neighboring_nodes[s][2] == 1)) { continue; }

			coord_neighbor.i = coord.i + neighboring_nodes[s][0];
			coord_neighbor.j = coord.j + neighboring_nodes[s][1];
			coord_neighbor.k = coord.k + neighboring_nodes[s][2];
			grid_pos_neighbor = (coord_neighbor.i * size_.j + coord_neighbor.j) * size_.k + coord_neighbor.k;

			if (layers_[grid_pos_neighbor] == 0) {
				if ((vp_[grid_pos_neighbor] == vp_[grid_pos]) &&
					(vs_[grid_pos_neighbor] == vs_[grid_pos]) &&
					(rho_[grid_pos_neighbor] == rho_[grid_pos])) {
					layers_[grid_pos_neighbor] = layer;
					nodes_to_view.push_back(coord_neighbor);
				}
			}
		}
		nodes_to_view.pop_front();
	}
}

void EnviromentModel::ConnectLayers() {
	if (layers_ == nullptr) {
		FindLayers();
	}

	Index coord_neighbor;
	long int grid_pos;
	long int grid_pos_neighbor;
	int layer;
	int layer_neighbor;
	for (int i = 0; i < size_.i; ++i) {
		for (int j = 0; j < size_.j; ++j) {
			for (int k = 0; k < size_.k; ++k) {
				grid_pos = (i * size_.j + j) * size_.k + k;
				for (int s = 0; s < count_neighbors; ++s) {
					if ((i == 0) && (neighboring_nodes[s][0] == -1)) { continue; }
					if ((j == 0) && (neighboring_nodes[s][1] == -1)) { continue; }
					if ((k == 0) && (neighboring_nodes[s][2] == -1)) { continue; }
					if ((i == size_.i - 1) && (neighboring_nodes[s][0] == 1)) { continue; }
					if ((j == size_.j - 1) && (neighboring_nodes[s][1] == 1)) { continue; }
					if ((k == size_.k - 1) && (neighboring_nodes[s][2] == 1)) { continue; }

					coord_neighbor.i = i + neighboring_nodes[s][0];
					coord_neighbor.j = j + neighboring_nodes[s][1];
					coord_neighbor.k = k + neighboring_nodes[s][2];
					grid_pos_neighbor = (coord_neighbor.i * size_.j + coord_neighbor.j) * size_.k + coord_neighbor.k;

					layer = layers_[grid_pos];
					layer_neighbor = layers_[grid_pos_neighbor];
					if (layer != layer_neighbor) {	
						bool is_refracted = false;
						if (vp_[grid_pos] < vp_[grid_pos_neighbor]) {
							is_refracted = true;		
						}
						connect_layers_[std::pair<int, int>(layer, layer_neighbor)] = is_refracted;					
					}
				}
			}
		}
	}
}



void EnviromentModel::FindBounds() {
	if (layers_ == nullptr) {
		FindLayers();
	}

	bool is_find_cb;
	long int grid_pos;
	long int grid_pos_neighbor;

	for (int i = 0; i < size_.i; ++i) {
		for (int j = 0; j < size_.j; ++j) {
			for (int k = 0; k < size_.k; ++k) {
				grid_pos = (i * size_.j + j) * size_.k + k;
				for (int s = 0; s < count_neighbors; ++s) {
					if ((i == 0) && (neighboring_nodes[s][0] == -1)) { continue; }
					if ((j == 0) && (neighboring_nodes[s][1] == -1)) { continue; }
					if ((k == 0) && (neighboring_nodes[s][2] == -1)) { continue; }
					if ((i == size_.i - 1) && (neighboring_nodes[s][0] == 1)) { continue; }
					if ((j == size_.j - 1) && (neighboring_nodes[s][1] == 1)) { continue; }
					if ((k == size_.k - 1) && (neighboring_nodes[s][2] == 1)) { continue; }

					grid_pos_neighbor = ((i + neighboring_nodes[s][0]) * size_.j + (j + neighboring_nodes[s][1])) * size_.k + (k + neighboring_nodes[s][2]);

					if (layers_[grid_pos] != layers_[grid_pos_neighbor]) {
						is_find_cb = false;
						for (auto it = bounds_.begin(); it != bounds_.end(); it++) {
							if (it->Layer() == layers_[grid_pos]) {
								it->add_node(Index(i, j, k));
								is_find_cb = true;
								break;
							}
						}

						if (!is_find_cb) {
							bounds_.push_back(ContactBoundary());
							bounds_.back().set_layer(layers_[grid_pos]);
							bounds_.back().add_node(Index(i, j, k));
						}
						break;
					}
				}
			}
		}
	}

	int count = bounds_.size();
	auto it = bounds_.begin();
	for (int i = 0; i < count; i++) {
		if (SplitBound(&*it)) {
			it = bounds_.erase(it);
		} else {
			it++;
		}
	}

	ConnectLayersWithBounds();
}

bool EnviromentModel::SplitBound(ContactBoundary* bound) {
	int layer = bound->Layer();
	std::list<Index> nodes_to_view;
	nodes_to_view.push_back(bound->Coords()[0]);

	std::map<Index, bool> nodes_viewed;
	for (long int i = 0; i < bound->Count(); ++i){
		nodes_viewed.insert(std::pair<Index, bool>(bound->Coords()[i], false));
	}

	std::list<Index> nodes_in_split_bound;

	Index coord;
	Index coord_neighbor;
	Index coord_nr_nr;
	bool is_bound_node;

	bool is_end = false;
	while (!is_end) {
		while (!nodes_to_view.empty()) {
			coord = nodes_to_view.front();
			for (int s = 0; s < count_neighbors; ++s) {
				if ((coord.i == 0) && (neighboring_nodes[s][0] == -1)) { continue; }
				if ((coord.j == 0) && (neighboring_nodes[s][1] == -1)) { continue; }
				if ((coord.k == 0) && (neighboring_nodes[s][2] == -1)) { continue; }
				if ((coord.i == size_.i - 1) && (neighboring_nodes[s][0] == 1)) { continue; }
				if ((coord.j == size_.j - 1) && (neighboring_nodes[s][1] == 1)) { continue; }
				if ((coord.k == size_.k - 1) && (neighboring_nodes[s][2] == 1)) { continue; }

				coord_neighbor.i = coord.i + neighboring_nodes[s][0];
				coord_neighbor.j = coord.j + neighboring_nodes[s][1];
				coord_neighbor.k = coord.k + neighboring_nodes[s][2];

				if (layers_[(coord_neighbor.i * size_.j + coord_neighbor.j) * size_.k + coord_neighbor.k] == layer) {
					is_bound_node = false;
					for (int ss = 0; ss < count_neighbors; ++ss) {
						if ((coord_neighbor.i == 0) && (neighboring_nodes[ss][0] == -1)) { continue; }
						if ((coord_neighbor.j == 0) && (neighboring_nodes[ss][1] == -1)) { continue; }
						if ((coord_neighbor.k == 0) && (neighboring_nodes[ss][2] == -1)) { continue; }
						if ((coord_neighbor.i == size_.i - 1) && (neighboring_nodes[ss][0] == 1)) { continue; }
						if ((coord_neighbor.j == size_.j - 1) && (neighboring_nodes[ss][1] == 1)) { continue; }
						if ((coord_neighbor.k == size_.k - 1) && (neighboring_nodes[ss][2] == 1)) { continue; }

						coord_nr_nr.i = coord_neighbor.i + neighboring_nodes[ss][0];
						coord_nr_nr.j = coord_neighbor.j + neighboring_nodes[ss][1];
						coord_nr_nr.k = coord_neighbor.k + neighboring_nodes[ss][2];

						if (layers_[(coord_nr_nr.i * size_.j + coord_nr_nr.j) * size_.k + coord_nr_nr.k] != layer) {
							is_bound_node = true;
							break;
						}
					}
					if (is_bound_node) {
						auto node = nodes_viewed.find(coord_neighbor);
						if (!node->second) {
							nodes_to_view.push_back(coord_neighbor);
							node->second = true;
						}
					}
				}
			}
			nodes_viewed.find(coord)->second = true;
			nodes_in_split_bound.push_back(coord);
			nodes_to_view.pop_front();
		}

		if (nodes_in_split_bound.size() == bound->Count()) {
			return false;
		}
		else {
			ContactBoundary* new_bound = new ContactBoundary;
			new_bound->set_layer(layer);
			for (auto it = nodes_in_split_bound.begin(); it != nodes_in_split_bound.end(); it++) {
				new_bound->add_node(*it);
			}
			bounds_.push_back(*new_bound);
			nodes_in_split_bound.clear();

			for (auto it  = nodes_viewed.begin(); it != nodes_viewed.end(); it++){
				if (!it->second) {
					nodes_to_view.push_back(it->first);
					it->second = true;
					break;
				}
			}
		}
		if (nodes_to_view.empty()) {
			is_end = true;
		}
	}
	return true;
}

void EnviromentModel::ConnectLayersWithBounds() {
	int layer_1, layer_2;
	Index coord;
	Index coord_neighbor;

	bool* is_viewed = new bool[count_layers_];
	for (int i = 0; i < count_layers_; i++) {
		is_viewed[i] = false;
	}

	for (auto i = bounds_.begin(); i != bounds_.end(); i++) {
		layer_1 = i->Layer();
		for (int j = 0; j < i->Count(); j++) {
			coord = i->Coords()[j];
			for (int s = 0; s < count_neighbors; ++s) {
				if ((coord.i == 0) && (neighboring_nodes[s][0] == -1)) { continue; }
				if ((coord.j == 0) && (neighboring_nodes[s][1] == -1)) { continue; }
				if ((coord.k == 0) && (neighboring_nodes[s][2] == -1)) { continue; }
				if ((coord.i == size_.i - 1) && (neighboring_nodes[s][0] == 1)) { continue; }
				if ((coord.j == size_.j - 1) && (neighboring_nodes[s][1] == 1)) { continue; }
				if ((coord.k == size_.k - 1) && (neighboring_nodes[s][2] == 1)) { continue; }

				coord_neighbor.i = coord.i + neighboring_nodes[s][0];
				coord_neighbor.j = coord.j + neighboring_nodes[s][1];
				coord_neighbor.k = coord.k + neighboring_nodes[s][2];

				layer_2 = layers_[(coord_neighbor.i * size_.j + coord_neighbor.j) * size_.k + coord_neighbor.k];

				if ((layer_1 != layer_2) && !is_viewed[layer_2 - 1]) {
					bound_on_connect_layers_[std::pair<int, int>(layer_1, layer_2)] = &*i;
					is_viewed[layer_2 - 1] = true;
				}
			}
		}
		for (int i = 0; i < count_layers_; i++) {
			is_viewed[i] = false;
		}
	}
}



int EnviromentModel::InitSize(const int& I, const int& J, const int& K) {
	if ((I > 0) && (J > 0) && (K > 0)) {
		size_ = Index(I, J, K);
		return 0;
	}
	return -1;
}

int EnviromentModel::InitEnviromentFromFile(const char* File_Vp, const char* File_Vs, const char* File_Rho) {
	return 0;
}

int EnviromentModel::InitEnviromentFromArray(double* vp, double* vs, double* rho) {
	if ((vp == nullptr) || (vs == nullptr) || (rho == nullptr)) {
		return -1;
	}

	vp_ = vp;
	vs_ = vs;
	rho_ = rho;
	return 0;
}

std::list<ContactBoundary>* EnviromentModel::Bounds() {
	if (bounds_.empty()) {
		FindBounds();
	}
	return &bounds_;
}

int* EnviromentModel::Layers()
{
	if (layers_ == nullptr)	{
		FindLayers();
	}
	return layers_;
}

double* EnviromentModel::Vp() {
	return this->vp_;
}

double* EnviromentModel::Vs() {
	return this->vs_;
}

double* EnviromentModel::Rho() {
	return rho_;
}

double EnviromentModel::h() {
	return 1. / (size_.i - 1);
}

Index EnviromentModel::Size() {
	return size_;
}

void EnviromentModel::Clear() {
	size_ = Index();
	vp_ = nullptr;
	vs_ = nullptr;
	rho_ = nullptr;
	bounds_.clear();
}