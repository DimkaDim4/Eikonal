#include "Containers.h"

int ContactBoundary::capacity_step = 100;

void ContactBoundary::add_node(const Index& coord) {
	if (contact_boundary_coords == nullptr) {
		contact_boundary_coords = new Index[capacity_step];
		contact_boundary_coords[0] = coord;
		count = 1;
		capacity += capacity_step;
	}
	else {
		if (count < capacity) {
			contact_boundary_coords[count] = coord;
			count++;
		}
		else {
			capacity += capacity_step;
			Index* temp = new Index[capacity];
			for (int i = 0; i < count; ++i) {
				temp[i] = contact_boundary_coords[i];
			}

			temp[count] = coord;
			++count;
			delete[] contact_boundary_coords;
			contact_boundary_coords = temp;
			temp = nullptr;
		}
	}
}

void ContactBoundary::remove_data() {
	delete[] contact_boundary_coords;
}

int ContactBoundary::Layer() {
	return layer;
}

void ContactBoundary::set_layer(const int &layer) {
	if (layer >= 0) {
		this->layer = layer;
	}
}
