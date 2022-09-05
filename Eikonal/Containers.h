#pragma once
#include <iostream>

struct Index {
	int i;
	int j;
	int k;

	Index() { i = 0; j = 0; k = 0; }
	Index(int _i, int _j, int _k) : i(_i), j(_j), k(_k) {}

	bool operator==(const Index &index) const {
		return !((i != index.i) || ((j != index.j)) || ((k != index.k)));
	}

	bool operator<(const Index& index) const {
		if (i < index.i) { return true; }
		if ((i == index.i) && (j < index.j)) { return true; }
		if ((i == index.i) && (j == index.j) && (k < index.k)) { return true; }
		return false;
	}
};

class ContactBoundary {
public:
	ContactBoundary() {}
	~ContactBoundary() { delete[] contact_boundary_coords; }

	void add_node(const Index& coord);
	void remove_data();

	long int contains(const Index &coord) {
		for (long int i = 0; i < count; i++) {
			if (contact_boundary_coords[i] == coord) {
				return i;
			}
		}
		return -1;
	}

	int Layer();
	void set_layer(const int &layer);

	Index* Coords() { return contact_boundary_coords; }
	int Count() { return count; }

	static int capacity_step;

private:
	Index* contact_boundary_coords = nullptr;
	int count = 0;
	int capacity = 0;
	int layer = 0;
};
