#include "Containers.h"

int WaveFront::capacity_step = 100;

void WaveFront::add_node(const Index& _coord, const double& _value) {
	if ((wave_front_coords == nullptr) && (travel_time_values == nullptr)) {
		wave_front_coords = new Index[capacity_step];
		travel_time_values = new double[capacity_step];
		wave_front_coords[0] = _coord;
		count = 1;
		capacity += capacity_step;
	}
	else {
		if (count < capacity) {
			wave_front_coords[count] = _coord;
			travel_time_values[count] = _value;
			count++;
		}
		else {
			capacity += capacity_step;
			Index* _temp_1 = new Index[capacity];
			double* _temp_2 = new double[capacity];
			for (int i = 0; i < count; ++i) {
				_temp_1[i] = wave_front_coords[i];
				_temp_2[i] = travel_time_values[i];
			}

			_temp_1[count] = _coord;
			_temp_2[count] = _value;
			++count;
			delete[] wave_front_coords;
			delete[] travel_time_values;
			wave_front_coords = _temp_1;
			travel_time_values = _temp_2;
			_temp_1 = nullptr;
			_temp_2 = nullptr;
		}
	}
}

void WaveFront::remove_data() {
	delete[] wave_front_coords;
	delete[] travel_time_values;
}