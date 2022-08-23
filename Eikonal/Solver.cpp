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
