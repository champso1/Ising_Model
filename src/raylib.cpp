#include "isingmodel.hpp"

#include "raylib.h"

namespace ising
{
	void IsingModel::show_simulation(
		double equil_temp,
		int win_w, int win_h,
		char const* win_title,
		int refresh_rate
	) {
		set_temp(equil_temp);
		
		auto n = _grid.size(); // for easier typing
		InitWindow(win_w, win_h, win_title);
		SetTargetFPS(refresh_rate);

		while (!WindowShouldClose()) {
			BeginDrawing();

			// want to be able to restart it if we press "r"
			if (IsKeyPressed(KEY_R))
				_grid.fill(_grid_init_flag);

			ClearBackground(BLACK);
			for (size_t i=0; i<n; i++) {
				for (size_t j=0; j<n; j++) {
					Color color = (_grid.at(i,j) > 0) ? WHITE : BLACK;
					DrawRectangle(j*(win_w/n), i*(win_h/n), (win_w/n), (win_h/n), color);
				}
			}
			update_grid();

			EndDrawing();
		}

		CloseWindow();
	}
}
