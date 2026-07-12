from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class LegacyEccentricityMapStyle:
    figure_size: tuple[float, float]
    figure_dpi: int
    x_limits: tuple[float, float]
    y_limits: tuple[float, float]
    x_grid_stop: float
    y_grid_stop: float
    top_ticks: tuple[float, ...]
    top_tick_labels: tuple[str, ...]
    panel_labels: tuple[str, ...]
    label_font_size: float
    tick_label_size: float
    top_tick_label_size: float
    major_tick_length: float
    major_tick_width: float
    minor_tick_length: float
    minor_tick_width: float
    top_minor_tick_width: float
    resonance_line_width: float
    colorbar_position: tuple[float, float, float, float]
    colorbar_tick_length: float
    subplot_wspace: float
    subplot_hspace: float
    subplot_right: float


LUNA_FIG02_STYLE = LegacyEccentricityMapStyle(figure_size=(14.0, 12.0), figure_dpi=300, x_limits=(3.5, 9.0), y_limits=(0.0, 0.3), x_grid_stop=9.01, y_grid_stop=0.31, top_ticks=(5.9678, 5.3222, 6.2119, 5.6769, 4.8769, 4.2948, 4.6070, 3.9278, 8.75), top_tick_labels=("8:1", "6:1", "9:1", "7:1", "5:1", "4:1", "9:2", "7:2", r"$\beta_\max$"), panel_labels=("a", "b", "c", "d", "e", "f"), label_font_size=40.0, tick_label_size=20.0, top_tick_label_size=10.0, major_tick_length=12.0, major_tick_width=4.0, minor_tick_length=6.0, minor_tick_width=2.0, top_minor_tick_width=4.0, resonance_line_width=2.0, colorbar_position=(0.92, 0.11, 0.015, 0.77), colorbar_tick_length=8.0, subplot_wspace=0.1, subplot_hspace=0.15, subplot_right=0.90)
LEGACY_UNSTABLE_FATES = ("a<0", "a>0.4R_H", "collision", "d<R_roche")
