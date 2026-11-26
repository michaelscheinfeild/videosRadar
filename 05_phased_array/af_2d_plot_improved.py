# af_2d_plot_improved.py
# AF(u,v) = 1/(NM) · Σₙ Σₘ wₙ wₘ · exp{ jk [ n·dx·(u−u₀) + m·dy·(v−v₀) ] }
# af_2d_plot.py
"""Animate the array factor of a rectangular N x M phased array in sine space."""

import os
import sys
from pathlib import Path

# Add native libraries when a Windows conda environment is not activated.
library_bin = Path(sys.prefix) / "Library" / "bin"
if sys.platform == "win32" and (Path(sys.prefix) / "conda-meta").is_dir() and library_bin.is_dir():
    os.environ["PATH"] = str(library_bin) + os.pathsep + os.environ.get("PATH", "")

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from scipy.constants import c, pi

# ---------------- Array parameters ----------------
F_0 = 10e9                      # center frequency (Hz)
WAVELENGTH = c / F_0            # m
K_0 = 2 * pi / WAVELENGTH       # wave number (rad/m)

N, M = 9, 9                     # elements along x and y
D_X = WAVELENGTH / 2            # x spacing (m)
D_Y = WAVELENGTH / 2            # y spacing (m)

THETA_0 = 20.0                   # steering angle from broadside (deg)
PHI_0 = 60.0                     # steering azimuth, from the x axis (deg)

WINDOW_X = np.ones(N)           # try np.hamming(N) to lower the sidelobes
WINDOW_Y = np.ones(M)           # try np.hamming(M)

# ---------------- Plot parameters ----------------
NPTS = 1001                     # odd, so u = 0 and v = 0 are sampled
Z_MIN, Z_MAX = -40, 0           # dB range
STRIDE = 5                      # mesh step in samples
GROW_FRAMES = 100               # frames used to grow the surface
TOTAL_FRAMES = 300              # growth + rotation
INTERVAL_MS = 30
DPI = 150
BACKGROUND = "#183340"


def af_1d(weights, d, U1, u0):
    """Normalized array factor of a uniformly spaced line of elements.

    weights: (n_elem,)   U1: (H, W)   returns: complex (H, W)
    """
    k_idx = np.arange(len(weights))
    phase = np.exp(1j * K_0 * d * k_idx[:, None, None] * (U1 - u0))  # (n_elem, H, W)
    return np.tensordot(weights, phase, axes=1) / weights.sum()


def compute_af_db(U, V, u0, v0):
    """Array factor in dB, clipped to [Z_MIN, Z_MAX]."""
    AF = af_1d(WINDOW_X, D_X, U, u0) * af_1d(WINDOW_Y, D_Y, V, v0)
    AF_db = 20 * np.log10(np.abs(AF) + 1e-12)   # field amplitude -> 20 log10
    return np.clip(AF_db, Z_MIN, Z_MAX)


def make_diagnostic_plots(u, U, V, Z, u0, v0):
    """Plot a v=0 cut and a top view of the array factor."""
    row = np.argmin(np.abs(u))
    cut = (
        af_1d(WINDOW_X, D_X, U[row:row + 1, :], u0)
        * af_1d(WINDOW_Y, D_Y, V[row:row + 1, :], v0)
    )
    magnitude = np.maximum(np.abs(cut[0]), 1e-12)
    floor = max(-25, Z_MIN)

    with plt.style.context("ggplot"):
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        axes[0].plot(u, 20 * np.log10(magnitude), label="20 log10 |AF| (field)")
        axes[0].plot(u, 10 * np.log10(magnitude), label="10 log10 |AF|")
        if np.allclose(WINDOW_X, WINDOW_X[0]):
            for k in range(1, len(WINDOW_X)):
                for sign in (-1, 1):
                    null_u = u0 + sign * k * WAVELENGTH / (len(WINDOW_X) * D_X)
                    if -1 <= null_u <= 1:
                        axes[0].axvline(null_u, color="gray", linestyle=":", linewidth=0.8)
        axes[0].set(
            title=f"1D cut at v = {u[row]:.3f}", xlabel="u", ylabel="dB",
            xlim=(-1, 1), ylim=(Z_MIN, Z_MAX),
        )
        axes[0].legend()

        cmap = cm.coolwarm.copy()
        cmap.set_bad("white")
        image = axes[1].pcolormesh(
            U, V, np.ma.masked_less_equal(Z, floor),
            cmap=cmap, vmin=floor, vmax=Z_MAX, shading="auto",
        )
        axes[1].add_patch(plt.Circle((0, 0), 1, fill=False, color="black", linestyle=":"))
        axes[1].set(
            title=f"Top view (white = below {floor:g} dB)",
            xlabel="u", ylabel="v", xlim=(-1, 1), ylim=(-1, 1),
        )
        axes[1].set_aspect("equal")
        fig.colorbar(image, ax=axes[1], label="dB")
        fig.tight_layout()
    return fig


def make_animation(U, V, Z):
    plt.style.use("ggplot")
    fig = plt.figure(figsize=(10, 10), facecolor=BACKGROUND)
    ax = fig.add_subplot(111, projection="3d", facecolor=BACKGROUND)
    ax.set(xlim=(-1, 1), ylim=(-1, 1), zlim=(Z_MIN, Z_MAX))
    ax.set_box_aspect(None, zoom=0.85)
    ax.set_axis_off()

    center = NPTS // 2
    surface = [None]

    def draw(half):
        """Replace the surface with a centered square of half-width `half` samples."""
        if surface[0] is not None:
            surface[0].remove()
        s = slice(center - half, center + half + 1)
        surface[0] = ax.plot_surface(
            U[s, s], V[s, s], Z[s, s],
            rstride=STRIDE, cstride=STRIDE,
            cmap=cm.coolwarm, vmin=Z_MIN, vmax=Z_MAX,   # fixed color scale
            antialiased=False,
        )

    def update(frame):
        if frame < GROW_FRAMES:
            # grow in whole strides so the mesh does not shimmer
            cells = max(1, round((frame + 1) / GROW_FRAMES * (center // STRIDE)))
            draw(cells * STRIDE)
        ax.view_init(elev=30, azim=-50 + 360 * frame / TOTAL_FRAMES)

    # init_func is a no-op so that frame 0 is not drawn twice
    ani = animation.FuncAnimation(
        fig, update, frames=TOTAL_FRAMES, interval=INTERVAL_MS, init_func=lambda: []
    )
    return fig, ani


def main():
    theta, phi = np.radians(THETA_0), np.radians(PHI_0)
    u0 = np.sin(theta) * np.cos(phi)
    v0 = np.sin(theta) * np.sin(phi)

    u = np.linspace(-1, 1, NPTS)
    U, V = np.meshgrid(u, u)
    Z = compute_af_db(U, V, u0, v0)

    output_dir = Path(__file__).resolve().parent / "static"
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_path = output_dir / "af_cuts.png"
    diagnostic_fig = make_diagnostic_plots(u, U, V, Z, u0, v0)
    diagnostic_fig.savefig(plot_path, dpi=DPI)
    plt.close(diagnostic_fig)
    print(f"Wrote {plot_path}", flush=True)

    fig, ani = make_animation(U, V, Z)

    output_path = output_dir / "af_2d_anim_rotating.mp4"
    ani.save(
        output_path,
        bitrate=-1,
        dpi=DPI,
        progress_callback=lambda i, total: print(f"Saved: {i + 1}/{total}"),
    )
    plt.close(fig)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()