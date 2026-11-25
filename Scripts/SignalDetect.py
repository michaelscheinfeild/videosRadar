# https://chatgpt.com/share/6b06a5d4-f2b8-83ed-b67f-5a11b2d9534a


import os
import sys
from pathlib import Path

# Add native libraries when a Windows conda environment is not activated.
library_bin = Path(sys.prefix) / "Library" / "bin"
if sys.platform == "win32" and (Path(sys.prefix) / "conda-meta").is_dir() and library_bin.is_dir():
    os.environ["PATH"] = str(library_bin) + os.pathsep + os.environ.get("PATH", "")

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve, find_peaks, spectrogram


def plot_rx(Fs, rx):
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(np.arange(len(rx)) / Fs * 1e6, np.abs(rx))
    ax.set_xlabel("Time (µs)")
    ax.set_ylabel("Amplitude")
    ax.set_title("Received Pulsed LFM Signal")
    ax.grid(True)
    #plt.show()


def plot_rx_fft(Fs, rx):
    # The received signal is complex, so include both positive and negative frequencies.
    frequencies = np.fft.fftshift(np.fft.fftfreq(len(rx), d=1 / Fs))
    magnitude = np.abs(np.fft.fftshift(np.fft.fft(rx))) / len(rx)

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(frequencies / 1e6, magnitude)
    ax.set_xlabel("Baseband frequency (MHz)")
    ax.set_ylabel("FFT magnitude / N")
    ax.set_title("Received Pulsed LFM Signal FFT")
    ax.grid(True)
    plt.show()


def plot_rx_spectogram(Fs, rx):
    # Overlapping windows show how the complex baseband frequencies change over time.
    frequencies, times, magnitude = spectrogram(
        rx, fs=Fs, window="hann", nperseg=256, noverlap=224,
        detrend=False, return_onesided=False, scaling="spectrum", mode="magnitude",
    )
    frequencies = np.fft.fftshift(frequencies)
    magnitude = np.fft.fftshift(magnitude, axes=0)
    peak = max(magnitude.max(), np.finfo(float).tiny)
    magnitude_db = 20 * np.log10(np.maximum(magnitude / peak, 1e-6))

    fig, ax = plt.subplots(figsize=(10, 4))
    image = ax.pcolormesh(
        times * 1e6, frequencies / 1e6, magnitude_db,
        shading="auto", vmin=-60, vmax=0,
    )
    ax.set_xlabel("Time (µs)")
    ax.set_ylabel("Baseband frequency (MHz)")
    ax.set_title("Received Pulsed LFM Signal Spectrogram")
    fig.colorbar(image, ax=ax, label="Magnitude (dB relative to peak)")
    fig.tight_layout()
    #plt.show()


def plot_compressed(ax, range_axis, range_profile, targets=None):
    """Plot the matched-filter pulse-compression range profile."""
    ax.plot(range_axis, range_profile)
    ax.set_xlim(0, 150)
    ax.set_title("Pulsed LFM: Matched Filter / Pulse Compression")
    ax.set_xlabel("Range (m)")
    ax.set_ylabel("Correlation magnitude")
    ax.grid(True)

    if targets is not None:
        for target in targets:
            ax.axvline(target, color="red", linestyle="--", alpha=0.5)

def plot_compressed_and_targets(range_axis, range_profile, targets=None):
    plt.figure()
    plt.plot(range_axis, range_profile)
    plt.xlabel("Range (m)")
    plt.ylabel("Correlation magnitude")
    plt.grid(True)
    if targets is not None:
        for target in targets:
            plt.axvline(target, color="red", linestyle="--", alpha=0.5)
        plt.title("Range Profile targets vs Estimated: " + ", ".join([f"{t:.1f} m" for t in targets]))
    plt.show()

# ============================================================
# 1. RADAR PARAMETERS
# ============================================================

c = 3e8
f0 = 9e9
B = 100e6
T = 20e-6
Fs = 200e6

S = B / T                    # Chirp slope (Hz/sec)
N = int(Fs * T)
t = np.arange(N) / Fs

targets = np.array([50.0, 80.0, 120.0])  # meters
amplitudes = np.array([1.0, 0.8, 0.7])

print("Chirp slope:", S / 1e12, "THz/sec")
print("Range resolution:", c / (2 * B), "m")

# Complex baseband LFM chirp.
# We do not directly sample the 9 GHz carrier.
tx = np.exp(1j * np.pi * S * t**2)


# ============================================================
# 2. PULSED LFM: SIMULATE RECEIVED ECHOES
# ============================================================

delays = 2 * targets / c
delay_samples = np.round(delays * Fs).astype(int)

rx_length = N + delay_samples.max() + 1
rx = np.zeros(rx_length, dtype=complex)

for delay, amplitude in zip(delay_samples, amplitudes):
    rx[delay:delay + N] += amplitude * tx


# plot signal


if 0:
    plot_rx(Fs, rx)
    plot_rx_fft(Fs, rx)
    plot_rx_spectogram(Fs, rx)
    plt.show()

# Matched filter: correlate received signal with transmitted chirp.
matched_filter = np.conj(tx[::-1])

compressed = fftconvolve(rx, matched_filter, mode="full")



# Correlation lag in samples.
lags = np.arange(len(compressed)) - (N - 1)

# Keep positive delays.
valid = lags >= 0
range_axis = lags[valid] / Fs * c / 2
range_profile = np.abs(compressed[valid])

#---------------------
if 0:
    plot_compressed_and_targets(range_axis, range_profile, targets=targets)
#---------------------

# Detect peaks.
peaks, _ = find_peaks(
    range_profile,
    distance=20,
    prominence=0.10 * np.max(range_profile)
)

print("\nPULSE COMPRESSION RESULTS")
for p in peaks:
    print(f"Estimated range: {range_axis[p]:.2f} m")


# ============================================================
# 3. FMCW: SIMULATE RECEIVED CHIRP ECHOES
# ============================================================

rx_fmcw = np.zeros(N, dtype=complex)

for R, amplitude in zip(targets, amplitudes):

    tau = 2 * R / c

    # Delayed chirp phase.
    delayed_phase = np.exp(
        1j * np.pi * S * (t - tau)**2
    )

    # Echo exists only after its propagation delay.
    valid_echo = (t >= tau) & (t < T + tau)

    rx_fmcw[valid_echo] += amplitude * delayed_phase[valid_echo]


# ============================================================
# 4. DECHIRPING: CREATE THE IF SIGNAL
# ============================================================

# Use the portion where all three echoes are present.
max_delay = delays.max()
valid_time = t >= max_delay

t_if = t[valid_time]

# Mixing TX with conjugated RX.
if_signal = tx[valid_time] * np.conj(rx_fmcw[valid_time])

# The beat frequencies should be S * 2R/c.
expected_beats = S * delays

print("\nEXPECTED IF BEAT FREQUENCIES")
for R, fb in zip(targets, expected_beats):
    print(f"Target {R:.1f} m: {fb / 1e6:.3f} MHz")


# ============================================================
# 5. FFT OF IF SIGNAL
# ============================================================

window = np.hanning(len(if_signal))

NFFT = 65536

spectrum = np.fft.fft(if_signal * window, n=NFFT)
frequencies = np.fft.fftfreq(NFFT, d=1/Fs)

positive = frequencies >= 0

frequencies = frequencies[positive]
spectrum = np.abs(spectrum[positive])

# Convert beat frequency to target range.
range_if = c * frequencies / (2 * S)

# Detect spectral peaks.
if_peaks, _ = find_peaks(
    spectrum,
    distance=100,
    prominence=0.10 * np.max(spectrum)
)

print("\nFMCW FFT RESULTS")
for p in if_peaks:
    print(
        f"Estimated range: {range_if[p]:.2f} m, "
        f"Beat frequency: {frequencies[p] / 1e6:.3f} MHz"
    )


# ============================================================
# 6. PLOTS
# ============================================================

fig, ax = plt.subplots(2, 1, figsize=(11, 9))

# Pulse compression range profile.
plot_compressed(ax[0], range_axis, range_profile, targets)

# FMCW range profile.
ax[1].plot(range_if, spectrum)
ax[1].set_xlim(0, 150)
ax[1].set_title("FMCW: IF Signal FFT / Range Profile")
ax[1].set_xlabel("Range (m)")
ax[1].set_ylabel("FFT magnitude")
ax[1].grid(True)

for R in targets:
    ax[1].axvline(R, color="red", linestyle="--", alpha=0.5)

plt.tight_layout()
plt.show()
