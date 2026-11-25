from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import chirp, spectrogram

# Generate a 10 ms chirp that rises linearly from 500 Hz to 3500 Hz.
fs = 20000
duration = 0.01
t = np.linspace(0, duration, int(fs * duration), endpoint=False)
x = chirp(t, f0=500, t1=duration, f1=3500, method="linear")

# Create the output folder next to this script.
output_dir = Path(__file__).resolve().with_name("Figures")
output_dir.mkdir(parents=True, exist_ok=True)

# Plot the pulse in time; its oscillations become faster as frequency rises.
fig, ax = plt.subplots(figsize=(8, 3))
ax.plot(t * 1e3, x)
ax.set_xlabel("Time (ms)")
ax.set_ylabel("Amplitude")
ax.set_title("Linear Chirp (0.5 - 3.5 kHz)")
ax.grid(True)
fig.tight_layout()
fig.savefig(output_dir / "SampleChirp.png", dpi=150)

# Compute the positive-frequency FFT of the entire pulse.
frequencies = np.fft.rfftfreq(len(x), d=1 / fs)
spectrum = np.abs(np.fft.rfft(x)) / len(x)

# Plot the frequency content across the pulse, without time information.
fig, ax = plt.subplots(figsize=(8, 3))
ax.plot(frequencies / 1e3, spectrum)
ax.set_xlabel("Frequency (kHz)")
ax.set_ylabel("FFT magnitude / N")
ax.set_title("Chirp FFT")
ax.set_xlim(0, 5)
ax.grid(True)
fig.tight_layout()
fig.savefig(output_dir / "SampleChirpFFT.png", dpi=150)

# Analyze overlapping 2 ms windows to see frequency change over time.
spectrogram_frequencies, spectrogram_times, magnitude = spectrogram(
    x, fs=fs, window="hann", nperseg=40, noverlap=36, nfft=512,
    scaling="spectrum", mode="magnitude",
)

# Express magnitude in dB relative to the strongest time-frequency bin.
magnitude_db = 20 * np.log10(np.maximum(magnitude / magnitude.max(), 1e-6))

# Plot time against frequency; the rising band shows the chirp sweep.
fig, ax = plt.subplots(figsize=(8, 3))
image = ax.pcolormesh(
    spectrogram_times * 1e3, spectrogram_frequencies / 1e3, magnitude_db,
    shading="auto", cmap="viridis", vmin=-60, vmax=0,
)
ax.set_xlabel("Time (ms)")
ax.set_ylabel("Frequency (kHz)")
ax.set_title("Chirp Spectrogram")
ax.set_ylim(0, 5)
fig.colorbar(image, ax=ax, label="Magnitude (dB relative to peak)")
fig.tight_layout()
fig.savefig(output_dir / "SampleChirpSpectrogram.png", dpi=150)

# Display all three figures until their windows are closed.
print(f"Saved figures to {output_dir}")
plt.show()