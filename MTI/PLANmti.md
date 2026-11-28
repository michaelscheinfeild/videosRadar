# Coherent MTI Radar Simulation — Implementation Plan

https://chatgpt.com/share/6b0ad5be-8138-83ed-858e-bbb08e4f3fb9

## 1. Project objective

Develop a Python simulation of a coherent pulsed MTI (Moving Target Indication) radar.

The simulation will contain five targets:

- Two stationary buildings.
- Three moving aircraft.
- Two aircraft approaching the radar.
- One aircraft moving away from the radar.

The main objective is to understand the complete radar processing chain:

1. Generate transmitted radar pulses.
2. Propagate signals to targets and generate echoes.
3. Receive echoes at different ranges.
4. Perform coherent I/Q phase detection.
5. Observe pulse-to-pulse phase changes.
6. Apply a one-PRI delay line.
7. Subtract consecutive pulses using an MTI canceller.
8. Calculate Doppler spectra.
9. Compare radar results before and after MTI.
10. Animate the target positions.

The simulation must provide mathematical explanations and plots at every important processing stage.

---

# 2. Radar configuration

Use the following initial radar parameters.

| Parameter | Value |
|---|---:|
| Carrier frequency | 10 GHz |
| Speed of light | 299792458 m/s |
| Wavelength | 0.02998 m |
| PRF | 5000 Hz |
| PRI | 0.2 ms |
| Number of pulses | 256 |
| Pulse width | 2 µs |
| Fast-time sampling frequency | 20 MHz |
| Maximum simulated range | 15 km |
| Signal representation | Complex baseband I/Q |

### Important design consideration: PRF

The maximum unambiguous Doppler frequency is:

$$
f_{d,\max}=\frac{PRF}{2}=2500\ Hz
$$

The maximum unambiguous range is:

$$
R_{\max}=\frac{c}{2PRF}
$$

$$
R_{\max}\approx 29.98\ km
$$

This PRF allows us to observe all three aircraft Doppler frequencies without Doppler aliasing.

The carrier wavelength is:

$$
\lambda=\frac{c}{f_c}
$$

$$
\lambda\approx 0.02998\ m
$$

---

# 3. Target configuration

Create five targets.

| ID | Target | Initial range | Radial velocity | Behavior |
|---|---|---:|---:|---|
| T1 | Building 1 | 3000 m | 0 m/s | Stationary |
| T2 | Building 2 | 6000 m | 0 m/s | Stationary |
| T3 | Aircraft 1 | 4000 m | -30 m/s | Approaching |
| T4 | Aircraft 2 | 8000 m | -12 m/s | Approaching |
| T5 | Aircraft 3 | 10000 m | +20 m/s | Receding |

Velocity convention:

- Negative radial velocity: approaching radar.
- Positive radial velocity: moving away from radar.
- Zero radial velocity: stationary.

### Expected Doppler frequencies

For a monostatic radar:

$$
f_d=-\frac{2v_r}{\lambda}
$$

The negative sign follows our selected velocity convention.

Expected values:

| Target | Radial velocity | Expected Doppler |
|---|---:|---:|
| Building 1 | 0 | 0 Hz |
| Building 2 | 0 | 0 Hz |
| Aircraft 1 | -30 m/s | +2001.4 Hz |
| Aircraft 2 | -12 m/s | +800.6 Hz |
| Aircraft 3 | +20 m/s | -1334.3 Hz |

These frequencies must be verified against the simulated Doppler spectra.

---

# 4. Mathematical signal model

## 4.1 Transmitted signal

Generate a rectangular pulse in complex baseband:

$$
s_T(t)=p(t)
$$

where:

$$
p(t)=
\begin{cases}
1 & 0\leq t<T_p\\
0 & \text{otherwise}
\end{cases}
$$

The actual RF signal is represented by:

$$
s_{RF}(t)=
\Re\{s_T(t)e^{j2\pi f_ct}\}
$$

Do not directly sample the 10 GHz carrier.

Use complex baseband to simulate the radar coherently.

## 4.2 Propagation delay

For target range $R_i$:

$$
\tau_i=\frac{2R_i}{c}
$$

The received echo is delayed by the round-trip propagation time.

The target range can change between pulses:

$$
R_i[m]=R_i[0]+v_{r,i}mT
$$

where:

- $m$: pulse index.
- $T$: PRI.
- $v_{r,i}$: radial velocity.

## 4.3 Doppler phase

The Doppler frequency is:

$$
f_{d,i}=-\frac{2v_{r,i}}{\lambda}
$$

The phase progression between pulses is:

$$
\phi_i[m]=\phi_{0,i}+2\pi f_{d,i}mT
$$

The complex echo is:

$$
x_i[m,n]=
A_i p(nT_s-\tau_i[m])
e^{j\phi_i[m]}
$$

where:

- $A_i$: target amplitude.
- $T_s$: fast-time sampling interval.
- $n$: fast-time sample index.

Use fractional-delay interpolation when necessary to avoid artificial range-bin jumps.

## 4.4 Total received signal

All five echoes are added:

$$
x[m,n]=\sum_{i=1}^{5}x_i[m,n]+w[m,n]
$$

where $w[m,n]$ is complex Gaussian noise.

Include configurable noise level and random seed.

---

# 5. Processing stages and required plots

## Step 1 — Generate transmitted pulse

Generate the transmitted rectangular pulse.

Required plots:

1. Pulse envelope versus fast time.
2. I component.
3. Q component.
4. Pulse train showing several consecutive PRIs.

Explain the difference between:

- Fast time.
- Slow time.
- Carrier frequency.
- Pulse repetition frequency.

Output:

`01_transmitted_pulse.png`

## Step 2 — Generate received echoes

Calculate the propagation delay for each target.

$$
\tau_i=\frac{2R_i}{c}
$$

Display the five received echoes.

Required plots:

1. Received amplitude versus fast time.
2. Received amplitude versus range in km.
3. Mark the five target ranges.
4. Show one transmitted pulse and its delayed echoes on a common time axis.

Output:

`02_received_echoes.png`

Explain why different ranges correspond to different echo arrival times.

## Step 3 — Build the received data matrix

Create a matrix:

$$
X\in\mathbb{C}^{N_p\times N_r}
$$

where:

- Rows represent pulses.
- Columns represent range samples.

For this simulation:

$$
X\in\mathbb{C}^{256\times N_r}
$$

Required plots:

1. Matrix magnitude as a heatmap.
2. Horizontal axis: range.
3. Vertical axis: pulse number.
4. Identify stationary and moving target locations.

Output:

`03_received_data_cube.png`

## Step 4 — Coherent phase detector

Use complex I/Q video:

$$
z[m,n]=I[m,n]+jQ[m,n]
$$

The phase is:

$$
\phi[m,n]=\operatorname{atan2}(Q[m,n],I[m,n])
$$

The magnitude is:

$$
A[m,n]=\sqrt{I^2[m,n]+Q^2[m,n]}
$$

Extract the samples corresponding to each target's range gate.

Required plots:

1. I versus pulse number.
2. Q versus pulse number.
3. Magnitude versus pulse number.
4. Phase versus pulse number.

Make separate plots for all five targets.

Expected behavior:

- Buildings: constant phase.
- Approaching aircraft: positive phase progression.
- Receding aircraft: negative phase progression.

Output:

`04_phase_detector_IQ.png`

## Step 5 — Pulse-to-pulse phase progression

Calculate unwrapped phase:

$$
\phi_u[m]=\operatorname{unwrap}(\angle z[m])
$$

The phase slope is:

$$
\frac{d\phi}{dt}=2\pi f_d
$$

Estimate Doppler from phase slope:

$$
\hat f_d=\frac{1}{2\pi}
\frac{d\phi_u}{dt}
$$

Compare estimated Doppler against the theoretical Doppler.

Required plots:

1. Unwrapped phase versus slow time.
2. Five targets on the same plot.
3. Separate plots for each target.
4. Estimated Doppler versus theoretical Doppler.

Output:

`05_phase_progression.png`

## Step 6 — Doppler FFT

For each target range gate, calculate the FFT along slow time.

$$
X_D[k,n]=
\sum_{m=0}^{N_p-1}
x[m,n]e^{-j2\pi km/N_p}
$$

Use a Hann window before the FFT.

Use zero padding for a smoother frequency display.

The Doppler frequency axis is:

$$
f_k=\frac{k}{N_pT}
$$

with proper FFT shifting.

Required plots:

1. Doppler spectrum for Building 1.
2. Doppler spectrum for Building 2.
3. Doppler spectrum for Aircraft 1.
4. Doppler spectrum for Aircraft 2.
5. Doppler spectrum for Aircraft 3.

Mark the theoretical Doppler frequency on every plot.

Output:

`06_doppler_spectra.png`

Expected behavior:

- Buildings have peaks at 0 Hz.
- Aircraft have peaks at their respective Doppler frequencies.

## Step 7 — One-PRI delay line

Implement the delay line:

$$
x_d[m,n]=x[m-1,n]
$$

The delayed signal must be aligned with the current pulse at the same range sample.

Required plots:

1. Current pulse video.
2. Previous pulse video.
3. Overlay of the two signals.
4. Show stationary and moving targets separately.

Output:

`07_delay_line.png`

Explain that the delay line stores the previous pulse.

## Step 8 — MTI canceller

Implement:

$$
\boxed{
y[m,n]=x[m,n]-x[m-1,n]
}
$$

This is the two-pulse MTI canceller.

For stationary clutter:

$$
x[m,n]=x[m-1,n]
$$

Therefore:

$$
y[m,n]=0
$$

For a moving target:

$$
x[m,n]\neq x[m-1,n]
$$

Therefore:

$$
y[m,n]\neq0
$$

Required plots:

1. MTI output magnitude for each target.
2. Before-MTI versus after-MTI comparison.
3. Stationary buildings should be suppressed.
4. Moving aircraft should remain visible.

Output:

`08_mti_output.png`

## Step 9 — Doppler spectrum after MTI

Calculate Doppler FFT again, this time using the MTI output.

$$
Y_D[k,n]=FFT_m\{y[m,n]\}
$$

Compare:

- Doppler spectrum before MTI.
- Doppler spectrum after MTI.

Output:

`09_doppler_before_after.png`

The stationary zero-Doppler clutter should be strongly suppressed.

The moving aircraft should remain.

## Step 10 — Range-Doppler map

Calculate a two-dimensional range-Doppler map.

$$
RD(f_d,R)=
\left|
FFT_{\text{slow-time}}\{X\}
\right|
$$

Produce two maps:

1. Before MTI.
2. After MTI.

Plot:

- X axis: range (km).
- Y axis: Doppler frequency (Hz).
- Color: normalized signal magnitude in dB.

Output:

`10_range_doppler_maps.png`

This is one of the most important outputs of the project.

It should clearly show the two stationary targets at zero Doppler before MTI and their suppression after MTI.

## Step 11 — MTI filter frequency response

The two-pulse canceller transfer function is:

$$
H(z)=1-z^{-1}
$$

Its frequency response:

$$
H(f_d)=1-e^{-j2\pi f_dT}
$$

Magnitude:

$$
\boxed{
|H(f_d)|=2|\sin(\pi f_dT)|
}
$$

Plot the frequency response.

Explain:

- Zero Doppler null.
- Doppler response.
- Blind speeds.

Blind Doppler frequencies:

$$
f_{d,\text{blind}}=kPRF
$$

Blind radial velocities:

$$
\boxed{
v_{\text{blind},k}=
\frac{k\lambda PRF}{2}
}
$$

Output:

`11_mti_frequency_response.png`

---

# 6. Target animation

Create an animated top-down radar scene.

Use Matplotlib animation and save as GIF.

### Target colors

| Object | Color | Marker |
|---|---|---|
| Radar | Black | Triangle |
| Building 1 | Blue | Square |
| Building 2 | Blue | Square |
| Aircraft 1 | Green | Circle |
| Aircraft 2 | Green | Circle |
| Aircraft 3 | Green | Circle |

The two buildings remain fixed.

The two approaching aircraft move toward the radar.

The receding aircraft moves away.

Show:

- Target names.
- Current target positions.
- Radar location.
- Elapsed simulation time.
- Motion trails.

Output:

`12_target_motion.gif`

The animation is a schematic illustration of target movement, not a geographically accurate radar display.

---

# 7. Interactive animation and visualization

If possible, implement a second visualization using Matplotlib with a time slider.

The user should be able to:

- Pause animation.
- Resume animation.
- Move forward and backward in time.
- Display target positions.
- Display the current radar pulse index.

Also consider a synchronized visualization:

Left panel:
- Moving targets.

Right panel:
- Current range-Doppler map.

This will connect physical target movement to radar signal processing.

---

# 8. Validation and expected results

The simulation must verify the following.

### Test 1 — Stationary buildings

For both buildings:

$$
f_d=0
$$

Expected:

- Phase remains constant.
- Doppler peak appears at zero.
- MTI output is approximately zero, apart from noise.

### Test 2 — Approaching aircraft

For Aircraft 1:

$$
v_r=-30\,m/s
$$

Expected:

$$
f_d\approx+2001.4\,Hz
$$

For Aircraft 2:

$$
v_r=-12\,m/s
$$

Expected:

$$
f_d\approx+800.6\,Hz
$$

Expected:

- Positive phase progression.
- Positive Doppler peak.
- Nonzero MTI output.

### Test 3 — Receding aircraft

For Aircraft 3:

$$
v_r=+20\,m/s
$$

Expected:

$$
f_d\approx-1334.3\,Hz
$$

Expected:

- Negative phase progression.
- Negative Doppler peak.
- Nonzero MTI output.

### Test 4 — MTI clutter cancellation

Calculate the suppression ratio:

$$
CR=20\log_{10}
\left(
\frac{A_{\text{before}}}{A_{\text{after}}}
\right)
$$

Evaluate this at the two building range gates, accounting for the noise floor.

### Test 5 — Doppler estimation accuracy

For each moving target:

$$
e_f=|\hat f_d-f_d|
$$

Compare estimated and theoretical Doppler.

Include an accuracy table in the final report.

---

# 9. Python implementation structure

Recommended project structure:

```text
mti_simulation/
│
├── plan.md
├── main.py
│
├── radar_config.py
├── targets.py
├── signal_generation.py
├── propagation.py
├── phase_detector.py
├── delay_line.py
├── mti_filter.py
├── doppler_processing.py
├── range_doppler.py
├── animation.py
├── visualization.py
│
└── outputs/
    ├── 01_transmitted_pulse.png
    ├── 02_received_echoes.png
    ├── 03_received_data_cube.png
    ├── 04_phase_detector_IQ.png
    ├── 05_phase_progression.png
    ├── 06_doppler_spectra.png
    ├── 07_delay_line.png
    ├── 08_mti_output.png
    ├── 09_doppler_before_after.png
    ├── 10_range_doppler_maps.png
    ├── 11_mti_frequency_response.png
    └── 12_target_motion.gif
```

### Main Python libraries

- NumPy: signal generation and mathematical operations.
- SciPy: signal processing and filtering.
- Matplotlib: plots and animation.
- Matplotlib.animation: GIF generation.

---

# 10. Final deliverables

The completed simulation must provide:

1. Fully executable Python project.
2. Five-target radar scenario.
3. Complex-baseband coherent signal simulation.
4. Mathematical equations documented in code.
5. Plots after every processing stage.
6. Doppler analysis before and after MTI.
7. Range-Doppler maps.
8. Animated target movement.
9. Validation table comparing theoretical and estimated Doppler.
10. README explaining how to run the simulation.

## 11. Important implementation rules

- Keep fast time and slow time separate.
- Use the same range gate when comparing consecutive pulses.
- Preserve complex I/Q information.
- Do not take absolute values before the MTI subtraction.
- Do not use magnitude-only signals for coherent phase cancellation.
- Apply the delay line along the pulse dimension, not the range dimension.
- Use a consistent Doppler sign convention.
- Ensure Doppler frequencies are within the unambiguous Doppler interval.
- Distinguish FFT zero padding from actual Doppler resolution.
- Include noise in the simulation but allow disabling it for mathematical verification.
- Test the ideal noiseless stationary-target cancellation separately.

**Primary learning objective:** Understand how a radar can distinguish stationary buildings from moving aircraft by measuring the phase difference between consecutive coherent radar echoes.