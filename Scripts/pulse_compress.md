# Pulse Compression and FMCW Range Processing

## Radar chirps, LFM signals, matched filtering, beat frequency, and FFT

This note collects the main ideas from our discussion. It compares
**pulsed LFM pulse compression** with **FMCW range processing**, derives
the equations, explains the Python signal models, and highlights
implementation details that can otherwise cause confusion.

------------------------------------------------------------------------

## 1. What is a chirp? What is LFM?

A **chirp** is a signal whose instantaneous frequency changes with time.
**LFM** means *Linear Frequency Modulation*: the instantaneous frequency
changes linearly.

For an up-chirp,

\[ f(t)=f_0+St, \]

where:

-   (f_0) is the starting RF frequency,
-   \(S\) is the chirp slope in Hz/s,
-   \(t\) is time measured from the beginning of the chirp.

The swept bandwidth is

\[ B=ST, `\qquad `{=tex}S=`\frac{B}{T}`{=tex}, \]

where (T) is the chirp duration.

A chirp is a general signal category; LFM is one particular type of
chirp. A nonlinear-FM chirp is also possible, but the equations in this
note concern LFM.

### Figure 1 --- Linear frequency sweep

``` text
Instantaneous
frequency
    ^
    |                         /  f(t) = f0 + S t
    |                      /
    |                   /
    |                /
 f0 +|-------------/
    |
    +---------------------------------> time
    0                                T

             swept bandwidth B = S T
```

------------------------------------------------------------------------

## 2. Why does an LFM signal contain t² in its phase?

A complex signal can be written as

\[ s(t)=A e\^{j`\phi`{=tex}(t)}, \]

where (`\phi`{=tex}(t)) is its phase in radians. Its instantaneous
frequency is

\[ f(t)=`\frac{1}{2\pi}`{=tex}`\frac{d\phi(t)}{dt}`{=tex}. \]

For an LFM up-chirp, choose the phase

\[ `\phi`{=tex}(t)=2`\pi `{=tex}f_0t+`\pi `{=tex}St\^2. \]

Differentiating gives

\[ f(t)=`\frac{1}{2\pi}`{=tex}(2`\pi `{=tex}f_0+2`\pi `{=tex}St)
=f_0+St. \]

That is the desired linear frequency sweep.

### Complex baseband representation

In our simulation, we do not directly generate or sample the 9 GHz RF
carrier. We represent the signal at complex baseband and omit the
carrier term (2`\pi `{=tex}f_0t):

\[ `\boxed{s_{\mathrm{TX}}(t)=e^{j\pi St^2}}`{=tex} \]

This has baseband instantaneous frequency

\[ `\boxed{f_{\mathrm{TX,bb}}(t)=St}`{=tex}. \]

The omitted RF carrier does not mean the physical radar has no carrier.
It means the fast carrier oscillation is represented by a complex
envelope instead of being sampled directly.

Python:

``` python
# Complex baseband LFM up-chirp
tx = np.exp(1j * np.pi * S * t**2)
```

The magnitude is constant:

\[ \|s\_{`\mathrm{TX}`{=tex}}(t)\|=1. \]

The quadratic term is in the **phase**, not the amplitude.

------------------------------------------------------------------------

## 3. Propagation delay: why does the received signal use (t - tau)?

For a stationary target at range (R), the transmitted waveform travels
to the target and back. The round-trip delay is

\[ `\boxed{\tau=\frac{2R}{c}}`{=tex}, \]

where (c`\approx 3`{=tex}`\times10`{=tex}\^8 `\mathrm{m/s}`{=tex}) is
the speed of light.

A delayed copy of the transmitted complex-baseband chirp is

\[ `\boxed{s_{\mathrm{RX}}(t)=A e^{j\pi S(t-\tau)^2}}`{=tex}, \]

where (A) represents the target echo amplitude (and can also be complex
if it includes an additional phase).

Python:

``` python
# Delayed LFM phase for a stationary target
delayed_phase = np.exp(
    1j * np.pi * S * (t - tau)**2
)
```

The two expressions have the same waveform shape. The second is shifted
in time:

-   Transmitted: (e^{j`\pi `{=tex}St^2})
-   Received: (e^{j`\pi `{=tex}S(t-`\tau`{=tex})^2})

The instantaneous baseband frequency of the received chirp is

\[ f\_{`\mathrm{RX,bb}`{=tex}}(t)
=`\frac{1}{2\pi}`{=tex}`\frac{d}{dt}`{=tex}`\left[\pi S(t-\tau)^2\right]`{=tex}=S(t-`\tau`{=tex}).
\]

Therefore, at the same observation time,

\[ f\_{`\mathrm{TX,bb}`{=tex}}(t)-f\_{`\mathrm{RX,bb}`{=tex}}(t)
=St-S(t-`\tau`{=tex}) =`\boxed{S\tau}`{=tex}. \]

### Figure 2 --- Transmit and delayed receive chirps

``` text
Frequency
   ^
   |                         TX: f_TX(t) = S t
   |                       /
   |                     /
   |                   /
   |                 /       vertical separation = S tau
   |               /       /
   |             /       /  RX: f_RX(t) = S(t - tau)
   |           /       /
   +--------------------------------------> time
                 same observation time

The received sweep is delayed by tau. For stationary-target LFM,
the two instantaneous frequencies differ by the constant S*tau.
```

**Important:** The delayed phase expression is a simplified
complex-baseband echo model. A fuller coherent RF-to-baseband model can
also include carrier propagation phase, target Doppler, antenna/channel
effects, and amplitude loss. The simplified model is useful for
understanding the basic beat-frequency relationship.

------------------------------------------------------------------------

## 4. Two ways to create a delayed echo in Python

### 4.1 Integer-sample delay by copying samples

``` python
N = len(tx)
delay_samples = int(round(tau * Fs))

rx_length = N + delay_samples
rx = np.zeros(rx_length, dtype=complex)

rx[delay_samples:delay_samples + N] += A * tx
```

This implements the discrete-time delay

\[ rx\[n\]=A,tx\[n-d\], \]

where (d) is an integer number of samples. The receive array is made
long enough to hold the entire delayed pulse.

For multiple targets:

``` python
rx = np.zeros(rx_length, dtype=complex)

for delay, amplitude in zip(delay_samples_all, amplitudes):
    rx[delay:delay + N] += amplitude * tx
```

The echoes are added together in one received array:

\[ rx\[n\]=`\sum`{=tex}\_i A_i,tx\[n-d_i\]. \]

If echoes overlap in time, their complex samples add.

**Limitation:** Rounding (`\tau `{=tex}F_s) to an integer loses
fractional-sample delay accuracy. A fractional-delay filter or
interpolation is needed when that error matters.

### 4.2 Evaluate the delayed phase at each sample time

``` python
tau = 2 * R / c
delayed_phase = np.exp(1j * np.pi * S * (t - tau)**2)
```

This evaluates the mathematical delayed waveform at the actual sample
times. It can represent a fractional-sample delay because `tau` need not
be an integer multiple of (1/F_s).

For a finite chirp that is nonzero only over (0`\leq `{=tex}t\<T), the
delayed echo is supported over

\[ `\tau`{=tex}`\leq `{=tex}t\<T+`\tau`{=tex}. \]

For example:

``` python
valid_echo = (t >= tau) & (t < T + tau)
rx[valid_echo] += A * delayed_phase[valid_echo]
```

The time array must extend far enough to include the echo tail up to
(T+`\tau`{=tex}). If `t` stops at (T), only the portion of the echo
before (T) is represented.

**Do not confuse the two methods:** integer-sample copying and
evaluating the delayed phase describe the same ideal delay, but have
different numerical accuracy and array-window requirements.

------------------------------------------------------------------------

## 5. Method A --- Pulsed LFM and matched-filter pulse compression

A pulsed radar transmits a finite-duration LFM pulse. A target returns a
delayed, attenuated copy. The receiver correlates the received signal
with a time-reversed complex-conjugated reference waveform
(equivalently, applies a matched filter).

One correlation convention is

\[ y(`\tau`{=tex})=`\int `{=tex}r(t)s\^\*(t-`\tau`{=tex}),dt. \]

The output becomes large when the reference waveform aligns with an
echo. The peak location estimates the delay:

\[ `\hat{\tau}`{=tex}=`\text{delay at correlation peak}`{=tex}. \]

Convert that delay to range:

\[ `\boxed{\hat R=\frac{c\hat\tau}{2}}`{=tex}. \]

### Figure 3 --- Pulsed LFM processing chain

``` mermaid
flowchart TD
    A["Transmit finite LFM pulse"] --> B["Target reflects pulse"]
    B --> C["Receive delayed echo(s)"]
    C --> D["Matched filter / correlation"]
    D --> E["Compressed-pulse peaks"]
    E --> F["Peak delay tau"]
    F --> G["Range R = c tau / 2"]
```

### Why is it called pulse compression?

The transmitted pulse can be relatively long, which provides energy, but
its bandwidth allows a narrow correlation peak. The matched filter
produces a narrow peak in delay: a long pulse is compressed in the
range-response domain.

The ideal range resolution for a sufficiently well-designed LFM pulse is
approximately

\[ `\boxed{\Delta R\approx\frac{c}{2B}}`{=tex}. \]

The time-bandwidth product is

\[ BT. \]

It describes how much the pulse can be compressed in time. An idealized
processing gain is often expressed as approximately
(10`\log`{=tex}\_{10}(BT)) dB, subject to definitions, windowing, and
implementation losses.

### Sidelobes

A matched-filter response generally has a mainlobe and sidelobes.
Sidelobes are responses around a strong target peak; they are not
automatically additional targets. Windowing or mismatched filtering can
reduce sidelobes, usually with a trade-off such as a wider mainlobe or
some SNR loss.

------------------------------------------------------------------------

## 6. Method B --- FMCW and IF (beat-frequency) processing

FMCW radar transmits a chirp ramp and receives delayed echoes while
using a coherent copy/reference of the transmitted chirp. The receiver
mixes (dechirps) the received signal with the transmit reference.

For one stationary target, use

\[ s\_{`\mathrm{TX}`{=tex}}(t)=e^{j`\pi `{=tex}St^2} \]

and

\[ s\_{`\mathrm{RX}`{=tex}}(t)=A e^{j`\pi `{=tex}S(t-`\tau`{=tex})^2}.
\]

One useful complex mixing convention is

\[ IF(t)=s\_{`\mathrm{TX}`{=tex}}(t)s\_{`\mathrm{RX}`{=tex}}\^\*(t). \]

Substitute:

\[ IF(t) =Ae^{j`\pi `{=tex}St^2}e^{-j`\pi `{=tex}S(t-`\tau`{=tex})^2}.
\]

Expand the square:

\[ (t-`\tau`{=tex})^2=t^2-2t`\tau`{=tex}+`\tau`{=tex}\^2. \]

Then

\[
```{=tex}
\begin{aligned}
IF(t)
&=A e^{j\pi S[t^2-(t^2-2t\tau+\tau^2)]}\\
&=A e^{j2\pi S\tau t}e^{-j\pi S\tau^2}.
\end{aligned}
```
\]

Define

\[ C=Ae^{-j`\pi `{=tex}S`\tau`{=tex}^2}. \]

For a stationary target, (A) and (`\tau`{=tex}) are constant during this
observation, so (C) is constant with respect to fast time (t). Thus

\[ `\boxed{IF(t)=C e^{j2\pi f_b t}}`{=tex} \]

with

\[ `\boxed{f_b=S\tau}`{=tex}. \]

The factor (e^{-j`\pi `{=tex}S`\tau`{=tex}^2}) is a constant phase
factor in this stationary-target model. Its magnitude is one, so it
changes phase but not amplitude. It does not change the beat frequency.

### Sign convention

If instead you calculate

\[ s\_{`\mathrm{RX}`{=tex}}(t)s\_{`\mathrm{TX}`{=tex}}\^\*(t), \]

the beat tone has the opposite frequency sign. This is not a physical
disagreement; it is a mixing convention. Be consistent about the sign
and whether you inspect positive or negative FFT frequencies.

For real-valued mixing, the product contains sum- and
difference-frequency terms; low-pass filtering retains the
difference-frequency component. Complex mixing with a conjugate gives
the signed difference-frequency component directly.

### Figure 4 --- FMCW processing chain

``` mermaid
flowchart TD
    A["Transmit LFM ramp"] --> B["Target reflection: delayed chirp"]
    A --> D["Coherent transmit reference"]
    B --> E["Mixer / dechirp: TX times conjugate RX"]
    D --> E
    E --> F["IF signal: beat tone(s)"]
    F --> G["Window and FFT over fast time"]
    G --> H["Detect frequency peaks"]
    H --> I["Convert beat frequency to range"]
```

------------------------------------------------------------------------

## 7. What do we do with the beat frequency?

The beat frequency is converted into range.

Since

\[ f_b=S`\tau`{=tex} `\quad`{=tex}`\text{and}`{=tex}`\quad`{=tex}
`\tau`{=tex}=`\frac{2R}{c}`{=tex}, \]

we get

\[ f_b=`\frac{2SR}{c}`{=tex}. \]

Rearrange:

\[ `\boxed{R=\frac{c f_b}{2S}}`{=tex}. \]

For multiple stationary targets, the IF signal is a sum of beat tones
(under the simplified model):

\[ IF(t)`\approx`{=tex}`\sum`{=tex}*i C_i e\^{j2`\pi `{=tex}f*{b,i}t}.
\]

An FFT over the samples within a chirp reveals the beat-frequency peaks.
Each peak maps to a range.

### Range-processing steps

1.  Collect complex IF samples during one chirp/ramp.
2.  Optionally remove DC or leakage components.
3.  Apply a window such as Hann to reduce FFT sidelobes.
4.  Compute the FFT over fast time.
5.  Inspect the appropriate positive/negative frequency half according
    to the mixing convention.
6.  Detect valid peaks, rejecting noise and sidelobes.
7.  Convert each detected beat frequency to range using (R=cf_b/(2S)).

The FFT is not itself a range measurement. It estimates the frequencies
in the IF signal; the radar slope converts those frequencies into range.

------------------------------------------------------------------------

## 8. Numerical example used in our discussion

Parameters:

  Parameter                                               Value
  ---------------------------- --------------------------------
  Speed of light (c)                 (3`\times10`{=tex}\^8) m/s
  Carrier example (f_0)                                   9 GHz
  Bandwidth (B)                                         100 MHz
  Chirp duration (T)                                      20 µs
  Slope (S=B/T)                  (5`\times10`{=tex}\^{12}) Hz/s
  Sampling frequency (F_s)                              200 MHz
  Samples per chirp (N=F_sT)                               4000

The complex-baseband simulation does not directly sample the 9 GHz
carrier.

For a stationary target:

\[ `\tau`{=tex}=`\frac{2R}{c}`{=tex},`\qquad `{=tex}f_b=S`\tau`{=tex}.
\]

    Target range   Round-trip delay (`\tau`{=tex})   Beat frequency (f_b)
  -------------- --------------------------------- ----------------------
            50 m                         0.3333 µs              1.667 MHz
            80 m                         0.5333 µs              2.667 MHz
           120 m                         0.8000 µs              4.000 MHz

For example, the 50 m target:

\[ R=`\frac{(3\times10^8)(1.667\times10^6)}`{=tex}
{2(5`\times10`{=tex}\^{12})} `\approx 50`{=tex} `\mathrm{m}`{=tex}. \]

The ideal bandwidth-limited range resolution is

\[
`\Delta `{=tex}R`\approx`{=tex}`\frac{3\times10^8}{2(100\times10^6)}`{=tex}
=1.5 `\mathrm{m}`{=tex}. \]

------------------------------------------------------------------------

## 9. Comparing the two methods

  -----------------------------------------------------------------------
  Aspect                  Pulsed LFM pulse        FMCW range processing
                          compression             
  ----------------------- ----------------------- -----------------------
  Transmitted waveform    Finite LFM pulse        Repeated/continuous
                                                  chirp ramps

  Received echo           Time delay              Beat frequency after
  information                                     dechirping

  Main operation          Matched filtering /     Mixing with reference,
                          correlation             then FFT

  Range estimate          (R=c`\tau`{=tex}/2)     (R=cf_b/(2S))

  What the peak           Correlation peak at     IF spectrum peak at
  represents              echo delay              beat frequency

  Sidelobes               Matched-filter response FFT/window spectral
                          sidelobes               sidelobes

  Ideal bandwidth-limited Approximately (c/(2B))  Approximately (c/(2B)),
  resolution                                      under standard FMCW
                                                  assumptions
  -----------------------------------------------------------------------

The two methods use different receiver processing, but both exploit the
delay of a bandwidth-(B) waveform. They are not identical systems: their
timing, transmit duty cycle, receiver architecture, leakage, Doppler
behavior, and practical constraints differ.

### Figure 5 --- The conceptual difference

``` text
PULSED LFM
Received echo ──> matched filter ──> peak position in DELAY ──> range

FMCW
Received echo + TX reference ──> dechirp ──> IF tone
                                      │
                                      v
                                  FFT peak in FREQUENCY
                                      │
                                      v
                                    range
```

------------------------------------------------------------------------

## 10. Important practical issues and common mistakes

### A. Do not confuse carrier frequency with beat frequency

The carrier example is (f_0=9) GHz. The beat frequencies in our example
are only about 1.7--4 MHz. The complex-baseband model avoids sampling
the 9 GHz carrier directly.

### B. The beat frequency depends on the chirp slope

\[ f_b=S`\tau`{=tex}. \]

For the same target range, a larger slope produces a larger beat
frequency. Range conversion must use the actual calibrated slope.

### C. The delayed echo must fit in the receive-time array

For a finite pulse starting at time zero and lasting (T), the echo
occupies (\[`\tau`{=tex},T+`\tau`{=tex})). A receive array that ends at
(T) truncates the echo tail.

### D. Integer sample delays are quantized

Using `round(tau * Fs)` approximates the delay to the nearest sample. At
(F_s=200) MHz, one sample is 5 ns, corresponding to a two-way range
increment of

\[ `\frac{c}{2F_s}`{=tex}=0.75 `\mathrm{m}`{=tex}. \]

This is the delay-grid spacing of that simple integer-delay simulation,
**not** the fundamental waveform range resolution. Fractional-delay
modeling and peak interpolation can provide finer delay estimates;
actual resolvability still depends strongly on bandwidth and SNR.

### E. FFT bin spacing is not exactly the same as range resolution

For (N) samples at sampling rate (F_s), the unpadded FFT bin spacing is

\[ `\Delta `{=tex}f=`\frac{F_s}{N}`{=tex}. \]

The corresponding range-bin spacing for the simple FMCW mapping is

\[ `\Delta `{=tex}R\_{`\mathrm{bin}`{=tex}}
=`\frac{c\Delta f}{2S}`{=tex}. \]

With (N=F_sT), this becomes (c/(2ST)=c/(2B)). Zero-padding makes the
displayed FFT grid denser, but does not create new physical resolving
power.

### F. Windowing changes the peak shape

A Hann window reduces sidelobes but broadens the mainlobe and changes
peak amplitude. Compare methods using consistent normalization and known
windowing.

### G. Moving targets need Doppler treatment

The simple derivation assumes a stationary target, so (`\tau`{=tex}) is
constant during a chirp. A moving target introduces Doppler and
potentially time-varying delay. In practical FMCW radar, range and
Doppler are commonly separated using processing across fast time and
across successive chirps (slow time).

### H. A single chirp and a sequence of chirps answer different questions

-   FFT across samples **within one chirp**: range profile.
-   FFT across repeated chirps at a fixed range bin: Doppler/velocity
    information.
-   Two-dimensional processing: range-Doppler map.

------------------------------------------------------------------------

## 11. Compact Python reference: FMCW IF and range FFT

This is a minimal stationary-target, complex-baseband illustration. It
assumes the observation time array includes the desired echo interval
and that the reference chirp is available for dechirping.

``` python
import numpy as np
import matplotlib.pyplot as plt

c = 3e8
B = 100e6
T = 20e-6
S = B / T
Fs = 200e6

# Include a receive interval long enough for the echoes.
ranges = np.array([50.0, 80.0, 120.0])
amplitudes = np.array([1.0, 0.8, 0.7])
taus = 2 * ranges / c

Nrx = int(np.ceil((T + taus.max()) * Fs)) + 1
t_rx = np.arange(Nrx) / Fs

# Reference chirp (mathematical phase reference).
tx_ref = np.exp(1j * np.pi * S * t_rx**2)

# Sum delayed echoes. Each echo exists only over [tau, T + tau).
rx = np.zeros(Nrx, dtype=complex)
for tau, amp in zip(taus, amplitudes):
    valid = (t_rx >= tau) & (t_rx < T + tau)
    rx[valid] += amp * np.exp(
        1j * np.pi * S * (t_rx[valid] - tau)**2
    )

# This simplified illustration uses the reference phase to dechirp.
# In a practical ramp system, use the correctly timed TX/LO reference
# and process the intended overlap interval.
if_signal = tx_ref * np.conj(rx)

# Select an interval where all three echoes are present.
valid_all = (t_rx >= taus.max()) & (t_rx < T)
if_data = if_signal[valid_all]

# Window and FFT.
N = len(if_data)
window = np.hanning(N)
spectrum = np.fft.fft(if_data * window)
freq = np.fft.fftfreq(N, d=1/Fs)

# This convention produces positive beat frequencies for the model above.
positive = freq >= 0
freq_pos = freq[positive]
mag_pos = np.abs(spectrum[positive])

plt.plot(freq_pos / 1e6, mag_pos)
plt.xlabel("Beat frequency (MHz)")
plt.ylabel("FFT magnitude")
plt.title("FMCW IF spectrum")
plt.grid(True)
plt.show()

# Example: convert a detected beat frequency to range.
fb_detected = 1.667e6
range_est = c * fb_detected / (2 * S)
print(f"Estimated range: {range_est:.2f} m")
```

**Modeling caution:** The code above illustrates the equations, but a
production FMCW simulator should carefully define the transmitted ramp
timing, receiver observation interval, reference/LO phase, chirp
repetition, and Doppler. In particular, do not interpret samples outside
the valid TX/RX overlap as ordinary target beat samples.

------------------------------------------------------------------------

## 12. Final mental model

Remember these three relationships:

1.  **LFM chirp:** quadratic phase creates a linear frequency sweep.

    \[ s(t)=e^{j`\pi `{=tex}St^2}
    `\quad`{=tex}`\Longrightarrow`{=tex}`\quad`{=tex} f(t)=St. \]

2.  **Propagation delay:** the received echo is a delayed copy.

    \[ s\_{`\mathrm{RX}`{=tex}}(t)=A
    e^{j`\pi `{=tex}S(t-`\tau`{=tex})^2},
    `\qquad `{=tex}`\tau`{=tex}=`\frac{2R}{c}`{=tex}. \]

3.  **Range extraction:**

    -   Pulsed LFM: matched-filter peak delay (`\tau`{=tex}), then
        (R=c`\tau`{=tex}/2).
    -   FMCW: dechirp, measure (f_b=S`\tau`{=tex}), then (R=cf_b/(2S)).

The central physical quantity is the round-trip propagation delay. Pulse
compression estimates it from a correlation peak; FMCW estimates it
indirectly from the beat frequency.
