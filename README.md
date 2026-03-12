# Pole Tracking (PT) Toolbox for Time-Varying Oscillatory Signal Analysis

This repository contains a MATLAB implementation of a **pole tracking framework for time-varying autoregressive (TV-AR) models** applied to oscillatory signal analysis.

The toolbox implements a **Kalman-filter based pole tracking algorithm with GA optimization**, allowing estimation of **time-varying spectral peaks** in signals such as EEG.

The approach was introduced in:

**Reference**

Kostoglou, K., & Müller-Putz, G. R. (2024).  
*Motor-Related EEG Analysis Using a Pole Tracking Approach.*  
IEEE Transactions on Neural Systems and Rehabilitation Engineering, 32, 3837-3847.

---

# Structure of the MATLAB Toolbox

```
/ [root]
├── code
│   ├── MAIN.m
│   ├── PT_optimize.m
│   ├── SIM_PT.m
│   ├── PT.m
│   ├── GAPT.m
│   ├── extractpoles.m
│   ├── generate_signals.m
```

---

# Preprocessing Recommendations (Important)

When applying pole tracking to **EEG or neural signals**, proper preprocessing is essential.

## High-pass filtering

If the goal is to track **alpha or beta rhythms**, the signals should first be **high-pass filtered at 5 Hz or higher**.

EEG typically contains strong low-frequency components:

• delta activity  
• slow drifts  
• baseline shifts  
• motion artifacts  

If these components remain in the signal, the pole tracker may allocate poles to them.  
In that case the algorithm may end up tracking **delta activity instead of alpha or beta rhythms**.

Recommended preprocessing example:

```
High-pass filter : 5 Hz
Low-pass filter  : 40 Hz
```

---

## Sampling rate and resampling

Pole tracking generally performs better at **moderate sampling rates** rather than extremely high ones.

At very high sampling rates (e.g. 500–1000 Hz):

• the pole angle corresponding to an oscillation becomes very small  
• the estimation problem becomes poorly conditioned  
• frequency estimates may drift toward lower frequencies  

For rhythms such as alpha and beta, recommended sampling rates are:

```
64 Hz
80 Hz
100 Hz
128 Hz
200 Hz
```

Typical preprocessing pipeline:

```
1. High-pass filter (≈5 Hz)
2. Optional low-pass filter
3. Resample / downsample the signal
4. Run pole tracking
```

Example:

```
Original EEG sampling rate : 1000 Hz
High-pass filter           : 5 Hz
Low-pass filter            : 40 Hz
Resample to                : 64 Hz
```

Similar preprocessing strategies were used in simulation and EEG analyses where signals were band-limited and resampled before tracking oscillatory dynamics. :contentReference[oaicite:1]{index=1}

---

# Code Description

## MAIN.m

Example script demonstrating the full pole tracking pipeline.

Steps performed:

1. Generate synthetic signals
2. Select training data
3. Optimize model parameters using a Genetic Algorithm
4. Run recursive pole tracking
5. Extract pole features
6. Plot estimated trajectories

---

## PT_optimize.m

Optimizes pole tracking parameters using a **Genetic Algorithm (GA)**.

The optimizer searches for parameters minimizing the **normalized prediction error** of the pole tracking model.

Inputs

```
ytr   : training signals [M x N x trials]
pmax  : number of oscillatory pole pairs
pmax2 : number of real poles
Fs    : sampling frequency
```

Outputs

```
lam    : optimized parameter vector
err_gm : final optimization cost
```

Parameter vector structure

```
lam =
[R2
 R1
 oscillatory_radii
 real_poles
 oscillatory_angles
 P0]
```

---

## SIM_PT.m

Recursive pole tracking algorithm.

The algorithm tracks pole parameters sample-by-sample using a Kalman filter update scheme.

Outputs include:

• pole radii  
• pole angles  
• prediction error  
• covariance norm  

The instantaneous frequency is obtained from the pole angle:

```
f = abs(theta) * Fs / (2*pi)
```

---

## PT.m

Cost function used during GA optimization.

Evaluates the normalized prediction error between the signal and the model output.

Output

```
J = normalized prediction error
```

---

## GAPT.m

Wrapper that applies the PT cost function across all signals and trials.

Used internally during GA optimization.

---

## extractpoles.m

Runs pole tracking on full signals and extracts pole features.

Feature layout:

```
[oscillatory radii
 oscillatory frequencies
 real poles
 covariance norm]
```

---

## generate_signals.m

Utility function generating synthetic signals with time-varying oscillations.

Possible dynamics include:

• linear frequency changes  
• sinusoidal modulation  
• abrupt transitions  
• mixed trajectories  

Outputs include the ground-truth frequency and amplitude trajectories for validation.

---

# Relation to the EKF Oscillator Tracking Repository

This repository is closely related to the **EKF Oscillator Tracking** repository.

However, the two methods differ conceptually.

## Pole Tracking (this repository)

• Tracks **multiple poles simultaneously**  
• Uses a **time-varying autoregressive cascade representation**  
• Focuses on **tracking spectral peak locations**

Important interpretation:

The **pole magnitude does NOT directly represent the true signal amplitude**.

Instead, the pole radius reflects **how strongly the tracked oscillatory frequency is represented in the signal at that time point**.

Thus, the magnitude of the pole mainly indicates **spectral prominence**, not instantaneous amplitude.

The framework can track **multiple oscillatory components simultaneously** by increasing the number of poles.

---

## EKF Oscillator Tracking

The EKF oscillator framework models the signal as a **time-varying damped harmonic oscillator** with explicit magnitude and frequency states. :contentReference[oaicite:2]{index=2}

In that formulation:

• the state vector contains oscillator states and parameters  
• the **instantaneous magnitude is directly estimated** from the oscillator states  

Magnitude is computed as

```
A(n) = sqrt(x1(n)^2 + x2(n)^2)
```

meaning that the EKF explicitly tracks **true amplitude fluctuations** in addition to frequency. :contentReference[oaicite:3]{index=3}

Key characteristics:

• tracks **one oscillator at a time**  
• uses a **nonlinear state-space model**  
• requires bandpass filtering to isolate the target rhythm  

---

## Key differences

| Feature | Pole Tracking | EKF Oscillator Tracking |
|-------|------|------|
Tracks multiple oscillations | ✓ | ✗ |
Tracks instantaneous amplitude | ✗ | ✓ |
Tracks spectral peak location | ✓ | ✓ |
Model type | TV-AR pole representation | nonlinear oscillator state-space |
Number of oscillators | multiple | single |

In practice:

• **Pole tracking** is useful for detecting and tracking **multiple spectral peaks simultaneously**.

• **EKF oscillator tracking** is better suited when the goal is to estimate **instantaneous amplitude and frequency of a specific rhythm**.

---

# Summary

This toolbox provides a framework for

• optimizing pole tracking models  
• tracking time-varying oscillatory poles  
• extracting frequency trajectories from signals  
• testing algorithms on synthetic or real data such as EEG  

For neural data, proper preprocessing is essential.  
High-pass filtering (≈5 Hz) and resampling to moderate sampling rates (e.g., 64–200 Hz) are recommended to ensure the tracker captures relevant oscillatory activity rather than slow drifts.
