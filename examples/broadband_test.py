import numpy as np
import matplotlib.pyplot as plt
import json
import os
import python_utils

class BroadbandTest:
    
    def __init__(self, fs, f0, amp, source_len, t):
        self.fs = fs
        self.f0 = f0
        self.amp = amp
        self.source_len = source_len
        self.t = t
        self.n = len(t)

    def hanning_window_source(self):

        # duration of the source in seconds: number of cycles / frequency
        T = self.source_len / self.f0
        window = np.zeros_like(self.t)
        mask = (self.t >= 0) & (self.t <= T)
        # standard Hann/Hanning window on [0, T]
        window[mask] = 0.5 * (1.0 - np.cos(2.0 * np.pi * self.t[mask] / T))

        pressure_time_series = self.amp * np.sin(2.0 * np.pi * self.f0 * self.t) * window

        # zero outside the window (mask already enforces this, but keep for clarity)
        pressure_time_series[~mask] = 0.0

        self.pressure_time_series = pressure_time_series

        return pressure_time_series
    
    def source_pulse_spectrum(self):
        # For a real-valued time series use rfft to get the one-sided spectrum
        source_spectrum = np.fft.rfft(self.pressure_time_series)
        frequencies = np.fft.rfftfreq(self.n, 1.0 / self.fs)

        return source_spectrum, frequencies

fs = 400        # Sampling frequency (Hz)
t = np.arange(0.0, 1.0, 1.0/fs)      # Time vector from 0 to 1 second
f0 = 50.0      # Center frequency of the Gaussian pulse (Hz)
amp = 1.0      # Amplitude of the pulse
source_len = 4

bbt = BroadbandTest(fs, f0, amp, source_len, t)

source_time_series = bbt.hanning_window_source()
source_spectrum, frequencies = bbt.source_pulse_spectrum()

plt.figure()
plt.subplot(2,1,1)
plt.plot(t, np.real(source_time_series), label='Real')
plt.plot(t, np.imag(source_time_series), label='Imaginary')
plt.plot(t, np.abs(source_time_series), label='Magnitude', color='black', lw=1)
plt.legend()
plt.title('Hanning windowed source time series')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.grid()

plt.subplot(2,1,2)
plt.plot(frequencies, np.real(source_spectrum), label='Real')
plt.plot(frequencies, np.imag(source_spectrum), label='Imaginary')
plt.plot(frequencies, np.abs(source_spectrum), label='Magnitude', color='black', lw=1)
plt.legend()
plt.title('Source spectrum (one-sided rfft)')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude')
# plt.xlim(0, fs/2)
plt.grid()
plt.tight_layout()


z = np.linspace(0.0, 5000.0, 50)
munk_ssp = python_utils.munk(z)
munk_ssp_3d = np.tile(munk_ssp, (2, 2, 1))
munk_ssp_3d_flat = munk_ssp_3d.flatten(order='C')
jsonfile = "examples/testbb.json"
h5file = "examples/testbb.h5"


env_bb = {
    "ssp": {
        "x_ssp_m": [0.0, 50000.0],
        "y_ssp_m": [0.0, 50000.0],
        "z_ssp_m": list(z),
        "c_m_s": list(munk_ssp_3d_flat)
    },
    "bathymetry": {
        "x_bty_m": [0.0, 50000.0],
        "y_bty_m": [0.0, 50000.0],
        "z_bty_m": np.array([[5000.0, 5000.0], [5000.0, 5000.0]]).flatten(order='C').tolist(),
        "density_g_cm3": 1.5,
        "c_bty_m_s": 1600.0,
        "attenuation_db": 0.0
    },
    "source": {
        "position": [0.0, 0.0, 1000.0],
        "freq_hz": frequencies.tolist(),
        "launch_elev_deg": np.linspace(-10.0, 10.0, 1000).tolist(),
        "launch_azim_deg": np.linspace(-0.05, 0.05, 3).tolist()
    },
    "receivers": {
        "x_rcvr_m": [0.0],
        "y_rcvr_m": [50000.0],
        "z_rcvr_m": [1000.0]
    },
    "beam": {
        "step_m": 15.0,
        "max_steps": 60_000,
        "max_range_m": 60_000.0
    }
}

with open(jsonfile, "w") as f:
    json.dump(env_bb, f, indent=2)

# run rtrs with the JSON-like input (keeps original behavior of the example)
os.system(f"./target/release/rtrs  {jsonfile}")

freq, x_m, y_m, z_m, pressure = python_utils.load_cmpx_pressure(h5file)
pressure = np.reshape(pressure, (len(freq), len(x_m), len(y_m), len(z_m)))

# Multiply frequency-domain pressure by source spectrum. Ensure frequency axes align.
# source_spectrum is one-sided (rfft) with frequencies = np.fft.rfftfreq(n, 1/fs)
if len(freq) != len(source_spectrum):
    # interpolate complex spectrum onto the frequency grid from the file
    src_freq = frequencies
    # interpolate real and imag parts separately
    real_interp = np.interp(freq, src_freq, source_spectrum.real)
    imag_interp = np.interp(freq, src_freq, source_spectrum.imag)
    source_spectrum_used = real_interp + 1j * imag_interp
else:
    source_spectrum_used = source_spectrum

# vectorized multiply: expand dims to broadcast across (freq, x, y, z)
pressure_scaled = pressure * source_spectrum_used[:, None, None, None]

# inverse transform back to time domain using irfft. Specify n to get original time length.
pressure_time = np.fft.irfft(pressure_scaled, axis=0, n=bbt.n)

plt.figure()
plt.subplot(2,1,1)
# choose last x,y and middle z index similar to original
plt.plot(t, np.real(pressure_time[:, -1, -1, len(z_m)//2]), label='Real')
plt.legend()
plt.title('Pressure time series at receiver (30 km, 20 m depth)')
plt.xlabel('Time (s)')
plt.ylabel('Pressure')
plt.grid() 
plt.subplot(2,1,2)
# plot the frequency-domain pressure at the same receiver
# use freq (from file) as x axis
plt.plot(freq, np.real(pressure[:, -1, -1, len(z_m)//2]), label='Real')
plt.legend()
plt.title('Pressure spectrum at receiver (30 km, 20 m depth)')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Pressure')
# plt.xlim(0, fs/2)
plt.grid()
plt.tight_layout()

# show plots
plt.show()
