import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ====================== PARAMETERS (normalized units) ======================
c = 1.0          # substrate propagation limit
d = 1.0          # wire diameter

L = 10.0
Nx = 900
Nr = 100
dx = L / (Nx - 1)
dr = d / (2 * (Nr - 1))
dt = 0.45 * dx / c

T_sim = 45.0
Nt = int(T_sim / dt)

# DC inlet: the source ramps up once and then stays on.
source_on = 2.0
source_ramp = 1.0
source_level = 1.0
source_gain = 1.6
source_cells = 6

# Wire bundle closure: many microchannels can locally jam while the aggregate stays smooth.
n_channels = 48
transverse_mix = 0.3
channel_seed = 7

# Transport closure: demand can exceed throughput, but actual flow is capped.
mobility = 2.8
drive_speed = 0.72 * c
backflow_speed = 0.75 * c
base_capacity = 0.56 * c
sat_threshold = 0.18
sat_steepness = 10.0
gate_tau = 0.18
overflow_gain = 3.0
sat_decay = 2.6
sat_pressure_gain = 0.4

print("DC inlet enabled: no temporal driver modulation.")
print(f"Probe FFT bin width: {1 / T_sim:.3f}")

# ====================== FIELD SETUP ======================
rng = np.random.default_rng(channel_seed)

r = np.linspace(0, d / 2, Nr)
tube_radius = 0.46 * d
tube_profile = np.zeros(Nr)
inside = r < tube_radius
tube_profile[inside] = 1.0 - (r[inside] / tube_radius) ** 4
tube_profile = np.clip(tube_profile, 0.0, None)

capacity_noise = rng.normal(size=(n_channels, Nx - 1))
threshold_noise = rng.normal(size=(n_channels, Nx - 1))
for _ in range(8):
    capacity_noise[:, 1:-1] = (
        0.2 * capacity_noise[:, :-2]
        + 0.6 * capacity_noise[:, 1:-1]
        + 0.2 * capacity_noise[:, 2:]
    )
    threshold_noise[:, 1:-1] = (
        0.2 * threshold_noise[:, :-2]
        + 0.6 * threshold_noise[:, 1:-1]
        + 0.2 * threshold_noise[:, 2:]
    )

channel_capacity = base_capacity * np.clip(
    1.0 + 0.28 * capacity_noise,
    0.4,
    1.7,
)
channel_threshold = sat_threshold * np.clip(
    1.0 + 0.35 * threshold_noise,
    0.35,
    1.9,
)
channel_drive = drive_speed * np.clip(
    1.0 + 0.08 * rng.normal(size=(n_channels, 1)),
    0.85,
    1.15,
)
channel_source = source_level * np.clip(
    1.0 + 0.04 * rng.normal(size=(n_channels, 1)),
    0.9,
    1.1,
)

signal = np.zeros((n_channels, Nx))        # deficit magnitude per pathway bundle
sat_face = np.zeros((n_channels, Nx - 1))  # local saturation backlog on each axial face
gate = np.ones((n_channels, Nx - 1))       # effective openness of each axial face

probe_z = Nx // 3
micro_probe = n_channels // 2
probe_series = []
micro_probe_series = []
spread_series = []
sat_series = []
history = []


def smooth_step(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3.0 - 2.0 * x)


def inlet_level(t):
    return source_level * smooth_step((t - source_on) / source_ramp)


def node_saturation(face_values):
    node_values = np.zeros((face_values.shape[0], face_values.shape[1] + 1))
    node_values[:, :-1] += 0.5 * face_values
    node_values[:, 1:] += 0.5 * face_values
    return node_values


def upwind_flux(speed, quantity):
    left = quantity[:, :-1]
    right = quantity[:, 1:]
    return np.where(speed >= 0.0, speed * left, speed * right)


def mix_channels(values, mix_rate):
    lap = np.zeros_like(values)
    lap[1:-1, :] = values[:-2, :] - 2.0 * values[1:-1, :] + values[2:, :]
    lap[0, :] = values[1, :] - values[0, :]
    lap[-1, :] = values[-2, :] - values[-1, :]
    return values + dt * mix_rate * lap


def step(signal, sat_face, gate, t):
    node_sat = node_saturation(sat_face)
    pressure = signal + sat_pressure_gain * node_sat

    pressure_drop = pressure[:, :-1] - pressure[:, 1:]
    dc_speed = channel_drive * smooth_step((t - source_on) / source_ramp)
    gradient_speed = backflow_speed * np.tanh(mobility * pressure_drop / max(source_level, 1e-6))
    demand_speed = dc_speed + gradient_speed

    gate_arg = np.clip(sat_steepness * (sat_face - channel_threshold), -60.0, 60.0)
    gate_target = 1.0 / (1.0 + np.exp(gate_arg))
    gate_next = gate + dt * (gate_target - gate) / gate_tau
    gate_next = np.clip(gate_next, 0.02, 1.0)

    speed_cap = channel_capacity * gate_next
    actual_speed = np.clip(demand_speed, -speed_cap, speed_cap)

    actual_flux = upwind_flux(actual_speed, signal)
    boundary_flux = np.zeros((n_channels, Nx + 1))
    boundary_flux[:, 1:-1] = actual_flux
    boundary_flux[:, -1] = channel_capacity[:, -1] * gate_next[:, -1] * signal[:, -1]

    signal_next = signal - (dt / dx) * (boundary_flux[:, 1:] - boundary_flux[:, :-1])
    inlet_target = channel_source * smooth_step((t - source_on) / source_ramp)
    signal_next[:, :source_cells] += dt * source_gain * (inlet_target - signal_next[:, :source_cells])
    signal_next = mix_channels(signal_next, transverse_mix)
    signal_next = np.clip(signal_next, 0.0, None)

    local_load = 0.5 * (signal[:, :-1] + signal[:, 1:])
    overflow = np.maximum(0.0, np.abs(demand_speed) - speed_cap) * local_load
    sat_next = sat_face + dt * (overflow_gain * overflow - sat_decay * sat_face)
    sat_next = np.clip(sat_next, 0.0, 4.0)

    return signal_next, sat_next, gate_next


def field_from_signal(signal):
    aggregate_signal = np.mean(signal, axis=0)
    return -tube_profile[:, None] * aggregate_signal[None, :]


def dominant_frequency(series):
    ac = series[fft_start:] - np.mean(series[fft_start:])
    freqs = np.fft.rfftfreq(ac.size, d=dt)
    fft = np.abs(np.fft.rfft(ac))

    if fft.size > 1 and np.any(fft[1:] > 0.0):
        peak_idx = np.argmax(fft[1:]) + 1
        return freqs[peak_idx]
    return None


# ====================== RUN SIMULATION ======================
history_stride = max(1, Nt // 220)
arrival_time = None
arrival_threshold = 0.12 * source_level

for n in range(Nt):
    t = n * dt
    signal, sat_face, gate = step(signal, sat_face, gate, t)

    node_sat = node_saturation(sat_face)
    aggregate_probe = np.mean(signal[:, probe_z])
    probe_series.append(aggregate_probe)
    micro_probe_series.append(signal[micro_probe, probe_z])
    spread_series.append(np.std(signal[:, probe_z]))
    sat_series.append(np.mean(node_sat[:, probe_z]))

    if arrival_time is None and aggregate_probe > arrival_threshold:
        arrival_time = t

    if n % history_stride == 0:
        history.append(field_from_signal(signal))

print("Simulation finished.")
if arrival_time is None:
    print("Probe did not receive the DC front inside the simulated time window.")
else:
    print(f"DC front arrival at probe: t = {arrival_time:.3f}")

# ====================== EMERGENT MODULATION DIAGNOSTIC ======================
probe_series = np.array(probe_series)
micro_probe_series = np.array(micro_probe_series)
spread_series = np.array(spread_series)
sat_series = np.array(sat_series)

fft_start = Nt // 3
aggregate_f = dominant_frequency(probe_series)
micro_f = dominant_frequency(micro_probe_series)

if aggregate_f is not None:
    print(f"Dominant aggregate probe modulation: {aggregate_f:.3f}")
else:
    print("No resolved aggregate modulation above DC at the probe.")

if micro_f is not None:
    print(f"Representative microchannel modulation: {micro_f:.3f}")
else:
    print("No resolved microchannel modulation above DC at the probe.")

print(f"Peak mean saturation backlog at probe: {np.max(sat_series):.3f}")
print(f"Peak cross-channel spread at probe: {np.max(spread_series):.3f}")
print(f"Minimum gate opening across channels: {np.min(gate):.3f}")

# ====================== ANIMATION ======================
fig, ax = plt.subplots(figsize=(9, 5))
im = ax.imshow(
    history[0],
    cmap="RdBu_r",
    extent=[0, L, 0, d / 2],
    aspect="auto",
    vmin=-1.35 * source_level,
    vmax=0.25 * source_level,
)
ax.set_xlabel("z along wire (normalized)")
ax.set_ylabel("r (radial)")
ax.set_title("Aptik Wire - bundled pathways with smooth aggregate flow")
plt.colorbar(im, label="Scalar field A")


def animate(frame):
    im.set_array(history[frame])
    return [im]


ani = FuncAnimation(fig, animate, frames=len(history), interval=40, blit=True)
plt.show()