import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# If Spyder stays inline/frozen, uncomment the next two lines and rerun.
# import matplotlib
# matplotlib.use("QtAgg")

NX, NY = 240, 240
DX = 1.0 / NX
DY = 1.0 / NY

DT = 0.006
SUBSTEPS_PER_FRAME = 4
C = 0.32
ATTRACTION = 0.95
AMBIENT = 1.0
VOID_GEOMETRY_GAIN = 1.2
VOID_SHELL_GAIN = 0.9
VOID_SHELL_BAND = 0.06
NUMERICAL_VISCOSITY = 0.0
NOISE = 0.0005

x = np.linspace(0.0, 1.0, NX)
y = np.linspace(0.0, 1.0, NY)
X, Y = np.meshgrid(x, y)


def ddx(a):
    ap = np.pad(a, ((0, 0), (1, 1)), mode="edge")
    return (ap[:, 2:] - ap[:, :-2]) / (2.0 * DX)


def ddy(a):
    ap = np.pad(a, ((1, 1), (0, 0)), mode="edge")
    return (ap[2:, :] - ap[:-2, :]) / (2.0 * DY)


def laplacian(a):
    ap = np.pad(a, ((1, 1), (1, 1)), mode="edge")
    return (
        (ap[1:-1, 2:] - 2.0 * a + ap[1:-1, :-2]) / (DX * DX)
        + (ap[2:, 1:-1] - 2.0 * a + ap[:-2, 1:-1]) / (DY * DY)
    )


VOID_X = 0.52
VOID_Y = 0.50
VOID_R = 0.055

rx = X - VOID_X
ry = Y - VOID_Y
dist = np.sqrt(rx * rx + ry * ry)
void_mask = dist <= VOID_R


def void_boundary_response():
    safe_dist = np.maximum(dist, VOID_R + 1.0e-6)
    ratio = np.clip(VOID_R / safe_dist, 0.0, 1.0)
    missing_fraction = np.arcsin(ratio) / np.pi
    shell_distance = (safe_dist - VOID_R) / VOID_SHELL_BAND
    shell_kernel = np.exp(-(shell_distance * shell_distance))

    outward_norm = np.maximum(dist, 1.0e-6)
    outward_x = rx / outward_norm
    outward_y = ry / outward_norm

    tension_boost = VOID_SHELL_GAIN * missing_fraction * shell_kernel
    fx = VOID_GEOMETRY_GAIN * missing_fraction * shell_kernel * outward_x
    fy = VOID_GEOMETRY_GAIN * missing_fraction * shell_kernel * outward_y

    tension_boost[void_mask] = 0.0
    fx[void_mask] = 0.0
    fy[void_mask] = 0.0
    return tension_boost, fx, fy


void_tension_boost, void_fx, void_fy = void_boundary_response()

rho = AMBIENT * np.ones((NY, NX), dtype=np.float64)
vx = np.zeros_like(rho)
vy = np.zeros_like(rho)

rng = np.random.default_rng(7)
rho += NOISE * rng.standard_normal(rho.shape)
rho = np.clip(rho, 0.0, None)


def enforce_constraints():
    global rho, vx, vy
    rho[void_mask] = 0.0
    vx[void_mask] = 0.0
    vy[void_mask] = 0.0

    # Hold the far field at ambient to represent the remote substrate tugging back.
    rho[:, 0] = AMBIENT
    rho[:, -1] = AMBIENT
    rho[0, :] = AMBIENT
    rho[-1, :] = AMBIENT

    vx[:, 0] = 0.0
    vx[:, -1] = 0.0
    vy[0, :] = 0.0
    vy[-1, :] = 0.0
    rho = np.clip(rho, 0.0, None)


enforce_constraints()
last_speed = np.zeros_like(rho)


def step():
    global rho, vx, vy, last_speed

    gx = ddx(rho)
    gy = ddy(rho)

    vx += (ATTRACTION * gx + void_fx * np.maximum(rho, 0.0)) * DT
    vy += (ATTRACTION * gy + void_fy * np.maximum(rho, 0.0)) * DT

    speed = np.sqrt(vx * vx + vy * vy)
    over = speed > C
    if np.any(over):
        scale = C / speed[over]
        vx[over] *= scale
        vy[over] *= scale
        speed[over] = C

    flux_x = rho * vx
    flux_y = rho * vy
    div_flux = ddx(flux_x) + ddy(flux_y)

    target_rho = AMBIENT + void_tension_boost
    rho += -DT * div_flux + 0.08 * (target_rho - rho) * DT

    if NUMERICAL_VISCOSITY > 0.0:
        rho += NUMERICAL_VISCOSITY * DT * laplacian(rho)

    enforce_constraints()
    last_speed = speed


fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax_rho, ax_speed = axes

im_rho = ax_rho.imshow(
    rho,
    origin="lower",
    extent=(0, 1, 0, 1),
    cmap="inferno",
    vmin=0.0,
    vmax=1.5,
    interpolation="nearest",
)
ax_rho.set_title("Void Boundary Tension Response")
ax_rho.set_xlabel("x")
ax_rho.set_ylabel("y")
fig.colorbar(im_rho, ax=ax_rho, fraction=0.046, pad=0.04)

im_speed = ax_speed.imshow(
    last_speed,
    origin="lower",
    extent=(0, 1, 0, 1),
    cmap="viridis",
    vmin=0.0,
    vmax=0.05,
    interpolation="nearest",
)
ax_speed.set_title("Local Flow Speed")
ax_speed.set_xlabel("x")
ax_speed.set_ylabel("y")
fig.colorbar(im_speed, ax=ax_speed, fraction=0.046, pad=0.04)

frame_counter = {"n": 0}


def update(_frame):
    for _ in range(SUBSTEPS_PER_FRAME):
        step()

    frame_counter["n"] += 1
    im_rho.set_data(rho)
    im_speed.set_data(last_speed)
    im_speed.set_clim(0.0, max(0.01, float(np.percentile(last_speed, 99.5))))
    ax_rho.set_title(f"Void Boundary Tension Response | frame={frame_counter['n']}")
    return [im_rho, im_speed]


ani = FuncAnimation(fig, update, interval=30, blit=False, cache_frame_data=False)
plt.tight_layout()
plt.show()