"""Two separately built GRASS modules in one simulation.

`dirt` (DIRT's Python module) runs a DEM bed of 200 glass spheres. `grass`
runs a tiny Python "fluid code": an upward flow through a vertical column
whose interstitial velocity rises where particles crowd it. A grass parent
App couples them every 50 DEM steps:

    DEM particles --(grass.Map, conservative)--> solid fraction per fluid cell
    fluid velocity --(grass.Map, interpolating)--> velocity at each particle

and a Python system *inside* the DIRT App applies a hindered-Stokes drag
force in DIRT's PostForce phase at every DEM step.

dirt.so and grass.so are different shared libraries, each with its own
compiled copy of GRASS; the DIRT App joins the grass parent through the
stable C ABI (grass_pyo3::abi), and data crosses as Python objects.

Run: PYTHONPATH=<dir with grass.abi3.so and dirt.abi3.so> python fluidized_bed.py
"""

import time
from dataclasses import dataclass, field

import numpy as np

import dirt
import grass

BOX = (0.04, 0.04, 0.08)
NZ = 16
MU = 0.05           # Pa s, a viscous liquid (50x water)
RHO_F = 1000.0      # kg/m3
N_RZ = 4.65         # Richardson-Zaki exponent for the hindrance factor
DEM_DT = 1e-5
EXCHANGE_EVERY = 50  # DEM steps per coupling window


# ---- resources -------------------------------------------------------------

@dataclass
class ParticleState:        # published by the DEM App for the parent
    x: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    radius: np.ndarray = field(default_factory=lambda: np.zeros(0))


@dataclass
class FluidAtParticles:     # written by the parent, read by the DEM drag system
    u: np.ndarray = field(default_factory=lambda: np.zeros(0))
    porosity: np.ndarray = field(default_factory=lambda: np.zeros(0))


@dataclass
class Column:               # the fluid code's state
    U: float                # superficial (inlet) velocity, m/s
    phi: np.ndarray         # solid fraction per cell
    u: np.ndarray           # interstitial velocity per cell


# ---- the DEM code (dirt module) ----------------------------------------------

def dem_app(steps: int) -> dirt.App:
    cfg = open("/home/claude/dirt/examples/hello_bed/config.toml").read()
    cfg = cfg.replace("steps = 200000", f"steps = {steps}").replace("thermo = 2000", "thermo = 100000000")
    cfg = cfg.replace("interval = 2000", "interval = 100000000")
    app = dirt.App(config=cfg, output_dir="/tmp/claude-0/pyg/bedout")
    for p in (dirt.plugins.Core(), dirt.plugins.GranularDefault(),
              dirt.plugins.Gravity(), dirt.plugins.Wall()):
        app.add_plugin(p)
    app.insert(ParticleState())
    app.insert(FluidAtParticles())

    def drag(p: dirt.Particles, fluid: FluidAtParticles) -> None:
        if len(fluid.u) != len(p.radius):
            return  # before the first exchange
        eps = np.clip(fluid.porosity, 0.3, 1.0)
        # hindered Stokes drag (Richardson-Zaki form) plus buoyancy
        k = 6 * np.pi * MU * p.radius * eps ** (-N_RZ + 1)
        p.f_ext[:, 2] = k * (fluid.u - p.v[:, 2]) + RHO_F * 9.81 * p.mass / 2500.0
        p.f_ext[:, :2] = -k[:, None] * p.v[:, :2]

    def publish(p: dirt.Particles, out: ParticleState) -> None:
        out.x, out.radius = p.x, p.radius

    app.add_system(drag, dirt.phases.PostForce)  # inside DIRT's step
    app.add_system(publish)                      # after DIRT's step
    return app


# ---- the fluid code (grass module) -----------------------------------------

def fluid_app(U: float) -> grass.App:
    app = grass.App()
    app.insert(Column(U=U, phi=np.zeros(NZ), u=np.full(NZ, U)))

    def continuity(c: Column) -> None:
        c.u = c.U / np.clip(1.0 - c.phi, 0.3, 1.0)

    app.add_system(continuity)
    return app


# ---- the parent (grass module) ---------------------------------------------

class Dem(grass.AppView):
    pass


class Fluid(grass.AppView):
    pass


grid = grass.UniformGrid((0, 0, 0), BOX, (1, 1, NZ), basis="trilinear")
cell_volume = BOX[0] * BOX[1] * BOX[2] / NZ


@dataclass
class BedHistory:
    t: list = field(default_factory=list)
    height: list = field(default_factory=list)


def exchange(sim: grass.Sim, dem: Dem, fluid: Fluid, hist: BedHistory) -> None:
    sim.tick(Dem, EXCHANGE_EVERY)
    ps = dem[ParticleState]
    cloud = grass.PointCloud(ps.x)
    volume = 4 / 3 * np.pi * ps.radius**3
    phi = grass.Map(cloud, grid, kind="extensive")(volume) / cell_volume
    col = fluid[Column]
    col.phi = phi
    sim.tick(Fluid)
    to_particles = grass.Map(grid, cloud, kind="intensive")
    at = dem[FluidAtParticles]
    at.u = to_particles(col.u)
    at.porosity = 1.0 - to_particles(col.phi)
    hist.t.append(sim.step * EXCHANGE_EVERY * DEM_DT)
    hist.height.append(float(np.percentile(ps.x[:, 2], 90)))


def run(U: float, seconds: float = 1.6):
    steps = int(seconds / DEM_DT)
    world = grass.App(config=f"[run]\nsteps = {steps // EXCHANGE_EVERY}\n")
    world.add_plugin(grass.plugins.Run())
    world.add_app(Dem, dem_app(steps))   # a dirt.App, through the C ABI
    world.add_app(Fluid, fluid_app(U))   # a grass.App
    world.insert(BedHistory())
    world.add_system(exchange)
    t0 = time.perf_counter()
    world.run()
    return world[BedHistory], time.perf_counter() - t0


if __name__ == "__main__":
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.2))
    for U in (0.0, 0.03, 0.06, 0.09):
        h, wall = run(U)
        print(f"U = {U*1000:4.1f} mm/s: bed height (90th pct) {h.height[-1]*1000:5.1f} mm "
              f"after {h.t[-1]:.2f} s  [{wall:.1f} s wall]")
        ax.plot(h.t, np.array(h.height) * 1000, label=f"U = {U*1000:.0f} mm/s")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("bed height, 90th percentile [mm]")
    ax.set_title("DIRT bed (dirt module) + Python fluid code (grass module)")
    ax.legend()
    fig.savefig("fluidized_bed.png", dpi=120, bbox_inches="tight")
    print("wrote fluidized_bed.png")
