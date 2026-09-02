# tensorSTHE-ThermalShell

OpenFOAM conjugate-heat-transfer case for a shell-and-tube heat exchanger (STHE). Tube-side and shell-side fluids are meshed as two regions. The metal wall is **not** a third mesh: it is a thin thermal shell on the fluid–fluid interface.

Solver: `chtMultiRegionFoam` (OpenFOAM v2412 / v2406 dictionaries).

| Stream | Region | Path |
|--------|--------|------|
| Tube side (hot) | `hot_fluid` | inlet at \(x = 0\), outlet at \(x = 2.4\,\mathrm{m}\) |
| Shell side (cold) | `cold_fluid` | inlet nozzle at \(z = -0.3\,\mathrm{m}\), outlet at \(z = +0.3\,\mathrm{m}\) |

Default duty: both streams \(2\,\mathrm{kg/s}\) of liquid water, \(T_\mathrm{hot} = 353\,\mathrm{K}\), \(T_\mathrm{cold} = 300\,\mathrm{K}\).

## How the case works

### Two fluids, one thin wall

`solid.stl` is the fluid–fluid enclosure, snapped as a **baffle**. After `topoSet` / `subsetMesh`, coincident baffle faces are merged back to internal faces (`mergeOrSplitBaffles`), then `splitMeshRegions` splits the mesh into `hot_fluid` and `cold_fluid`. Those coincident faces (tube walls) become a pair of `mappedWall` patches:

- `hot_fluid_to_cold_fluid`
- `cold_fluid_to_hot_fluid`

Heat crosses the interface with `compressible::turbulentTemperatureRadCoupledMixed` using turbulent `kappaEff` on both fluids (no extra `thicknessLayers` by default — strongest coupling). Optional metal resistance can be added in `0.orig/include/thermalShellCoupled`.

Remaining single-sided faces from `solid.stl` (shell outer wall, heads) stay as the `solid` wall patch and are adiabatic. No prism / inflation layers.

```
  inlet_hot (x=0)                          outlet_hot (x=2.4 m)
       |                                          |
       v                                          ^
  +----|---------------- tube bundle -------------|----+
  |    |              hot_fluid                   |    |
  |    +-------------- thermal shell -------------+    |
  |                   cold_fluid                       |
  |         ^                               |          |
  +---------|-------------------------------|----------+
      inlet_cold (z=-0.3)              outlet_cold (z=+0.3)
```

### Mesh pipeline (`./buildMesh`)

1. `blockMesh` — background hex box that fully contains the STLs.
2. `surfaceFeatureExtract` — feature edges from the STLs in `constant/triSurface/`.
3. `snappyHexMesh -overwrite` — snap to `solid.stl` as baffle `faceZone` (no prism layers).
4. `topoSet` — walk from a seed in each fluid; baffles stop the walk so the split follows `solid.stl`.
5. `subsetMesh keepCells` — drop exterior cells; expose leftover faces on `solid`.
6. `mergeOrSplitBaffles -overwrite` — turn coincident `solid` / `solid_slave` pairs into internal faces.
7. `splitMeshRegions -cellZonesOnly -overwrite` — create the two region meshes and `mappedWall` coupling.
8. `restore0Dir` — copy `0.orig/` → `0/` (fields are already set per region; there is no `changeDictionary`).

Snappy seed points must sit **inside** each fluid and **off cell faces**:

- hot (tube header): `(0.11 0.01 0.01)`
- cold (shell near inlet nozzle): `(1.80 0.01 -0.22)`

The same points are in `system/snappyHexMeshDict` (`locationsInMesh`) and `system/topoSetDict` (`insidePoints`).

### Physics

- Incompressible liquid (`heRhoThermo` + `rhoConst`) with constant \(\mu\), \(C_p\), Pr. Both regions default to water at ~300 K.
- RAS `kEpsilon`, radiation off, gravity `(0 -9.81 0)`.
- Transient PIMPLE (`nOuterCorrectors 8`, `nNonOrthogonalCorrectors 2` so mapped \(T\) is iterated each step), adjustable \(\Delta t\) with `maxCo 0.5`. `endTime` is in `system/controlDict`.
- Outlet area-average \(T\) and interface `wallHeatFlux` are logged from `controlDict` function objects.

## Directory map

```
constant/triSurface/           STLs in metres — what snappy actually reads
                               (solid.stl, inlet/outlet STLs)
constant/<region>/             thermo, turbulence, radiation, region polyMesh
0.orig/<region>/               initial/boundary fields (edit these)
0.orig/include/                shared thermal-shell BC
system/include/caseSettings    mass flow, T, p, inlet turbulence
system/blockMeshDict
system/snappyHexMeshDict
system/topoSetDict
system/controlDict
system/<region>/               fvSchemes, fvSolution, decomposeParDict
buildMesh  Run  Continue  reconstruct  Allclean
```

`0/` is generated from `0.orig/` and is gitignored. Edit `0.orig/`, then copy again (or re-run `buildMesh`).

## Run

Needs an OpenFOAM environment (this repo was run with the `opencfd/openfoam-run` image, v2412).

```bash
./buildMesh          # mesh + restore 0/ from 0.orig/
./Run                # decompose both regions, then chtMultiRegionFoam -parallel
./Continue           # resume from the latest processor time (does not decompose)
./reconstruct        # assemble missing processor times (both regions, lockstep)
./Allclean           # wipe mesh, 0/, logs, processors
```

`./Continue` sets `startFrom latestTime` and launches `mpirun` only. Raise `endTime` in `system/controlDict` first if the previous run already reached it. Do not run `./Run` to resume — it re-decomposes from `0/` and wipes processor times.

`./Run` decomposes **each region** (`decomposePar -region …`). Bare `decomposePar` looks for a default mesh and fails.

Processor count is `numberOfSubdomains` in **all three** files (keep them equal):

- `system/decomposeParDict` — `mpirun -n` in `./Run`
- `system/hot_fluid/decomposeParDict`
- `system/cold_fluid/decomposeParDict`

Current setting: **8** subdomains, method `hierarchical` with `n (2 2 2)` so both regions use the same xyz cuts. Independent `scotch` partitions put mapped twin faces on different ranks and the wall looks adiabatic. `./Run` overwrites previous `log.decomposePar.*`.

`./reconstruct` reconstructs **both regions at each time** (`reconstructPar -allRegions`). Doing one region’s full series first leaves the other empty if you Ctrl-C. A time dir that exists but is missing `T`/`U`/`phi` (interrupted write) is retried. Progress is printed to the terminal and appended to `log.reconstructPar`.

Serial solve (after `./buildMesh`):

```bash
chtMultiRegionFoam
```

## Edit boundary conditions

### Operating conditions (usual first stop)

`system/include/caseSettings` is included by the field files:

```
mdotCold        2.0;       // kg/s
mdotHot         2.0;       // kg/s
T_cold          300;       // K
T_hot           353;       // K
pRef            1e5;       // Pa
rhoInlet        997;       // kg/m3  (must match constant/<region> rho)
kInlet          0.01;
epsilonInlet    0.01;
```

Inlets use `flowRateInletVelocity` with `massFlowRate`; outlets use `inletOutlet` (velocity) and `fixedValue` `p_rgh`. Temperature at inlets is `fixedValue`.

After changing `0.orig/` on an existing mesh:

```bash
rm -rf 0
restore0Dir          # or: cp -r 0.orig 0
```

If processors already exist, re-decompose (or copy fields into `processor*/0/`).

### Thermal coupling (mapped walls)

`0.orig/include/thermalShellCoupled` is the mapped mixed-T coupling (`compressible::turbulentTemperatureRadCoupledMixed`).

By default there are **no** `thicknessLayers` / `kappaLayers`: heat transfers through turbulent `kappaEff` on both fluids (strongest coupling). To add a thin metal wall, put the **same** lists in the shared include on **both** mapped walls:

```
thicknessLayers (0.002);   // m
kappaLayers     (16);      // W/m/K  (stainless steel)
```

Leftover `solid` patches (unpaired baffle faces after `mergeOrSplitBaffles`) stay `zeroGradient` — they are not mapped. Heat crosses only `hot_fluid_to_cold_fluid` / `cold_fluid_to_hot_fluid`. Outlets use `zeroGradient` on \(T\).

### Per-patch field files

| Field | File | Inlet | Wall / mappedWall |
|-------|------|-------|-------------------|
| \(U\) | `0.orig/<region>/U` | `flowRateInletVelocity` | `noSlip` |
| \(T\) | `0.orig/<region>/T` | `fixedValue` | coupled on mapped walls; `zeroGradient` on outlets and `solid` |
| \(p_\mathrm{rgh}\) | `0.orig/<region>/p_rgh` | `fixedFluxPressure` | `fixedFluxPressure` |
| \(k\) | `…/k` | 5% intensity | `kqRWallFunction` |
| \(\varepsilon\) | `…/epsilon` | mixing length | `epsilonWallFunction` |

Mixing lengths (≈ 7% of nozzle diameter) live in the epsilon files, not in `caseSettings`:

- cold: `0.011` m (~160 mm nozzle)
- hot: `0.014` m (~200 mm nozzle)

Name `T` patches **exactly** for the mapped walls (`hot_fluid_to_cold_fluid`, `cold_fluid_to_hot_fluid`). Do not use a catch-all `".*"` on `T`: it can override the mapped-wall entry and leave the interface `zeroGradient` (adiabatic). Other fields may still use regex (`".*"`). Named inlets/outlets must match the STL names (`inlet_cold`, `outlet_hot`, …).

### Fluid properties

`constant/hot_fluid/thermophysicalProperties` and `constant/cold_fluid/thermophysicalProperties`. If you change density, update `rhoInlet` in `caseSettings` as well.

Turbulence model: `constant/<region>/turbulenceProperties`. Radiation is off in `radiationProperties`.

### Time control

`system/controlDict`: `endTime`, `writeInterval`, `maxCo`, `maxDeltaT`.

## Edit geometry

### STL layout and units

| File | Role |
|------|------|
| `solid.stl` | Enclosure baffle (tube bundle + shell). Splits hot/cold fluids. |
| `inlet_hot.stl` / `outlet_hot.stl` | Tube-side openings (patches). |
| `inlet_cold.stl` / `outlet_cold.stl` | Shell-side nozzle openings (patches). |

**Units:** `constant/triSurface/*.stl` are in **metres**. Snappy, feature extract, and seed points all use this.

`./buildMesh` does **not** copy or scale CAD into `triSurface/`. After a CAD change:

1. Export `solid.stl` and the four opening STLs.
2. Scale millimetres → metres (divide coordinates by 1000) into `constant/triSurface/`.
3. Keep the same file names, or update `snappyHexMeshDict` and `surfaceFeatureExtractDict` together.

Inlets/outlets must lie on `solid.stl` (same plane, same outline) so snappy can cut patches out of the wall.

### Bounding box

`system/blockMeshDict` must enclose the **metre** STLs with a small margin. Current box:

```
x: -0.05 .. 2.45
y: -0.30 .. 0.30
z: -0.35 .. 0.35
```

`(100 40 50)` is the background resolution. Increase it (or snappy surface levels) for a finer mesh.

### Seed points after a geometry change

If tubes, nozzles, or the shell move, the old seeds may land in the wrong fluid or in a solid. Update **both**:

- `system/snappyHexMeshDict` → `locationsInMesh`
- `system/topoSetDict` → `insidePoints` for `hot_fluid` and `cold_fluid`

Pick a point clearly inside the tube header for hot, and clearly inside the shell (e.g. in the inlet nozzle) for cold. If `splitMeshRegions` does not write both `constant/hot_fluid/polyMesh/points` and `constant/cold_fluid/polyMesh/points`, `buildMesh` exits and you should inspect `log.topoSet` / `log.splitMeshRegions`.

### Refinement

`system/snappyHexMeshDict`:

- `features` / `refinementSurfaces` levels (currently 2 on `solid` and inlets/outlets)
- `maxGlobalCells`
- `addLayers false` (no prism / inflation layers)

Feature-edge extraction: `system/surfaceFeatureExtractDict` (`includedAngle 150`).

### Rebuild

```bash
./Allclean
# refresh constant/triSurface/*.stl if geometry changed
./buildMesh
./Run
```

`./Allclean` removes region meshes and `0/` but not the STLs in `constant/triSurface/`.

## Patches after a successful mesh

| Region | Patch | Type |
|--------|--------|------|
| `cold_fluid` | `inlet_cold`, `outlet_cold` | `patch` |
| `cold_fluid` | `solid` | `wall` (outer shell / heads) |
| `cold_fluid` | `cold_fluid_to_hot_fluid` | `mappedWall` → hot |
| `hot_fluid` | `inlet_hot`, `outlet_hot` | `patch` |
| `hot_fluid` | `solid` | `wall` (headers / unused STL faces) |
| `hot_fluid` | `hot_fluid_to_cold_fluid` | `mappedWall` → cold |

## Notes

- There is no `solid` region in `constant/regionProperties`. Adding a meshed metal region would be a different case (resolved CHT), not this thermal-shell setup.
- `0.orig` is the source of truth for fields. Do not rely on editing live `0/` or `processor*/0/` across remeshes.
- Logs (`log.*`), `processor*`, `polyMesh`, and time directories are gitignored.
