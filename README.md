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

`solid.stl` is snapped as a **baffle**. After `topoSet` / `subsetMesh`, coincident baffle faces are merged back to internal faces (`mergeOrSplitBaffles`), then `splitMeshRegions` splits the mesh into `hot_fluid` and `cold_fluid`. Those coincident faces become a pair of `mappedWall` patches:

- `hot_fluid_to_cold_fluid`
- `cold_fluid_to_hot_fluid`

Heat crosses the interface with
`compressible::turbulentTemperatureRadCoupledMixed` plus a 1-D conduction layer (`thicknessLayers`, `kappaLayers`). That is the “thermal shell”: wall thickness and conductivity set the thermal resistance without resolving a solid mesh.

Remaining single-sided faces from the STL (shell outer wall, heads) stay as the `solid` wall patch and are adiabatic.

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
3. `snappyHexMesh -overwrite` — snap to geometry; `solid.stl` is a baffle `faceZone`.
4. `topoSet` — walk from a seed in each fluid (baffles stop the walk, so the split follows the STL).
5. `subsetMesh keepCells` — drop exterior cells; expose the STL as patch `solid`.
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
- Transient PIMPLE (`nOuterCorrectors 5` so mapped T is iterated each step), adjustable \(\Delta t\) with `maxCo 0.5`. `endTime` is in `system/controlDict`.
- Outlet area-average \(T\) and interface `wallHeatFlux` are logged from `controlDict` function objects.

## Directory map

```
geom/                          CAD export (millimetres) — source STLs
constant/triSurface/           STLs in metres — what snappy actually reads
constant/<region>/             thermo, turbulence, radiation, region polyMesh
0.orig/<region>/               initial/boundary fields (edit these)
0.orig/include/                shared thermal-shell BC
system/include/caseSettings    mass flow, T, p, inlet turbulence
system/blockMeshDict
system/snappyHexMeshDict
system/topoSetDict
system/controlDict
system/<region>/               fvSchemes, fvSolution, decomposeParDict
buildMesh  Run  reconstruct  Allclean
```

`0/` is generated from `0.orig/` and is gitignored. Edit `0.orig/`, then copy again (or re-run `buildMesh`).

## Run

Needs an OpenFOAM environment (this repo was run with the `opencfd/openfoam-run` image, v2412).

```bash
./buildMesh          # mesh + restore 0/ from 0.orig/
./Run                # decompose both regions, then chtMultiRegionFoam -parallel
./reconstruct        # assemble processor time directories (use after the run)
./Allclean           # wipe mesh, 0/, logs, processors
```

`./Run` decomposes **each region** (`decomposePar -region …`). Bare `decomposePar` looks for a default mesh and fails.

Processor count is `numberOfSubdomains` in **all three** files (keep them equal):

- `system/decomposeParDict` — `mpirun -n` in `./Run`
- `system/hot_fluid/decomposeParDict`
- `system/cold_fluid/decomposeParDict`

Current setting: **8** subdomains, method `hierarchical` with `n (2 2 2)` so both regions use the same xyz cuts. Independent `scotch` partitions put mapped twin faces on different ranks and the wall looks adiabatic. `./Run` overwrites previous `log.decomposePar.*`.

`./reconstruct` must also use `-region` (the script does this).

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

### Thermal shell (wall thickness and conductivity)

`0.orig/include/thermalShellCoupled` is the mapped mixed-T coupling (`compressible::turbulentTemperatureRadCoupledMixed`).

Put `thicknessLayers` / `kappaLayers` on **one region only** (here `0.orig/hot_fluid/T`). The same lists on both fluids stack two shells in series and cut the flux in half.

```
thicknessLayers (0.002);   // m
kappaLayers     (16);      // W/m/K  (stainless steel)
```

Several layers are allowed, e.g. metal + fouling: `thicknessLayers (0.002 5e-5); kappaLayers (16 0.5);`.

### Per-patch field files

| Field | File | Inlet | Wall / mappedWall |
|-------|------|-------|-------------------|
| \(U\) | `0.orig/<region>/U` | `flowRateInletVelocity` | `noSlip` |
| \(T\) | `0.orig/<region>/T` | `fixedValue` | shell couple on `.*_to_.*`; `zeroGradient` on `solid.*` |
| \(p_\mathrm{rgh}\) | `0.orig/<region>/p_rgh` | `fixedFluxPressure` | `fixedFluxPressure` |
| \(k\) | `…/k` | 5% intensity | `kqRWallFunction` |
| \(\varepsilon\) | `…/epsilon` | mixing length | `epsilonWallFunction` |

Mixing lengths (≈ 7% of nozzle diameter) live in the epsilon files, not in `caseSettings`:

- cold: `0.011` m (~160 mm nozzle)
- hot: `0.014` m (~200 mm nozzle)

Name `T` patches **exactly** (`hot_fluid_to_cold_fluid`, `cold_fluid_to_hot_fluid`). Do not use a catch-all `".*"` on `T`: it can override the mapped-wall entry and leave the interface `zeroGradient` (adiabatic). Other fields may still use regex (`"solid.*"`, `".*"`). Named inlets/outlets must match the STL names (`inlet_cold`, `outlet_hot`, …).

### Fluid properties

`constant/hot_fluid/thermophysicalProperties` and `constant/cold_fluid/thermophysicalProperties`. If you change density, update `rhoInlet` in `caseSettings` as well.

Turbulence model: `constant/<region>/turbulenceProperties`. Radiation is off in `radiationProperties`.

### Time control

`system/controlDict`: `endTime`, `writeInterval`, `maxCo`, `maxDeltaT`.

## Edit geometry

### STL layout and units

| File | Role |
|------|------|
| `solid.stl` | Wetted envelope: shell outer, tubes, headers. Snapped as a baffle. |
| `inlet_hot.stl` / `outlet_hot.stl` | Tube-side openings (patches). |
| `inlet_cold.stl` / `outlet_cold.stl` | Shell-side nozzle openings (patches). |

**Two copies, two unit systems:**

- `geom/*.stl` — millimetre CAD export (envelope about \(2400 \times 500 \times 600\,\mathrm{mm}\)).
- `constant/triSurface/*.stl` — **metres**. Snappy, feature extract, and seed points all use this.

`./buildMesh` does **not** copy or scale `geom/` → `triSurface/`. After a CAD change:

1. Export the five named STLs.
2. Scale millimetres → metres (divide coordinates by 1000) into `constant/triSurface/`.
3. Keep the same file names, or update `snappyHexMeshDict` and `surfaceFeatureExtractDict` together.

Inlets/outlets must lie on the `solid.stl` surface (same plane, same outline) so snappy can cut patches out of the wall.

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

- `features` / `refinementSurfaces` levels (currently 1–2 on `solid`, 2 on inlets/outlets)
- `maxGlobalCells`
- `addLayers` is off

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
