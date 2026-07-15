# Gas Leak Blast Damage Simulator — Project Summary

*Last updated: 2026-07-08 (rev 14)*

A physics-based training simulator for engineers and emergency responders to
evaluate structural damage and human casualties from indoor gas-leak explosions.
Built in Python with a Plotly Dash interactive web UI.

---

## Quick Start

```bash
cd blast_sim
source venv/bin/activate        # activate virtual environment
python main.py                  # start the server
# → open http://127.0.0.1:8050
```

---

## Architecture

```
blast_sim/
├── core/
│   ├── blast.py          # Blast physics (Kinney-Graham + Friedlander)
│   ├── geometry.py       # Parametric building mesh + SDOF panel properties
│   ├── materials.py      # Material library (concrete, brick, steel, glass)
│   ├── simulation.py     # SDOF panel solver + shear-building FEM
│   ├── casualties.py     # Room-level casualty model (Probit)
│   ├── persons.py        # Individual person injury assessment
│   └── rescue.py         # Dijkstra rescue-route pathfinding through room graph
├── viz/
│   └── dashboard.py      # Plotly Dash interactive UI (~1 600 lines)
├── assets/
│   └── style.css         # Custom dark-theme stylesheet
├── main.py               # Entry point
├── SUMMARY.md            # This file
└── requirements.txt
```

**Total project code: ~2 500 lines (excluding venv)**

---

## Physics Layers

### 1 · Blast Source — `core/blast.py`

| Item | Detail |
|---|---|
| Model | Kinney & Graham (1985) empirical scaling laws |
| Inputs | Charge mass (kg), explosive type, standoff distance (m) |
| Outputs | Peak side-on overpressure `Pso` (kPa), positive-phase duration `td` (ms), specific impulse `Is` (kPa·ms) |
| Waveform | Modified Friedlander exponential; waveform parameter `b` solved numerically from impulse constraint |
| Reflected pressure | Rankine-Hugoniot normal reflection, with cosine-law angle-of-incidence correction |
| Burst types | **Surface** (gas-leak ground burst — TNT × **1.8** for hemispherical confinement, UFC 3-340-02), **free-air** |
| Arrival time | Ambient speed of sound (340 m/s) |

**Impulse formula corrected (rev 13):** the Kinney-Graham specific-impulse denominator takes the
**cube root** of `1 + (Z/1.55)³`, not the square root. The previous square-root implementation
under-predicted impulse by up to ~2.5× at Z = 10 m/kg^⅓ (far field), which suppressed
throw-velocity (tertiary) injuries for distant people. The leading constant is now written
explicitly as 6.7 kPa·ms (= 0.067 bar·ms).

**Surface-burst factor corrected (rev 13):** was 2.0 (ideal rigid-ground reflection); now **1.8**
per UFC 3-340-02 practice — part of the blast energy is lost to ground shock and cratering.
The factor is exported as `SURFACE_BURST_FACTOR` and shared with the dashboard (no more
hardcoded 2.0 in the UI).

#### TNT Equivalency — `TNT_EQUIVALENCY` (UFC 3-340-02)

The charge mass is converted to an effective TNT mass before applying Kinney-Graham scaling:

| Explosive | Overpressure factor | Impulse factor |
|---|---|---|
| TNT | 1.00 | 1.00 |
| C4 (Composition C-4) | 1.34 | 1.19 |
| ANFO | 0.82 | 0.82 |
| PETN | 1.27 | 1.00 |
| Semtex 1A | 1.25 | 1.00 |
| RDX | 1.14 | 1.09 |
| HMX | 1.14 | 1.02 |
| Ammonium Nitrate (pure) | 0.42 | 0.45 |
| Methane / Natural Gas (VCE) | 0.48 | 0.48 |
| Propane / LPG (VCE) | 0.44 | 0.44 |

**Gas VCE entries added (rev 13):** since the simulator's primary use case is gas-leak
explosions, the explosive dropdown now includes vapour-cloud-explosion equivalencies. The
factors are energy-based — `yield efficiency (~4 %) × ΔHc / ΔH_TNT` — landing inside the
0.2–0.5 band typical of TNO Multi-Energy / Baker-Strehlow assessments (methane
ΔHc ≈ 50 MJ/kg → 0.48; propane ≈ 46 MJ/kg → 0.44). The charge-mass input is interpreted as
the mass of gas in the flammable cloud. Note the waveform is still a TNT-shaped Friedlander
pulse; the longer positive-phase duration of a real VCE is not reproduced (see Limitations).

The overpressure factor is used for structural response (`tnt_equivalent()` → `W_eff`).
The impulse factor is available in the dict but currently informational only (Kinney-Graham
impulse is derived consistently from `W_eff`). The sidebar displays `W_TNT` live as the user
changes explosive type or charge mass.

### 2 · Building Geometry — `core/geometry.py`

Two building generators are available, both producing the same `Panel / Column / Room` objects with full SDOF properties:

#### `create_building()` — Parametric grid generator
Configurable parameters:

- Width, depth, number of floors, floor height
- Wall and floor thickness
- Column size
- Rooms per floor (X × Y grid)
- Window fraction of exterior wall area
- Wall and floor material

#### `create_building_from_blueprint(blueprint)` — JSON-driven generator
Accepts a blueprint dict (see **Blueprint Format** below) that specifies
arbitrary room layouts per floor. Rooms are axis-aligned rectangles; any number
of rooms can be placed per floor and rooms can vary between floors.

Supports **variable floor heights** via the `floor_heights` list in the
building config — each entry overrides the floor height for that storey
independently. Floor `z_bot` values are accumulated correctly:
`z_bot[f] = Σ floor_heights[0..f-1]`.

Both generators produce:
- **Exterior wall panels** (solid + window sub-panels — *windows currently
  disabled, see the `WINDOWS_ENABLED` flag below*)
- **Interior partition walls** (room subdivisions)
- **Floor slabs** and **roof**
- **Columns** at every grid intersection

**Windows temporarily disabled (rev 14):** the `WINDOWS_ENABLED = False` flag at the top of
`viz/dashboard.py` forces `window_frac = 0` for both the parametric generator and blueprint
imports, removes the window-fraction control from the sidebar, and drops the window entries
from the 3D legend and filter. All window code paths (geometry, shielding, colouring) are
intact — set the flag to `True` to restore them once the glazing damage model is recalibrated.

Each panel is assigned SDOF structural properties at creation time using
**Timoshenko plate theory** (simply-supported boundary conditions with a
partial-fixity correction factor):

```
Flexural rigidity:  D = E·t³ / 12(1 − ν²)
Equivalent stiffness: Ke = BC_factor × D·A / (α·a⁴)
Equivalent mass:    Me = KLM × ρ·A·t
Yield displacement: wy = α·Py·a⁴ / (D · BC_factor)
Failure displacement: wf = wy × ductility_ratio (μ)
```

#### Blueprint Format — Native (blast-sim)
```json
{
  "building": {
    "n_floors": 2,
    "floor_height": 3.0,
    "floor_heights": [4.27, 3.05, 3.05],
    "wall_thickness": 0.20,
    "floor_thickness": 0.25,
    "col_size": 0.40,
    "wall_material": "concrete",
    "window_frac": 0.30,
    "total_occupants": 0
  },
  "floors": [
    {
      "rooms": [
        {"x_min": 0, "x_max": 5,  "y_min": 0, "y_max": 4, "name": "Room A"},
        {"x_min": 5, "x_max": 12, "y_min": 0, "y_max": 4, "name": "Room B"}
      ]
    }
  ]
}
```

`floor_heights` overrides `floor_height` per storey when provided.
Room coordinates start at (0, 0) and must tile the footprint without gaps.

#### Blueprint Format — OpenStudio FloorSpaceJS (auto-detected)
The **FloorSpaceJS** format exported by OpenStudio's floor-plan editor is also
accepted directly. It is detected automatically by the presence of a `stories`
top-level key. The converter (`parse_floorspacejs` in `dashboard.py`):

1. Iterates each `story` and reads its `geometry` (vertices, edges, faces)
2. Reconstructs each `space` polygon by walking the face edge list
3. Computes the axis-aligned bounding box of each polygon (rooms are
   approximated as rectangles for the blast simulation)
4. Converts coordinates from feet to metres (× 0.3048)
5. Maps each story's `floor_to_ceiling_height` to a per-floor height entry

Example FloorSpaceJS output: a 4-storey, 200 ft × 65 ft building with 18–20
spaces per floor is converted to a blueprint with
`floor_heights = [4.267, 3.048, 3.048, 3.048]` metres and 78 rooms total.

### 3 · Materials — `core/materials.py`

| Material | E (GPa) | ρ (kg/m³) | fy (MPa) | μ |
|---|---|---|---|---|
| Reinforced Concrete | 30 | 2 400 | 25 | 4.0 |
| Brick Masonry | 6 | 1 800 | **0.7** | 1.5 |
| Steel Plate | 200 | 7 850 | 250 | 12.0 |
| Annealed Glass | 70 | 2 500 | 30 | 1.0 |
| Structural Steel | 200 | 7 850 | 345 | 15.0 |

**Brick masonry `fy` corrected (rev 12):** The previous value of 6 MPa was 6–15× too high relative to the
BS EN 1996-1-1 design value of 0.4–1.0 MPa for unreinforced masonry in flexure. The updated value of
0.7 MPa (mid-range of the code band) gives a yield displacement ≈ 8× smaller and a failure displacement
≈ 12× smaller than before — brick walls now fail at much lower blast pressures, consistent with observed
masonry fragility in blast events.

### 4 · Structural Simulation — `core/simulation.py`

Two parallel models run simultaneously:

#### Panel SDOF (per-panel dynamics)
Each wall/floor panel solved as an independent single-degree-of-freedom system:

```
Me · ÿ + C · ẏ + K · y = P(t) · A
  where  Me = KLM · ρ·A·t  (KLM = 0.55, set once in geometry.py)
```

**Double-counted load factor removed (rev 13):** the Biggs equivalent SDOF system is
`KM·M·ÿ + KL·K·y = KL·F(t)`; dividing through by KL gives `KLM·M·ÿ + K·y = F(t)` with
`KLM = KM/KL`. Since `Me` already includes KLM, the load must be the **unfactored** `P(t)·A`.
The previous code additionally multiplied the load by `KL = 0.53`, under-loading every panel
by ~47% — panels now yield/fail at the correct (lower) blast pressures.

- **Three-phase exact integration (rev 14):**
  1. *Before arrival* (t < toa) — force and response are identically zero; skipped.
  2. *During the pulse* (toa → toa + td) — `solve_ivp` (RK45) with
     `max_step = td/20` and the Friedlander force evaluated directly (no
     interpolant), so the pulse can never fall between solver steps.
  3. *After the pulse* — closed-form damped free vibration from the
     end-of-pulse state; the peak (which usually occurs in the first free
     half-cycle) is sampled analytically at 256 points over one natural period.

  This replaced a single full-window `solve_ivp` call with unbounded
  `max_step`: RK45 grew its step size over the quiet pre-arrival phase and
  routinely leapt over the entire millisecond-scale pulse, so whether a panel
  was loaded at all depended on where the adaptive steps landed. The visible
  symptom was *inverted* damage — walls 5 m from the blast showing zero
  displacement under 1 000 kPa while windows 12 m away failed. (The rev-12
  "fine sub-grid" only refined the pressure interpolant and reported peak;
  it never constrained the integrator's own sampling, so it did **not**
  actually fix near-field under-damage as previously claimed.) The new scheme
  is also ~10× faster since the ODE solver only runs over the pulse window.
- **Sheltered-side factor**: `facing = dot(−r/R, n)`. If `facing < −0.1` the
  panel's back face is towards the blast; it receives 10% of the direct load.
- **Damage index** = peak displacement / yield displacement
  - DI ≥ 1 → yielded (structural damage)
  - DI ≥ μ → failed (element removed)

#### Shear-Building FEM (global lateral response)
N-story lateral model with inter-story stiffness from column `12EI/h³`:

```
[M]{ÿ} + [C]{ẏ} + [K]{y} = {F_blast(t)}
```

- Rayleigh damping (5% at 1st and 3rd modes)
- Solved with the same RK45 integrator, with `max_step = dt` (0.5 ms) so the
  adaptive stepper cannot leap over the blast pulse (rev 14)
- Outputs inter-story drift ratio per floor
- Limit states: Life Safety > 2.5%, Collapse prevention > 5.0%

#### Interior Blast Propagation
- Failed exterior panels → interior room pressure = 50% of exterior `Pso` (incident side-on overpressure, not reflected pressure)
- Intact exterior panels → 5% leakage (cracks, gaps)
- Failed interior walls → 80% pressure transmission to adjacent room

### 5 · Room Casualties — `core/casualties.py`

Baker (1983) Probit function for lung-haemorrhage fatality:

```
Y = −77.1 + 6.91 · ln(Pso_Pa)
P(fatality) = Φ(Y − 5)
```

Threshold-based linear-ramp models:

| Category | 0 % threshold | 100 % threshold | Formula |
|---|---|---|---|
| Severe injury | 35 kPa | 200 kPa | `(Pso − 35) / 165` |
| Minor injury  | 7 kPa  | 70 kPa  | `(Pso − 7) / 63`   |

Casualty categories are mutually exclusive (fatal > severe > minor):
`p_severe` has `p_fatal` subtracted; `p_minor` has both `p_fatal` and `p_severe` subtracted.
The metric tiles show **expected counts** (probability × room occupants), not raw percentages.

### 6 · Individual Person Injuries — `core/persons.py`

For each placed person:

**Shielding** — Möller-Trumbore ray-panel intersection test traces a ray from
the blast source to the person. Each intercepting wall attenuates the pressure
by a material-dependent transmission factor:

| Material | Transmission |
|---|---|
| Reinforced Concrete | 5% |
| Brick Masonry | 8% |
| Steel Plate | 3% |
| Structural Steel | 3% |
| Annealed Glass (intact) | 35% |
| Unknown material | 10% |
| Any failed panel | 50% (breach opening) |

**Transmission factors recalibrated (rev 13):** the previous values (concrete 12%, brick 22%,
glass 70%, failed 88%) were physically too high and inconsistent with the interior-propagation
model in `simulation.py`. An intact heavy wall passes only leakage pressure through cracks and
joints, so concrete is now **5%** — matching the 5% intact-exterior-wall leakage factor used
for room pressures — with the other materials scaled by areal mass / stiffness (steel is
airtight and heavier; brick has more joints; 6 mm glazing flexes and re-radiates). A failed
panel is now **50%**, matching the room-propagation breach factor: the wave diffracts and
expands through the opening rather than passing unimpeded.

**Co-located wall + window grouping (rev 13):** wall panels and their window sub-panels share
the same quad in the geometry model, so a single ray used to hit both and attenuate twice
(e.g. 5% × 35% through one physical wall). Hits are now grouped by panel centre and the most
transmissive path through each physical wall is used — the dominant leakage path is the
glazing, not the series product.

**Vertical blast paths (rev 13):** floor and roof slabs are now included in the ray test
(they were previously skipped), so a person on the floor above an interior blast is correctly
shielded by the intact concrete slab instead of receiving the full free-field pressure.

**Frontal area (rev 12):** `_FRONTAL_AREA = 0.6 m²` — representative standing-adult frontal area per UFC 3-340-02 / the implicit
assumptions register. The previous value of 0.7 m² over-estimated throw velocity by ~17%.

**Injury types computed:**

| Injury | Model |
|---|---|
| Eardrum rupture | Baker (1983) Probit: `Y = −15.6 + 1.93·ln(Pso_Pa)` |
| Pulmonary barotrauma | Baker (1983) Probit: `Y = −77.1 + 6.91·ln(Pso_Pa)` |
| Abdominal / GI injury | Linear threshold (150–350 kPa) |
| Blast traumatic brain injury | Linear threshold (80–430 kPa) |
| Whole-body displacement | Impulse × 0.6 m² ÷ body mass → throw velocity → impact severity |
| Fragment / debris | Distance to nearest failed panel |

**Overall severity:** Uninjured / Minor / Moderate / Severe / Fatal

---

## User Interface

### Layout
```
┌─ Sidebar (260 px) ─┬──── 3D View + Charts ────┬─ Right Panel (245 px) ─┐
│  [▶ Run Simulation] │  Interactive 3D building  │  [Injuries][Rescue Plan]│
│  ─────────────────  │  (click floor to place    │  ──────────────────────│
│  Building controls  │   👤/🚪 in 3D)            │  Injuries tab:         │
│  Blueprint JSON ↑   ├──────────────────────────┤   per-person injury    │
│  Blast controls     │  4-panel results chart    │   assessment cards     │
│  People section     │                           │  Rescue Plan tab:      │
│  Place-mode toggle  │                           │   priority-ordered     │
│  Floor plan (2D)    │                           │   rescue route cards   │
│  Exits section      │                           │                        │
│  (scrollable area)  │                           │                        │
└─────────────────────┴──────────────────────────┴────────────────────────┘
```

The sidebar scrolls as a whole (`overflow-y:auto` on `.sidebar`). The **Run Simulation button uses `position:sticky; top:0`** so it stays visible at the top while the user scrolls through the controls below it. This avoids all flex-height complexity that previously caused the button to disappear or the controls to be un-scrollable.

All numeric parameters (building dimensions, blast position, TNT mass, etc.) use **`dcc.Input(type='number', debounce=True)`** boxes instead of sliders. Values commit on Enter or focus-out. The floor-plan click grid and the two dropdowns (wall material, burst type) are unchanged.

### Controls
- **Building**: width, depth, floors, floor height (all number inputs), wall material (dropdown), room grid X/Y. *(The window-fraction input was removed in rev 14 — windows are disabled via `WINDOWS_ENABLED`.)*
- **Blueprint JSON** (below building controls):
  - **Drag-and-drop or file-select** — accepts any `.json` file in native
    blast-sim format or OpenStudio FloorSpaceJS format. Format is
    auto-detected. On upload the file is parsed, the `bld-floors` counter is
    set to the blueprint's floor count, the floor-selector is updated, and
    the **3D structural model is rendered immediately** (no need to run the
    simulation first). Status line shows filename, floor count, and room count.
  - **Clear Blueprint** — removes the loaded blueprint and reverts to the
    parametric grid generator for the next simulation run.
- **Blast source**: Explosive type (dropdown — TNT, C4, ANFO, PETN, Semtex 1A, RDX, HMX, Ammonium Nitrate),
  charge mass (kg), X/Y/Z position (all number inputs), burst type (dropdown: surface / free-air).
  A live label below the explosive dropdown shows the effective TNT equivalent mass (`W_TNT = mass × factor`)
  and the overpressure factor so the user understands the conversion at a glance.
- **People**:
  - **Floor selector** — auto-updates to match the current number of floors (reads `n_floors` from blueprint when loaded)
  - **Click floor plan** (in Person mode) — drops a person at the clicked grid position on the 2D plan (0.5 m resolution). Blueprint rooms are drawn as blue-outlined shapes; parametric grid rooms are grey.
  - **Click 3D model** (in Person mode) — drops a person at the clicked position. An invisible 1 m click grid at standing height (floor_bot + 1.2 m) on every floor level provides hover targets; clicking any point on it places a person at that (x, y, z).
  - **＋ Add** — places a person at (3.0, 3.0) on the selected floor
  - **Clear** — removes all people
  - **×** — removes one person from the list
  - Each person in the list shows name, and X / Y / Z in metres
- **Place mode toggle** (above the floor plan): switches between `👤 Person` and `🚪 Exit` modes — controls what *both* floor-plan clicks and 3D model clicks create
- **Exits & Entrances**:
  - **Click floor plan** (in Exit mode) — drops an exit at the clicked position on the 2D plan
  - **Click 3D model** (in Exit mode) — drops an exit at the clicked floor position; floor index is read from the `customdata` attached to each click-grid point
  - **＋ Add Exit** — places an exit at (0.0, 5.0) on the selected floor
  - **Clear** — removes all exits
  - **×** — removes one exit from the list
  - Each exit shows name and X / Y / floor in metres

### Outputs
- **3D view**: each panel type has its own colour family so damage is immediately readable at a glance:

  | Element | Intact colour | Yielded (DI ≥ 1) | Failed |
  |---|---|---|---|
  | Exterior wall | Green | Orange | Deep red |
  | Interior wall | Teal | Amber | Deep red |
  | Window (glass) — *disabled rev 14* | Steel blue (22% opacity) | Pale orange (55% opacity) | Deep red (55% opacity) |
  | Floor / Roof slab | Slate blue | Amber | Deep red |
  | Column | Dark blue-grey | — | Red line |

  A colour legend is displayed in the top-left corner of the 3D view; the window legend entry appears only when window panels exist. Hovering a panel shows its type, material, dimensions, thickness, and DI. People are shown as coloured spheres (colour = injury severity).

  A **filter dropdown** (inside the 3D view, below the legend) lets you isolate damage to specific element types without any page reload — all client-side via Plotly `updatemenus`:

  | Filter option | Elements shown |
  |---|---|
  | All elements | Everything (reset — added rev 14) |
  | Exterior walls | Exterior wall panels only |
  | Interior walls | Interior partition panels only |
  | All walls | Exterior + interior walls (shown while windows are disabled) |
  | Windows / Walls & windows | Only offered when window panels exist |
  | Floors & roof | Floor slabs and roof only |

  Columns, people, and blast source marker are always visible regardless of the filter.
- **Metric tiles**: four live indicators after each simulation run:
  - **Fatal** — count of placed people whose `overall_severity == 'Fatal'`
  - **Severe Inj.** — count of placed people with `overall_severity` of `'Severe'` or `'Moderate'`
  - **Minor Inj.** — count of placed people with `overall_severity == 'Minor'`
  - **Stability %** — overall structural health index: `(panel_health × 0.6 + drift_health × 0.4) × 100`, where `panel_health = 1 − failed_panels / total_panels` and `drift_health = max(0, 1 − max_story_drift / 5%)`. Green > 70 %, amber 40–70 %, red < 40 %.
- **Results charts**: room overpressures, per-room casualties (stacked bar: fatal / severe / minor), panel damage index histogram, story drift ratios
- **People injury cards** (right panel — Injuries tab): per-person card showing:
  - Raw pressure, effective pressure after wall shielding, number of walls intercepted
  - Eardrum rupture probability (Baker 1983)
  - Pulmonary barotrauma severity and probability
  - Abdominal / GI injury risk
  - Blast TBI risk
  - Throw velocity and tertiary impact severity
  - Debris / fragment risk from nearest failed panel
  - Overall severity (Uninjured / Minor / Moderate / Severe / Fatal) and survival probability
- **Rescue Plan** (right panel — Rescue Plan tab): priority-ordered list of rescue routes computed after simulation. For each person:
  - **Priority level** (1 Immediate → 5 Expectant) colour-coded by urgency
  - Target exit name and total route cost
  - Room-by-room path description (e.g., `R0(F1) → R1(F1)`)
  - Hazard warning if the route passes through damaged or failed panels
  - **Floor plan overlay**: dotted coloured route lines drawn from each person through room centroids to their nearest exit (current floor only)
  - Rescue routes update automatically when exits are placed / removed or when the simulation is re-run; no separate refresh needed

---

## References

| Source | Used for |
|---|---|
| Kinney & Graham (1985) *Explosive Shocks in Air* | Blast scaling laws (overpressure, duration, impulse) |
| Baker et al. (1983) *Explosion Hazards and Evaluation* | Lung-haemorrhage and eardrum Probit functions |
| UFC 3-340-02 (2008) *Structures to Resist the Effects of Accidental Explosions* | SDOF transformation factors (KL, KM, KLM), limit states, TNT equivalency factors, person frontal area (0.6 m²) |
| BS EN 1996-1-1 *Design of Masonry Structures* | Brick masonry flexural tensile strength (0.4–1.0 MPa) |
| Timoshenko & Woinowsky-Krieger (1959) *Theory of Plates and Shells* | Plate deflection and moment coefficients (α, β tables) |
| Möller & Trumbore (1997) | Ray-triangle intersection algorithm for shielding |

---

## Known Limitations & Simplifications

1. **No gas dispersion model — partially addressed (rev 13)** — the explosive dropdown now includes **Methane / Natural Gas (VCE, ×0.48)** and **Propane / LPG (VCE, ×0.44)** equivalencies inside the recommended 0.2–0.5 band, so gas-leak scenarios no longer need a manual correction factor. The source is still a point charge with a TNT-shaped Friedlander waveform, however: the lower-peak / longer-duration signature of a real volume deflagration (and where the gas actually accumulated) is not reproduced. A full Baker-Strehlow-Tang or Multi-Energy source model remains on the roadmap.
2. **Plane-stress shell only** — panels modelled as thin plates; no membrane-bending coupling or geometric nonlinearity.
3. **Independent panels** — no load redistribution between panels after failure; in reality failed panels transfer load to edges and adjacent structure.
4. **Linear elastic material** — plasticity is approximated by a ductility-based failure criterion, not a full elasto-plastic constitutive model.
5. **Horizontal blast path only — fixed (rev 13)** — the shielding ray test now includes floor and roof slabs, so vertical blast propagation between storeys is attenuated by the slab material (5% through intact concrete, 50% through a failed slab). Note the *room-pressure* propagation model (`_interior_pressure`) still only propagates through failed interior *walls*, not failed slabs — room-level casualties on other floors remain driven by each room's own exterior loading.
6. **No fire / thermal effects** — gas explosions often trigger secondary fires; burn injuries are not modelled.
7. **No progressive collapse** — once panels fail they are removed, but the load path redistribution and potential pancake collapse are not simulated.
8. **People position is fixed** — persons cannot be moved after placement except by deleting and re-adding; there is no drag-to-reposition interaction on the floor plan.
9. **Floor plan click resolution is 0.5 m** — the click grid uses 0.5 m spacing so the placed position snaps to the nearest grid point, not the exact pixel clicked. This is accurate enough for training purposes.
10. **Floor plan click uses invisible marker grid** — the floor plan detects clicks via a dense 0.5 m invisible `Scatter` marker grid rather than a polygon, to work around a Plotly limitation where polygon `clickData` returns the nearest vertex rather than the actual click coordinates.
11. **Floor plan room shapes use `layer='below'`** — Plotly layout shapes default to `layer='above'`, which paints the dark room-fill rectangles on top of person marker dots. All room shapes are explicitly set to `layer='below'` so trace markers (people dots, blast star) render above the grid lines.
12. **Floor plan aspect ratio is free** — `scaleanchor='y'` was removed from the floor-plan x-axis. The locked 1:1 scale forced the plot area taller than its 200 px box for wide buildings, clipping the bottom. Without the lock the building always fits; x and y axes simply stretch to fill the available space independently.
13. **Floor plan redraws on floor change** — `floor-selector` is now an `Input` (not `State`) of `update_people_ui`, so switching floors immediately redraws the plan showing only the people on that floor.
14. **`build_results_figure` tuple bug fixed** — Casualty bar-chart loop `for key, name, color in [...]` was given 2-element tuples; each tuple now has all three values `(key, label, color)`.
15. **Sidebar scrolling** — The sidebar uses CSS Grid (`grid-template-rows: auto 1fr`). The button row takes its natural height; the controls row takes the rest with `overflow-y: auto; min-height: 0` (class `sidebar-scroll`). CSS Grid is more reliable than flexbox for this layout because grid rows do not have the implicit `min-height: auto` growth that flex items do.
16. **Floor plan axis** — Range is fixed to the building footprint `[−0.5, width+0.5]` × `[−0.5, depth+0.5]`. Scroll-zoom (`scrollZoom: true`) lets the user zoom out to see the blast source if it is outside the building. The blast marker is visible in the floor plan only when it lies within the building footprint; otherwise it is clipped by the axis range.
17. **People-store delete-on-re-render bug fixed** — When `update_people_ui` rewrites `people-list-sidebar.children` (e.g. on floor change), Dash treats the `×` delete buttons as *new* components and fires `manage_people` with `n_clicks = 0`. The old code deleted the person even at `n_clicks = 0`. Fix: `ctx.triggered[0]['value'] < 1 → return no_update`. Only a genuine user click (`n_clicks ≥ 1`) now removes a person. All other no-op paths also return `no_update` to prevent any stale-store overwrite.
18. **UI zoom hack removed (rev 14)** — the former `html { zoom: 0.85 }` density hack interacted badly with the `calc(100vh - …)` panel heights: Chrome does not rescale `vh` units under `zoom`, so the whole app fell ~15% short of the window and left a dead band at the bottom. The zoom was removed; the layout now genuinely fills the viewport and the small font sizes render at their intended size.
19. **Blast source in 3D view** — `build_3d_figure` accepts an optional `blast_pos` tuple. After running the simulation the blast source is shown as a red diamond marker with a "💥" label in the 3D structural model.
20. **Panel Response tab removed** — The right-panel `dbc.Tabs` wrapper and Panel Response tab (pressure/displacement time history of the most-damaged panel) have been removed. Error tracebacks from the simulation callback are printed to the server console instead of the UI.
21. **Metric tiles show actual person counts** — Fatal, Severe Inj., and Minor Inj. tiles show the count of individually placed people at that severity level, not statistical room-occupancy percentages.
22. **3D view — per-type colour coding and filter dropdown** — Each panel type (exterior wall, interior wall, window, floor/roof) now uses a distinct colour family so damage to different structural elements is immediately distinguishable. Windows use low opacity (22%) when intact and higher opacity (55%) when failed. A legend is shown in the top-left of the 3D view. Hover text now includes material name, panel dimensions, thickness, and DI.
23. **Exits & rescue routes** — `core/rescue.py` implements Dijkstra pathfinding through a room-adjacency graph. Edge costs: intact wall = 1.0, damaged (DI ≥ 0.5) = 2.0, failed wall = 5.0. Exits are placed by the user on the floor plan (in Exit mode). After simulation the Rescue Plan tab shows priority-ordered routes (1 Immediate → 5 Expectant) and dotted route lines are drawn on the floor plan. Room data (`x_min`, `x_max`, `y_min`, `y_max`, `floor_idx`, `z_bot`, `z_top`) is serialised into `sim-store` so routes can be re-computed whenever exits or people change without re-running the simulation. Panel `room_inside` / `room_outside` IDs are also stored so the graph builder can identify which rooms each interior wall connects.
24. **Interior-wall sheltered-side check** — each interior partition is assigned a single outward normal (+x or +y). The sheltered-side heuristic (10% load when the blast is behind the normal direction) is therefore orientation-dependent; the same interior wall receives full direct load from one room and 10% from the other. Interior loads are predominantly driven by the interior-pressure propagation model rather than direct SDOF integration, so the impact on final room pressures is small, but the SDOF damage index for individual interior panels can be asymmetric.
25. **Sheltered-side factor** *(historical — superseded by the rev-14 integrator, note 35)* — `facing_factor` (0.1 for back-facing panels, 1.0 otherwise) is applied to the panel's peak reflected pressure before integration. The rev-12 fine/coarse interpolant branching this note originally described no longer exists; the force is now evaluated analytically inside the pulse window.
26. **Blueprint rooms approximated as axis-aligned bounding boxes** — `create_building_from_blueprint` and the FloorSpaceJS converter both represent rooms as rectangles (`x_min`, `x_max`, `y_min`, `y_max`). Non-rectangular or L-shaped spaces (common in FloorSpaceJS plans) are approximated by their bounding box. This over-estimates the room area and can create small geometry mismatches at re-entrant corners, but does not affect blast load calculations significantly since overpressure attenuates smoothly with distance.
27. **FloorSpaceJS cell-to-room lookup uses first-match** — the interior grid-cell → room mapping iterates rooms in definition order and returns the first match. When two bounding boxes overlap (e.g. a wide corridor and a narrow room sharing an edge), the first room wins and the overlapping cell is assigned to it. This can produce extra or missing interior partition walls in overlap regions, but is acceptable for training-level fidelity.
28. **3D click-to-place uses invisible Scatter3d grid** — clicking on panel meshes (Mesh3d traces) in Plotly 3D does not return reliable x/y/z coordinates. Instead, an invisible dense 1 m Scatter3d grid is added at standing height (floor_bot + 1.2 m) on each floor. Clicks on this grid are identified by the presence of `customdata` in the click event. Clicking directly on a wall surface (no customdata) is ignored. This means the floor area must be unobstructed for a click to register — if the user clicks on a panel face rather than the floor space, the placement is silently ignored and they should click again in open floor area.
29. **Blueprint-driven 3D preview uses `allow_duplicate`** — `manage_blueprint` and `run_sim` both update `graph-3d`. Dash 4.x requires the non-primary callback to set `allow_duplicate=True` on that output. `run_sim` holds the primary output; `manage_blueprint` uses `allow_duplicate=True`. If both fire in the same update cycle (e.g. uploading a blueprint and immediately clicking Run), `run_sim`'s result overwrites the preview — the correct behaviour.
30. **Kinney-Graham impulse cube-root fix (rev 13)** — `specific_impulse` used `sqrt(1+(Z/1.55)³)` where the published formula has a cube root, under-predicting far-field impulse by up to ~2.5× at Z = 10. Affects throw-velocity (tertiary) injury for distant people; overpressure-driven injuries were unaffected (different formula).
31. **SDOF load factor double-count fix (rev 13)** — `_integrate_panel` multiplied the blast force by `KL = 0.53` even though the equivalent mass already includes `KLM = KM/KL`. In the Biggs formulation the load on the KLM-transformed system is unfactored, so every panel was under-loaded by ~47%. Expect noticeably more panel damage at the same charge mass than in rev ≤ 12 runs.
32. **Shielding transmission factors recalibrated (rev 13)** — intact-wall transmissions lowered to leakage-level values consistent with the interior-propagation model (concrete 12%→5%, brick 22%→8%, steel 10%→3%, glass 70%→35%) and failed-panel transmission lowered 88%→50% (a breach diffracts the wave; it is not free field). Wall+window sub-panels sharing one quad are grouped so a single physical wall attenuates once (most-transmissive path), not twice in series.
33. **Surface-burst factor 2.0 → 1.8 (rev 13)** — ideal rigid-ground doubling replaced by the UFC 3-340-02 value of 1.8 accounting for ground-shock/cratering energy loss. Exported as `SURFACE_BURST_FACTOR`; the dashboard's W_TNT display and store round-trip use the constant instead of a hardcoded 2.0.
34. **Column-check dead code removed (rev 13)** — `_check_columns` computed an unused `trib_area` accumulator (loop with an unconditional `break`); the actual demand has always used the `floor_height × min(5·b, 3 m)` tributary strip. The dead loop and unused imports were deleted; results are unchanged.
35. **Adaptive-integrator pulse skipping fixed (rev 14)** — the root cause of "panels near the blast show less damage than panels far away": `solve_ivp` with unbounded `max_step` grew its step size over the quiet pre-arrival phase and jumped clean over the millisecond-scale blast pulse, so panel loading was effectively a lottery (a wall 5 m away received 1 093 kPa and computed 0.000 mm displacement — 0 of its 8 solver steps landed inside the 2.95 ms pulse). The SDOF solver now integrates in three exact phases (zero / forced-with-`max_step=td/20` / analytic free vibration) and the shear-building FEM solver got `max_step = dt`. See § 4.
36. **Windows temporarily disabled (rev 14)** — `WINDOWS_ENABLED = False` in `viz/dashboard.py` removes glass sub-panels from the parametric generator, blueprint imports (both formats), the sidebar control, the 3D legend, and the filter dropdown. All code paths remain; flip the flag to restore. Rationale: annealed glass fails at a few kPa, so window failures dominated the damage view and (before note 35 was fixed) amplified the confusing far-field failures.
37. **Custom stylesheet was never being served — fixed (rev 14)** — `dash.Dash(__name__)` resolves the assets folder relative to the module that creates the app (`viz/`), which has no `assets/` directory, so `/assets/style.css` returned an error and the entire custom stylesheet silently never loaded; the UI was rendering with the CYBORG CDN theme plus inline styles only. `create_app` now passes `assets_folder` pointing at the repo-root `assets/`.
38. **UI refresh (rev 14)** — section headers with accent bars and icons; the place-mode toggle is a proper segmented control; secondary buttons share `.btn-ghost` CSS classes instead of repeated inline styles; metric tiles have severity-coloured top borders; right-panel tabs restyled; friendly empty-state placeholders (3D pane, Injuries, Rescue Plan); Dash 4 `dash-dropdown-*` components themed (the legacy `.Select-*` CSS rules targeted markup that Dash 4 no longer renders and were dead).
39. **`.claude/launch.json` added (rev 14)** — defines the `blast-sim` dev-server launch configuration (`./venv/bin/python main.py`, port 8050) used for in-IDE preview.

---

## Suggested Next Steps

### Short term (1–2 weeks)

- [ ] **Drag-to-reposition people on floor plan**
  Implement Plotly `relayoutData` or a click-then-move interaction so users can
  drag existing person markers rather than delete-and-re-add. Store selected
  person ID in a `dcc.Store` and update coordinates on the second click.

- [ ] **Editable person name and body mass**
  Add an inline edit mode in the people list — clicking a name opens a small
  form with name and mass (kg) fields. Body mass affects the throw-velocity
  tertiary injury calculation.

- [ ] **Vapour cloud explosion (VCE) source model**
  Replace the TNT point-charge with a Baker-Strehlow-Tang (BST) or Multi-Energy
  method model. Accept gas type, volume, and confinement factor as inputs and
  compute the appropriate TNT equivalency and overpressure-distance curve.

- [ ] **Scenario save / load**
  Serialize the full scenario (building config + blast + people) to a JSON file
  so users can save, reload, and compare runs. Add "Export Results" button that
  writes a PDF or CSV summary report.

- [ ] **Animation / time slider**
  Add a `dcc.Slider` that scrubs through the simulation time. Update the 3D
  view to show panel displacement at each time step (stored from the SDOF
  solution), making the propagation of the blast wave visible.

### Medium term (1–2 months)

- [ ] **Progressive collapse model**
  After panel failure, redistribute loads to adjacent elements. Use an
  iterative static analysis (or continue the dynamic simulation with updated
  stiffness matrix) to check for cascade failures — particularly important for
  multi-storey buildings.

- [ ] **Proper nonlinear FEM**
  Replace the per-panel SDOF model with a full nonlinear finite element solver
  (e.g., integrate with **OpenSeesPy** or implement a corotational shell
  element). This would capture load redistribution, membrane action, and
  post-yield behaviour.

- [ ] **CFD gas dispersion pre-processor**
  Add a pre-simulation step using a simplified Gaussian plume or zone model to
  determine where gas accumulates before ignition. The gas concentration map
  then drives which rooms have internal blast sources, not just the exterior.

- [ ] **Multi-floor person tracking**
  Allow people to be placed on different floors simultaneously and show all
  floors in the 3D view. Add a "floor plan stack" view that displays all floors
  as layered 2D plans.

- [ ] **Secondary fire / burn injuries**
  After the blast wave, model fire spread using a simple room-to-room
  propagation model. Add burn injury categories (1st/2nd/3rd degree) based on
  thermal flux and exposure duration.

### Long term (3–6 months)

- [ ] **Building database / library**
  Pre-load standard building archetypes (residential apartment block,
  commercial office, industrial warehouse, hospital) with calibrated structural
  properties so users can select a type rather than configure everything manually.

- [ ] **Probabilistic / Monte-Carlo mode**
  Run hundreds of scenarios with randomly varied TNT equivalent, blast location,
  and occupant positions to produce probability distributions of casualties and
  damage — more useful for risk assessment than a single deterministic run.

- [ ] **Real-time multi-user mode**
  Deploy as a hosted web application (e.g., on AWS or Azure) so multiple
  trainees can run and compare scenarios simultaneously during a training session.
  Add a session management panel for an instructor to control scenarios.

- [ ] **Calibration and validation**
  Compare simulation outputs against published experimental data (e.g., DTRA
  blast test series, BlastEM database) and adjust model parameters. Produce a
  validation report with error bounds.

- [ ] **Export to standard formats**
  Write a converter that exports the building mesh and blast loads to LS-DYNA
  keyword format (`.k` file) or Abaqus input format so the simplified model can
  be used as a starting point for a high-fidelity commercial solver run.
