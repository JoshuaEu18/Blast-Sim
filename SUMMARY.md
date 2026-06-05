# Gas Leak Blast Damage Simulator — Project Summary

*Last updated: 2026-06-05 (rev 8)*

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
│   └── dashboard.py      # Plotly Dash interactive UI (~950 lines)
├── assets/
│   └── style.css         # Custom dark-theme stylesheet
├── main.py               # Entry point
├── SUMMARY.md            # This file
└── requirements.txt
```

**Total project code: ~1 850 lines (excluding venv)**

---

## Physics Layers

### 1 · Blast Source — `core/blast.py`

| Item | Detail |
|---|---|
| Model | Kinney & Graham (1985) empirical scaling laws |
| Inputs | TNT equivalent mass (kg), standoff distance (m) |
| Outputs | Peak side-on overpressure `Pso` (kPa), positive-phase duration `td` (ms), specific impulse `Is` (kPa·ms) |
| Waveform | Modified Friedlander exponential; waveform parameter `b` solved numerically from impulse constraint |
| Reflected pressure | Rankine-Hugoniot normal reflection, with cosine-law angle-of-incidence correction |
| Burst types | **Surface** (gas-leak ground burst — TNT doubled for hemispherical confinement), **free-air** |
| Arrival time | Ambient speed of sound (340 m/s) |

### 2 · Building Geometry — `core/geometry.py`

Parametric generator. Configurable parameters:

- Width, depth, number of floors, floor height
- Wall and floor thickness
- Column size
- Rooms per floor (X × Y grid)
- Window fraction of exterior wall area
- Wall and floor material

Generates:
- **Exterior wall panels** (solid + window sub-panels)
- **Interior partition walls** (room subdivisions)
- **Floor slabs** and **roof**
- **Columns** at every grid intersection

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

### 3 · Materials — `core/materials.py`

| Material | E (GPa) | ρ (kg/m³) | fy (MPa) | μ |
|---|---|---|---|---|
| Reinforced Concrete | 30 | 2 400 | 25 | 4.0 |
| Brick Masonry | 6 | 1 800 | 6 | 1.5 |
| Steel Plate | 200 | 7 850 | 250 | 12.0 |
| Annealed Glass | 70 | 2 500 | 30 | 1.0 |
| Structural Steel | 200 | 7 850 | 345 | 15.0 |

### 4 · Structural Simulation — `core/simulation.py`

Two parallel models run simultaneously:

#### Panel SDOF (per-panel dynamics)
Each wall/floor panel solved as an independent single-degree-of-freedom system:

```
Me · ÿ + C · ẏ + K · y = KL · P(t) · A
  where  Me = KLM · ρ·A·t  (KLM = 0.55, set once in geometry.py)
         KL = 0.53          (load transformation factor, fixed-fixed)
```

- Integrated with `scipy.integrate.solve_ivp` (RK45)
- Blast pressure interpolated from the Friedlander waveform. The coarse output
  time step (dt = 0.5 ms) under-samples near-field pulses: at R < ~3 m the
  positive-phase duration td can be shorter than dt, so the peak is missed
  entirely on the coarse grid. When td < 10 · dt the simulation inserts a fine
  sub-grid at td/20 resolution over the positive-phase window before building
  the interpolant, ensuring near-field panels are not under-damaged relative to
  far-field panels.
- Panels sheltered from the blast receive 10% of the direct load
- **Damage index** = peak displacement / yield displacement
  - DI ≥ 1 → yielded (structural damage)
  - DI ≥ μ → failed (element removed)

#### Shear-Building FEM (global lateral response)
N-story lateral model with inter-story stiffness from column `12EI/h³`:

```
[M]{ÿ} + [C]{ẏ} + [K]{y} = {F_blast(t)}
```

- Rayleigh damping (5% at 1st and 3rd modes)
- Solved with the same RK45 integrator
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
| Reinforced Concrete | 12% |
| Brick Masonry | 22% |
| Steel Plate | 10% |
| Structural Steel | 10% |
| Annealed Glass (intact) | 70% |
| Any failed panel | 88% (open gap) |

**Injury types computed:**

| Injury | Model |
|---|---|
| Eardrum rupture | Baker (1983) Probit: `Y = −15.6 + 1.93·ln(Pso_Pa)` |
| Pulmonary barotrauma | Baker (1983) Probit: `Y = −77.1 + 6.91·ln(Pso_Pa)` |
| Abdominal / GI injury | Linear threshold (150–350 kPa) |
| Blast traumatic brain injury | Linear threshold (80–430 kPa) |
| Whole-body displacement | Impulse ÷ specific body mass → throw velocity → impact severity |
| Fragment / debris | Distance to nearest failed panel |

**Overall severity:** Uninjured / Minor / Moderate / Severe / Fatal

---

## User Interface

### Layout
```
┌─ Sidebar (260 px) ─┬──── 3D View + Charts ────┬─ Right Panel (245 px) ─┐
│  [▶ Run Simulation] │  Interactive 3D building  │  [Injuries][Rescue Plan]│
│  ─────────────────  │  (colour-coded damage)    │  ──────────────────────│
│  Building controls  ├──────────────────────────┤  Injuries tab:         │
│  Blast controls     │  4-panel results chart    │   per-person injury    │
│  People section     │                           │   assessment cards     │
│  Place-mode toggle  │                           │  Rescue Plan tab:      │
│  Floor plan         │                           │   priority-ordered     │
│  Exits section      │                           │   rescue route cards   │
│  (scrollable area)  │                           │                        │
└─────────────────────┴──────────────────────────┴────────────────────────┘
```

The sidebar scrolls as a whole (`overflow-y:auto` on `.sidebar`). The **Run Simulation button uses `position:sticky; top:0`** so it stays visible at the top while the user scrolls through the controls below it. This avoids all flex-height complexity that previously caused the button to disappear or the controls to be un-scrollable.

All numeric parameters (building dimensions, blast position, TNT mass, etc.) use **`dcc.Input(type='number', debounce=True)`** boxes instead of sliders. Values commit on Enter or focus-out. The floor-plan click grid and the two dropdowns (wall material, burst type) are unchanged.

### Controls
- **Building**: width, depth, floors, floor height (all number inputs), wall material (dropdown), room grid X/Y, window fraction
- **Blast source**: TNT equivalent, X/Y/Z position (all number inputs), burst type (dropdown: surface / free-air)
- **People**:
  - **Floor selector** — auto-updates to match the current number of floors
  - **Click floor plan** (in Person mode) — drops a person at the clicked grid position (0.5 m resolution)
  - **＋ Add** — places a person at (3.0, 3.0) on the selected floor
  - **Clear** — removes all people
  - **×** — removes one person from the list
  - Each person in the list shows name, and X / Y / Z in metres
- **Place mode toggle** (above the floor plan): switches between `👤 Person` and `🚪 Exit` modes to control what a floor-plan click creates
- **Exits & Entrances**:
  - **Click floor plan** (in Exit mode) — drops an exit at the clicked position
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
  | Window (glass) | Steel blue (22% opacity) | Pale orange (55% opacity) | Deep red (55% opacity) |
  | Floor / Roof slab | Slate blue | Amber | Deep red |
  | Column | Dark blue-grey | — | Red line |

  Windows are nearly transparent when intact so you can see inside; they become opaque and red when they fail, making blast breaches obvious. A colour legend is displayed in the top-left corner of the 3D view. Hovering a panel shows its type, material, dimensions, thickness, and DI. People are shown as coloured spheres (colour = injury severity).

  A **filter dropdown** (inside the 3D view, below the legend) lets you isolate damage to specific element types without any page reload — all client-side via Plotly `updatemenus`:

  | Filter option | Elements shown |
  |---|---|
  | Exterior walls | Exterior wall panels only |
  | Interior walls | Interior partition panels only |
  | Windows | Glass window panels only |
  | Floors & roof | Floor slabs and roof only |
  | Walls & windows | Exterior walls + interior walls + windows |

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
| UFC 3-340-02 (2008) *Structures to Resist the Effects of Accidental Explosions* | SDOF transformation factors (KL, KM, KLM), limit states |
| Timoshenko & Woinowsky-Krieger (1959) *Theory of Plates and Shells* | Plate deflection and moment coefficients (α, β tables) |
| Möller & Trumbore (1997) | Ray-triangle intersection algorithm for shielding |

---

## Known Limitations & Simplifications

1. **No gas dispersion model** — blast source is a point charge (TNT equivalent), not a volume deflagration. A real gas-leak explosion has a lower peak pressure but longer duration than TNT; a correction factor (typically 0.2–0.5 TNT equivalent for vapour cloud explosions) should be applied.
2. **Plane-stress shell only** — panels modelled as thin plates; no membrane-bending coupling or geometric nonlinearity.
3. **Independent panels** — no load redistribution between panels after failure; in reality failed panels transfer load to edges and adjacent structure.
4. **Linear elastic material** — plasticity is approximated by a ductility-based failure criterion, not a full elasto-plastic constitutive model.
5. **Horizontal blast path only** — the shielding model skips floor/roof panels, so vertical blast propagation (e.g., blast entering through a roof) is not captured.
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
18. **UI zoom** — `html { zoom: 0.85 }` in the CSS scales the entire interface to 85%, giving more visible space for all panels (Chrome / Edge / Safari; Firefox ignores `zoom` but the layout still functions).
19. **Blast source in 3D view** — `build_3d_figure` accepts an optional `blast_pos` tuple. After running the simulation the blast source is shown as a red diamond marker with a "💥" label in the 3D structural model.
20. **Panel Response tab removed** — The right-panel `dbc.Tabs` wrapper and Panel Response tab (pressure/displacement time history of the most-damaged panel) have been removed. Error tracebacks from the simulation callback are printed to the server console instead of the UI.
21. **Metric tiles show actual person counts** — Fatal, Severe Inj., and Minor Inj. tiles show the count of individually placed people at that severity level, not statistical room-occupancy percentages.
22. **3D view — per-type colour coding and filter dropdown** — Each panel type (exterior wall, interior wall, window, floor/roof) now uses a distinct colour family so damage to different structural elements is immediately distinguishable. Windows use low opacity (22%) when intact and higher opacity (55%) when failed. A legend is shown in the top-left of the 3D view. Hover text now includes material name, panel dimensions, thickness, and DI.
24. **Interior-wall sheltered-side check** — each interior partition is assigned a single outward normal (+x or +y). The sheltered-side heuristic (10% load when the blast is behind the normal direction) is therefore orientation-dependent; the same interior wall receives full direct load from one room and 10% from the other. Interior loads are predominantly driven by the interior-pressure propagation model rather than direct SDOF integration, so the impact on final room pressures is small, but the SDOF damage index for individual interior panels can be asymmetric.
23. **Exits & rescue routes** — `core/rescue.py` implements Dijkstra pathfinding through a room-adjacency graph. Edge costs: intact wall = 1.0, damaged (DI ≥ 0.5) = 2.0, failed wall = 5.0. Exits are placed by the user on the floor plan (in Exit mode). After simulation the Rescue Plan tab shows priority-ordered routes (1 Immediate → 5 Expectant) and dotted route lines are drawn on the floor plan. Room data (`x_min`, `x_max`, `y_min`, `y_max`, `floor_idx`, `z_bot`, `z_top`) is serialised into `sim-store` so routes can be re-computed whenever exits or people change without re-running the simulation. Panel `room_inside` / `room_outside` IDs are also stored so the graph builder can identify which rooms each interior wall connects.

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
