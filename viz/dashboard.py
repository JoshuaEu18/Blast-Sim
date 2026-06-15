"""
Blast Damage Simulator – Plotly Dash UI.
Run via:  python main.py  →  http://127.0.0.1:8050
"""

import base64
import json
import traceback
from types import SimpleNamespace

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import dash
from dash import dcc, html, Input, Output, State, ALL, ctx, no_update
import dash_bootstrap_components as dbc
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.materials  import MATERIALS
from core.geometry   import create_building, create_building_from_blueprint
from core.blast      import BlastSource, TNT_EQUIVALENCY
from core.simulation import run_simulation
from core.persons    import Person, assess_person_injuries
from core.rescue     import Exit, RescueRoute, compute_rescue_routes, RESCUE_PRIORITY, PRIORITY_LABEL

# ── FloorSpaceJS → blast-sim blueprint converter ──────────────────────────

def parse_floorspacejs(data: dict) -> dict:
    """
    Convert an OpenStudio FloorSpaceJS JSON (from OpenStudio's floorspace.js editor)
    into the blast-sim blueprint dict format.

    FloorSpaceJS coordinates are in real-world units (feet by convention).
    This function converts them to metres.
    Each space polygon is approximated by its axis-aligned bounding box.
    """
    FT_TO_M = 0.3048
    stories = data.get('stories', [])

    floors = []
    for story in stories:
        fh_raw = story.get('floor_to_ceiling_height') or 10.0
        fh_m = round(float(fh_raw) * FT_TO_M, 3)

        geom       = story.get('geometry', {})
        vertices_d = {v['id']: v for v in geom.get('vertices', [])}
        edges_d    = {e['id']: e for e in geom.get('edges', [])}
        faces_d    = {f['id']: f for f in geom.get('faces', [])}

        rooms = []
        for space in story.get('spaces', []):
            face = faces_d.get(space.get('face_id', ''))
            if not face:
                continue
            poly_x, poly_y = [], []
            for eid, order in zip(face['edge_ids'], face['edge_order']):
                e = edges_d.get(eid)
                if not e:
                    continue
                v = vertices_d.get(e['vertex_ids'][int(order)])
                if v:
                    poly_x.append(float(v['x']))
                    poly_y.append(float(v['y']))
            if not poly_x:
                continue
            rooms.append({
                'x_min': round(min(poly_x) * FT_TO_M, 3),
                'x_max': round(max(poly_x) * FT_TO_M, 3),
                'y_min': round(min(poly_y) * FT_TO_M, 3),
                'y_max': round(max(poly_y) * FT_TO_M, 3),
                'name':  space.get('name', f'Space_{space["id"]}'),
            })

        if rooms:
            floors.append({'rooms': rooms, '_floor_height_m': fh_m})

    if not floors:
        raise ValueError('No stories with rooms found in FloorSpaceJS data.')

    floor_heights = [f['_floor_height_m'] for f in floors]
    default_fh    = floor_heights[0]

    return {
        'building': {
            'n_floors':      len(floors),
            'floor_height':  round(default_fh, 3),
            'floor_heights': floor_heights,
            'wall_material': 'concrete',
            'window_frac':   0.25,
            'total_occupants': 0,
        },
        'floors': [{'rooms': f['rooms']} for f in floors],
        '_source': 'floorspacejs',
    }


def _auto_parse_blueprint(raw: dict) -> dict:
    """Detect format and return a blast-sim blueprint dict."""
    if 'stories' in raw:
        return parse_floorspacejs(raw)
    # Already in blast-sim format (has 'building' + 'floors')
    return raw


# ── Shared theme ───────────────────────────────────────────────────────────
_BG     = '#0d1117'
_BG2    = '#161b22'
_GRID   = '#21262d'
_TEXT   = '#e6edf3'
_MUTED  = '#8b949e'
_RED    = '#f85149'
_GREEN  = '#3fb950'
_ORANGE = '#d29922'
_BLUE   = '#388bfd'
_TEAL   = '#39d353'

_SEV_COLOR = {
    'Uninjured': _GREEN,
    'Minor':     '#d29922',
    'Moderate':  '#f0883e',
    'Severe':    _RED,
    'Fatal':     '#9b1c1c',
}
_SEV_BG = {
    'Uninjured': '#071a07',
    'Minor':     '#1a140a',
    'Moderate':  '#1a0e07',
    'Severe':    '#1a0707',
    'Fatal':     '#120404',
}

# Rescue-priority colours (1=highest urgency)
_PRI_COLOR = {
    1: _RED,
    2: '#f0883e',
    3: _ORANGE,
    4: _GREEN,
    5: _MUTED,
}


def _layout(**kwargs):
    base = dict(
        paper_bgcolor=_BG, plot_bgcolor=_BG2,
        font=dict(color=_TEXT, family='ui-monospace,monospace', size=11),
    )
    return {**base, **kwargs}


# ── 3-D building figure ────────────────────────────────────────────────────

def _damage_color(di, failed, panel_type='ext_wall'):
    """Return an RGB colour string that encodes both element type and damage state."""
    if failed or di >= 1.0:
        frac = min((di - 1.0) / 3.0, 1.0) if di > 1.0 else 0.0
        return f'rgb({int(200 + frac*55)},25,25)'

    if panel_type == 'window':
        r = int(60  + di * 160); g = int(160 - di * 80); b = int(230 - di * 110)
        return f'rgb({r},{g},{b})'
    elif panel_type in ('floor', 'roof'):
        r = int(70  + di * 150); g = int(110 + di * 40); b = int(160 - di * 140)
        return f'rgb({r},{g},{b})'
    elif panel_type == 'int_wall':
        r = int(40  + di * 180); g = int(160 - di * 60); b = int(140 - di * 120)
        return f'rgb({r},{g},{b})'
    else:  # ext_wall
        return f'rgb({int(di*220)},{int((1-di*0.65)*185)},40)'


_PANEL_TYPE_LABEL = {
    'ext_wall': 'Exterior Wall',
    'int_wall': 'Interior Wall',
    'window':   'Window (Glass)',
    'floor':    'Floor Slab',
    'roof':     'Roof Slab',
}

_LEGEND_ENTRIES = [
    ('Ext. Wall — intact',    'rgb(0,185,40)'),
    ('Int. Wall — intact',    'rgb(40,160,140)'),
    ('Window — intact',       'rgb(60,160,230)'),
    ('Floor / Roof — intact', 'rgb(70,110,160)'),
    ('Any — yielded (DI ≥ 1)', _ORANGE),
    ('Any — failed',          'rgb(200,25,25)'),
]


def build_3d_figure(panels, columns, result=None, people_injuries=None, blast_pos=None):
    """Build the 3D structural model figure."""
    traces       = []
    trace_groups = []

    # ── Panels ────────────────────────────────────────────────────────────────
    for panel in panels:
        di    = panel.damage_index if result else 0.0
        fail  = panel.failed       if result else False
        ptype = panel.panel_type
        color = _damage_color(di, fail, ptype)

        if ptype == 'window':
            opacity = 0.55 if (fail or di >= 1.0) else 0.22
        elif ptype in ('floor', 'roof'):
            opacity = 0.45
        elif ptype == 'int_wall':
            opacity = 0.55
        else:
            opacity = 0.82

        c      = panel.corners
        tlabel = _PANEL_TYPE_LABEL.get(ptype, ptype)
        status = ('<b style="color:#f85149">FAILED</b>' if fail
                  else '<span style="color:#3fb950">Intact</span>')
        lbl = (f"<b>{tlabel}</b>  ·  Floor {panel.floor_idx + 1}<br>"
               f"Material : {panel.material.name}<br>"
               f"Size     : {panel.width:.1f} × {panel.height:.1f} m  "
               f"| t = {panel.thickness * 1000:.0f} mm<br>"
               f"DI = {di:.2f}  {status}")

        traces.append(go.Mesh3d(
            x=c[:,0], y=c[:,1], z=c[:,2],
            i=[0,0], j=[1,2], k=[2,3],
            color=color, opacity=opacity, flatshading=True,
            hovertext=lbl, hoverinfo='text', showscale=False,
            showlegend=False,
            lighting=dict(ambient=0.7, diffuse=0.6, specular=0.1),
        ))
        trace_groups.append(ptype)

    # ── Legend dummy traces ───────────────────────────────────────────────────
    for leg_name, leg_col in _LEGEND_ENTRIES:
        traces.append(go.Scatter3d(
            x=[None], y=[None], z=[None],
            mode='markers',
            marker=dict(size=8, color=leg_col),
            name=leg_name,
            showlegend=True,
        ))
        trace_groups.append('legend')

    # ── Columns ───────────────────────────────────────────────────────────────
    for col in columns:
        ccolor = _RED if (result and col.failed) else '#3d4f6b'
        traces.append(go.Scatter3d(
            x=[col.base[0], col.top[0]], y=[col.base[1], col.top[1]],
            z=[col.base[2], col.top[2]],
            mode='lines', line=dict(color=ccolor, width=5),
            hoverinfo='skip', showlegend=False,
        ))
        trace_groups.append('column')

    # ── People markers ────────────────────────────────────────────────────────
    if people_injuries:
        for inj in people_injuries:
            px, py, pz = inj.position
            sev_col = _SEV_COLOR.get(inj.overall_severity, _MUTED)
            lbl = (f"<b>{inj.name}</b><br>"
                   f"Severity: <b>{inj.overall_severity}</b><br>"
                   f"Pressure: {inj.effective_pressure_kPa:.0f} kPa "
                   f"(shielded {int((1-inj.shielding_factor)*100)}%)<br>"
                   f"Standoff: {inj.standoff_m:.1f} m")
            traces.append(go.Scatter3d(
                x=[px], y=[py], z=[pz], mode='markers+text',
                text=[inj.name], textposition='top center',
                textfont=dict(size=9, color=sev_col),
                marker=dict(size=8, color=sev_col, symbol='circle',
                            line=dict(color='#fff', width=1)),
                hovertext=lbl, hoverinfo='text', showlegend=False,
            ))
            trace_groups.append('person')

    # ── Blast source marker ───────────────────────────────────────────────────
    if blast_pos is not None:
        bx3, by3, bz3 = float(blast_pos[0]), float(blast_pos[1]), float(blast_pos[2])
        lbl3 = f'<b>💥 Blast source</b><br>({bx3:.1f}, {by3:.1f}, {bz3:.1f}) m'
        traces.append(go.Scatter3d(
            x=[bx3], y=[by3], z=[bz3],
            mode='markers+text',
            text=['💥'], textposition='top center',
            textfont=dict(size=11, color=_RED),
            marker=dict(size=10, color=_RED, symbol='diamond',
                        line=dict(color='#ff8080', width=2)),
            hovertext=lbl3, hoverinfo='text', showlegend=False,
        ))
        trace_groups.append('blast')

    # ── Invisible 3D click grid (one layer per floor at standing height) ─────────
    if panels:
        xs_all = [float(c[0]) for p in panels for c in p.corners]
        ys_all = [float(c[1]) for p in panels for c in p.corners]
        x_lo, x_hi = min(xs_all), max(xs_all)
        y_lo, y_hi = min(ys_all), max(ys_all)

        floor_z_bots: dict = {}
        for p in panels:
            z_b = float(min(float(c[2]) for c in p.corners))
            floor_z_bots[p.floor_idx] = min(floor_z_bots.get(p.floor_idx, 1e9), z_b)

        step = 1.0
        gx = np.arange(x_lo + 0.5, x_hi, step)
        gy = np.arange(y_lo + 0.5, y_hi, step)
        if len(gx) and len(gy):
            xx, yy = np.meshgrid(gx, gy)
            flat_x = xx.flatten(); flat_y = yy.flatten()
            for fi, z_b in sorted(floor_z_bots.items()):
                z_s = z_b + 1.2   # standing height ≈ 1.2 m
                traces.append(go.Scatter3d(
                    x=flat_x, y=flat_y, z=np.full(len(flat_x), z_s),
                    mode='markers',
                    marker=dict(size=8, color='rgba(0,0,0,0)', opacity=0),
                    hovertemplate=(
                        f'Floor {fi + 1}  ·  (%{{x:.1f}} m, %{{y:.1f}} m)<br>'
                        '<i>Click to place person / exit</i><extra></extra>'
                    ),
                    customdata=np.column_stack([
                        np.full(len(flat_x), fi),
                        np.full(len(flat_x), z_s),
                    ]),
                    showlegend=False,
                ))
                trace_groups.append('click_grid')

    # ── Build visibility arrays for the filter dropdown ───────────────────────
    _ALWAYS = {'legend', 'column', 'person', 'blast', 'click_grid'}

    def _vis(show_panel_types):
        return [True if (g in _ALWAYS or g in show_panel_types) else False
                for g in trace_groups]

    _FILTER_BUTTONS = [
        ('Exterior walls',  {'ext_wall'}),
        ('Interior walls',  {'int_wall'}),
        ('Windows',         {'window'}),
        ('Floors & roof',   {'floor', 'roof'}),
        ('Walls & windows', {'ext_wall', 'int_wall', 'window'}),
    ]

    dropdown_buttons = [
        dict(label=lbl, method='restyle',
             args=[{'visible': _vis(groups)}])
        for lbl, groups in _FILTER_BUTTONS
    ]

    fig = go.Figure(data=traces)
    fig.update_layout(
        **_layout(margin=dict(l=0, r=0, t=0, b=0)),
        scene=dict(
            xaxis=dict(title='X (m)', gridcolor=_GRID, zerolinecolor=_GRID,
                       backgroundcolor=_BG, color=_MUTED),
            yaxis=dict(title='Y (m)', gridcolor=_GRID, zerolinecolor=_GRID,
                       backgroundcolor=_BG, color=_MUTED),
            zaxis=dict(title='Z (m)', gridcolor=_GRID, zerolinecolor=_GRID,
                       backgroundcolor=_BG, color=_MUTED),
            aspectmode='data', bgcolor=_BG,
        ),
        showlegend=True,
        legend=dict(
            x=0.01, y=0.99,
            xanchor='left', yanchor='top',
            bgcolor='rgba(22,27,34,0.85)',
            bordercolor=_GRID, borderwidth=1,
            font=dict(size=10, color=_TEXT),
            itemsizing='constant',
        ),
        updatemenus=[dict(
            type='dropdown',
            buttons=dropdown_buttons,
            direction='down',
            showactive=True,
            x=0.01, y=0.78,
            xanchor='left', yanchor='top',
            bgcolor='#161b22', bordercolor='#30363d',
            font=dict(color=_TEXT, size=10),
            pad=dict(t=2, b=2, l=4, r=4),
        )],
    )
    return fig


# ── Floor plan (top-down 2D view) ─────────────────────────────────────────

def build_floor_plan(width, depth, n_rooms_x, n_rooms_y, floor_idx, floor_height,
                     blast_pos, people, people_injuries=None,
                     exits=None, place_mode='person',
                     rescue_routes=None, rooms_data=None,
                     blueprint_rooms=None):
    """
    Top-down 2D floor plan.

    Parameters:
      exits          : list of exit dicts from exits-store
      place_mode     : 'person' | 'exit' — changes click hint text
      rescue_routes  : list of route dicts (from rescue-routes-store) for path lines
      rooms_data     : list of room dicts from sim-store, used for route waypoints
      blueprint_rooms: list of room dicts from blueprint JSON for the current floor;
                       when provided, draws actual room shapes instead of a uniform grid
    """
    shapes = []
    if blueprint_rooms:
        for room in blueprint_rooms:
            shapes.append(dict(
                type='rect',
                x0=float(room['x_min']), x1=float(room['x_max']),
                y0=float(room['y_min']), y1=float(room['y_max']),
                line=dict(color='#388bfd', width=1.5),
                fillcolor='#161b22',
                layer='below',
            ))
    else:
        xs = np.linspace(0, width, n_rooms_x + 1)
        ys = np.linspace(0, depth, n_rooms_y + 1)
        for xi in range(n_rooms_x):
            for yi in range(n_rooms_y):
                shapes.append(dict(
                    type='rect',
                    x0=float(xs[xi]), x1=float(xs[xi+1]),
                    y0=float(ys[yi]), y1=float(ys[yi+1]),
                    line=dict(color='#3d444d', width=1),
                    fillcolor='#161b22',
                    layer='below',
                ))

    traces = []

    # Dense invisible marker grid — gives accurate click coordinates.
    step = 0.5
    gx = np.arange(step / 2, width,  step)
    gy = np.arange(step / 2, depth, step)
    xx, yy = np.meshgrid(gx, gy)
    click_hint = 'person' if place_mode == 'person' else 'exit'
    traces.append(go.Scatter(
        x=xx.flatten().tolist(),
        y=yy.flatten().tolist(),
        mode='markers',
        marker=dict(size=14, color='rgba(0,0,0,0)', opacity=0),
        hovertemplate=f'(%{{x:.1f}} m, %{{y:.1f}} m)<br><i>Click to place {click_hint}</i><extra></extra>',
        showlegend=False,
    ))

    floor_idx_int = int(floor_idx)
    z_bot = floor_idx_int * float(floor_height)
    z_top = z_bot + float(floor_height)

    # ── Rescue route lines ────────────────────────────────────────────────────
    if rescue_routes and rooms_data:
        room_map  = {r['id']: r for r in rooms_data}
        exit_map  = {e['id']: e for e in (exits or [])}

        for route in rescue_routes:
            if not route.get('reachable', True):
                continue
            path = route.get('path_room_ids', [])
            if not path:
                continue

            # Only draw for routes whose first room is on the current floor
            first_room = room_map.get(path[0])
            if not first_room or first_room.get('floor_idx') != floor_idx_int:
                continue

            wpts_x, wpts_y = [], []

            # Person start position
            pid = route.get('person_id')
            inj_map = {inj.person_id: inj for inj in (people_injuries or [])}
            inj = inj_map.get(pid)
            if inj and z_bot <= inj.position[2] < z_top:
                wpts_x.append(inj.position[0])
                wpts_y.append(inj.position[1])

            # Room centroids along path (current floor only)
            for rid in path:
                r = room_map.get(rid)
                if r and r.get('floor_idx') == floor_idx_int:
                    wpts_x.append((r['x_min'] + r['x_max']) / 2.0)
                    wpts_y.append((r['y_min'] + r['y_max']) / 2.0)

            # Exit end position
            eid = route.get('exit_id', -1)
            ex_obj = exit_map.get(eid)
            if ex_obj and ex_obj.get('floor_idx') == floor_idx_int:
                wpts_x.append(ex_obj['x'])
                wpts_y.append(ex_obj['y'])

            if len(wpts_x) >= 2:
                sev_col = _SEV_COLOR.get(route.get('person_severity', 'Uninjured'), _MUTED)
                traces.append(go.Scatter(
                    x=wpts_x, y=wpts_y,
                    mode='lines+markers',
                    line=dict(color=sev_col, width=2, dash='dot'),
                    marker=dict(size=4, color=sev_col, opacity=0.7),
                    hoverinfo='skip',
                    showlegend=False,
                    opacity=0.65,
                ))

    # ── People on this floor ──────────────────────────────────────────────────
    inj_map = {inj.person_id: inj for inj in (people_injuries or [])}
    fp_x, fp_y, fp_txt, fp_col = [], [], [], []
    fp_people_on_floor = []
    for p in (people or []):
        if z_bot <= float(p['z']) < z_top:
            inj = inj_map.get(p['id'])
            col = _SEV_COLOR.get(inj.overall_severity, _GREEN) if inj else _GREEN
            fp_x.append(p['x']); fp_y.append(p['y'])
            fp_txt.append(p['name']); fp_col.append(col)
            fp_people_on_floor.append(p)

    if fp_x:
        hover = [f'<b>{n}</b><br>X={x:.1f} m  Y={y:.1f} m  Z={float(p["z"]):.1f} m'
                 for n, x, y, p in zip(fp_txt, fp_x, fp_y, fp_people_on_floor)]
        traces.append(go.Scatter(
            x=fp_x, y=fp_y, mode='markers+text',
            text=fp_txt, textposition='top center',
            textfont=dict(size=9, color=_TEXT, family='ui-monospace,monospace'),
            marker=dict(size=12, color=fp_col, symbol='circle',
                        line=dict(color='#ffffff', width=1.5)),
            hovertext=hover, hoverinfo='text',
            showlegend=False,
        ))

    # ── Exits on this floor ───────────────────────────────────────────────────
    ex_x, ex_y, ex_txt = [], [], []
    for e in (exits or []):
        if e.get('floor_idx', 0) == floor_idx_int:
            ex_x.append(e['x']); ex_y.append(e['y'])
            ex_txt.append(e['name'])
    if ex_x:
        hover_ex = [f'<b>🚪 {n}</b><br>X={x:.1f} m  Y={y:.1f} m  Floor {floor_idx_int+1}'
                    for n, x, y in zip(ex_txt, ex_x, ex_y)]
        traces.append(go.Scatter(
            x=ex_x, y=ex_y, mode='markers+text',
            text=ex_txt, textposition='top center',
            textfont=dict(size=9, color=_TEAL, family='ui-monospace,monospace'),
            marker=dict(size=13, color=_TEAL, symbol='square',
                        line=dict(color='#ffffff', width=1.5)),
            hovertext=hover_ex, hoverinfo='text',
            showlegend=False,
        ))

    # ── Blast source ──────────────────────────────────────────────────────────
    bx, by = float(blast_pos[0]), float(blast_pos[1])
    traces.append(go.Scatter(
        x=[bx], y=[by], mode='markers+text',
        text=['💥 blast'], textposition='top center',
        textfont=dict(size=8, color=_RED),
        marker=dict(size=14, color=_RED, symbol='star',
                    line=dict(color='#ff8080', width=1.5)),
        hovertemplate=f'Blast source<br>({bx:.1f} m, {by:.1f} m)<extra></extra>',
        showlegend=False,
    ))

    if blueprint_rooms:
        all_x = [float(r['x_min']) for r in blueprint_rooms] + [float(r['x_max']) for r in blueprint_rooms]
        all_y = [float(r['y_min']) for r in blueprint_rooms] + [float(r['y_max']) for r in blueprint_rooms]
        x_range = [min(all_x) - 0.5, max(all_x) + 0.5]
        y_range = [min(all_y) - 0.5, max(all_y) + 0.5]
    else:
        x_range = [-0.5, width + 0.5]
        y_range = [-0.5, depth + 0.5]

    fig = go.Figure(data=traces)
    fig.update_layout(
        **_layout(margin=dict(l=4, r=4, t=4, b=4)),
        shapes=shapes,
        xaxis=dict(range=x_range, showgrid=False, color=_MUTED, zeroline=False),
        yaxis=dict(range=y_range, showgrid=False, color=_MUTED, zeroline=False),
        showlegend=False,
        clickmode='event+select',
        dragmode='pan',
        height=200,
    )
    return fig


# ── Results charts ─────────────────────────────────────────────────────────

def build_results_figure(rooms, result):
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=('Room Peak Overpressure (kPa)', 'Casualties by Room',
                        'Panel Damage Index', 'Story Drift Ratio (%)'),
        vertical_spacing=0.20, horizontal_spacing=0.10,
    )
    room_ids  = [f'R{r.id}' for r in rooms]
    room_pres = [result.room_peak_pressure.get(r.id, 0.0) for r in rooms]

    fig.add_trace(go.Bar(
        x=room_ids, y=room_pres, name='kPa',
        marker=dict(color=[_RED if p>35 else _GREEN for p in room_pres], line=dict(width=0)),
        hovertemplate='%{x}: %{y:.1f} kPa<extra></extra>',
    ), row=1, col=1)

    for key, name, color in [('fatal','Fatal',_RED),('severe','Severe',_ORANGE),('minor','Minor',_BLUE)]:
        vals = [result.casualties.get(r.id,{}).get(key, 0.0) for r in rooms]
        fig.add_trace(go.Bar(x=room_ids, y=vals, name=name.title(),
                             marker=dict(color=color, line=dict(width=0))), row=1, col=2)

    dis = [r.damage_index for r in result.panel_results.values()]
    fig.add_trace(go.Histogram(x=dis, nbinsx=22,
                               marker=dict(color=_BLUE, line=dict(width=0))), row=2, col=1)
    fig.add_vline(x=1.0, line=dict(dash='dot', color=_ORANGE, width=1.5),
                  annotation=dict(text='Yield', font=dict(color=_ORANGE, size=9)), row=2, col=1)

    drift_pct = result.building.story_drift * 100.0
    fig.add_trace(go.Bar(
        x=[f'F{i+1}' for i in range(len(drift_pct))], y=drift_pct,
        marker=dict(color=[_RED if d>2.5 else (_ORANGE if d>1.0 else _GREEN) for d in drift_pct],
                    line=dict(width=0)),
        hovertemplate='%{x}: %{y:.2f}%<extra></extra>',
    ), row=2, col=2)
    for y, lbl, col in [(2.5,'Life Safety',_ORANGE),(5.0,'Collapse',_RED)]:
        fig.add_hline(y=y, line=dict(dash='dot', color=col, width=1.5),
                      annotation=dict(text=lbl, font=dict(color=col, size=9)), row=2, col=2)

    fig.update_layout(**_layout(margin=dict(l=10,r=10,t=38,b=8)),
                      barmode='stack', showlegend=True,
                      legend=dict(bgcolor=_BG2, bordercolor=_GRID, borderwidth=1,
                                  font=dict(size=10), x=1.01))
    for ann in fig.layout.annotations:
        ann.font = dict(color=_MUTED, size=11)
    fig.update_xaxes(showgrid=False, color=_MUTED, tickfont=dict(size=10))
    fig.update_yaxes(gridcolor=_GRID, color=_MUTED, tickfont=dict(size=10))
    return fig


def _empty_fig(msg=''):
    fig = go.Figure()
    fig.update_layout(**_layout(margin=dict(l=8,r=8,t=8,b=8)),
                      annotations=[dict(text=msg, x=0.5, y=0.5, xref='paper',
                                        yref='paper', showarrow=False,
                                        font=dict(color=_MUTED, size=12))])
    return fig


# ── People injury cards ────────────────────────────────────────────────────

def _mini_stat(label, value):
    return html.Div([
        html.Span(label + ':', style={'color': _MUTED, 'fontSize': '10px', 'marginRight': '4px'}),
        html.Span(value, style={'color': _TEXT, 'fontSize': '10px', 'fontFamily': 'monospace'}),
    ], style={'marginBottom': '2px'})


def render_injury_card(inj):
    sev   = inj.overall_severity
    col   = _SEV_COLOR.get(sev, _MUTED)
    bg    = _SEV_BG.get(sev, _BG2)
    _rank = {'None':0,'Low Risk':0,'Low':0,'Mild':1,'Minor':1,'Moderate':2,'Severe':3,'Fatal':4}

    rows = []
    for detail in inj.injuries:
        if detail.probability < 0.03 and detail.severity in ('None', 'Low Risk', 'Low'):
            continue
        d_col = (_rank.get(detail.severity, 0) and
                 [_MUTED, _ORANGE, '#f0883e', _RED, '#9b1c1c'][min(_rank.get(detail.severity,0),4)])
        d_col = d_col or _MUTED
        rows.append(html.Div([
            html.Div([
                html.Span(detail.name, style={'fontWeight':'600','fontSize':'11px'}),
                html.Span(f' {detail.probability*100:.0f}%',
                          style={'fontSize':'10px','color':_MUTED}),
                html.Span(f'  {detail.severity}',
                          style={'fontSize':'10px','color':d_col,'fontWeight':'600'}),
            ]),
            html.Div(detail.mechanism,
                     style={'fontSize':'9px','color':_MUTED,'marginLeft':'8px',
                            'marginBottom':'4px','lineHeight':'1.3'}),
        ]))

    shielding_pct = int((1 - inj.shielding_factor) * 100)

    return html.Div([
        html.Div([
            html.Span(f'👤  {inj.name}', style={'fontWeight':'700','fontSize':'12px'}),
            html.Span(sev, style={'float':'right','fontSize':'11px',
                                   'fontWeight':'700','color':col}),
        ], style={'marginBottom':'6px'}),
        html.Div([
            _mini_stat('Distance', f'{inj.standoff_m:.1f} m'),
            _mini_stat('Raw pressure', f'{inj.raw_pressure_kPa:.0f} kPa'),
            _mini_stat('After shielding',
                       f'{inj.effective_pressure_kPa:.0f} kPa  '
                       f'({shielding_pct}% blocked, {inj.walls_between} wall{"s" if inj.walls_between!=1 else ""})'),
            _mini_stat('Throw velocity', f'{inj.throw_velocity_ms:.1f} m/s'),
            _mini_stat('Survival prob.', f'{inj.survival_probability*100:.0f}%'),
        ], style={'marginBottom':'8px'}),
        html.Div(rows or [html.Span('No significant injuries',
                                     style={'fontSize':'11px','color':_GREEN})]),
    ], style={
        'background': bg,
        'border': f'1px solid {col}40',
        'borderLeft': f'3px solid {col}',
        'borderRadius': '8px',
        'padding': '10px 10px 6px',
        'marginBottom': '8px',
    })


# ── People list (sidebar) ──────────────────────────────────────────────────

def render_people_list(people):
    if not people:
        return html.Div('No people placed yet.\nClick the floor plan or press ＋ Add.',
                        style={'fontSize':'10px','color':_MUTED,'padding':'4px 0',
                               'whiteSpace':'pre-line'})
    rows = []
    for p in people:
        rows.append(html.Div([
            html.Div([
                html.Div(p['name'],
                         style={'fontSize':'11px','fontWeight':'700','color':_TEXT,
                                'lineHeight':'1.2'}),
                html.Div(f"X {p['x']:.1f} m  ·  Y {p['y']:.1f} m  ·  Z {p['z']:.1f} m",
                         style={'fontSize':'9px','color':_MUTED,
                                'fontFamily':'ui-monospace,monospace','marginTop':'1px'}),
            ], style={'flex':'1','minWidth':0}),
            html.Button('×',
                id={'type':'person-del','index': p['id']},
                n_clicks=0,
                style={'background':'none','border':'none','color':'#666',
                       'cursor':'pointer','fontSize':'16px','padding':'0 4px',
                       'lineHeight':'1','flexShrink':'0'}),
        ], style={'display':'flex','alignItems':'center','padding':'5px 0',
                  'borderBottom':f'1px solid {_GRID}'}))
    return html.Div(rows)


# ── Exits list (sidebar) ───────────────────────────────────────────────────

def render_exits_list(exits):
    if not exits:
        return html.Div('No exits placed yet.\nSwitch to Exit mode and click the plan.',
                        style={'fontSize':'10px','color':_MUTED,'padding':'4px 0',
                               'whiteSpace':'pre-line'})
    rows = []
    for e in exits:
        rows.append(html.Div([
            html.Div([
                html.Div([
                    html.Span('🚪 ', style={'fontSize':'11px'}),
                    html.Span(e['name'],
                              style={'fontSize':'11px','fontWeight':'700','color':_TEAL}),
                ], style={'lineHeight':'1.2'}),
                html.Div(f"X {e['x']:.1f} m  ·  Y {e['y']:.1f} m  ·  F{e['floor_idx']+1}",
                         style={'fontSize':'9px','color':_MUTED,
                                'fontFamily':'ui-monospace,monospace','marginTop':'1px'}),
            ], style={'flex':'1','minWidth':0}),
            html.Button('×',
                id={'type':'exit-del','index': e['id']},
                n_clicks=0,
                style={'background':'none','border':'none','color':'#666',
                       'cursor':'pointer','fontSize':'16px','padding':'0 4px',
                       'lineHeight':'1','flexShrink':'0'}),
        ], style={'display':'flex','alignItems':'center','padding':'5px 0',
                  'borderBottom':f'1px solid {_GRID}'}))
    return html.Div(rows)


# ── Rescue route card ──────────────────────────────────────────────────────

def render_rescue_card(route: dict):
    """Render one rescue route card from a dict (serialised RescueRoute)."""
    pri     = route.get('rescue_priority', 4)
    sev     = route.get('person_severity', 'Uninjured')
    col     = _PRI_COLOR.get(pri, _MUTED)
    sev_col = _SEV_COLOR.get(sev, _MUTED)
    lbl     = PRIORITY_LABEL.get(pri, str(pri))
    name    = route.get('person_name', '?')
    reachable = route.get('reachable', True)
    exit_name = route.get('exit_name', '?')
    cost      = route.get('path_cost', 0.0)
    desc      = route.get('path_description', '')

    if reachable:
        route_body = html.Div([
            html.Div([
                html.Span('🚪 ', style={'color': _TEAL}),
                html.Span(exit_name,
                          style={'color': _TEAL, 'fontWeight': '600', 'fontSize': '11px'}),
                html.Span(f'  cost {cost:.1f}',
                          style={'color': _MUTED, 'fontSize': '9px', 'marginLeft': '4px'}),
            ], style={'marginBottom': '3px'}),
            html.Div(desc, style={'fontSize': '9px', 'color': _MUTED, 'lineHeight': '1.4'}),
        ])
    else:
        route_body = html.Div([
            html.Span('⚠  ', style={'color': _RED}),
            html.Span(desc, style={'color': _RED, 'fontSize': '10px'}),
        ])

    return html.Div([
        html.Div([
            html.Span(f'👤  {name}',
                      style={'fontWeight': '700', 'fontSize': '12px', 'color': _TEXT}),
            html.Span(f'P{pri} · {lbl}',
                      style={'float': 'right', 'fontSize': '10px',
                             'color': col, 'fontWeight': '700'}),
        ], style={'marginBottom': '4px', 'overflow': 'hidden'}),
        html.Div([
            html.Span(sev, style={'color': sev_col, 'fontSize': '10px', 'fontWeight': '600'}),
        ], style={'marginBottom': '5px'}),
        route_body,
    ], style={
        'background': _BG2,
        'border': f'1px solid {col}40',
        'borderLeft': f'3px solid {col}',
        'borderRadius': '8px',
        'padding': '8px 10px',
        'marginBottom': '6px',
    })


# ── Sidebar helpers ────────────────────────────────────────────────────────

def _numbox(id_, label, min_, max_, step, val):
    return html.Div([
        html.Div(label, style={'color': _MUTED, 'fontSize': '10px', 'marginBottom': '2px'}),
        dcc.Input(
            id=id_, type='number', value=val, min=min_, max=max_, step=step,
            debounce=True,
            style={
                'width': '100%', 'background': '#21262d',
                'border': '1px solid #30363d', 'borderRadius': '5px',
                'color': _TEXT, 'padding': '3px 6px',
                'fontSize': '11px', 'outline': 'none', 'boxSizing': 'border-box',
            }
        ),
    ], style={'marginBottom': '5px'})


def _dropdown(id_, options, val):
    return dcc.Dropdown(options=options, value=val, id=id_, clearable=False,
                        style={'backgroundColor':'#21262d','color':'#e6edf3',
                               'border':'1px solid #30363d','borderRadius':'6px',
                               'fontSize':'12px'})


def _tile(value, label, col):
    return html.Div([
        html.Div(value, className='metric-val', style={'color': col}),
        html.Div(label, className='metric-lbl'),
    ], className='metric-tile')


# ── App factory ────────────────────────────────────────────────────────────

def create_app():
    app = dash.Dash(__name__, external_stylesheets=[dbc.themes.CYBORG],
                    title='Blast Damage Simulator',
                    suppress_callback_exceptions=True)

    # ── Place-mode toggle ────────────────────────────────────────────────────
    place_mode_row = html.Div([
        html.Span('Click places:',
                  style={'color': _MUTED, 'fontSize': '10px', 'marginRight': '8px',
                         'lineHeight': '22px'}),
        dcc.RadioItems(
            id='place-mode',
            options=[
                {'label': '👤 Person', 'value': 'person'},
                {'label': '🚪 Exit',   'value': 'exit'},
            ],
            value='person',
            inline=True,
            inputStyle={'marginRight': '3px', 'cursor': 'pointer'},
            labelStyle={'color': _TEXT, 'fontSize': '11px',
                        'marginRight': '12px', 'cursor': 'pointer'},
        ),
    ], style={'display': 'flex', 'alignItems': 'center', 'marginBottom': '6px',
              'padding': '4px 6px', 'background': '#21262d',
              'borderRadius': '6px', 'border': f'1px solid {_GRID}'})

    sidebar = html.Div([
      html.Div([
        html.Button('▶  Run Simulation', id='btn-run', n_clicks=0),
        html.Div(id='run-status', className='status-idle',
                 style={'textAlign':'center','marginTop':'6px'}),
      ], style={'padding':'10px 12px','borderBottom':'1px solid #21262d',
                'background':'#161b22'}),
      html.Div([
        html.Div('Building', className='section-label accent-blue'),
        dbc.Row([
            dbc.Col(_numbox('bld-width',  'Width (m)',    5,  30, 1,   12), width=6),
            dbc.Col(_numbox('bld-depth',  'Depth (m)',    5,  20, 1,   10), width=6),
        ], className='g-1'),
        dbc.Row([
            dbc.Col(_numbox('bld-floors', 'Floors',       1,   6, 1,    3), width=6),
            dbc.Col(_numbox('bld-fh',     'Floor ht (m)', 2.5, 5, 0.5, 3.0), width=6),
        ], className='g-1'),
        html.Div('Wall material',
                 style={'color':_MUTED,'fontSize':'11px','marginBottom':'3px','marginTop':'2px'}),
        _dropdown('bld-mat',
                  [{'label':'Reinforced Concrete','value':'concrete'},
                   {'label':'Brick Masonry',       'value':'brick'},
                   {'label':'Steel Plate',         'value':'steel_plate'}],
                  'concrete'),
        html.Div(style={'height':'2px'}),
        dbc.Row([
            dbc.Col(_numbox('bld-rx',  'Rooms X',       1, 4, 1, 2), width=6),
            dbc.Col(_numbox('bld-ry',  'Rooms Y',       1, 4, 1, 2), width=6),
        ], className='g-1'),
        _numbox('bld-wf', 'Window frac.', 0, 0.7, 0.05, 0.30),

        # ── Blueprint upload ───────────────────────────────────────────────────
        html.Div('Blueprint JSON', className='section-label accent-blue'),
        dcc.Upload(
            id='blueprint-upload',
            children=html.Div([
                html.Span('Drag & drop or '),
                html.A('select a JSON file', style={'color': _BLUE, 'cursor': 'pointer'}),
            ], style={'fontSize': '11px', 'color': _MUTED, 'textAlign': 'center',
                      'padding': '8px 4px'}),
            style={
                'border': f'1px dashed {_GRID}',
                'borderRadius': '6px',
                'background': '#0d1117',
                'marginBottom': '4px',
                'cursor': 'pointer',
            },
            multiple=False,
        ),
        html.Div(id='blueprint-status',
                 style={'fontSize': '10px', 'color': _MUTED, 'marginBottom': '4px',
                        'minHeight': '14px'}),
        html.Button('Clear Blueprint', id='btn-clear-blueprint', n_clicks=0,
                    style={'background': 'none', 'border': f'1px solid {_GRID}',
                           'color': _MUTED, 'borderRadius': '5px', 'padding': '3px 8px',
                           'fontSize': '10px', 'cursor': 'pointer', 'marginBottom': '6px'}),

        html.Div('Blast Source', className='section-label accent-red'),
        html.Div('Explosive type',
                 style={'color':_MUTED,'fontSize':'11px','marginBottom':'3px'}),
        _dropdown('blast-explosive',
                  [{'label': f'{v[0]}  (×{v[1]:.2f} TNT)', 'value': k}
                   for k, v in TNT_EQUIVALENCY.items()],
                  'tnt'),
        html.Div(id='explosive-label',
                 style={'fontSize':'9px','color':_MUTED,'marginBottom':'4px',
                        'minHeight':'12px'}),
        _numbox('blast-mass', 'Charge mass (kg)', 1, 500, 1, 25),
        dbc.Row([
            dbc.Col(_numbox('blast-x', 'X (m)', -10, 30, 0.5,  6.0), width=6),
            dbc.Col(_numbox('blast-y', 'Y (m)', -10, 20, 0.5, -3.0), width=6),
        ], className='g-1'),
        dbc.Row([
            dbc.Col(_numbox('blast-z', 'Z (m)', 0, 5, 0.25, 0.5), width=6),
            dbc.Col([
                html.Div('Burst type',
                         style={'color':_MUTED,'fontSize':'11px','marginBottom':'3px'}),
                _dropdown('blast-type',
                          [{'label':'Surface — gas leak','value':'surface'},
                           {'label':'Free-air',           'value':'free_air'}],
                          'surface'),
            ], width=6),
        ], className='g-1'),

        # ── People ────────────────────────────────────────────────────────────
        html.Div('People', className='section-label'),
        dbc.Row([
            dbc.Col([
                html.Div('Floor', style={'fontSize':'10px','color':_MUTED,'marginBottom':'2px'}),
                dcc.Dropdown(id='floor-selector',
                             options=[{'label':f'Floor {i+1}','value':i} for i in range(6)],
                             value=0, clearable=False,
                             style={'backgroundColor':'#21262d','fontSize':'12px',
                                    'border':'1px solid #30363d','borderRadius':'6px'}),
            ], width=5),
            dbc.Col(html.Div([
                html.Button('＋ Add', id='btn-add-person', n_clicks=0,
                            style={'background':'#21262d','border':'1px solid #30363d',
                                   'color':_TEXT,'borderRadius':'5px','padding':'4px 8px',
                                   'fontSize':'11px','cursor':'pointer','marginRight':'4px'}),
                html.Button('Clear', id='btn-clear-people', n_clicks=0,
                            style={'background':'none','border':'1px solid #f8514940',
                                   'color':_RED,'borderRadius':'5px','padding':'4px 8px',
                                   'fontSize':'11px','cursor':'pointer'}),
            ], style={'paddingTop':'16px'}), width=7),
        ], className='g-1', style={'marginBottom':'6px'}),

        # ── Place mode + floor plan ───────────────────────────────────────────
        place_mode_row,
        dcc.Graph(id='floor-plan', config={'displayModeBar': False, 'scrollZoom': True},
                  style={'borderRadius':'8px','overflow':'hidden',
                         'border':f'1px solid {_GRID}'}),
        html.Div('Click plan · 👤 person  🚪 exit  💥 blast',
                 style={'fontSize':'9px','color':_MUTED,'textAlign':'center','marginTop':'3px'}),
        html.Div(id='people-list-sidebar', style={'marginTop':'6px'}),

        # ── Exits & Entrances ─────────────────────────────────────────────────
        html.Div('Exits & Entrances', className='section-label'),
        dbc.Row([
            dbc.Col(html.Button('＋ Add Exit', id='btn-add-exit', n_clicks=0,
                                style={'background':'#21262d',
                                       'border':f'1px solid {_TEAL}60',
                                       'color':_TEAL,'borderRadius':'5px',
                                       'padding':'4px 8px','fontSize':'11px',
                                       'cursor':'pointer','width':'100%'}), width=6),
            dbc.Col(html.Button('Clear', id='btn-clear-exits', n_clicks=0,
                                style={'background':'none','border':'1px solid #f8514940',
                                       'color':_RED,'borderRadius':'5px','padding':'4px 8px',
                                       'fontSize':'11px','cursor':'pointer','width':'100%'}),
                    width=6),
        ], className='g-1', style={'marginBottom':'6px'}),
        html.Div(id='exits-list-sidebar', style={'marginTop':'2px'}),

      ], className='sidebar-scroll', style={'padding':'10px 10px'}),
    ], className='sidebar')

    metrics_row = html.Div([
        dbc.Row([
            dbc.Col(_tile('—', 'Fatal',         _RED),    width=3),
            dbc.Col(_tile('—', 'Severe Inj.',   _ORANGE), width=3),
            dbc.Col(_tile('—', 'Minor Inj.',    _BLUE),   width=3),
            dbc.Col(_tile('—', 'Stability',    _MUTED),   width=3),
        ], id='metrics-row', className='g-2'),
    ], style={'padding':'10px 10px 6px'})

    # ── Right panel: Injuries + Rescue Plan tabs ──────────────────────────────
    right_panel = html.Div([
        dbc.Tabs([
            dbc.Tab(
                label='Injuries',
                tab_id='tab-injuries',
                children=html.Div(
                    id='people-injuries-panel',
                    children=html.Div(
                        'Run simulation to see injury assessment.',
                        style={'fontSize':'11px','color':_MUTED,'padding':'10px'},
                    ),
                    style={'overflowY':'auto','maxHeight':'calc(90vh - 95px)','padding':'8px'},
                ),
            ),
            dbc.Tab(
                label='Rescue Plan',
                tab_id='tab-rescue',
                children=html.Div(
                    id='rescue-plan-panel',
                    children=html.Div(
                        'Place exits and run simulation to generate rescue routes.',
                        style={'fontSize':'11px','color':_MUTED,'padding':'10px'},
                    ),
                    style={'overflowY':'auto','maxHeight':'calc(90vh - 95px)','padding':'8px'},
                ),
            ),
        ], id='right-tabs', active_tab='tab-injuries',
           style={'borderBottom': f'1px solid {_GRID}'}),
    ], className='panel-card')

    app.layout = html.Div([
        html.Div([
            html.Span('💥', style={'fontSize':'18px'}),
            html.H4('Gas Leak Blast Damage Simulator'),
            html.Span('SDOF + Shear-Building FEM  ·  Kinney-Graham blast model  ·  Baker (1983) injury model',
                      style={'color':_MUTED,'fontSize':'11px','marginLeft':'8px'}),
            html.Div(id='header-status', style={'marginLeft':'auto','fontSize':'11px','color':_MUTED}),
        ], className='sim-header'),

        html.Div([
            html.Div(sidebar, style={'width':'260px','flexShrink':'0',
                                     'height':'calc(100vh - 45px)'}),
            html.Div([
                metrics_row,
                html.Div([
                    html.Div([
                        html.Div([
                            html.Span('3D Structural Model', className='card-head',
                                      style={'display':'inline'}),
                            html.Span('  ·  click floor to place 👤 / 🚪',
                                      style={'fontSize':'9px','color':_MUTED,
                                             'marginLeft':'8px'}),
                        ], style={'display':'flex','alignItems':'center',
                                  'padding':'4px 10px 2px'}),
                        dcc.Graph(id='graph-3d',
                                  style={'flex':'1','minHeight':'0'},
                                  config={'displayModeBar':True,
                                          'modeBarButtonsToRemove':['toImage']}),
                    ], className='panel-card',
                       style={'display':'flex','flexDirection':'column',
                              'flex':'1','minHeight':'0'}),
                    # Hidden — keeps run_sim callback outputs valid
                    dcc.Graph(id='graph-results', style={'display':'none'},
                              config={'displayModeBar':False}),
                ], style={'padding':'0 8px 8px','display':'flex',
                          'flexDirection':'column','flex':'1','minHeight':'0'}),
            ], style={'flex':'1','minWidth':'0','overflow':'hidden',
                      'height':'calc(100vh - 45px)','display':'flex',
                      'flexDirection':'column'}),
            html.Div(right_panel,
                     style={'width':'245px','flexShrink':'0',
                            'padding':'8px 8px 8px 0',
                            'overflowY':'auto',
                            'height':'calc(100vh - 45px)'}),
        ], style={'display':'flex','height':'calc(100vh - 45px)','overflow':'hidden'}),

        dcc.Store(id='sim-store'),
        dcc.Store(id='people-store',        data=[]),
        dcc.Store(id='exits-store',         data=[]),
        dcc.Store(id='rescue-routes-store', data=[]),
        dcc.Store(id='blueprint-store',     data=None),
    ])

    # ── Callbacks ──────────────────────────────────────────────────────────

    @app.callback(
        Output('floor-selector', 'options'),
        Output('floor-selector', 'value'),
        Input('bld-floors',      'value'),
        Input('blueprint-store', 'data'),
        State('floor-selector',  'value'),
    )
    def sync_floor_selector(n_floors, blueprint, current):
        if blueprint:
            bldg = blueprint.get('building', {})
            n = int(bldg.get('n_floors', n_floors or 3))
        else:
            n = int(n_floors or 3)
        opts = [{'label': f'Floor {i + 1}', 'value': i} for i in range(n)]
        val  = min(int(current or 0), n - 1)
        return opts, val

    # ── Blueprint management ─────────────────────────────────────────────────

    @app.callback(
        Output('blueprint-store',  'data'),
        Output('blueprint-status', 'children'),
        Output('blueprint-status', 'style'),
        Output('bld-floors',       'value'),
        Output('graph-3d',         'figure',   allow_duplicate=True),
        Input('blueprint-upload',     'contents'),
        Input('btn-clear-blueprint',  'n_clicks'),
        State('blueprint-upload',     'filename'),
        State('bld-floors',           'value'),
        prevent_initial_call=True,
    )
    def manage_blueprint(contents, _clear, filename, current_floors):
        tid = ctx.triggered_id
        base_style = {'fontSize': '10px', 'marginBottom': '4px', 'minHeight': '14px'}
        empty_3d   = _empty_fig('Upload a blueprint or run simulation.')

        if tid == 'btn-clear-blueprint':
            return None, 'No blueprint loaded.', {**base_style, 'color': _MUTED}, current_floors, empty_3d

        if not contents:
            return no_update, no_update, no_update, no_update, no_update

        try:
            _, content_string = contents.split(',', 1)
            decoded = base64.b64decode(content_string).decode('utf-8')
            raw = json.loads(decoded)

            # Auto-detect and convert format
            bp = _auto_parse_blueprint(raw)

            n_f         = int(bp['building']['n_floors'])
            total_rooms = sum(len(fd.get('rooms', [])) for fd in bp.get('floors', []))
            source_tag  = '  ·  FloorSpaceJS' if bp.get('_source') == 'floorspacejs' else ''
            status      = f'✓  {filename}  ·  {n_f} floor(s)  ·  {total_rooms} room(s){source_tag}'

            # Build 3D preview (no simulation)
            try:
                panels, columns, _ = create_building_from_blueprint(bp)
                fig3d = build_3d_figure(panels, columns,
                                        result=None, people_injuries=None, blast_pos=None)
            except Exception as e3d:
                fig3d = _empty_fig(f'3D preview error: {e3d}')

            return bp, status, {**base_style, 'color': _GREEN}, n_f, fig3d

        except Exception as e:
            msg = str(e)[:120]
            return None, f'⚠  {msg}', {**base_style, 'color': _RED}, no_update, no_update

    # ── People placement ────────────────────────────────────────────────────

    @app.callback(
        Output('people-store', 'data'),
        Input('floor-plan',       'clickData'),
        Input('graph-3d',         'clickData'),
        Input('btn-add-person',   'n_clicks'),
        Input('btn-clear-people', 'n_clicks'),
        Input({'type':'person-del','index':ALL}, 'n_clicks'),
        State('people-store',  'data'),
        State('floor-selector','value'),
        State('bld-fh',        'value'),
        State('place-mode',    'value'),
        prevent_initial_call=True,
    )
    def manage_people(click_2d, click_3d, _add, _clear, _del_ns,
                      store, floor_idx, fh, place_mode):
        store = store or []
        tid   = ctx.triggered_id

        if tid == 'floor-plan':
            if (place_mode or 'person') != 'person':
                return no_update
            if not click_2d:
                return no_update
            pts = click_2d.get('points', [])
            if not pts:
                return no_update
            x   = round(float(pts[0].get('x', 3.0)), 1)
            y   = round(float(pts[0].get('y', 3.0)), 1)
            z   = round(int(floor_idx or 0) * float(fh or 3.0) + 1.0, 1)
            nid = max((p['id'] for p in store), default=-1) + 1
            return store + [{'id':nid,'name':f'P{nid+1}','x':x,'y':y,'z':z,'mass_kg':70}]

        if tid == 'graph-3d':
            if (place_mode or 'person') != 'person':
                return no_update
            if not click_3d:
                return no_update
            pts = click_3d.get('points', [])
            if not pts:
                return no_update
            pt = pts[0]
            # Only accept clicks on our invisible grid (customdata present)
            if 'customdata' not in pt:
                return no_update
            x   = round(float(pt.get('x', 3.0)), 1)
            y   = round(float(pt.get('y', 3.0)), 1)
            z   = round(float(pt.get('z', 1.2)), 1)
            nid = max((p['id'] for p in store), default=-1) + 1
            return store + [{'id':nid,'name':f'P{nid+1}','x':x,'y':y,'z':z,'mass_kg':70}]

        if tid == 'btn-add-person':
            z   = round(int(floor_idx or 0) * float(fh or 3.0) + 1.0, 1)
            nid = max((p['id'] for p in store), default=-1) + 1
            return store + [{'id':nid,'name':f'P{nid+1}','x':3.0,'y':3.0,'z':z,'mass_kg':70}]

        if tid == 'btn-clear-people':
            return []

        if isinstance(tid, dict) and tid.get('type') == 'person-del':
            triggered_nclicks = (ctx.triggered[0].get('value') or 0) if ctx.triggered else 0
            if triggered_nclicks < 1:
                return no_update
            return [p for p in store if p['id'] != tid['index']]

        return no_update

    # ── Exit placement ──────────────────────────────────────────────────────

    @app.callback(
        Output('exits-store', 'data'),
        Input('floor-plan',      'clickData'),
        Input('graph-3d',        'clickData'),
        Input('btn-add-exit',    'n_clicks'),
        Input('btn-clear-exits', 'n_clicks'),
        Input({'type':'exit-del','index':ALL}, 'n_clicks'),
        State('exits-store',   'data'),
        State('floor-selector','value'),
        State('bld-fh',        'value'),
        State('place-mode',    'value'),
        prevent_initial_call=True,
    )
    def manage_exits(click_2d, click_3d, _add, _clear, _del_ns,
                     store, floor_idx, fh, place_mode):
        store = store or []
        tid   = ctx.triggered_id

        if tid == 'floor-plan':
            if (place_mode or 'person') != 'exit':
                return no_update
            if not click_2d:
                return no_update
            pts = click_2d.get('points', [])
            if not pts:
                return no_update
            x   = round(float(pts[0].get('x', 6.0)), 1)
            y   = round(float(pts[0].get('y', 0.0)), 1)
            fi  = int(floor_idx or 0)
            z   = round(fi * float(fh or 3.0), 1)
            nid = max((e['id'] for e in store), default=-1) + 1
            return store + [{'id':nid,'name':f'E{nid+1}','x':x,'y':y,'floor_idx':fi,'z':z}]

        if tid == 'graph-3d':
            if (place_mode or 'person') != 'exit':
                return no_update
            if not click_3d:
                return no_update
            pts = click_3d.get('points', [])
            if not pts:
                return no_update
            pt = pts[0]
            if 'customdata' not in pt:
                return no_update
            x   = round(float(pt.get('x', 6.0)), 1)
            y   = round(float(pt.get('y', 0.0)), 1)
            cd  = pt.get('customdata', [0, 0.0])
            fi  = int(cd[0]) if cd else 0
            z   = round(float(cd[1]) - 1.2, 1) if cd else 0.0   # floor z_bot
            nid = max((e['id'] for e in store), default=-1) + 1
            return store + [{'id':nid,'name':f'E{nid+1}','x':x,'y':y,'floor_idx':fi,'z':z}]

        if tid == 'btn-add-exit':
            fi  = int(floor_idx or 0)
            z   = round(fi * float(fh or 3.0), 1)
            nid = max((e['id'] for e in store), default=-1) + 1
            return store + [{'id':nid,'name':f'E{nid+1}','x':0.0,'y':5.0,'floor_idx':fi,'z':z}]

        if tid == 'btn-clear-exits':
            return []

        if isinstance(tid, dict) and tid.get('type') == 'exit-del':
            triggered_nclicks = (ctx.triggered[0].get('value') or 0) if ctx.triggered else 0
            if triggered_nclicks < 1:
                return no_update
            return [e for e in store if e['id'] != tid['index']]

        return no_update

    # ── Floor-plan + sidebar lists ──────────────────────────────────────────

    @app.callback(
        Output('people-list-sidebar', 'children'),
        Output('floor-plan',          'figure'),
        Output('exits-list-sidebar',  'children'),
        Input('people-store',        'data'),
        Input('sim-store',           'data'),
        Input('floor-selector',      'value'),
        Input('exits-store',         'data'),
        Input('rescue-routes-store', 'data'),
        Input('blueprint-store',     'data'),
        State('bld-width',  'value'), State('bld-depth',  'value'),
        State('bld-floors', 'value'), State('bld-fh',     'value'),
        State('bld-rx',     'value'), State('bld-ry',     'value'),
        State('blast-x',    'value'), State('blast-y',    'value'),
        State('blast-z',    'value'),
        State('place-mode', 'value'),
        prevent_initial_call=False,
    )
    def update_plan_ui(people, sim_store, floor_idx, exits_data, rescue_routes,
                       blueprint, width, depth, n_floors, fh, n_rx, n_ry, bx, by, bz,
                       place_mode):
        people       = people       or []
        exits_data   = exits_data   or []
        rescue_routes = rescue_routes or []

        people_injuries = _get_people_injuries(people, sim_store)
        rooms_data      = (sim_store or {}).get('rooms', [])

        # When a blueprint is loaded, derive dimensions from it
        blueprint_rooms = None
        if blueprint:
            bldg = blueprint.get('building', {})
            fh   = float(bldg.get('floor_height', fh or 3.0))
            fi   = int(floor_idx or 0)
            floor_defs = blueprint.get('floors', [])
            if fi < len(floor_defs):
                blueprint_rooms = floor_defs[fi].get('rooms', [])
            elif floor_defs:
                blueprint_rooms = floor_defs[-1].get('rooms', [])
            if blueprint_rooms:
                all_x = [r['x_max'] for r in blueprint_rooms]
                all_y = [r['y_max'] for r in blueprint_rooms]
                width = max(all_x)
                depth = max(all_y)

        blast_pos = np.array([float(bx or 6), float(by or -3), float(bz or 0.5)])
        fp = build_floor_plan(
            width=float(width or 12), depth=float(depth or 10),
            n_rooms_x=int(n_rx or 2), n_rooms_y=int(n_ry or 2),
            floor_idx=int(floor_idx or 0), floor_height=float(fh or 3),
            blast_pos=blast_pos, people=people,
            people_injuries=people_injuries,
            exits=exits_data,
            place_mode=place_mode or 'person',
            rescue_routes=rescue_routes,
            rooms_data=rooms_data,
            blueprint_rooms=blueprint_rooms,
        )
        return render_people_list(people), fp, render_exits_list(exits_data)

    # ── Injury panel ────────────────────────────────────────────────────────

    @app.callback(
        Output('people-injuries-panel', 'children'),
        Input('sim-store',    'data'),
        Input('people-store', 'data'),
        prevent_initial_call=True,
    )
    def update_injury_panel(sim_store, people):
        people = people or []
        if not people:
            return html.Div('Place people on the floor plan, then run the simulation.',
                            style={'fontSize':'11px','color':_MUTED,'padding':'10px'})
        if not sim_store or 'panel_states' not in sim_store:
            return html.Div('Run the simulation to see injury assessments.',
                            style={'fontSize':'11px','color':_MUTED,'padding':'10px'})
        injuries = _get_people_injuries(people, sim_store)
        if not injuries:
            return html.Div('No injuries to display.',
                            style={'fontSize':'11px','color':_MUTED})
        return [render_injury_card(inj) for inj in injuries]

    # ── Rescue plan panel ────────────────────────────────────────────────────

    @app.callback(
        Output('rescue-plan-panel',    'children'),
        Output('rescue-routes-store',  'data'),
        Input('sim-store',    'data'),
        Input('exits-store',  'data'),
        Input('people-store', 'data'),
        prevent_initial_call=True,
    )
    def update_rescue_panel(sim_store, exits_data, people_data):
        exits_data  = exits_data  or []
        people_data = people_data or []

        placeholder = lambda msg: (
            html.Div(msg, style={'fontSize':'11px','color':_MUTED,'padding':'10px',
                                 'whiteSpace':'pre-line','lineHeight':'1.6'}),
            [],
        )

        if not exits_data:
            return placeholder(
                '🚪  No exits placed.\n\n'
                'Switch to Exit mode (toggle above the floor plan) '
                'and click the plan to mark exit / entrance locations.',
            )
        if not people_data:
            return placeholder(
                '👤  No people placed.\n\n'
                'Switch to Person mode and click the floor plan to place people.',
            )
        if not sim_store or 'panel_states' not in sim_store:
            return placeholder(
                '▶  Run the simulation to generate rescue routes.',
            )

        rooms_data = sim_store.get('rooms', [])
        if not rooms_data:
            return placeholder(
                'Building room data unavailable — please re-run the simulation.',
            )

        people_injuries = _get_people_injuries(people_data, sim_store)
        if not people_injuries:
            return placeholder('No injury data available.')

        # Build panel-damage dict for pathfinding
        panel_damage = {}
        for pid_str, pd in (sim_store.get('panels') or {}).items():
            panel_damage[int(pid_str)] = {
                'damage_index': pd.get('damage_index', 0.0),
                'failed':       pd.get('failed', False),
            }

        exits = [Exit(**e) for e in exits_data]

        routes = compute_rescue_routes(
            people_injuries=people_injuries,
            exits=exits,
            rooms_data=rooms_data,
            panel_states=sim_store.get('panel_states', []),
            panel_damage=panel_damage,
        )

        if not routes:
            return placeholder('No rescue routes could be computed.')

        # Serialise for rescue-routes-store (used by floor-plan overlay)
        routes_serial = [
            {
                'person_id':       r.person_id,
                'person_name':     r.person_name,
                'person_severity': r.person_severity,
                'rescue_priority': r.rescue_priority,
                'exit_id':         r.exit_id,
                'exit_name':       r.exit_name,
                'path_room_ids':   r.path_room_ids,
                'path_cost':       r.path_cost,
                'room_path_names': r.room_path_names,
                'path_description': r.path_description,
                'reachable':       r.reachable,
            }
            for r in routes
        ]

        children = [
            html.Div([
                html.Span('🏃  Rescue priority order',
                          style={'fontSize':'11px','color':_MUTED,
                                 'fontWeight':'600'}),
                html.Span(f'  ({len(routes)} people)',
                          style={'fontSize':'10px','color':_MUTED}),
            ], style={'marginBottom':'8px','padding':'0 2px'}),
        ]
        for r_dict in routes_serial:
            children.append(render_rescue_card(r_dict))

        return children, routes_serial

    # ── Explosive-label live update ──────────────────────────────────────────

    @app.callback(
        Output('explosive-label', 'children'),
        Input('blast-explosive',  'value'),
        Input('blast-mass',       'value'),
    )
    def update_explosive_label(exp_key, mass):
        exp_key = exp_key or 'tnt'
        mass    = float(mass or 25)
        name, op_f, imp_f = TNT_EQUIVALENCY.get(exp_key, ('TNT', 1.0, 1.0))
        w_tnt = mass * op_f
        return f'→ W_TNT = {w_tnt:.1f} kg  (×{op_f:.2f} overpressure)'

    # ── Main simulation run ──────────────────────────────────────────────────

    @app.callback(
        Output('graph-3d',      'figure'),
        Output('graph-results', 'figure'),
        Output('metrics-row',   'children'),
        Output('run-status',    'children'),
        Output('run-status',    'className'),
        Output('header-status', 'children'),
        Output('sim-store',     'data'),
        Input('btn-run', 'n_clicks'),
        State('bld-width',  'value'), State('bld-depth',  'value'),
        State('bld-floors', 'value'), State('bld-fh',     'value'),
        State('bld-mat',    'value'), State('bld-rx',     'value'),
        State('bld-ry',     'value'), State('bld-wf',     'value'),
        State('blast-mass',      'value'),
        State('blast-explosive', 'value'),
        State('blast-x',    'value'), State('blast-y',    'value'),
        State('blast-z',    'value'), State('blast-type', 'value'),
        State('people-store',    'data'),
        State('blueprint-store', 'data'),
        prevent_initial_call=True,
    )
    def run_sim(n_clicks,
                width, depth, n_floors, fh, mat_key, n_rx, n_ry, wfrac,
                mass, explosive_key, bx, by, bz, btype, people_data, blueprint):
        try:
            if blueprint:
                # Re-parse in case the store holds raw FloorSpaceJS
                bp = _auto_parse_blueprint(blueprint)
                panels, columns, rooms = create_building_from_blueprint(bp)
                bldg      = bp.get('building', {})
                n_floors  = int(bldg.get('n_floors', n_floors or 3))
                fh        = float(bldg.get('floor_height', fh or 3.0))
                floor_defs = bp.get('floors', [])
                all_rooms  = [r for fd in floor_defs for r in fd.get('rooms', [])]
                width  = max((r['x_max'] for r in all_rooms), default=float(width or 12))
                depth  = max((r['y_max'] for r in all_rooms), default=float(depth or 10))
            else:
                wall_mat = MATERIALS.get(mat_key or 'concrete')
                panels, columns, rooms = create_building(
                    width=float(width or 12), depth=float(depth or 10),
                    n_floors=int(n_floors or 3), floor_height=float(fh or 3),
                    wall_mat=wall_mat,
                    n_rooms_x=int(n_rx or 2), n_rooms_y=int(n_ry or 2),
                    window_frac=float(wfrac or 0.3),
                    total_occupants=0,
                )
            exp_key = explosive_key or 'tnt'
            blast = BlastSource(
                x=float(bx or 6), y=float(by or -3), z=float(bz or 0.5),
                tnt_kg=float(mass or 25), burst_type=btype or 'surface',
                explosive_type=exp_key,
            )
            result = run_simulation(
                blast=blast, panels=panels, columns=columns, rooms=rooms,
                n_floors=int(n_floors or 3), floor_height=float(fh or 3),
                building_width=float(width or 12), building_depth=float(depth or 10),
            )

            # Serialise to store
            store = {
                'panels': {},
                'panel_states': [],
                'blast': {'x': float(bx or 6), 'y': float(by or -3),
                          'z': float(bz or 0.5), 'W_eff': blast.W_eff,
                          'burst_type': btype or 'surface'},
                'config': {'width': float(width or 12), 'depth': float(depth or 10),
                           'n_floors': int(n_floors or 3), 'floor_height': float(fh or 3)},
                # Room data needed for rescue pathfinding
                'rooms': [
                    {
                        'id':       r.id,
                        'floor_idx': r.floor_idx,
                        'x_min':    r.x_min, 'x_max': r.x_max,
                        'y_min':    r.y_min, 'y_max': r.y_max,
                        'z_bot':    r.z_bot, 'z_top': r.z_top,
                    }
                    for r in rooms
                ],
            }
            for p_id, pr in result.panel_results.items():
                panel = panels[p_id]
                store['panels'][str(p_id)] = {
                    'time':          pr.time.tolist(),
                    'displacement':  pr.displacement.tolist(),
                    'pressure':      pr.blast_pressure.tolist(),
                    'yield_disp_mm': panel.yield_disp * 1000.0,
                    'ductility':     panel.material.ductility,
                    'damage_index':  pr.damage_index,
                    'failed':        pr.failed,
                }
            for panel in panels:
                store['panel_states'].append({
                    'id':            panel.id,
                    'center':        panel.center.tolist(),
                    'normal':        panel.normal.tolist(),
                    'corners':       panel.corners.tolist(),
                    'failed':        panel.failed,
                    'panel_type':    panel.panel_type,
                    'material_name': panel.material.name,
                    'room_inside':   panel.room_inside,
                    'room_outside':  panel.room_outside,
                })

            # People injuries for 3D view
            people_injuries = _get_people_injuries(people_data or [], store)

            fig_3d      = build_3d_figure(panels, columns, result, people_injuries,
                                          blast_pos=(float(bx or 6), float(by or -3),
                                                     float(bz or 0.5)))
            fig_results = build_results_figure(rooms, result)

            # Casualty counts from individually placed people
            n_fatal  = sum(1 for inj in people_injuries if inj.overall_severity == 'Fatal')
            n_severe = sum(1 for inj in people_injuries
                           if inj.overall_severity in ('Severe', 'Moderate'))
            n_minor  = sum(1 for inj in people_injuries if inj.overall_severity == 'Minor')

            # Structural stability
            n_total      = max(len(panels), 1)
            n_failed_p   = sum(1 for p in panels if p.failed)
            max_drift    = float(max(result.building.story_drift)) \
                           if len(result.building.story_drift) else 0.0
            panel_health = 1.0 - n_failed_p / n_total
            drift_health = max(0.0, 1.0 - max_drift / 0.05)
            stability    = (panel_health * 0.6 + drift_health * 0.4) * 100
            stab_color   = (_GREEN  if stability > 70
                            else _ORANGE if stability > 40
                            else _RED)

            metrics_children = dbc.Row([
                dbc.Col(_tile(str(n_fatal),  'Fatal',       _RED),       width=3),
                dbc.Col(_tile(str(n_severe), 'Severe Inj.', _ORANGE),    width=3),
                dbc.Col(_tile(str(n_minor),  'Minor Inj.',  _BLUE),      width=3),
                dbc.Col(_tile(f"{stability:.0f}%", 'Stability', stab_color), width=3),
            ], className='g-2').children

            exp_name = TNT_EQUIVALENCY.get(exp_key, ('TNT',))[0]
            status_txt = f'✓  {len(panels)} panels · {len(rooms)} rooms · done'
            header_txt = (f'{float(mass or 25):.0f} kg {exp_name}  ·  '
                          f'W_TNT={blast.W_eff / (2.0 if (btype or "surface") == "surface" else 1.0):.1f} kg  ·  '
                          f'({float(bx):.1f}, {float(by):.1f}, {float(bz):.1f}) m')

            return (fig_3d, fig_results,
                    metrics_children, status_txt, 'status-ok', header_txt, store)

        except Exception:
            print(traceback.format_exc())
            ef  = _empty_fig()
            mtiles = dbc.Row([
                dbc.Col(_tile('Err', 'Fatal',       _RED),   width=3),
                dbc.Col(_tile('Err', 'Severe Inj.', _ORANGE),width=3),
                dbc.Col(_tile('Err', 'Minor Inj.',  _BLUE),  width=3),
                dbc.Col(_tile('Err', 'Stability',   _MUTED), width=3),
            ], className='g-2').children
            return (ef, ef, mtiles, '⚠  Error', 'status-error', 'Error', {})

    return app


# ── Shared helper: reconstruct persons + compute injuries from store ────────

def _get_people_injuries(people_data, sim_store):
    if not people_data or not sim_store or 'panel_states' not in sim_store:
        return []

    panels_proxy = []
    pr_proxy = {}
    for ps in sim_store['panel_states']:
        mat = SimpleNamespace(name=ps['material_name'])
        p   = SimpleNamespace(
            id=ps['id'],
            center=np.array(ps['center']),
            normal=np.array(ps['normal']),
            corners=np.array(ps['corners']),
            failed=ps['failed'],
            panel_type=ps['panel_type'],
            material=mat,
        )
        panels_proxy.append(p)
        pr_proxy[p.id] = SimpleNamespace(failed=p.failed)

    bi = sim_store['blast']
    blast = BlastSource(x=bi['x'], y=bi['y'], z=bi['z'],
                        tnt_kg=bi['W_eff'] / 2.0,
                        burst_type=bi['burst_type'])
    blast.W_eff = bi['W_eff']

    injuries = []
    for pd in people_data:
        person = Person(id=pd['id'], name=pd['name'],
                        x=pd['x'], y=pd['y'], z=pd['z'],
                        mass_kg=pd.get('mass_kg', 70))
        inj = assess_person_injuries(person, blast, panels_proxy, pr_proxy)
        injuries.append(inj)
    return injuries
