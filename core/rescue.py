"""
Rescue route computation for blast simulations.

Builds a room-adjacency graph and uses Dijkstra's algorithm to find
the optimal evacuation route from each person to the nearest exit.

Edge costs (traversal difficulty):
  intact wall (DI < 0.5)  : 1.0
  damaged wall (DI >= 0.5) : 2.0   — debris slows movement
  failed wall              : 5.0   — structural hazard / major debris
  (exterior panels are not room-to-room edges)

Rescue priority (lower = rescued first):
  1 Severe    — immediate life threat
  2 Moderate  — urgent
  3 Minor     — delayed
  4 Uninjured — self-evacuation capable
  5 Fatal     — expectant / body recovery
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import heapq


# ---------------------------------------------------------------------------
# Priority mapping
# ---------------------------------------------------------------------------

RESCUE_PRIORITY: Dict[str, int] = {
    'Severe':    1,
    'Moderate':  2,
    'Minor':     3,
    'Uninjured': 4,
    'Fatal':     5,
}

PRIORITY_LABEL: Dict[int, str] = {
    1: 'Immediate',
    2: 'Urgent',
    3: 'Delayed',
    4: 'Minor',
    5: 'Expectant',
}

_COST_INTACT  = 1.0
_COST_DAMAGED = 2.0
_COST_FAILED  = 5.0


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Exit:
    id: int
    name: str
    x: float
    y: float
    floor_idx: int
    z: float


@dataclass
class RescueRoute:
    person_id: int
    person_name: str
    person_severity: str
    rescue_priority: int          # 1 = highest

    exit_id: int
    exit_name: str

    path_room_ids: List[int]
    path_cost: float
    room_path_names: List[str]
    path_description: str

    reachable: bool = True


# ---------------------------------------------------------------------------
# Graph helpers
# ---------------------------------------------------------------------------

def _find_room(x: float, y: float, floor_idx: int, rooms_data: list) -> Optional[int]:
    """
    Return the id of the room containing (x, y) on floor_idx.
    Falls back to nearest room centroid if the point is outside all rooms.
    """
    # Exact containment
    for r in rooms_data:
        if r['floor_idx'] != floor_idx:
            continue
        if r['x_min'] <= x <= r['x_max'] and r['y_min'] <= y <= r['y_max']:
            return r['id']
    # Nearest centroid fallback (handles points on perimeter / outside building)
    best_id, best_d = None, float('inf')
    for r in rooms_data:
        if r['floor_idx'] != floor_idx:
            continue
        cx = (r['x_min'] + r['x_max']) / 2.0
        cy = (r['y_min'] + r['y_max']) / 2.0
        d  = (x - cx) ** 2 + (y - cy) ** 2
        if d < best_d:
            best_d, best_id = d, r['id']
    return best_id


def _build_graph(
    rooms_data: list,
    panel_states: list,
    panel_damage: Dict[int, dict],
) -> Dict[int, List[Tuple[int, float]]]:
    """
    Build undirected adjacency list from interior wall panels.

    rooms_data   : list of dicts with keys id, floor_idx, x_min/max, y_min/max
    panel_states : list of dicts with keys id, panel_type, room_inside, room_outside, failed
    panel_damage : {panel_id -> {damage_index, failed}}
    """
    adj: Dict[int, List[Tuple[int, float]]] = {r['id']: [] for r in rooms_data}

    for ps in panel_states:
        if ps.get('panel_type') != 'int_wall':
            continue
        ri = ps.get('room_inside',  -1)
        ro = ps.get('room_outside', -1)
        if ri < 0 or ro < 0:
            continue
        if ri not in adj or ro not in adj:
            continue

        pd      = panel_damage.get(ps['id'], {})
        failed  = ps.get('failed', False) or pd.get('failed', False)
        di      = pd.get('damage_index', 0.0)

        if failed:
            cost = _COST_FAILED
        elif di >= 0.5:
            cost = _COST_DAMAGED
        else:
            cost = _COST_INTACT

        adj[ri].append((ro, cost))
        adj[ro].append((ri, cost))

    return adj


def _dijkstra(
    adj: Dict[int, List[Tuple[int, float]]],
    start: Optional[int],
    goals: set,
) -> Tuple[Optional[int], float, List[int]]:
    """
    Single-source shortest path from *start* to any node in *goals*.
    Returns (reached_goal, total_cost, path_as_list_of_room_ids).
    """
    if start is None or not goals:
        return None, float('inf'), []

    dist: Dict[int, float] = {start: 0.0}
    prev: Dict[int, Optional[int]] = {start: None}
    pq: List[Tuple[float, int]] = [(0.0, start)]

    while pq:
        d, u = heapq.heappop(pq)
        if d > dist.get(u, float('inf')):
            continue
        if u in goals:
            path: List[int] = []
            node: Optional[int] = u
            while node is not None:
                path.append(node)
                node = prev.get(node)
            path.reverse()
            return u, d, path
        for v, w in adj.get(u, []):
            nd = d + w
            if nd < dist.get(v, float('inf')):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))

    return None, float('inf'), []


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_rescue_routes(
    people_injuries,          # list of PersonInjury (from persons.py)
    exits: List[Exit],
    rooms_data: list,         # list of room dicts (from sim-store)
    panel_states: list,       # list of panel dicts (from sim-store)
    panel_damage: Dict[int, dict],  # {panel_id: {damage_index, failed}}
) -> List[RescueRoute]:
    """
    Compute optimal rescue route for every person.

    Returns a list of RescueRoute objects sorted by rescue_priority (1 first).
    """
    if not people_injuries:
        return []
    if not exits:
        return [
            RescueRoute(
                person_id=inj.person_id, person_name=inj.name,
                person_severity=inj.overall_severity,
                rescue_priority=RESCUE_PRIORITY.get(inj.overall_severity, 4),
                exit_id=-1, exit_name='No exits defined',
                path_room_ids=[], path_cost=float('inf'),
                room_path_names=[], path_description='No exits placed on this floor.',
                reachable=False,
            )
            for inj in people_injuries
        ]

    adj = _build_graph(rooms_data, panel_states, panel_damage)

    # Map each exit to its room id
    exit_room: Dict[int, Optional[int]] = {}
    for e in exits:
        exit_room[e.id] = _find_room(e.x, e.y, e.floor_idx, rooms_data)

    # Goal set: room ids that contain at least one exit
    goal_rooms = {rid for rid in exit_room.values() if rid is not None}

    # Reverse map: room_id -> list of exit_ids in that room
    room_exits: Dict[int, List[int]] = {}
    for eid, rid in exit_room.items():
        if rid is not None:
            room_exits.setdefault(rid, []).append(eid)

    exit_by_id = {e.id: e for e in exits}
    room_map   = {r['id']: r for r in rooms_data}

    routes: List[RescueRoute] = []

    for inj in people_injuries:
        px, py, pz = inj.position

        # Determine which floor the person is on
        person_floor = 0
        for r in rooms_data:
            if r['z_bot'] <= pz < r['z_top']:
                person_floor = r['floor_idx']
                break

        person_room = _find_room(px, py, person_floor, rooms_data)

        if not goal_rooms:
            routes.append(RescueRoute(
                person_id=inj.person_id, person_name=inj.name,
                person_severity=inj.overall_severity,
                rescue_priority=RESCUE_PRIORITY.get(inj.overall_severity, 4),
                exit_id=-1, exit_name='No exits reachable',
                path_room_ids=[], path_cost=float('inf'),
                room_path_names=[],
                path_description='No exits placed on any accessible floor.',
                reachable=False,
            ))
            continue

        goal_exit_room, cost, path = _dijkstra(adj, person_room, goal_rooms)

        if goal_exit_room is None:
            routes.append(RescueRoute(
                person_id=inj.person_id, person_name=inj.name,
                person_severity=inj.overall_severity,
                rescue_priority=RESCUE_PRIORITY.get(inj.overall_severity, 4),
                exit_id=-1, exit_name='Blocked',
                path_room_ids=[], path_cost=float('inf'),
                room_path_names=[],
                path_description='No passable route — structural collapse may block all paths.',
                reachable=False,
            ))
        else:
            # Pick first exit in the goal room
            eids = room_exits.get(goal_exit_room, [])
            chosen = exit_by_id.get(eids[0]) if eids else None

            path_names = []
            for rid in path:
                r = room_map.get(rid)
                if r:
                    path_names.append(f'R{rid}(F{r["floor_idx"]+1})')

            desc_parts = [' → '.join(path_names) if path_names else '(same room as exit)']
            if cost >= _COST_FAILED:
                desc_parts.append('⚠️ Route through failed structure — extreme hazard')
            elif cost >= _COST_DAMAGED:
                desc_parts.append('⚡ Damaged panels on route — debris risk')

            routes.append(RescueRoute(
                person_id=inj.person_id, person_name=inj.name,
                person_severity=inj.overall_severity,
                rescue_priority=RESCUE_PRIORITY.get(inj.overall_severity, 4),
                exit_id=chosen.id if chosen else -1,
                exit_name=chosen.name if chosen else 'Exit',
                path_room_ids=path,
                path_cost=cost,
                room_path_names=path_names,
                path_description=' · '.join(desc_parts),
                reachable=True,
            ))

    routes.sort(key=lambda r: r.rescue_priority)
    return routes
