from dataclasses import dataclass


@dataclass
class Material:
    name: str
    E: float          # Young's modulus [Pa]
    nu: float         # Poisson's ratio
    rho: float        # density [kg/m³]
    fy: float         # yield / tensile strength [Pa]
    fu: float         # ultimate strength [Pa]
    ductility: float  # displacement ductility ratio μ (max / yield disp)
    damping_ratio: float = 0.05
    color: str = '#888888'


CONCRETE = Material(
    name='Reinforced Concrete',
    E=30e9, nu=0.2, rho=2400,
    fy=25e6, fu=35e6,
    ductility=4.0, damping_ratio=0.05,
    color='#B0B0B0',
)

BRICK = Material(
    name='Brick Masonry',
    E=6e9, nu=0.2, rho=1800,
    fy=6e6, fu=8e6,
    ductility=1.5, damping_ratio=0.07,
    color='#C46B40',
)

STEEL_PLATE = Material(
    name='Steel Plate',
    E=200e9, nu=0.3, rho=7850,
    fy=250e6, fu=400e6,
    ductility=12.0, damping_ratio=0.02,
    color='#607080',
)

GLASS = Material(
    name='Annealed Glass',
    E=70e9, nu=0.23, rho=2500,
    fy=30e6, fu=30e6,
    ductility=1.0, damping_ratio=0.01,
    color='#AADDFF',
)

STEEL_STRUCTURAL = Material(
    name='Structural Steel',
    E=200e9, nu=0.3, rho=7850,
    fy=345e6, fu=450e6,
    ductility=15.0, damping_ratio=0.02,
    color='#505060',
)

MATERIALS = {
    'concrete': CONCRETE,
    'brick': BRICK,
    'steel_plate': STEEL_PLATE,
    'glass': GLASS,
    'steel_structural': STEEL_STRUCTURAL,
}
