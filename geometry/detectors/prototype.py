"""The MuTel prototype: two pairs of super layers and a lead wall.

Four super layers of CMS-type cells (4.2 x 100 cm, 23 per layer, 4 layers per
super layer) in two identical pairs:

    z = +27 cm   pair B: SL1 phi, SL2 eta, back to back (no gap inside a pair)
    z =   0 cm   a 5 cm lead wall, 17.5 cm of air on each side
    z = -27 cm   pair A: SL1 phi, SL2 eta, identical to pair B

The 40 cm accessible gap, the lead included, is only between the pairs. Every
station is 14 cm thick (two super layers of 7), so the whole detector is
100 x 100 x 68 cm in the world.

Each pair is one station of two super layers, because the simulator decodes
at most three super layers per chamber. The chamber keys (0, 1, 1) and
(0, 1, 2) are what the Station_W0_Sec1_St<n> names carry, so the existing
StationSD and DriftCellSD decode the prototype without a new C++ class.
"""

from mpldts.geometry.specs import Box, CellGrid, LayerSpec, Pose, StationSpec, SuperLayerSpec

NAME = "MuTel"

#: One cell: drift gap along z, length along y, pitch along x.
PITCH = 4.2
LENGTH = 100.0
THICKNESS = 1.3
COUNT = 23
LAYERS = 4
#: Distance between layer centres: the cell plus the 0.6 cm between layers.
LAYER_PITCH = 1.9
#: Air between the cells and the chamber envelope.
MARGIN = 0.1

#: Extent of the four layers of one super layer: 7.0 cm, no air inside.
LAYER_SPAN = (LAYERS - 1) * LAYER_PITCH + THICKNESS
SUPER_LAYER = LAYER_SPAN
#: A pair: two super layers back to back, plus the envelope margin.
STATION = 2 * SUPER_LAYER + 2 * MARGIN
#: The row of cells, and the pair footprint (the eta cells are 100 cm long).
WIDTH = COUNT * PITCH
AREA = max(WIDTH, LENGTH) + 2 * MARGIN

#: The accessible gap between the pairs, and the lead wall in its middle.
PAIR_GAP = 40.0
LEAD = 5.0

MATERIAL = "GasMixture"
PANEL = "G4_Pb"

#: The two pairs, named one after the other from the positive z side.
PAIRS = (("MuTel_B", +1), ("MuTel_A", -1))
#: The chamber keys the C++ decoders read (wheel, sector, station).
KEYS = {"MuTel_A": (0, 1, 1), "MuTel_B": (0, 1, 2)}


def _layers():
    """Four layers of COUNT cells, span across x, length along y."""
    grid = CellGrid(
        pitch=PITCH,
        count=COUNT,
        shape=Box(PITCH, LENGTH, THICKNESS),
        first_id=1,
        material=MATERIAL,
    )
    return [
        LayerSpec(number=index, pose=Pose.translated(0.0, 0.0, (index - 2.5) * LAYER_PITCH),
                  cells=grid)
        for index in range(1, LAYERS + 1)
    ]


def station(label="MuTel", z=0.0):
    """One pair: SL1 with the cells along x (phi), SL2 turned 90 degrees (eta)."""
    super_layers = [
        SuperLayerSpec(
            number=1,
            pose=Pose.translated(0.0, 0.0, -SUPER_LAYER / 2.0),
            layers=_layers(),
        ),
        SuperLayerSpec(
            number=2,
            pose=Pose.rotated_z(90.0, (0.0, 0.0, SUPER_LAYER / 2.0)),
            layers=_layers(),
        ),
    ]
    return StationSpec(
        name=label,
        pose=Pose.translated(0.0, 0.0, z),
        super_layers=super_layers,
        shape=Box(AREA, AREA, STATION),
        material="G4_AIR",
    )


def lead_wall(z=0.0):
    """The 5 cm lead slab in the middle of the gap: a passive station solid."""
    return StationSpec(
        name="LeadWall",
        pose=Pose.translated(0.0, 0.0, z),
        super_layers=[],
        shape=Box(AREA, AREA, LEAD),
        material=PANEL,
    )


def chambers():
    """[(label, key, StationSpec)]: the two pairs and the wall between them."""
    half = PAIR_GAP / 2.0 + STATION / 2.0
    result = [(label, KEYS[label], station(label, sign * half)) for label, sign in PAIRS]
    result.append(("LeadWall", None, lead_wall(0.0)))
    return result


def families():
    """The volume families the concentrator styles and overlap-checks."""
    return ["Station", "DriftCell", "LeadWall"]
