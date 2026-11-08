"""T-005 smoke test: build one CMS station and check the .tg output.

Run against the mplDTs checkout:

    cd /home/destrada/Geant4/DTG4Simulator/geometry
    PYTHONPATH=/home/destrada/mplDTs ../env/bin/python smoke_tg.py

The station file is written by the mplDTs exporter, so the check is also the
downstream contract: the names the sensitive detectors parse, the copy numbers
they decode the cell from, and the SL2 rotation.
"""

import re
import sys
import tempfile
from pathlib import Path

GEOMETRY = Path(__file__).resolve().parent
sys.path.insert(0, str(GEOMETRY))

STATION = "Station_W0_Sec1_St1"
CELL_VOLUMES = [
    "DriftCell_W0_Sec1_St1_SL1_49504948",
    "DriftCell_W0_Sec1_St1_SL2_57585756",
    "DriftCell_W0_Sec1_St1_SL3_49504948",
]
EXPECTED_CELLS = 620
#: the 90 degree super layer 2 rotation, by value: the exporter names it
SL2_ROTATION = "0.000000 1.000000 0.000000 -1.000000 0.000000 0.000000 0.000000 0.000000 1.000000"


def main():
    import mpldts
    from generate_geometry_data import generate_station_geometry_ascii

    print("mplDTs:", mpldts.__file__)

    with tempfile.TemporaryDirectory() as tmp:
        generate_station_geometry_ascii(
            stations_list=[(0, 1, 1)],
            output_dir=tmp,
            concentrator_template=str(GEOMETRY / "geometry_concentrator_template"),
            concentrator_output=str(Path(tmp) / "geometry_concentrator.tg"),
            ref_frame="Station",
        )
        path = Path(tmp) / "station_W0_Sec1_St1.tg"
        assert path.exists(), "station file was not generated"
        text = path.read_text()

    assert text.strip(), "station file is empty"

    # the station solid and volume the yoke and StationSD refer to
    assert ":VOLU %s Solid_%s G4_AIR" % (STATION, STATION) in text
    assert ":PLACE %s 1 world R0" % STATION in text

    # one cell volume per super layer, with the encoding DriftCellSD parses
    for volume in CELL_VOLUMES:
        assert ":VOLU %s Solid_%s GasMixture" % (volume, volume) in text, volume

    place = [line for line in text.splitlines() if line.startswith(":PLACE")]
    assert len(place) == 1 + EXPECTED_CELLS, "expected %d :PLACE lines, got %d" % (
        1 + EXPECTED_CELLS,
        len(place),
    )

    # the SL2 rotation matrix is the 90 degree one
    assert SL2_ROTATION in text, "the SL2 rotation matrix is missing or changed"

    # every cell volume numbers its copies from 1: the decoder maps copyNo to
    # (layer, wire) inside the super layer
    copies = {}
    for line in text.splitlines():
        match = re.match(r":PLACE (\S+) (\d+) ", line)
        if match and match.group(1).startswith("DriftCell"):
            copies.setdefault(match.group(1), []).append(int(match.group(2)))
    assert sorted(copies) == sorted(CELL_VOLUMES), sorted(copies)
    for volume, numbers in copies.items():
        assert numbers == list(range(1, len(numbers) + 1)), volume

    print(
        "ok: %d placements, %d cell volumes, SL2 rotation present"
        % (len(place), len(copies))
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
