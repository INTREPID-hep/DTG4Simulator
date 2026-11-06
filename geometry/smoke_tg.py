"""T-005 smoke test: build one CMS station and check the .tg output.

Run against the mplDTs checkout:

    cd /home/destrada/Geant4/DTG4Simulator/geometry
    PYTHONPATH=/home/destrada/mplDTs ../env/bin/python smoke_tg.py

Checks: the file is non-empty, the cell volume and every placement exist, and
the SL2 rotation matrix is the expected one.
"""

import sys
import tempfile
from pathlib import Path

GEOMETRY = Path(__file__).resolve().parent
sys.path.insert(0, str(GEOMETRY))

EXPECTED_CELLS = 620
SL2_ROTATION = (
    ":ROTM RM_SL2_W0_Sec1_St1 0.000000 1.000000 0.000000 "
    "-1.000000 0.000000 0.000000 0.000000 0.000000 1.000000"
)


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
        )
        path = Path(tmp) / "station_W0_Sec1_St1.tg"
        assert path.exists(), "station file was not generated"
        text = path.read_text()

    assert text.strip(), "station file is empty"
    assert ":VOLU " in text, "no :VOLU (cell volume) defined"

    place = [line for line in text.splitlines() if line.startswith(":PLACE")]
    assert len(place) == 1 + EXPECTED_CELLS, "expected %d :PLACE lines, got %d" % (
        1 + EXPECTED_CELLS,
        len(place),
    )
    assert SL2_ROTATION in text, "the SL2 rotation matrix is missing or changed"

    print("ok: %d placements, SL2 rotation present" % len(place))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
