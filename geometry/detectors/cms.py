"""The CMS drift tube detector: 36 chambers, their honeycomb and the iron yoke.

The chambers come from the CMS XML through mplDTs' source; the yoke cannot be
data yet, because a wedge with the chambers subtracted needs boolean solids, so
it is written here as ASCII text (T-060 replaces it with a CSG backend).

The build is the development subset: three wheels, three sectors, four station
types. The whole detector is the same comprehension over every sector, 1...14.
"""

import numpy as np

from mpldts.geometry.sources.cms_xml import CMSXMLSource

NAME = "CMS"

#: The development subset (36 chambers, 25896 cells, a few seconds).
WHEELS = (-1, 0, 1)
SECTORS = (1, 2, 12)
STATIONS = (1, 2, 3, 4)


def chambers():
    """[(label, key, StationSpec)] of the CMS chambers, in build order."""
    source = CMSXMLSource()
    return [
        ("W%d_Sec%d_St%d" % (wheel, sector, station), (wheel, sector, station),
         source.station(wheel, sector, station))
        for wheel in WHEELS
        for sector in SECTORS
        for station in STATIONS
    ]


def families():
    """The volume families the concentrator styles and overlap-checks."""
    return ["Station", "Yoke", "Honeycomb", "DriftCell"]


def _rotation(rotation, name):
    """A ROTM line; -0.0 would print as "-0.000000" and break a byte comparison."""
    values = np.asarray(rotation, dtype=float).flatten()
    values = np.where(values == 0.0, 0.0, values)
    return ":ROTM %s %s\n" % (name, " ".join("%.6f" % value for value in values))


def _sector_bounds(specs):
    """The iron wedge that wraps one sector: radial and z extent of its chambers.

    The radial extent is the chamber thickness and its z extent the chamber
    length, because a CMS chamber lies with its length along the beam.
    """
    r_min, r_max, z_min, z_max = 1e9, -1e9, 1e9, -1e9
    for spec in specs:
        centre = spec.pose.translation
        _, length, thickness = spec.shape.aabb()
        radius = float(np.hypot(centre[0], centre[1]))
        r_min = min(r_min, radius - thickness / 2.0)
        r_max = max(r_max, radius + thickness / 2.0)
        z_min = min(z_min, centre[2] - length / 2.0)
        z_max = max(z_max, centre[2] + length / 2.0)
    # the inner radius leaves room for the coil, as the hand-written yoke did
    return r_min - 100.0, r_max + 5.0, z_min - 5.0, z_max + 5.0, (z_min + z_max) / 2.0


def yoke(chambers, path):
    """Write the fused iron yoke with every chamber subtracted, or None."""
    sectors = {}
    for _, key, spec in chambers:
        sectors.setdefault((key[0], key[1]), []).append(spec)
    if not sectors:
        return None

    lines = ["// Unified Yoke Assembly with Station Holes\n\n"]
    lines.append("// Base Iron Wedges (TUBS) for all sectors\n")
    wedges = {}
    for (wheel, sector), specs in sectors.items():
        r_min, r_max, z_min, z_max, z_centre = _sector_bounds(specs)
        centre_deg = (sector - 1) * 30.0
        name = "IronTrap_W%d_Sec%d" % (wheel, sector)
        lines.append("// Wheel %d Sector %d\n" % (wheel, sector))
        lines.append(
            ":SOLID %s TUBS %.2f*cm %.2f*cm %.2f*cm %.2f*deg 30.00*deg\n"
            % (name, r_min, r_max, (z_max - z_min) / 2.0, centre_deg - 15.0)
        )
        wedges[(wheel, sector)] = (name, z_centre)

    lines.append("\n// Union all sector wedges into single volume\n")
    keys = sorted(wedges)
    reference_name, reference_z = wedges[keys[0]]
    current = reference_name
    for index, key in enumerate(keys[1:], start=1):
        name, z_centre = wedges[key]
        union = "Solid_Union_%d" % index
        lines.append(
            ":SOLID %s UNION %s %s R0 0*cm 0*cm %.6f*cm\n"
            % (union, current, name, z_centre - reference_z)
        )
        current = union

    lines.append("\n// Subtract all station volumes to create holes\n")
    count = 0
    for _, key, spec in chambers:
        wheel, sector, station = key
        centre = spec.pose.translation
        count += 1
        rotation = "RM_Yoke_Station_%d_%d_%d" % (wheel, sector, station)
        lines.append(_rotation(spec.pose.rotation, rotation))
        lines.append(
            ":SOLID Solid_Yoke_Sub%d SUBTRACTION %s Solid_Station_W%d_Sec%d_St%d %s "
            "%.6f*cm %.6f*cm %.6f*cm\n"
            % (count, current, wheel, sector, station, rotation,
               centre[0], centre[1], centre[2] - reference_z)
        )
        current = "Solid_Yoke_Sub%d" % count

    lines.append("\n// Final yoke volume and placement\n")
    lines.append(":VOLU Yoke_Fused %s G4_Fe\n" % current)
    lines.append(":PLACE Yoke_Fused 1 world R0 0*cm 0*cm %.6f*cm\n" % reference_z)
    path.write_text("".join(lines))
    return path.name
