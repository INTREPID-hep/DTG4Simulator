"""
CMS Drift Tube Geometry Generator for Geant4

This module generates Geant4 ASCII text geometry files for CMS drift tube (DT)
stations: drift cells and honeycomb layers, plus the iron yoke structures.

A station file is the mplDTs node tree written by the .tg exporter; the yoke is
built here because it needs boolean solids (UNION / SUBTRACTION).
"""

from dataclasses import replace
from pathlib import Path

import numpy as np

from mpldts.export import TgSink, export_volumes, iter_volume_placements
from mpldts.geometry import Station
from mpldts.geometry.model import build_station
from mpldts.geometry.presets import cms
from mpldts.geometry.sources.cms_xml import CMSXMLSource
from mpldts.geometry.specs import Pose

#: the CMS chambers: the XML is parsed once for the whole run
SOURCE = CMSXMLSource()


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def extract_rotation(matrix):
    """
    Extract 3x3 rotation matrix from 4x4 transformation matrix.
    
    Args:
        matrix: 4x4 transformation matrix
        
    Returns:
        Flattened rotation matrix as 9 values (XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ)
    """
    return matrix[:3, :3].flatten()


def format_rotation_matrix(rot, name):
    """
    Format rotation matrix for Geant4 geometry file.
    
    Args:
        rot: Flattened 9-element rotation matrix
        name: Name for the rotation matrix
        
    Returns:
        Formatted ROTM line
    """
    return (f":ROTM {name} "
            f"{rot[0]:.6f} {rot[1]:.6f} {rot[2]:.6f} "
            f"{rot[3]:.6f} {rot[4]:.6f} {rot[5]:.6f} "
            f"{rot[6]:.6f} {rot[7]:.6f} {rot[8]:.6f}\n")


# =============================================================================
# STATION GEOMETRY
# =============================================================================

def _volume_names(wheel, sector, station_num):
    """
    Names of a station .tg, as the simulator parses them back.

    StationSD reads the chamber from "Station_W<w>_Sec<s>_St<n>", DriftCellSD
    reads the super layer and the cells-per-layer encoding from
    "DriftCell_W<w>_Sec<s>_St<n>_SL<n>_<encoding>" and maps the copy number of
    that volume onto (layer, wire). The scheme lives here, not in mplDTs: the
    exporter takes it as a hook (T-005, T-057).
    """

    def name_of(node):
        if node.kind == "station":
            return "Station_W%d_Sec%d_St%d" % (wheel, sector, station_num)
        if node.kind == "attachment":
            return "Honeycomb_W%d_Sec%d_St%d" % (wheel, sector, station_num)
        if node.kind == "cell":
            super_layer = node.parent.parent
            encoding = "".join("%02d" % len(layer.children) for layer in super_layer.children)
            return "DriftCell_W%d_Sec%d_St%d_%s_%s" % (
                wheel, sector, station_num, super_layer.name, encoding
            )
        return None

    return name_of


def station_node(wheel, sector, station_num, ref_frame="CMS", include_honeycomb=False):
    """
    mplDTs node tree of one chamber, ready for the exporter.

    Args:
        wheel, sector, station_num: Station identifiers
        ref_frame: "CMS" places the station in the detector, "Station" writes it
            at the world origin with no rotation (single chamber files)
        include_honeycomb: add the aluminium panel between SL2 and SL1

    Returns:
        The station Node
    """
    spec = SOURCE.station(wheel, sector, station_num)
    if include_honeycomb:
        spec = cms.with_attachments(spec, SOURCE.record(wheel, sector, station_num))
    if ref_frame == "Station":
        spec = replace(spec, pose=Pose.identity())
    return build_station(spec)


def export_station(wheel, sector, station_num, path, ref_frame="CMS", include_honeycomb=False):
    """
    Write one station .tg through the mplDTs exporter.

    Returns:
        Number of drift cells placed in the file
    """
    node = station_node(wheel, sector, station_num, ref_frame, include_honeycomb)
    export_volumes(node, TgSink(str(path), name_of=_volume_names(wheel, sector, station_num)))
    return sum(1 for p in iter_volume_placements(node) if p.node.kind == "cell")


# =============================================================================
# STATION COMPONENT GENERATORS
# =============================================================================


# =============================================================================
# YOKE GENERATION
# =============================================================================

def _calculate_sector_bounds(wheel, sector, stations_in_sector):
    """
    Calculate bounding cylinder for all stations in a sector.
    
    Args:
        wheel: Wheel number
        sector: Sector number
        stations_in_sector: List of station numbers
        
    Returns:
        Tuple (r_min, r_max, z_min, z_max, z_center, valid) or None if invalid
    """
    r_min, r_max = 1e9, -1e9
    z_min, z_max = 1e9, -1e9
    valid = False
    
    for st_num in stations_in_sector:
        try:
            st = Station(wheel=wheel, sector=sector, station=st_num)
            
            # Radial extent
            center = st.global_center
            r = np.sqrt(center[0]**2 + center[1]**2)
            h = st.bounds[1]
            r_min = min(r_min, r - h / 2)
            r_max = max(r_max, r + h / 2)
            
            # Z extent
            z, l = center[2], st.bounds[2]
            z_min = min(z_min, z - l / 2)
            z_max = max(z_max, z + l / 2)
            
            valid = True
        except:
            continue
    
    if not valid:
        return None
    
    # Add safety margins (r_min - 100 cm to compensate the presence of the coil)
    return (r_min - 100.0, r_max + 5.0, z_min - 5.0, z_max + 5.0, 
            (z_max + z_min) / 2.0, valid)


def generate_yoke(stations_list, sector_stations):
    """
    Generate unified yoke assembly in a single file using boolean operations.
    
    Creates yoke.tg containing:
    1. Base TUBS solids for each sector
    2. Union operations to combine all sectors
    3. Subtraction operations to create holes for stations
    4. Final yoke volume and placement
    
    Args:
        stations_list: List of (wheel, sector, station) tuples
        sector_stations: Dict {(wheel, sector): [station_nums]}
        
    Returns:
        Filename of the generated yoke file ("yoke.tg")
    """
    filename = "yoke.tg"
    filepath = Path(filename)
    
    geometry = "// Unified Yoke Assembly with Station Holes\n\n"
    
    # 1. Generate base TUBS solids for each sector
    geometry += "// Base Iron Wedges (TUBS) for all sectors\n"
    base_yoke_info = {}  # {(wheel, sector): (solid_name, z_center)}
    
    for (wheel, sector), stations_in_sector in sector_stations.items():
        bounds = _calculate_sector_bounds(wheel, sector, stations_in_sector)
        if bounds is None:
            continue
        
        r_min, r_max, z_min, z_max, z_center, _ = bounds
        
        # Calculate TUBS parameters
        z_half = (z_max - z_min) / 2.0
        
        # Phi angle for sector (30° per sector)
        sector_center_deg = (sector - 1) * 30.0
        phi_start = sector_center_deg - 15.0
        phi_delta = 30.0
        
        solid_name = f"IronTrap_W{wheel}_Sec{sector}"
        
        geometry += f"// Wheel {wheel} Sector {sector}\n"
        geometry += (f":SOLID {solid_name} TUBS "
                    f"{r_min:.2f}*cm {r_max:.2f}*cm {z_half:.2f}*cm "
                    f"{phi_start:.2f}*deg {phi_delta:.2f}*deg\n")
        
        base_yoke_info[(wheel, sector)] = (solid_name, z_center)
    
    if not base_yoke_info:
        return ""
    
    geometry += "\n"
    
    # 2. Union all base solids
    geometry += "// Union all sector wedges into single volume\n"
    sorted_keys = sorted(base_yoke_info.keys())
    ref_key = sorted_keys[0]
    ref_solid_name, ref_z_center = base_yoke_info[ref_key]
    current_solid_name = ref_solid_name
    
    for i, key in enumerate(sorted_keys[1:], start=1):
        next_solid, next_z = base_yoke_info[key]
        union_name = f"Solid_Union_{i}"
        z_shift = next_z - ref_z_center
        
        geometry += (f":SOLID {union_name} UNION {current_solid_name} {next_solid} R0 "
                    f"0*cm 0*cm {z_shift:.6f}*cm\n")
        current_solid_name = union_name
    
    geometry += "\n"
    
    # 3. Subtract all station volumes to create holes
    geometry += "// Subtract all station volumes to create holes\n"
    subtraction_count = 0
    
    for w, s, st_num in stations_list:
        try:
            st = Station(wheel=w, sector=s, station=st_num)
        except:
            continue
        
        # The solid comes from the station .tg, named Solid_Station_... there.
        # The rotation is defined here under its own name: the station file
        # defines the one it uses for its own placement, and both are the same
        # Station -> CMS matrix.
        st_solid_name = f"Solid_Station_W{w}_Sec{s}_St{st_num}"
        rot_name = f"RM_Yoke_Station_{w}_{s}_{st_num}"
        st_center = st.global_center
        geometry += format_rotation_matrix(
            extract_rotation(
                st.transformer.get_transformation(from_frame="Station", to_frame="CMS")
            ),
            rot_name
        )
        
        # Position relative to reference frame
        rel_z = st_center[2] - ref_z_center
        
        subtraction_count += 1
        next_solid_name = f"Solid_Yoke_Sub{subtraction_count}"
        
        geometry += (f":SOLID {next_solid_name} SUBTRACTION {current_solid_name} {st_solid_name} {rot_name} "
                    f"{st_center[0]:.6f}*cm {st_center[1]:.6f}*cm {rel_z:.6f}*cm\n")
        
        current_solid_name = next_solid_name
    
    # 4. Create final volume and place in world
    geometry += f"\n// Final yoke volume and placement\n"
    geometry += f":VOLU Yoke_Fused {current_solid_name} G4_Fe\n"
    geometry += f":PLACE Yoke_Fused 1 world R0 0*cm 0*cm {ref_z_center:.6f}*cm\n"
    
    filepath.write_text(geometry)
    
    return filename


# =============================================================================
# MAIN GENERATION FUNCTIONS
# =============================================================================

def generate_station_geometry_ascii(
    stations_list,
    output_dir, 
    concentrator_template,
    concentrator_output,
    include_yoke=False,
    include_honeycomb=False,
    ref_frame="CMS"
):
    """
    Generate Geant4 ASCII text geometry for multiple DT stations.
    Each station is saved in a separate file, and a concentrator file includes them all.
    
    Args:
        stations_list: List of (wheel, sector, station_num) tuples
        output_dir: Directory where individual station files will be saved
        concentrator_template: Path to concentrator template file
        concentrator_output: Output concentrator geometry file
        include_yoke: If True, generate sector-wide iron yokes with holes
        include_honeycomb: If True, add aluminum honeycomb layer between SL2 and SL1
        ref_frame: reference frame for the station placement, "CMS" or "Station"
    """
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Generate individual station files
    station_files = []
    
    # Group stations by sector for yoke generation
    sector_stations = {} # (wheel, sector) -> [station_nums]
    total_cells = 0
    
    for wheel, sector, station_num in stations_list:
        # Add to sector group
        key = (wheel, sector)
        if key not in sector_stations:
            sector_stations[key] = []
        sector_stations[key].append(station_num)
        
        # Create filename for this station
        station_filename = f"station_W{wheel}_Sec{sector}_St{station_num}.tg"
        station_files.append(station_filename)
        
        # Write the station: the mplDTs node tree through the .tg exporter
        cells = export_station(
            wheel, sector, station_num, Path(output_dir) / station_filename,
            ref_frame, include_honeycomb
        )
        total_cells += cells
        
        print(f"  ✓ Generated: {station_filename} ({cells} cells)")
    
    # Generate Fused Yoke if requested
    yoke_file = ""
    if include_yoke:
        print("\nGenerating Fused Yoke...")
        yoke_file = generate_yoke(stations_list, sector_stations)
        if yoke_file:
            print(f"  ✓ Generated: {yoke_file}")

    # Create concentrator file with #include directives
    _create_concentrator_file(
        station_files, yoke_file, concentrator_template, concentrator_output,
        len(stations_list), total_cells
    )


def _create_concentrator_file(station_files, yoke_file, template_path, output_path,
                              n_stations, n_cells):
    """
    Create the main concentrator geometry file that includes all station files.
    
    Args:
        station_files: List of station filenames
        yoke_file: Filename of the yoke geometry file
        template_path: Path to concentrator template file
        output_path: Output path for concentrator file
        n_stations: Number of stations included
        n_cells: Total number of drift cells placed
    """
    include_lines = "\n".join([f"#include geometry/stations/{fname}" for fname in station_files])
    if yoke_file:
        include_lines += f"\n#include geometry/{yoke_file}"
    
    # Read concentrator template
    with open(template_path, 'r') as f:
        concentrator_content = f.read()
    
    # Replace {STATIONS} placeholder with include directives
    concentrator_content = concentrator_content.replace('{STATIONS}', include_lines)
    
    # Write concentrator file
    with open(output_path, 'w') as f:
        f.write(concentrator_content)
    
    # Print summary statistics
    print(f"\n✓ Concentrator file: {output_path}")
    print(f"  Stations: {n_stations}")
    print(f"  Total cells: {n_cells}")


if __name__ == "__main__":
    # The development subset: three wheels, three sectors, four station types
    # (36 chambers, 25896 cells, ~10 s). The whole detector is the same loops
    # over every sector, 1...14.
    stations_to_generate = [
        (wheel, sector, station)
        for wheel in (-1, 0, 1)
        for sector in (1, 2, 12)
        for station in (1, 2, 3, 4)
    ]

    print("Generating DT geometry files...")
    print(f"Stations to generate: {len(stations_to_generate)}")

    # The concentrator template asks for Yoke*, Honeycomb* and DriftCell*
    # volumes and it places the chambers in CMS, so all three switches have to
    # agree with it: without the yoke or the honeycomb Geant4 aborts at setup,
    # and in the "Station" frame every chamber sits on the world origin.
    generate_station_geometry_ascii(
        stations_list=stations_to_generate,
        output_dir='stations',
        concentrator_template='geometry_concentrator_template',
        concentrator_output='geometry_concentrator.tg',
        include_yoke=True,
        include_honeycomb=True,
        ref_frame="CMS",
    )
