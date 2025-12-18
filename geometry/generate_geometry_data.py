"""
CMS Drift Tube Geometry Generator for Geant4

This module generates Geant4 ASCII text geometry files for CMS drift tube (DT) stations,
including drift cells, honeycomb layers, and iron yoke structures.
"""

from mpldts.geometry import Station
import numpy as np
from pathlib import Path


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
# STATION COMPONENT GENERATORS
# =============================================================================
def create_station_volume(station, wheel, sector, station_num):
    """
    Create the station volume definition.
    
    Args:
        station: Station object from mplDTs
        wheel, sector, station_num: Station identifiers
        
    Returns:
        String with station volume and placement definitions
    """
    geometry = ""
    station_bounds = station.bounds
    station_name = f"Station_W{wheel}_Sec{sector}_St{station_num}"
    solid_name = f"Solid_{station_name}"
    
    # Define station solid and volume (BOX: half-widths in X, Y, Z)
    geometry += f":SOLID {solid_name} BOX "
    geometry += f"{station_bounds[0]/2:.6f}*cm {station_bounds[2]/2:.6f}*cm {station_bounds[1]/2:.6f}*cm\n"
    geometry += f":VOLU {station_name} {solid_name} G4_AIR\n"
    
    # Get station transformation and rotation
    station_transform = station.transformer.get_transformation(from_frame="Station", to_frame="CMS")
    station_rot = extract_rotation(station_transform)
    station_center = station.global_center
    
    # Create rotation matrix
    rot_name = f"RM_Station_{wheel}_{sector}_{station_num}"
    geometry += format_rotation_matrix(station_rot, rot_name)
    
    # Place station in world
    geometry += (f":PLACE {station_name} 1 world {rot_name} "
                 f"{station_center[0]:.6f}*cm {station_center[1]:.6f}*cm {station_center[2]:.6f}*cm\n\n")
    
    return geometry, station_name, station_bounds


def create_superlayer_cells(station, wheel, sector, station_num, station_name):
    """
    Create drift cell volumes for all superlayers in the station.
    
    Defines ONE logical volume per superlayer and places multiple copies.
    The volume name encodes the number of cells per layer for identification.
    
    Args:
        station: Station object from mplDTs
        wheel, sector, station_num: Station identifiers
        station_name: Name of the station volume
        
    Returns:
        Geometry string with cell volume definitions and placements
    """
    geometry = ""
    
    for sl in station.super_layers:
        geometry += f"// SuperLayer {sl.number} - Drift Cell Volume\n"
        
        # Count cells per layer and encode in volume name
        # Format: 2 digits per layer count, e.g., [60, 60, 59, 59] -> "60605959"
        cells_per_layer = [len(layer.cells) for layer in sl.layers]
        layer_encoding = "".join([f"{count:02d}" for count in cells_per_layer])
        
        # Define drift cell volume using first cell as reference
        if not sl.layers or not sl.layers[0].cells:
            continue
            
        first_cell = sl.layers[0].cells[0]
        cell_volume_name = f"DriftCell_W{wheel}_Sec{sector}_St{station_num}_SL{sl.number}_{layer_encoding}"
        cell_bounds = first_cell.bounds
        
        geometry += f":VOLU {cell_volume_name} BOX "
        geometry += f"{cell_bounds[0]/2:.6f}*cm {cell_bounds[2]/2:.6f}*cm {cell_bounds[1]/2:.6f}*cm GasMixture\n\n"
        
        # SL2 needs its own rotation matrix (rotated 90° around Z)
        if sl.number == 2:
            sl_transform = sl.transformer.get_transformation(from_frame="SuperLayer", to_frame="Station")
            sl_rot = extract_rotation(sl_transform)
            rotation_to_use = f"RM_SL2_W{wheel}_Sec{sector}_St{station_num}"
            geometry += format_rotation_matrix(sl_rot, rotation_to_use)
        else:
            # SL1 and SL3 use identity rotation
            rotation_to_use = "R0"
            sl_transform = None
        
        # Place multiple copies of the cell volume
        geometry += f"// Placing {sum(cells_per_layer)} cells ({' '.join([str(c) for c in cells_per_layer])} per layer)\n"
        copy_number = 1
        
        for layer in sl.layers:
            for cell in layer.cells:
                # Get cell position relative to station
                cell_center = cell.local_center
                if sl.number == 2 and sl_transform is not None:
                    cell_center = np.dot(sl_transform[:3, :3], cell.local_center)
                
                # Place copy of the cell volume in station
                geometry += f":PLACE {cell_volume_name} {copy_number} {station_name} {rotation_to_use} "
                geometry += f"{cell_center[0]:.6f}*cm {cell_center[1]:.6f}*cm {cell_center[2]:.6f}*cm\n"
                copy_number += 1
        geometry += "\n"
    
    return geometry


def create_honeycomb(station, wheel, sector, station_num, station_name, station_bounds):
    """
    Create aluminum honeycomb layer between SL2 and SL1 (or SL3 for MB4).

    Args:
        station: Station object from mplDTs
        wheel, sector, station_num: Station identifiers
        station_name: Name of the station volume
        station_bounds: Station dimensions [width, height, length]
        
    Returns:
        Geometry string with honeycomb volume and placement
    """
    geometry = ""
    honeycomb_name = f"Honeycomb_W{wheel}_Sec{sector}_St{station_num}"
    
    # Get superlayer objects to calculate actual gap
    sl1 = [sl for sl in station.super_layers if sl.number == 1]
    sl2 = [sl for sl in station.super_layers if sl.number == 2]
    
    # MB4 stations don't have SL2, use SL3 instead
    if not sl2:
        sl2 = [sl for sl in station.super_layers if sl.number == 3]
    
    if not sl1 or not sl2:
        return ""
    
    sl1, sl2 = sl1[0], sl2[0]
    
    # Calculate honeycomb thickness from actual gap
    # Z coordinate is the height/stacking direction in local station frame
    sl1_z_min = sl1.local_cords_at_min[2]  # Bottom of SL1
    sl2_z_max = sl2.local_cords_at_min[2] + sl2.bounds[1]  # Top of SL2
    honeycomb_thickness = sl1_z_min - sl2_z_max  # Gap between them (~12.8 cm)
    
    # Honeycomb volume
    geometry += f"// Aluminum Honeycomb (between SL2 and SL1)\n"
    geometry += f":VOLU {honeycomb_name} BOX "
    geometry += f"{station_bounds[0]/2 - 0.5:.6f}*cm {station_bounds[2]/2 - 0.5:.6f}*cm "
    geometry += f"{honeycomb_thickness/2 - 0.5:.6f}*cm G4_Al\n"
    
    # Position: midpoint between SL2 (top) and SL1 (bottom) in Z
    honeycomb_z = (sl2_z_max + sl1_z_min) / 2.0
    
    # Place honeycomb in station (no rotation needed)
    geometry += f":PLACE {honeycomb_name} 1 {station_name} R0 "
    geometry += f"0.0*cm 0.0*cm {honeycomb_z:.6f}*cm\n\n"
    
    return geometry


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
        
        st_solid_name = f"Solid_Station_W{w}_Sec{s}_St{st_num}"
        rot_name = f"RM_Station_{w}_{s}_{st_num}"
        st_center = st.global_center
        
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
    include_honeycomb=False
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
    """
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Generate individual station files
    station_files = []
    
    # Group stations by sector for yoke generation
    sector_stations = {} # (wheel, sector) -> [station_nums]
    
    for wheel, sector, station_num in stations_list:
        # Add to sector group
        key = (wheel, sector)
        if key not in sector_stations:
            sector_stations[key] = []
        sector_stations[key].append(station_num)
        
        station = Station(wheel=wheel, sector=sector, station=station_num)
        
        # Create filename for this station
        station_filename = f"station_W{wheel}_Sec{sector}_St{station_num}.tg"
        station_filepath = Path(output_dir) / station_filename
        station_files.append(station_filename)
        
        # Generate geometry for this station
        station_geometry = (f"// ===== Station: Wheel={wheel}, Sector={sector}, Station={station_num} =====\n"
                           f"// {station.name}\n\n")
        
        # Create station volume and placement
        station_vol, station_name, station_bounds = create_station_volume(
            station, wheel, sector, station_num
        )
        station_geometry += station_vol
        
        # NOTE: Yoke generation is now handled per-sector, not per-station
        
        # Add drift cells for all superlayers
        station_geometry += create_superlayer_cells(
            station, wheel, sector, station_num, station_name
        )
        
        # Add honeycomb layer if requested
        if include_honeycomb:
            station_geometry += create_honeycomb(
                station, wheel, sector, station_num, station_name, station_bounds
            )
        
        # Write this station's geometry to its own file
        station_filepath.write_text(station_geometry)
        
        print(f"  ✓ Generated: {station_filename}")
    
    # Generate Fused Yoke if requested
    yoke_file = ""
    if include_yoke:
        print("\nGenerating Fused Yoke...")
        yoke_file = generate_yoke(stations_list, sector_stations)
        if yoke_file:
            print(f"  ✓ Generated: {yoke_file}")

    # Create concentrator file with #include directives
    _create_concentrator_file(
        station_files, yoke_file, concentrator_template, concentrator_output, stations_list
    )


def _create_concentrator_file(station_files, yoke_file, template_path, output_path, stations_list):
    """
    Create the main concentrator geometry file that includes all station files.
    
    Args:
        station_files: List of station filenames
        yoke_file: Filename of the yoke geometry file
        template_path: Path to concentrator template file
        output_path: Output path for concentrator file
        stations_list: List of station tuples for statistics
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
    print(f"  Stations: {len(stations_list)}")
    
    # Count total cells
    total_cells = 0
    for wheel, sector, station_num in stations_list:
        station = Station(wheel=wheel, sector=sector, station=station_num)
        for sl in station.super_layers:
            for layer in sl.layers:
                total_cells += len(layer.cells)
    
    print(f"  Total cells: {total_cells}")


if __name__ == "__main__":
    # Define stations to include (wheel, sector, station)
    stations_to_generate = [
        (0, 12, 1),  # MB1
        (0, 12, 2),  # MB2
        (0, 12, 3),  # MB3
        (0, 12, 4),  # MB4
        (0, 1, 1),  # MB1
        (0, 1, 2),  # MB2
        (0, 1, 3),  # MB3
        (0, 1, 4),  # MB4
        (0, 2, 1),  # MB1
        (0, 2, 2),  # MB2
        (0, 2, 3),  # MB3
        (0, 2, 4),  # MB4
    ]
    
    print("Generating DT geometry files...")
    print(f"Stations to generate: {len(stations_to_generate)}")
    
    # Generate geometry
    generate_station_geometry_ascii(
        stations_list=stations_to_generate,
        output_dir='stations',
        concentrator_template='geometry_concentrator_template',
        concentrator_output='geometry_concentrator.tg',
        include_yoke=True,       # Include iron yoke
        include_honeycomb=True   # Include aluminum honeycomb between SL2 and SL1
    )
