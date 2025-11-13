from mpldts.geometry import Station
import numpy as np
from pathlib import Path


def calculate_yoke_thicknesses(wheel, sector):
    """
    Calculate yoke thickness for each station type by measuring the radial gap
    between consecutive stations.
    
    - MB1 yoke: uses MB1 station height (innermost station)
    - MB2 yoke: gap between MB1 and MB2
    - MB3 yoke: gap between MB2 and MB3
    - MB4 yoke: gap between MB3 and MB4
    
    Args:
        wheel: Wheel number
        sector: Sector number
        
    Returns:
        Dictionary mapping station_num -> yoke_thickness in cm
    """
    yoke_thickness = {}
    
    # Get all stations
    stations = {}
    for st_num in [1, 2, 3, 4]:
        try:
            station = Station(wheel=wheel, sector=sector, station=st_num)
            stations[st_num] = station
        except:
            continue
    
    # MB1 yoke uses the station height itself (innermost station)
    if 1 in stations:
        yoke_thickness[1] = stations[1].bounds[1]
    
    # Calculate gaps between consecutive stations for MB2, MB3, MB4
    for st_num in [1, 2, 3]:
        next_st = st_num + 1
        
        if st_num not in stations or next_st not in stations:
            continue
        
        station = stations[st_num]
        next_station = stations[next_st]
        
        # Outer edge of current station (radial distance + half height)
        center = station.global_center
        radial_dist = np.sqrt(center[0]**2 + center[1]**2)
        outer_edge = radial_dist + station.bounds[1]/2
        
        # Inner edge of next station (radial distance - half height)
        next_center = next_station.global_center
        next_radial_dist = np.sqrt(next_center[0]**2 + next_center[1]**2)
        next_inner_edge = next_radial_dist - next_station.bounds[1]/2
        
        # Gap between them is the yoke thickness for the next station
        gap = next_inner_edge - outer_edge - 5.0  # Subtract 5 cm safety margin
        yoke_thickness[next_st] = gap
    
    return yoke_thickness


# Calculate yoke thicknesses once (using a reference station)
# These values are consistent across all wheels/sectors for same station type
YOKE_THICKNESS = calculate_yoke_thicknesses(wheel=-1, sector=2)


def extract_rotation(matrix):
    """
    Extract rotation from 4x4 transformation matrix.
    
    Args:
        matrix: 4x4 transformation matrix
        
    Returns:
        rotation_flat: 9 values (XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ)
    """
    rotation = matrix[:3, :3]
    rotation_flat = rotation.flatten()
    return rotation_flat


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
    
    # Define station volume (BOX: half-widths in X, Y, Z) - merged syntax
    geometry += f":VOLU {station_name} BOX "
    geometry += f"{station_bounds[0]/2:.6f} {station_bounds[2]/2:.6f} {station_bounds[1]/2:.6f} G4_AIR\n"
    
    # Get station transformation and rotation
    station_transform = station.transformer.get_transformation(from_frame="Station", to_frame="CMS")
    station_rot = extract_rotation(station_transform)
    station_center = station.global_center
    
    # Create rotation matrix
    geometry += f":ROTM RM_Station_{wheel}_{sector}_{station_num} "
    geometry += f"{station_rot[0]:.6f} {station_rot[1]:.6f} {station_rot[2]:.6f} "
    geometry += f"{station_rot[3]:.6f} {station_rot[4]:.6f} {station_rot[5]:.6f} "
    geometry += f"{station_rot[6]:.6f} {station_rot[7]:.6f} {station_rot[8]:.6f}\n"
    
    # Place station in world
    geometry += f":PLACE {station_name} 1 world RM_Station_{wheel}_{sector}_{station_num} "
    geometry += f"{station_center[0]:.6f} {station_center[1]:.6f} {station_center[2]:.6f}\n\n"
    
    return geometry, station_name, station_bounds


def create_yoke(station, wheel, sector, station_num, station_bounds):
    """
    Create iron yoke volume positioned radially below the station.
    
    Yoke thickness fills the radial gap between consecutive stations:
    - MB1 yoke: uses MB1 station height (innermost, no station below)
    - MB2 yoke: fills gap between MB1 and MB2 (radially inward from MB2)
    - MB3 yoke: fills gap between MB2 and MB3 (radially inward from MB3)
    - MB4 yoke: fills gap between MB3 and MB4 (radially inward from MB4)
    
    Args:
        station: Station object from mplDTs
        wheel, sector, station_num: Station identifiers
        station_bounds: Station dimensions [width, height, length]
        
    Returns:
        String with yoke volume and placement definitions
    """
    geometry = ""
    yoke_name = f"Yoke_W{wheel}_Sec{sector}_St{station_num}"
    station_center = station.global_center
    
    # Get yoke thickness for this station type (gap to next station)
    yoke_thickness = YOKE_THICKNESS.get(station_num, 40.0)  # Default if not found
    
    # Yoke volume: merged solid+volume syntax
    geometry += f"// Iron Yoke (radially inward, thickness={yoke_thickness:.2f} cm)\n"
    geometry += f":VOLU {yoke_name} BOX "
    geometry += f"{station_bounds[0]/2:.6f} {station_bounds[2]/2:.6f} {yoke_thickness/2:.6f} G4_Fe\n"
    
    # Calculate yoke position: radially inward from station in XY plane
    station_pos_xy = np.array([station_center[0], station_center[1], 0.0])
    radial_distance_xy = np.linalg.norm(station_pos_xy)
    radial_unit_xy = station_pos_xy / radial_distance_xy
    
    # Move yoke inward by station_height/2 + yoke_thickness/2
    offset = (station_bounds[1]/2 + yoke_thickness/2)
    yoke_pos_xy = station_pos_xy - radial_unit_xy * offset
    yoke_pos = np.array([yoke_pos_xy[0], yoke_pos_xy[1], station_center[2]])
    
    # Place yoke with same rotation as station
    geometry += f":PLACE {yoke_name} 1 world RM_Station_{wheel}_{sector}_{station_num} "
    geometry += f"{yoke_pos[0]:.6f} {yoke_pos[1]:.6f} {yoke_pos[2]:.6f}\n\n"
    
    return geometry


def create_superlayer_cells(station, wheel, sector, station_num, station_name):
    """
    Create all drift cells for all superlayers in the station.
    Defines ONE logical volume per superlayer and places multiple copies.
    Volume name encodes cells-per-layer information.
    
    Args:
        station: Station object from mplDTs
        wheel, sector, station_num: Station identifiers
        station_name: Name of the station volume
        
    Returns:
        String with cell volume definitions and placements
    """
    geometry = ""
    
    for sl in station.super_layers:
        geometry += f"// SuperLayer {sl.number} - Drift Cell Volume\n"
        
        # Count cells per layer to encode in volume name
        cells_per_layer = [len(layer.cells) for layer in sl.layers]
        
        # Format: up to 2 digits per layer count (pad with leading zero if needed)
        # E.g., [60, 60, 59, 59] -> "60605959"
        layer_encoding = "".join([f"{count:02d}" for count in cells_per_layer])
        
        # Define ONE drift cell volume for this superlayer with encoded name
        first_cell = None
        if len(sl.layers) > 0 and len(sl.layers[0].cells) > 0:
            first_cell = sl.layers[0].cells[0]
        
        cell_volume_name = f"DriftCell_W{wheel}_Sec{sector}_St{station_num}_SL{sl.number}_{layer_encoding}"
        
        if first_cell:
            cell_bounds = first_cell.bounds
            geometry += f":VOLU {cell_volume_name} BOX "
            geometry += f"{cell_bounds[0]/2:.6f} {cell_bounds[2]/2:.6f} {cell_bounds[1]/2:.6f} GasMixture\n\n"
        
        # SL2 needs its own rotation matrix (rotated 90° in Z)
        if sl.number == 2:
            sl_transform = sl.transformer.get_transformation(from_frame="SuperLayer", to_frame="Station")
            sl_rot = extract_rotation(sl_transform)
            
            sl_rot_name = f"RM_SL2_W{wheel}_Sec{sector}_St{station_num}"
            geometry += f":ROTM {sl_rot_name} "
            geometry += f"{sl_rot[0]:.6f} {sl_rot[1]:.6f} {sl_rot[2]:.6f} "
            geometry += f"{sl_rot[3]:.6f} {sl_rot[4]:.6f} {sl_rot[5]:.6f} "
            geometry += f"{sl_rot[6]:.6f} {sl_rot[7]:.6f} {sl_rot[8]:.6f}\n"
            
            rotation_to_use = sl_rot_name
        else:
            # SL1 and SL3 use identity rotation
            rotation_to_use = "R0"
            sl_transform = None
        
        # Place multiple copies of the same logical volume for all cells
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
                geometry += f"{cell_center[0]:.6f} {cell_center[1]:.6f} {cell_center[2]:.6f}\n"
                copy_number += 1
        geometry += "\n"
    
    return geometry


def create_honeycomb(station, wheel, sector, station_num, station_name, station_bounds):
    """
    Create aluminum honeycomb layer between SL2 and SL1.

    Args:
        station: Station object from mplDTs
        wheel, sector, station_num: Station identifiers
        station_name: Name of the station volume
        station_bounds: Station dimensions [width, height, length]
        
    Returns:
        String with honeycomb volume and placement definitions
    """
    geometry = ""
    honeycomb_name = f"Honeycomb_W{wheel}_Sec{sector}_St{station_num}"
    
    # Get superlayer objects to calculate actual gap
    sl1 = [sl for sl in station.super_layers if sl.number == 1]
    sl2 = [sl for sl in station.super_layers if sl.number == 2]
    
    # Handle MB4 stations that don't have SL2 (use SL3 instead)
    if not sl2:
        sl2 = [sl for sl in station.super_layers if sl.number == 3]

    sl1 = sl1[0]
    sl2 = sl2[0]
    
    # Calculate honeycomb thickness from actual gap
    # Z coordinate is the height/stacking direction in local station frame
    sl1_z_min = sl1.local_cords_at_min[2]  # Bottom of SL1
    sl2_z_max = sl2.local_cords_at_min[2] + sl2.bounds[1]  # Top of SL2
    honeycomb_thickness = sl1_z_min - sl2_z_max  # Gap between them (~12.8 cm)
    
    # Honeycomb volume: merged syntax
    geometry += f"// Aluminum Honeycomb (between SL2 and SL1)\n"
    geometry += f":VOLU {honeycomb_name} BOX "
    geometry += f"{station_bounds[0]/2 - 0.5:.6f} {station_bounds[2]/2 - 0.5:.6f} {honeycomb_thickness/2 - 0.5:.6f} G4_Al\n"
    
    # Position: midpoint between SL2 (top) and SL1 (bottom) in Z
    honeycomb_z = (sl2_z_max + sl1_z_min) / 2.0
    honeycomb_pos = [0.0, 0.0, honeycomb_z]
    
    # Place honeycomb in station (no rotation needed)
    geometry += f":PLACE {honeycomb_name} 1 {station_name} R0 "
    geometry += f"{honeycomb_pos[0]:.6f} {honeycomb_pos[1]:.6f} {honeycomb_pos[2]:.6f}\n\n"
    
    return geometry


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
        include_yoke: If True, add iron yoke blocks radially below each station
        include_honeycomb: If True, add aluminum honeycomb layer between SL2 and SL1
    """
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Generate individual station files
    station_files = []
    
    for wheel, sector, station_num in stations_list:
        station = Station(wheel=wheel, sector=sector, station=station_num)
        
        # Create filename for this station
        station_filename = f"station_W{wheel}_Sec{sector}_St{station_num}.tg"
        station_filepath = f"{output_dir}/{station_filename}"
        station_files.append(station_filename)
        
        # Generate geometry for this station
        station_geometry = ""
        station_geometry += f"// ===== Station: Wheel={wheel}, Sector={sector}, Station={station_num} =====\n"
        station_geometry += f"// {station.name}\n\n"
        
        # Create station volume and placement
        station_vol, station_name, station_bounds = create_station_volume(
            station, wheel, sector, station_num
        )
        station_geometry += station_vol
        
        # Add yoke if requested
        if include_yoke:
            station_geometry += create_yoke(
                station, wheel, sector, station_num, station_bounds
            )
        
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
        with open(station_filepath, 'w') as f:
            f.write(station_geometry)
        
        print(f"  ✓ Generated: {station_filename}")
    
    # Create concentrator file with #include directives
    _create_concentrator_file(
        station_files, concentrator_template, concentrator_output, stations_list
    )


def _create_concentrator_file(station_files, template_path, output_path, stations_list):
    """
    Create the main concentrator geometry file that includes all station files.
    
    Args:
        station_files: List of station filenames
        template_path: Path to concentrator template file
        output_path: Output path for concentrator file
        stations_list: List of station tuples for statistics
    """
    include_lines = "\n".join([f"#include geometry/stations/{fname}" for fname in station_files])
    
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
        (-1, 1, 1),  # MB1
        (-1, 1, 2),  # MB2
        (-1, 1, 3),  # MB3
        (-1, 1, 4),  # MB4
        (-1, 2, 1),  # MB1
        (-1, 2, 2),  # MB2
        (-1, 2, 3),  # MB3
        (-1, 2, 4),  # MB4
        (-1, 3, 1),  # MB1
        (-1, 3, 2),  # MB2
        (-1, 3, 3),  # MB3
        (-1, 3, 4),  # MB4
    ]
    
    print("Generating DT geometry files...")
    print(f"Stations to generate: {len(stations_to_generate)}")
    
    # Generate geometry
    generate_station_geometry_ascii(
        stations_list=stations_to_generate,
        output_dir='stations',
        concentrator_template='geometry_concentrator_template',
        concentrator_output='geometry_concentrator.tg',
        include_yoke=True,       # Include iron yokes below stations
        include_honeycomb=True   # Include aluminum honeycomb between SL2 and SL1
    )
