from mpldts.geometry import Station
import numpy as np
from pathlib import Path

def extract_rotation(matrix):
    """
    Extract rotation from 4x4 transformation matrix
    Returns: rotation_flat
    - rotation_flat: 9 values (XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ)
    """
    rotation = matrix[:3, :3]
    rotation_flat = rotation.flatten()
    
    return rotation_flat

def generate_station_geometry_ascii(
    stations_list,
    output_dir, 
    concentrator_template,
    concentrator_output,
    include_yoke=False
):
    """
    Generate Geant4 ASCII text geometry for multiple DT stations
    Each station is saved in a separate file, and a concentrator file includes them all.
    
    Args:
        stations_list: List of (wheel, sector, station_num) tuples
        output_dir: Directory where individual station files will be saved
        concentrator_template: Path to concentrator template file
        concentrator_output: Output concentrator geometry file
        include_yoke: If True, add iron yoke blocks radially below each station
    """
    
    # Create output directory if it doesn't exist
    from pathlib import Path
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
        
        station_bounds = station.bounds
        station_name = f"Station_W{wheel}_Sec{sector}_St{station_num}"
        station_geometry += f":VOLU {station_name} BOX {station_bounds[0]/2:.6f} {station_bounds[2]/2:.6f} {station_bounds[1]/2:.6f} G4_AIR\n"
        
        # Get station transformation
        station_transform = station.transformer.get_transformation(from_frame="Station", to_frame="CMS")
        station_rot = extract_rotation(station_transform)
        station_center = station.global_center
        
        # Create rotation matrix for station
        station_geometry += f":ROTM RM_Station_{wheel}_{sector}_{station_num} "
        station_geometry += f"{station_rot[0]:.6f} {station_rot[1]:.6f} {station_rot[2]:.6f} "
        station_geometry += f"{station_rot[3]:.6f} {station_rot[4]:.6f} {station_rot[5]:.6f} "
        station_geometry += f"{station_rot[6]:.6f} {station_rot[7]:.6f} {station_rot[8]:.6f}\n"
        
        # Place station in world
        station_geometry += f":PLACE {station_name} 1 world RM_Station_{wheel}_{sector}_{station_num} "
        station_geometry += f"{station_center[0]:.6f} {station_center[1]:.6f} {station_center[2]:.6f}\n\n"
        
        # Add yoke if requested
        if include_yoke:
            yoke_name = f"Yoke_W{wheel}_Sec{sector}_St{station_num}"
            # Yoke has same dimensions as station
            station_geometry += f"// Iron Yoke (radially below station)\n"
            station_geometry += f":VOLU {yoke_name} BOX {station_bounds[0]/2:.6f} {station_bounds[2]/2:.6f} {station_bounds[1]/2:.6f} G4_Fe\n"
            
            # Calculate yoke position: radially inward from station
            # Direction vector from origin (0,0,0) to station center in XY plane only
            station_pos_xy = np.array([station_center[0], station_center[1], 0.0])
            radial_distance_xy = np.linalg.norm(station_pos_xy)
            
            # Unit vector pointing radially outward in XY plane
            radial_unit_xy = station_pos_xy / radial_distance_xy
            
            # Yoke thickness (same as station height)
            yoke_thickness = station_bounds[1]  # height dimension
            
            # Yoke center: move inward in XY plane by station_height/2 + yoke_height/2
            # Keep same Z coordinate as station
            offset = (station_bounds[1]/2 + yoke_thickness/2)
            yoke_pos_xy = station_pos_xy - radial_unit_xy * offset
            yoke_pos = np.array([yoke_pos_xy[0], yoke_pos_xy[1], station_center[2]])
            
            # Yoke uses same rotation as station
            station_geometry += f":PLACE {yoke_name} 1 world RM_Station_{wheel}_{sector}_{station_num} "
            station_geometry += f"{yoke_pos[0]:.6f} {yoke_pos[1]:.6f} {yoke_pos[2]:.6f}\n\n"
        
        # Process all superlayers (1, 2, 3)
        # Place cells directly in station (no superlayer/layer volumes)
        for sl in station.super_layers:
            if sl.number not in [1, 2, 3]:
                continue
            
            station_geometry += f"// SuperLayer {sl.number} - Cells\n"
            
            # Define drift cell solid for this superlayer (dimensions may vary between superlayers)
            first_cell = None
            if len(sl.layers) > 0 and len(sl.layers[0].cells) > 0:
                first_cell = sl.layers[0].cells[0]
            
            cell_solid_name = f"DriftCellSolid_W{wheel}_Sec{sector}_St{station_num}_SL{sl.number}"
            if first_cell:
                cell_bounds = first_cell.bounds
                station_geometry += f":SOLID {cell_solid_name} BOX {cell_bounds[0]/2:.6f} {cell_bounds[2]/2:.6f} {cell_bounds[1]/2:.6f}\n"
            
            # SL2 is rotated, so it needs its own rotation matrix
            if sl.number == 2:
                # Get superlayer transformation relative to station
                sl_transform = sl.transformer.get_transformation(from_frame="SuperLayer", to_frame="Station")
                sl_rot = extract_rotation(sl_transform)
                
                # Create rotation matrix for SL2
                sl_rot_name = f"RM_SL2_W{wheel}_Sec{sector}_St{station_num}"
                station_geometry += f":ROTM {sl_rot_name} "
                station_geometry += f"{sl_rot[0]:.6f} {sl_rot[1]:.6f} {sl_rot[2]:.6f} "
                station_geometry += f"{sl_rot[3]:.6f} {sl_rot[4]:.6f} {sl_rot[5]:.6f} "
                station_geometry += f"{sl_rot[6]:.6f} {sl_rot[7]:.6f} {sl_rot[8]:.6f}\n"
                
                rotation_to_use = sl_rot_name
            else:
                # SL1 and SL3 use identity rotation
                rotation_to_use = "R0"
            
            # Process each layer - place all cells
            for layer in sl.layers:
                for cell in layer.cells:
                    cell_name = f"Cell_W{wheel}_Sec{sector}_St{station_num}_SL{sl.number}_L{layer.number}_C{cell.number}"
                    
                    # Define volume using the superlayer-specific solid
                    station_geometry += f":VOLU {cell_name} {cell_solid_name} GasMixture\n"
                    
                    # Get cell position relative to station (local coordinates)
                    cell_center = cell.local_center

                    if sl.number == 2:
                        cell_center = np.dot(sl_transform[:3, :3], cell.local_center)

                    # Place cell in station with appropriate rotation
                    station_geometry += f":PLACE {cell_name} 1 {station_name} {rotation_to_use} "
                    station_geometry += f"{cell_center[0]:.6f} {cell_center[1]:.6f} {cell_center[2]:.6f}\n"
                    break # Only one cell to test
            station_geometry += "\n"
        
        # Write this station's geometry to its own file
        with open(station_filepath, 'w') as f:
            f.write(station_geometry)
        
        print(f"  ✓ Generated: {station_filename}")
    
    # Create concentrator file with #include directives
    include_lines = "\n".join([f"#include geometry/stations/{fname}" for fname in station_files])
    
    # Read concentrator template
    with open(concentrator_template, 'r') as f:
        concentrator_content = f.read()
    
    # Replace {STATIONS} placeholder with include directives
    concentrator_content = concentrator_content.replace('{STATIONS}', include_lines)
    
    # Write concentrator file
    with open(concentrator_output, 'w') as f:
        f.write(concentrator_content)
    
    print(f"\n✓ Concentrator file: {concentrator_output}")
    print(f"  Stations: {len(stations_list)}")
    
    # Count total cells
    total_cells = 0
    for wheel, sector, station_num in stations_list:
        station = Station(wheel=wheel, sector=sector, station=station_num)
        for sl in station.super_layers:
            if sl.number in [1, 2, 3]:
                for layer in sl.layers:
                    total_cells += len(layer.cells)
    
    print(f"  Total cells: {total_cells}")

if __name__ == "__main__":
    # Define stations to include (wheel, sector, station)
    stations_to_generate = [
        (-1, 2, 1),  # MB1
        (-1, 2, 2),  # MB2
        (-1, 2, 3),  # MB3
        (-1, 2, 4),  # MB4
    ]
    
    print("Generating DT geometry files...")
    print(f"Stations to generate: {len(stations_to_generate)}")
    
    # Generate geometry
    generate_station_geometry_ascii(
        stations_list=stations_to_generate,
        output_dir='stations',
        concentrator_template='geometry_concentrator_template',
        concentrator_output='geometry_concentrator.tg',
        include_yoke=True  # Set to True to include iron yokes
    )
