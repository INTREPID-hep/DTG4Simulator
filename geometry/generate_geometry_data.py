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

def generate_station_geometry_ascii(stations_list, template_file='templates/geometry_template.tg', 
                                  output_file='geometry/dt_geometry.tg'):
    """
    Generate Geant4 ASCII text geometry for multiple DT stations
    
    Args:
        stations_list: List of (wheel, sector, station_num) tuples
        template_file: Path to template file with fixed geometry parts
        output_file: Output geometry file
    """
    
    # Read template
    with open(template_file, 'r') as f:
        template_content = f.read()
    
    # Generate dynamic geometry for each station
    dynamic_geometry = ""
    
    for wheel, sector, station_num in stations_list:
        station = Station(wheel=wheel, sector=sector, station=station_num)
        
        dynamic_geometry += f"\n// ===== Station: Wheel={wheel}, Sector={sector}, Station={station_num} =====\n"
        dynamic_geometry += f"// {station.name}\n\n"
        
        station_bounds = station.bounds
        station_name = f"Station_W{wheel}_Sec{sector}_St{station_num}"
        dynamic_geometry += f":VOLU {station_name} BOX {station_bounds[0]/2:.6f} {station_bounds[2]/2:.6f} {station_bounds[1]/2:.6f} G4_AIR\n"
        
        # Get station transformation
        station_transform = station.transformer.get_transformation(from_frame="Station", to_frame="CMS")
        station_rot, _ = extract_rotation_translation(station_transform)
        station_trans = station.global_center
        
        # Create rotation matrix for station
        dynamic_geometry += f":ROTM RM_Station_{wheel}_{sector}_{station_num} "
        dynamic_geometry += f"{station_rot[0]:.6f} {station_rot[1]:.6f} {station_rot[2]:.6f} "
        dynamic_geometry += f"{station_rot[3]:.6f} {station_rot[4]:.6f} {station_rot[5]:.6f} "
        dynamic_geometry += f"{station_rot[6]:.6f} {station_rot[7]:.6f} {station_rot[8]:.6f}\n"
        
        # Place station in world
        dynamic_geometry += f":PLACE {station_name} 1 world RM_Station_{wheel}_{sector}_{station_num} "
        dynamic_geometry += f"{station_trans[0]:.6f} {station_trans[1]:.6f} {station_trans[2]:.6f}\n\n"
        
        # Process all superlayers (1, 2, 3)
        # Place cells directly in station (no superlayer/layer volumes)
        for sl in station.super_layers:
            if sl.number not in [1, 2, 3]:
                continue
            
            dynamic_geometry += f"// SuperLayer {sl.number} - Cells\n"
            
            # Define drift cell solid for this superlayer (dimensions may vary between superlayers)
            first_cell = None
            if len(sl.layers) > 0 and len(sl.layers[0].cells) > 0:
                first_cell = sl.layers[0].cells[0]
            
            cell_solid_name = f"DriftCellSolid_W{wheel}_Sec{sector}_St{station_num}_SL{sl.number}"
            if first_cell:
                cell_bounds = first_cell.bounds
                dynamic_geometry += f":SOLID {cell_solid_name} BOX {cell_bounds[0]/2:.6f} {cell_bounds[2]/2:.6f} {cell_bounds[1]/2:.6f}\n"
            
            # SL2 is rotated, so it needs its own rotation matrix
            if sl.number == 2:
                # Get superlayer transformation relative to station
                sl_transform = sl.transformer.get_transformation(from_frame="SuperLayer", to_frame="Station")
                sl_rot, _ = extract_rotation_translation(sl_transform)
                
                # Create rotation matrix for SL2
                sl_rot_name = f"RM_SL2_W{wheel}_Sec{sector}_St{station_num}"
                dynamic_geometry += f":ROTM {sl_rot_name} "
                dynamic_geometry += f"{sl_rot[0]:.6f} {sl_rot[1]:.6f} {sl_rot[2]:.6f} "
                dynamic_geometry += f"{sl_rot[3]:.6f} {sl_rot[4]:.6f} {sl_rot[5]:.6f} "
                dynamic_geometry += f"{sl_rot[6]:.6f} {sl_rot[7]:.6f} {sl_rot[8]:.6f}\n"
                
                rotation_to_use = sl_rot_name
            else:
                # SL1 and SL3 use identity rotation
                rotation_to_use = "R0"
            
            # Process each layer - only place FIRST cell for testing
            for layer in sl.layers:
                # Get first cell only
                for cell in layer.cells:
                    cell_name = f"Cell_W{wheel}_Sec{sector}_St{station_num}_SL{sl.number}_L{layer.number}_C{cell.number}"
                    
                    # Define volume using the station-specific solid
                    dynamic_geometry += f":VOLU {cell_name} {cell_solid_name} GasMixture\n"
                    
                    # Get cell position relative to station (local coordinates)
                    cell_center = cell.local_center

                    if sl.number == 2:
                        cell_center = np.dot(sl_transform[:3, :3], cell.local_center)

                    # Place cell in station with appropriate rotation
                    dynamic_geometry += f":PLACE {cell_name} 1 {station_name} {rotation_to_use} "
                    dynamic_geometry += f"{cell_center[0]:.6f} {cell_center[1]:.6f} {cell_center[2]:.6f}\n"
                    # break  # Only first cell
            dynamic_geometry += "\n"
    
    # Write combined output
    with open(output_file, 'w') as f:
        f.write(template_content)
        f.write(dynamic_geometry)
    
    print(f"✓ Geometry written to {output_file}")
    print(f"  Stations: {len(stations_list)}")
    
    # Count total cells
    total_cells = 0
    for wheel, sector, station_num in stations_list:
        station = Station(wheel=wheel, sector=sector, station=station_num)
        for sl in station.super_layers:
            if sl.number in [1, 2, 3]:
                total_cells += len(sl.layers.cells)
    
    print(f"  Total cells placed: {total_cells} (first cell of each layer)")

if __name__ == "__main__":
    # Define stations to include (wheel, sector, station)
    stations_to_generate = [
        (-1, 2, 2),  # MB1
        (-1, 2, 1),  # MB2
        (-1, 2, 3),  # MB3
        (-1, 2, 4),  # MB4
    ]
    
    # Generate geometry
    generate_station_geometry_ascii(
        stations_list=stations_to_generate,
        template_file='templates/geometry_template.txt',
        output_file='dt_geometry.txt'
    )
