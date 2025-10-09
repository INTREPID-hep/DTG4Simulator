#!/usr/bin/env python3
"""
Generate C++ geometry data header from mplDTs Python package
"""
import sys
sys.path.insert(0, '/home/destrada/mplDTs') # Adjust this path as needed

from mpldts.geometry import Station
from pathlib import Path

def extract_station_geometry(wheel, sector, station_num):
    """Extract geometry data from mplDTs Station"""
    try:
        station = Station(wheel=wheel, sector=sector, station=station_num)
        print(f"  Processing: {station.name}")
        
        station_data = {
            'wheel': wheel,
            'sector': sector,
            'station': station_num,
            'bounds': station.bounds,
            'center': station.local_center,
            'superlayers': []
        }
        
        for sl in station.super_layers:
            sl_data = {
                'number': sl.number,
                'bounds': sl.bounds,
                'center': sl.local_center,
                'rotation': 90.0 if sl.number == 2 else 0.0,  # SL2 is rotated
                'cells': []  # Flatten: all cells directly under SuperLayer
            }
            
            # Flatten: collect all cells from all layers
            for layer in sl.layers:
                for cell in layer.cells:
                    cell_data = {
                        'number': cell.number,
                        'layer': layer.number,  # Keep layer number for reference
                        'center': cell.local_center,
                        'bounds': cell.bounds
                    }
                    sl_data['cells'].append(cell_data)
            
            station_data['superlayers'].append(sl_data)
        
        return station_data
    except Exception as e:
        print(f"  WARNING: Could not load Wh={wheel} Sec={sector} St={station_num}: {e}")
        return None

def generate_geometry_data_header(stations_to_include, output_file, template_file):
    """
    Generate C++ header with geometry data initialization
    
    Args:
        stations_to_include: List of (wheel, sector, station) tuples
        output_file: Output file path
        template_file: Template file path
    """
    print("Extracting geometry data from mplDTs...")
    
    # Extract all station geometries
    all_stations = []
    for wheel, sector, station_num in stations_to_include:
        station_data = extract_station_geometry(wheel, sector, station_num)
        if station_data:
            all_stations.append(station_data)
    
    print(f"\nSuccessfully extracted {len(all_stations)} station geometries")
    
    # Generate station data code
    station_data_code = ""
    
    for st_data in all_stations:
        w, sec, st = st_data['wheel'], st_data['sector'], st_data['station']
        
        station_data_code += f"""
    // ===== Station: Wheel={w}, Sector={sec}, Station={st} =====
    stGeom = StationGeometry();
    stGeom.wheel = {w};
    stGeom.sector = {sec};
    stGeom.station = {st};
    stGeom.centerX = {st_data['center'][0]:.6f};
    stGeom.centerY = {st_data['center'][1]:.6f};
    stGeom.centerZ = {st_data['center'][2]:.6f};
    stGeom.width = {st_data['bounds'][0]:.6f};
    stGeom.height = {st_data['bounds'][1]:.6f};
    stGeom.length = {st_data['bounds'][2]:.6f};
"""
        
        # Generate superlayer data
        for sl_data in st_data['superlayers']:
            sl_num = sl_data['number']
            station_data_code += f"""
    // SuperLayer {sl_num}
    slGeom = SuperLayerGeometry();
    slGeom.slNumber = {sl_num};
    slGeom.centerX = {sl_data['center'][0]:.6f};
    slGeom.centerY = {sl_data['center'][1]:.6f};
    slGeom.centerZ = {sl_data['center'][2]:.6f};
    slGeom.width = {sl_data['bounds'][0]:.6f};
    slGeom.height = {sl_data['bounds'][1]:.6f};
    slGeom.length = {sl_data['bounds'][2]:.6f};
    slGeom.rotation = {sl_data['rotation']:.1f};
"""
            
            # Generate all cells directly under SuperLayer
            for cell_data in sl_data['cells']:
                cell_num = cell_data['number']
                layer_num = cell_data['layer']
                cx, cy, cz = cell_data['center']
                cw, ch, cl = cell_data['bounds']
                
                station_data_code += f"""    cellGeom.cellNumber = {cell_num};
    cellGeom.layerNumber = {layer_num};
    cellGeom.x = {cx:.6f}; cellGeom.y = {cy:.6f}; cellGeom.z = {cz:.6f};
    cellGeom.width = {cw:.6f}; cellGeom.height = {ch:.6f}; cellGeom.length = {cl:.6f};
    slGeom.cells.push_back(cellGeom);
"""
            
            station_data_code += f"""    stGeom.superlayers.push_back(slGeom);
"""
        
        station_data_code += f"""
    fGeometryDB[MakeKey({w}, {sec}, {st})] = stGeom;
"""
    
    # Read template file
    template_path = Path(template_file)
    if not template_path.exists():
        raise FileNotFoundError(f"Template file not found: {template_file}")
    
    with open(template_path, 'r') as f:
        template = f.read()
    
    # Replace placeholder with generated station data
    output_content = template.replace('{STATION_DATA}', station_data_code)
    
    # Write to file
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        f.write(output_content)
    
    print(f"\nGenerated: {output_path}")
    print(f"File size: {output_path.stat().st_size / 1024:.1f} KB")
    
    # Estimate memory usage
    n_stations = len(all_stations)
    n_cells_total = sum(
        len(sl['cells']) 
        for st in all_stations 
        for sl in st['superlayers']
    )
    
    # Rough estimate: ~100 bytes per cell structure
    estimated_memory_kb = (n_cells_total * 100) / 1024
    
    print(f"\nEstimated runtime memory usage: ~{estimated_memory_kb:.1f} KB")
    print(f"({n_stations} stations, ~{n_cells_total} total cells)")

if __name__ == "__main__":
    # Define which stations to include
    # You can expand this list as needed
    stations_to_include = [
        # Format: (wheel, sector, station)
        (-1, 1, 1),  # MB1
        (-1, 1, 2),  # MB2
        (-1, 1, 3),  # MB3
        # Add more as needed:
        # (-1, 1, 4),
        # (-1, 2, 1),
        # etc.
    ]
    
    template_file = "templates/DTGeometryData.hh.template"
    output_file = "include/DTGeometryData.hh"
    
    generate_geometry_data_header(stations_to_include, output_file, template_file)
    
    print("\n✓ Done! Include this header in DTGeometryBuilder.cc")
