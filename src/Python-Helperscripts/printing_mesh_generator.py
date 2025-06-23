# Disclaimer: Written by ChatGPT, not by me. 
import json
import numpy as np
from shapely.geometry import Polygon, LineString, MultiLineString
from shapely.ops import unary_union
import trimesh
import sys
import matplotlib.pyplot as plt

def load_fold(filepath, target_size=100.0):
    with open(filepath, 'r') as f:
        data = json.load(f)

    vertices = np.array(data['vertices_coords'])
    edges = np.array(data['edges_vertices'])
    assignments = data['edges_assignment']

    # Scale vertices so bounding box becomes target_size x target_size
    min_xy = vertices.min(axis=0)
    max_xy = vertices.max(axis=0)
    bbox_size = max_xy - min_xy

    scale_factor = target_size / max(bbox_size)
    scaled_vertices = (vertices - min_xy) * scale_factor

    return scaled_vertices, edges, assignments


def extract_geometry(vertices, edges, assignments):
    boundary_coords = []
    crease_edges = []

    for idx, (v1, v2) in enumerate(edges):
        if assignments[idx] == "B":
            boundary_coords.append((vertices[v1], vertices[v2]))
        elif assignments[idx] in ("M", "V"):
            crease_edges.append(LineString([vertices[v1], vertices[v2]]))

    # Reconstruct the outer boundary polygon directly
    # Start with all segments and stitch them together
    all_boundary_lines = [LineString(seg) for seg in boundary_coords]
    boundary_outline = unary_union(all_boundary_lines)

    # Let shapely try to polygonize it
    from shapely.ops import polygonize
    polygons = list(polygonize(boundary_outline))

    if len(polygons) == 0:
        raise ValueError("Could not reconstruct boundary polygon.")

    boundary_polygon = polygons[0]  # pick the biggest if multiple?

    creases = unary_union(crease_edges)
    return boundary_polygon, creases

def debug_plot(boundary, creases, cleared_polygon):
    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    if boundary:
        x, y = boundary.exterior.xy
        ax.plot(x, y, label='Boundary')

    if isinstance(creases, LineString):
        crease_list = [creases]
    elif isinstance(creases, MultiLineString):
        crease_list = list(creases.geoms)
    else:
        crease_list = []

    for c in crease_list:
        x, y = c.xy
        ax.plot(x, y, color='red', linewidth=1, label='Crease')


    if cleared_polygon:
        if hasattr(cleared_polygon, "geoms"):
            for poly in cleared_polygon.geoms:
                x, y = poly.exterior.xy
                ax.plot(x, y, color='green', linestyle='--', label='PLA area')
        else:
            x, y = cleared_polygon.exterior.xy
            ax.plot(x, y, color='green', linestyle='--', label='PLA area')

    ax.legend()
    plt.show()


def create_layer_from_polygon(polygon, z_min, z_max):
    height = z_max - z_min

    if polygon.is_empty:
        return None

    if polygon.geom_type == 'Polygon':
        polys = [polygon]
    elif polygon.geom_type == 'MultiPolygon':
        polys = list(polygon.geoms)
    else:
        raise ValueError(f"Unexpected geometry type: {polygon.geom_type}")

    meshes = []
    for poly in polys:
        try:
            mesh = trimesh.creation.extrude_polygon(poly, height)
            mesh.apply_translation([0, 0, z_min])
            meshes.append(mesh)
        except Exception as e:
            print(f"Failed to extrude polygon: {e}")

    if not meshes:
        return None

    return trimesh.util.concatenate(meshes)


def subtract_creases_from_boundary(boundary, creases, buffer_mm=1.0):
    crease_buffer = creases.buffer(buffer_mm, cap_style=2)
    cleared_polygon = boundary.difference(crease_buffer)
    return cleared_polygon

def save_mesh(mesh, filename):
    mesh.export(filename)

def main():
    # Modify this path to your actual .fold file
    fold_path = "/home/cwernke/BA/origami-sim/data/crease_patterns/defaultsquare.fold"
    vertices, edges, assignments = load_fold(fold_path, target_size=50)
    boundary, creases = extract_geometry(vertices, edges, assignments)
    

    # TPU Layer (middle)
    tpu_mesh = create_layer_from_polygon(boundary, 1.0, 1.8)
    save_mesh(tpu_mesh, "tpu_layer.obj")

    # PLA Layers
    cleared_polygon = subtract_creases_from_boundary(boundary, creases, buffer_mm=1.5)
    debug_plot(boundary, creases, cleared_polygon)
    if cleared_polygon.is_empty:
        print("Warning: PLA geometry is empty after crease subtraction.")
        return

    # Bottom PLA
    pla_bottom = create_layer_from_polygon(cleared_polygon, 0.0, 1.0)
    if pla_bottom is not None:
        save_mesh(pla_bottom, "pla_bottom.obj")
    else:
        print("Warning: PLA bottom layer extrusion failed.")

    # Top PLA
    pla_top = create_layer_from_polygon(cleared_polygon, 1.8, 2.8)
    if pla_top is not None:
        save_mesh(pla_top, "pla_top.obj")
    else:
        print("Warning: PLA top layer extrusion failed.")


if __name__ == "__main__":
    main()
