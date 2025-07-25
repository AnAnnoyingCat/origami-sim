import json
import matplotlib.pyplot as plt
import os

# Hardcoded path to your .fold file
FOLD_FILE_PATH = "/home/cwernke/origami-sim/src/Python-Helperscripts/3d_print_maker/Objects/Crease_Robot/crease_robot_with_motor_and_flaps.fold"

# Color map for edge assignments
ASSIGNMENT_COLORS = {
    "M": "red",
    "V": "blue",
    "B": "black",
    "A": "lightblue",
    "F": "lightblue"
}

def main():
    # Load the fold file
    with open(FOLD_FILE_PATH, 'r') as f:
        data = json.load(f)

    vertices = data['vertices_coords']
    edges = data['edges_vertices']
    assignments = data.get('edges_assignment', [""] * len(edges))  # fallback to blank if missing

    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    
    for i, ((v1_idx, v2_idx), assign) in enumerate(zip(edges, assignments)):
        x1, y1 = vertices[v1_idx]
        x2, y2 = vertices[v2_idx]
        color = ASSIGNMENT_COLORS.get(assign, "gray")  # default to gray if unknown
        ax.plot([x1, x2], [y1, y2], color=color)

        # Label the edge
        mid_x = (x1 + x2) / 2
        mid_y = (y1 + y2) / 2
        ax.text(mid_x, mid_y, str(i), color='red', fontsize=8, ha='center', va='center')

    # Plot vertices
    for x, y in vertices:
        ax.plot(x, y, 'bo', markersize=3)

    ax.set_title("Edges with Labels and Assignments")
    plt.axis('off')
    plt.tight_layout()

    # Save output
    os.makedirs("Out", exist_ok=True)
    base_name = os.path.basename(FOLD_FILE_PATH).rsplit('.', 1)[0]
    output_path = os.path.join("Out", f"{base_name}.png")
    plt.savefig(output_path, dpi=300)
    print(f"Saved plot to: {output_path}")

if __name__ == "__main__":
    main()
