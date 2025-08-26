import matplotlib.pyplot as plt
import os
# DISCLAIMER: Written by ChatGPT, not by me

def plot_strain(file_path):
    frames = []
    strains = []
    avg_strain = None

    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                # Extract average strain from the first line
                if 'Average strain' in line:
                    try:
                        avg_strain = float(line.strip().split(':')[-1])
                    except ValueError:
                        pass
                continue
            parts = line.strip().split()
            if len(parts) == 2:
                frames.append(int(parts[0]))
                strains.append(float(parts[1]) * 100)  # Convert to percent

    plt.figure(figsize=(6, 6))
    plt.plot(frames, strains, marker='o', linewidth=2, color='black', markersize=0)  # Thinner line

    # Label peak point
    max_strain = max(strains)
    max_index = strains.index(max_strain)
    plt.annotate(f'Peak: {max_strain:.2f}%', 
                 xy=(frames[max_index], max_strain), 
                 xytext=(frames[max_index], max_strain + 0.8),
                 arrowprops=dict(arrowstyle='->', lw=0.8),
                 fontsize=9)

    # Add average strain annotation
    if avg_strain is not None:
        plt.text(0.98, 0.02, f'Avg Strain: {avg_strain * 100:.2f}%', 
                 ha='right', va='bottom', 
                 transform=plt.gca().transAxes, fontsize=9, 
                 bbox=dict(facecolor='white', alpha=0.6, edgecolor='gray'))

    plt.xlabel('Frame')
    plt.ylabel('Strain (%)')
    plt.title('Strain Over Time')
    plt.grid(True)
    plt.tight_layout()
    
    save_path = os.path.splitext(file_path)[0] + '.png'
    plt.savefig(save_path, dpi=300)
    plt.show()

if __name__ == "__main__":
    file_path = "/home/cwernke/origami-sim/data/simulation_strains/hypar_strain.txt"
    plot_strain(file_path)
