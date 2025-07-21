import json
import math
import argparse

def generate_simulation_data(creases, start_time, total_time):
    simulation_data = {}
    for t in range(start_time, start_time + total_time + 1):
        time_entry = []
        for crease in creases:
            amplitude = (crease["max_angle"] - crease["min_angle"]) / 2
            offset = crease["min_angle"] + amplitude
            phase_radians = math.radians(crease.get("phase_offset", 0))  # Default 0°
            angle = offset + amplitude * math.sin(2 * math.pi * (t - start_time) / crease["period"] + phase_radians)
            time_entry.append([crease["name"], round(angle, 3)])
        simulation_data[str(t)] = time_entry
    return simulation_data

def main():
    parser = argparse.ArgumentParser(description="Generate JSON input for origami fold oscillation.")
    parser.add_argument("--start", type=int, default=15, help="Start time (default: 0)")
    parser.add_argument("--duration", type=int, default=2000, help="Duration in time units (default: 1000)")
    parser.add_argument("--output", type=str, default="origami_oscillation_input.json", help="Output file name")

    args = parser.parse_args()

    # Define your crease waveforms with int names and optional phase offsets in degrees
    creases = [
        {"name": 32, "min_angle": -60, "max_angle": -5, "period": 200, "phase_offset": 0},
        {"name": 34, "min_angle": -60, "max_angle": -5, "period": 200, "phase_offset": 180},
    ]

    data = generate_simulation_data(creases, args.start, args.duration)

    with open(args.output, "w") as f:
        json.dump(data, f, indent=2)

    print(f"JSON file '{args.output}' generated successfully from t={args.start} to t={args.start + args.duration}.")

if __name__ == "__main__":
    main()
