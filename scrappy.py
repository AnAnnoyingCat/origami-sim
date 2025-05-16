# Set how many temporary variables you want to declare
N = 260  # Change to however many you need

# Generate and print the line
print('double ' + ', '.join(f't{i}' for i in range(1, N + 1)) + ' = 0.0;')