import xml.etree.ElementTree as ET
import math

# Load the XML file and read it line by line
xml_file = "/path/to/nan-errors-logs/iterationy-replicax-statez-state.xml"  # Change filename as needed
forces = []

with open(xml_file, "r") as f:
    particle_counter = 0
    for line in f:
        if "<Force" in line:  # Check if the line contains a force entry
            particle_counter += 1
            try:
                element = ET.fromstring(line.strip())  # Parse the XML element
                x, y, z = float(element.get("x")), float(element.get("y")), float(element.get("z"))
                magnitude = math.sqrt(x**2 + y**2 + z**2)
                forces.append((particle_counter, x, y, z, magnitude))
            except ET.ParseError:
                pass  # Ignore malformed lines (should not happen if XML is valid)

# Compute the average
force_average = sum([f[4] for f in forces])/len(forces)
print(f"Average force experienced by a particle: {round(force_average, 2)}\n")

# Sort by magnitude in descending order
forces.sort(key=lambda f: f[4], reverse=True)

# Top N forces (change N as needed)
N = 50
print(f"Top {N} forces by magnitude:")
for i, (particle_count, x, y, z, mag) in enumerate(forces[:N], 1):
    print(f"{i}: Particle {particle_count} -> Force: {mag:.4f}")

