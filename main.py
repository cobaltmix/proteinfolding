from vpython import sphere, vector, color, scene, pi, cos, sin, rate

# Set up the scene
scene.background = color.white

# Parameters
lipid_radius = 0.05  # Radius of lipid head spheres
sphere_radius = 2  # Radius of the entire spherical membrane
num_layers = 180  # Number of layers or "latitude bands" of lipids
num_lipids_per_layer = 180  # Number of lipids per layer or "longitude lines"

# Create the spherical lipid bilayer
for layer in range(num_layers):
    theta = pi * layer / (num_layers - 1)  # Latitude angle
    for i in range(num_lipids_per_layer):
        phi = 2 * pi * i / num_lipids_per_layer  # Longitude angle
        
        # Calculate x, y, z coordinates on the sphere surface
        x = sphere_radius * sin(theta) * cos(phi)
        y = sphere_radius * sin(theta) * sin(phi)
        z = sphere_radius * cos(theta)
        
        # Place spheres on the top hemisphere
        sphere(pos=vector(x, y, z), radius=lipid_radius, color=color.orange)

# Adjust the camera view
scene.camera.pos = vector(0, 0, sphere_radius*3)
scene.camera.axis = vector(0, 0, -sphere_radius)

while True:
    rate(10)