from vpython import *
import random

# Define properties for all 20 standard amino acids
amino_acids_properties = {
    "ala": vec(0.8, 0.8, 0.8),  # Alanine
    "arg": vec(0, 1, 1),        # Arginine
    "asn": vec(0.5, 0.5, 1),    # Asparagine
    "asp": vec(1, 0, 0),        # Aspartic Acid
    "cys": vec(1, 1, 0),        # Cysteine
    "gln": vec(0.5, 0.5, 0.5),  # Glutamine
    "glu": vec(1, 0.5, 0),      # Glutamic Acid
    "gly": vec(0.5, 0.5, 0.5),  # Glycine
    "his": vec(0.5, 0.5, 0.5),  # Histidine
    "ile": vec(0.5, 0.5, 0.5),  # Isoleucine
    "leu": vec(0.5, 0.5, 0.5),  # Leucine
    "lys": vec(0.5, 0.5, 0.5),  # Lysine
    "met": vec(0.5, 0.5, 0.5),  # Methionine
    "phe": vec(0.5, 0.5, 0.5),  # Phenylalanine
    "pro": vec(0.5, 0.5, 0.5),  # Proline
    "ser": vec(0.5, 0.5, 0.5),  # Serine
    "thr": vec(0.5, 0.5, 0.5),  # Threonine
    "trp": vec(0.5, 0.5, 0.5),  # Tryptophan
    "tyr": vec(0.5, 0.5, 0.5),  # Tyrosine
    "val": vec(0.5, 0.5, 0.5)   # Valine
}

# Define hydrophilic and hydrophobic amino acids
hydrophilic_amino_acids = {"arg", "asn", "asp", "gln", "glu", "his", "lys", "ser", "thr", "tyr"}
hydrophobic_amino_acids = {"ala", "cys", "gly", "ile", "leu", "met", "phe", "pro", "trp", "val"}

# Map single-letter codes to three-letter amino acid codes
single_to_three_letter = {
    'A': 'ala', 'R': 'arg', 'N': 'asn', 'D': 'asp', 'C': 'cys', 'E': 'glu', 'Q': 'gln', 'G': 'gly',
    'H': 'his', 'I': 'ile', 'L': 'leu', 'K': 'lys', 'M': 'met', 'F': 'phe', 'P': 'pro', 'S': 'ser',
    'T': 'thr', 'W': 'trp', 'Y': 'tyr', 'V': 'val'
}

# Convert a large string of single-letter amino acids into a list of three-letter codes
def convert_string_to_amino_acid_list(amino_acid_string):
    amino_acid_list = []
    for letter in amino_acid_string:
        amino_acid = single_to_three_letter[letter]
        amino_acid_list.append(amino_acid)
    return amino_acid_list

# Example input string
amino_acid_string = "HSGRVHYEAALQEIDSDFDDGIIKYTYQLG"

# Convert the string to a list of amino acids
amino_acid_chain = convert_string_to_amino_acid_list(amino_acid_string)

# Function to initialize the amino acid chain
def initialize_amino_chain(amino_acid_chain):
    amino_pos_list = []
    amino_spheres = []
    for counter, amino in enumerate(amino_acid_chain):
        amino_color = amino_acids_properties[amino]
        amino_init = vec(3 * counter, counter^2, counter^3)
        amino_pos_list.append(amino_init)
        amino_sphere = sphere(pos=amino_init, radius=1, color=amino_color)
        amino_spheres.append(amino_sphere)
    return amino_pos_list, amino_spheres

# Initialize the amino acid chain and get positions
amino_pos_list, amino_spheres = initialize_amino_chain(amino_acid_chain)

# Function to initialize bonds between amino acids
def initialize_bonds(amino_pos_list):
    bonds = []
    for i in range(len(amino_pos_list) - 1):
        start = amino_pos_list[i]
        end = amino_pos_list[i + 1]
        bond = cylinder(pos=start, axis=end - start, color=color.white, radius=0.1)
        bonds.append(bond)
    return bonds

# Initialize bonds
bonds = initialize_bonds(amino_pos_list)

# Function to calculate the center of mass of the molecule
def calculate_center_of_mass(amino_spheres):
    total_mass = len(amino_spheres)
    center_of_mass = vec(0, 0, 0)
    for sphere in amino_spheres:
        center_of_mass += sphere.pos
    center_of_mass /= total_mass
    return center_of_mass

# Function to apply forces based on hydrophilic and hydrophobic properties
def apply_forces(amino_acid_chain, amino_spheres, bonds):
    center_of_mass = calculate_center_of_mass(amino_spheres)
    correction_force = -center_of_mass * 0.1  # Correction force to keep the molecule centered
    for i, amino in enumerate(amino_acid_chain):
        amino_sphere = amino_spheres[i]
        initial_pos = amino_sphere.pos
        if amino in hydrophilic_amino_acids:
            for j in range(len(bonds)):
                bonds[j].pos = amino_spheres[j].pos
                bonds[j].axis = amino_spheres[j + 1].pos - amino_spheres[j].pos
                # Ensure bonds do not stretch more than 12 units or compress less than 8 units
                if mag(bonds[j].axis) > 3 or mag(bonds[j].axis) < 4:
                    bonds[j].axis = norm(bonds[j].axis) * 4
                    amino_spheres[j + 1].pos = bonds[j].pos + bonds[j].axis
            # Move hydrophilic amino acids away from the center
            force = norm(amino_sphere.pos - center_of_mass)
            amino_sphere.pos += force
        elif amino in hydrophobic_amino_acids:
            # Move hydrophobic amino acids towards the center
            for j in range(len(bonds)):
                bonds[j].pos = amino_spheres[j].pos
                bonds[j].axis = amino_spheres[j + 1].pos - amino_spheres[j].pos
                # Ensure bonds do not stretch more than 12 units or compress less than 8 units
                if mag(bonds[j].axis) > 3 or mag(bonds[j].axis) < 4:
                    bonds[j].axis = norm(bonds[j].axis) * 4
                    amino_spheres[j + 1].pos = bonds[j].pos + bonds[j].axis
            # Move hydrophobic amino acids towards the center
            force = norm(center_of_mass - amino_sphere.pos)
            amino_sphere.pos += force
        # Apply correction force to keep the molecule centered
        amino_sphere.pos += correction_force
    # Apply repulsive force if amino acids are within 1 unit of each other
    for i in range(len(amino_spheres)):
        for j in range(i + 1, len(amino_spheres)):
            distance = mag(amino_spheres[i].pos - amino_spheres[j].pos)
            min_distance = amino_spheres[i].radius + amino_spheres[j].radius
            if distance < min_distance:
                repulsive_force = norm(amino_spheres[i].pos - amino_spheres[j].pos) * (min_distance - distance)
                amino_spheres[i].pos += repulsive_force
                amino_spheres[j].pos -= repulsive_force

# Animation loop
while True:
    rate(20)
    apply_forces(amino_acid_chain, amino_spheres, bonds)