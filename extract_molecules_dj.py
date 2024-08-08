import sys
import numpy as np
import os

class MoleculeExtractor:
    def __init__(self, gro_file, output_file):
        self.gro_file = gro_file
        self.output_file = output_file
        self.atoms = []
        self.box_dimensions = []

    def read_gro_file(self):
        if not os.path.exists(self.gro_file):
            raise FileNotFoundError(f"GRO file not found: {self.gro_file}")
        
        with open(self.gro_file, 'r') as f:
            lines = f.readlines()
        
        self.box_dimensions = np.array(lines[-1].strip().split(), dtype=float)
        for line in lines[2:-1]:  # Skip header, number of atoms, and box dimensions lines
            resid = int(line[0:5].strip())
            resname = line[5:10].strip()
            atomname = line[10:15].strip()
            atomnum = int(line[15:20].strip())
            x, y, z = map(float, (line[20:28], line[28:36], line[36:44]))
            self.atoms.append((resid, resname, atomname, atomnum, x, y, z))

    def get_distance(self, coord1, coord2):
        delta = np.abs(coord1 - coord2)
        delta = np.where(delta > 0.5 * self.box_dimensions, self.box_dimensions - delta, delta)
        return np.sqrt((delta ** 2).sum(axis=-1))

    def extract_molecules(self, source_resid, cutoff_distance):
        source_coords = np.array([atom[4:] for atom in self.atoms if atom[0] == source_resid])
        if not len(source_coords):
            raise ValueError(f"Source residue ID {source_resid} not found in the .gro file.")
        
        extracted_resids = set()
        source_atoms = set()
        
        for atom in self.atoms:
            resid, x, y, z = atom[0], atom[4], atom[5], atom[6]
            if resid == source_resid:
                source_atoms.add((x, y, z))
                continue  # Skip adding atoms from the source molecule
            
            dist = self.get_distance(np.array([x, y, z]), source_coords)
            if np.any(dist < cutoff_distance):
                extracted_resids.add(resid)
        
        with open(self.output_file, 'w') as f:
            f.write('Extracted molecules\n')
            f.write(f'{len([atom for atom in self.atoms if atom[0] in extracted_resids])}\n')
            for atom in self.atoms:
                if atom[0] in extracted_resids:
                    f.write(f'{atom[0]:5d}{atom[1]:>5}{atom[2]:>5}{atom[3]:5d}{atom[4]:8.3f}{atom[5]:8.3f}{atom[6]:8.3f}\n')
            f.write(f'{self.box_dimensions[0]:10.5f}{self.box_dimensions[1]:10.5f}{self.box_dimensions[2]:10.5f}\n')

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 extract_molecules_dj.py <gro_file> <output_file> <source_resid> <cutoff_distance>")
        sys.exit(1)
    
    gro_file = sys.argv[1]
    output_file = sys.argv[2]
    source_resid = int(sys.argv[3])
    cutoff_distance = float(sys.argv[4])
    
    extractor = MoleculeExtractor(gro_file, output_file)
    extractor.read_gro_file()
    extractor.extract_molecules(source_resid, cutoff_distance)
    print(f"Molecules extracted and saved to {output_file}")

if __name__ == "__main__":
    main()