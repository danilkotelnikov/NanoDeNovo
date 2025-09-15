import os
import tempfile
from openmm import app, unit, LocalEnergyMinimizer, VerletIntegrator
from openmm.app import PDBFile, Modeller
from pdbfixer.pdbfixer import PDBFixer

# Get input and output directories interactively
input_dir = input("Enter the path to the input directory containing PDB files: ").strip()
output_dir = input("Enter the path to the output directory: ").strip()
os.makedirs(output_dir, exist_ok=True)

# Get list of PDB files in input directory
pdb_files = [f for f in os.listdir(input_dir) if f.endswith('.pdb')]

for pdb_file in pdb_files:
    input_path = os.path.join(input_dir, pdb_file)
    complex_name = os.path.splitext(pdb_file)[0]
    output_pdb = os.path.join(output_dir, f"{complex_name}_em.pdb")
    log_path = os.path.join(output_dir, f"{complex_name}_em.log")

    with open(log_path, 'w') as log_file:
        log_file.write(f"Minimization log for {complex_name}\n")
        log_file.write("===================================\n")

        try:
            # Use PDBFixer to fix terminal groups
            log_file.write("Using PDBFixer to add missing terminal groups\n")
            fixer = PDBFixer(filename=input_path)
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            # Write the fixed structure to a temporary file
            import uuid
            fixed_pdb_path = os.path.join('/tmp', f'fixed_{complex_name}_{uuid.uuid4()}.pdb')
            with open(fixed_pdb_path, 'w') as output_file:
                app.PDBFile.writeFile(fixer.topology, fixer.positions, output_file)
            
            # Load the fixed PDB file
            pdb = app.PDBFile(fixed_pdb_path)
            os.unlink(fixed_pdb_path)  # Clean up temporary file
            
            # Prepare forcefield and modeller
            forcefield = app.ForceField('amber14/protein.ff14SB.xml')
            modeller = Modeller(pdb.topology, pdb.positions)
            
            # Add hydrogens and extra particles
            modeller.addHydrogens(forcefield)
            modeller.addExtraParticles(forcefield)
            
            # Create system and simulation
            system = forcefield.createSystem(modeller.topology)
            integrator = VerletIntegrator(0.001 * unit.femtoseconds)
            simulation = app.Simulation(modeller.topology, system, integrator)
            simulation.context.setPositions(modeller.positions)

            # Get initial energy
            state = simulation.context.getState(getEnergy=True)
            initial_energy = state.getPotentialEnergy()
            log_file.write(f"Initial potential energy: {initial_energy}\n")
            print(f"\nProcessing {complex_name}")
            print(f"Initial energy: {initial_energy}")

            # Minimize in chunks to track progress
            # Minimize with convergence checking
            max_iterations = 1000
            energy_tolerance = 10.0  # kJ/mol - energy change threshold for convergence
            previous_energy = initial_energy
            total_steps = 0
            
            for iteration in range(max_iterations):
                # Run minimization for a chunk of steps
                LocalEnergyMinimizer.minimize(simulation.context, maxIterations=100)
                total_steps += 100
                
                # Get current energy
                current_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
                energy_change = abs((current_energy - previous_energy) / unit.kilojoule_per_mole)
                
                # Log progress
                log_file.write(f"After {total_steps} steps: Energy = {current_energy}, Change = {energy_change} kJ/mol\n")
                print(f"Steps: {total_steps} \t Current energy: {current_energy} \t Change: {energy_change} kJ/mol")
                
                # Check for convergence
                if energy_change < energy_tolerance:
                    log_file.write(f"Minimization converged (energy change < {energy_tolerance} kJ/mol)\n")
                    print(f"Minimization converged (energy change < {energy_tolerance} kJ/mol)")
                    break
                    
                previous_energy = current_energy
            # Get final energy and positions
            final_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
            log_file.write(f"\nMinimization completed in {total_steps} steps\n")
            log_file.write(f"Final potential energy: {final_energy}\n")
            print(f"Minimization completed in {total_steps} steps")
            print(f"Final energy: {final_energy}")

            # Save minimized structure
            positions = simulation.context.getState(getPositions=True).getPositions()
            with open(output_pdb, 'w') as output_file:
                PDBFile.writeFile(simulation.topology, positions, output_file)
            
            print(f"Saved minimized structure to {output_pdb}")

        except Exception as e:
            log_file.write(f"\nError processing {complex_name}: {str(e)}\n")
            print(f"\nError processing {complex_name}: {str(e)}")
            continue

print("\nAll complexes processed successfully!")
