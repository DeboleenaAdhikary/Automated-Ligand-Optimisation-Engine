#!/usr/bin/env python3
"""
Enzyme-Ligand Processor with Active Site and Substitution Point Selection
Includes BRICS fragmentation and manual atom marking
"""

import os
import sys
import subprocess
import time
import json
from pathlib import Path
import re
import csv
from flask import Flask, request, jsonify, render_template, send_from_directory, abort
from werkzeug.utils import secure_filename
import tempfile
import json
import threading
import uuid
import io
from contextlib import redirect_stdout, redirect_stderr
try:
    from openbabel import openbabel as ob
    from Bio import PDB
    from Bio.PDB import PDBIO, Select
    from rdkit import Chem
    from rdkit.Chem import BRICS, Descriptors, AllChem, Draw
    from rdkit import RDLogger
except ImportError as e:
    print(f"Missing package: {e}")
    print("\nInstall required packages:")
    print("pip install openbabel-wheel biopython rdkit pillow --break-system-packages")
    sys.exit(1)

# ── Silence noisy warnings from OpenBabel (UFFTYPER etc.) and RDKit ──────────
# OpenBabel: only show actual errors, not warnings/info (suppresses UFFTYPER
# "Unrecognized atom type: *_" messages for wildcard/dummy atoms)
ob.obErrorLog.SetOutputLevel(ob.obError)
# RDKit: suppress rdApp.warning / rdApp.error chatter (e.g. sanitization notes)
RDLogger.DisableLog('rdApp.*')


class ChainAProteinSelect(Select):
    """Select only chain A and protein residues (no heteroatoms)."""
    
    def __init__(self):
        self.removed_hetero = []
        self.removed_chains = []
    
    def accept_chain(self, chain):
        if chain.id == 'A':
            return 1
        else:
            self.removed_chains.append(chain.id)
            return 0
    
    def accept_residue(self, residue):
        hetero_flag = residue.id[0]
        if hetero_flag == ' ':
            return 1
        else:
            self.removed_hetero.append({
                'name': residue.resname,
                'id': residue.id[1],
                'chain': residue.parent.id
            })
            return 0


class SubstitutionJob:
    """Represents a single substitution job with marked positions."""
    
    def __init__(self, job_id, smiles_with_attachment, attachment_atoms, active_site_atoms):
        self.job_id = job_id
        self.smiles_with_attachment = smiles_with_attachment
        self.attachment_atoms = attachment_atoms  # List of atom indices with [*]
        self.active_site_atoms = active_site_atoms
        self.mol = None
        self.mw_range = None  # (min_mw, max_mw) tuple
        self.substituted_products = []  # List of (name, smiles) tuples
        self.validate()
    
    def validate(self):
        """Validate the SMILES and check for conflicts."""
        # Try to parse SMILES
        self.mol = Chem.MolFromSmiles(self.smiles_with_attachment)
        if self.mol is None:
            raise ValueError(f"Invalid SMILES for Job {self.job_id}: {self.smiles_with_attachment}")
        
        # Find attachment points
        # IMPORTANT: Only check [*] (isotope 0) for fragment attachment
        # Ignore [16*] or other isotope dummies - those are placeholders, not substitution points!
        dummy_atoms = []
        attachment_carbons = []
        
        for atom in self.mol.GetAtoms():
            if atom.GetAtomicNum() == 0:  # Dummy atom
                # ONLY process [*] with isotope 0 (fragment attachment)
                # Skip [16*], [2H*], etc. (isotope != 0)
                if atom.GetIsotope() == 0:
                    dummy_idx = atom.GetIdx()
                    dummy_atoms.append(dummy_idx)
                    
                    # Find the carbon this dummy is attached to
                    neighbors = list(atom.GetNeighbors())
                    if neighbors:
                        carbon_idx = neighbors[0].GetIdx()
                        attachment_carbons.append(carbon_idx)
        
        if not dummy_atoms:
            raise ValueError(f"No attachment points [*] found in Job {self.job_id}")
        
        self.attachment_atoms = dummy_atoms
        self.attachment_carbons = attachment_carbons
        
        # Check if the CARBONS being substituted overlap with active sites
        # (Not the dummy atoms themselves, which can have any index)
        conflicts = set(attachment_carbons) & set(self.active_site_atoms)
        if conflicts:
            raise ValueError(
                f"Job {self.job_id}: Substitution carbons {conflicts} are active site atoms. "
                f"Active sites cannot be substitution points!"
            )
        
        print(f"✓ Job {self.job_id} validated successfully")
        print(f"  Fragment attachment dummy at: {self.attachment_atoms}")
        print(f"  Substituting carbons at: {attachment_carbons}")
        print(f"  Active sites: {self.active_site_atoms}")
    
    def set_mw_range(self, min_mw, max_mw):
        """Set molecular weight range for fragment selection."""
        self.mw_range = (min_mw, max_mw)
        print(f"✓ MW range set: {min_mw}-{max_mw} Da")
    
    def attach_fragment(self, fragment_smiles):
        """
        Attach a fragment to the scaffold at [*] position.
        Simple string replacement method: removes [*] and inserts fragment.
        
        Example:
            Scaffold: O[C@@H](c1ccc([*])cc1)C(=O)OC
            Fragment: [*]NCC
            Result:   O[C@@H](c1ccc(NCC)cc1)C(=O)OC
        """
        try:
            # Remove [*] from fragment
            fragment_core = fragment_smiles.replace('[*]', '')
            
            # Find position of [*] in scaffold and replace with fragment
            star_pos = self.smiles_with_attachment.find('[*]')
            if star_pos == -1:
                return None
            
            # Build combined SMILES: before + fragment + after
            before = self.smiles_with_attachment[:star_pos]
            after = self.smiles_with_attachment[star_pos + 3:]  # +3 for len('[*]')
            combined_smiles = before + fragment_core + after
            
            # Validate with RDKit
            combined_mol = Chem.MolFromSmiles(combined_smiles)
            if combined_mol is None:
                return None
            
            # Return canonical SMILES
            return Chem.MolToSmiles(combined_mol, canonical=True)
            
        except Exception:
            return None
    
    def __repr__(self):
        mw_info = f", MW={self.mw_range[0]}-{self.mw_range[1]}Da" if self.mw_range else ""
        return f"Job{self.job_id}(SMILES={self.smiles_with_attachment}, attachments={self.attachment_atoms}{mw_info})"


class EnzymeLigandProcessor:
    """Main processing class with active site marking and substitution points."""
    
    def __init__(self, enzyme_file=None, molecule_file=None, output_dir=None, 
                 use_original_as_scaffold=None, fragment_file=None,
                 substitution_jobs=None, interactive=True):
        self.interactive = interactive
        
        # Set default fragment library if not provided
        if fragment_file is None:
            _base = "/Users/deboleenaadhikary/Downloads/Fragment_lib/350_frag_sorted.csv"
            _fallback = "/Users/deboleenaadhikary/Downloads/Fragment_lib/350_frag_cleaned.txt"
            fragment_file = _base if os.path.exists(_base) else _fallback

        # Collect all inputs upfront if interactive
        if interactive and (enzyme_file is None or molecule_file is None):
            inputs = self._collect_all_inputs_upfront(
                enzyme_file, molecule_file, output_dir,
                use_original_as_scaffold, fragment_file, substitution_jobs
            )
            self.enzyme_file = os.path.abspath(inputs['enzyme_file'])
            self.molecule_file = os.path.abspath(inputs['molecule_file'])
            self.output_dir = inputs['output_dir']
            self.use_original_as_scaffold = inputs['use_original_as_scaffold']
            self.fragment_file = inputs['fragment_file']
            self.substitution_jobs_config = inputs['substitution_jobs']
        else:
            self.enzyme_file = os.path.abspath(enzyme_file) if enzyme_file else None
            self.molecule_file = os.path.abspath(molecule_file) if molecule_file else None
            self.output_dir = output_dir or '.'
            self.use_original_as_scaffold = use_original_as_scaffold
            self.fragment_file = fragment_file
            self.substitution_jobs_config = substitution_jobs or []
        
        # Molecule data
        self.original_smiles = None
        self.scaffold_smiles = None
        self.scaffold_mol = None
        self.is_scaffold_original = None
        self.brics_fragments = []
        
        # Active site and substitution data
        self.active_site_atoms = []
        self.substitution_jobs = []  # List of SubstitutionJob objects
        
        # Fragment library
        self.fragment_library = []
        
        # Enzyme data
        self.enzyme_apo_file = None
        self.chains_info = None
        self.removed_hetero = []
        self.removed_chains = []
        
        if self.enzyme_file and self.molecule_file:
            self._validate_inputs()
        os.makedirs(self.output_dir, exist_ok=True)
    
    def _collect_all_inputs_upfront(self, enzyme_file, molecule_file, output_dir,
                                     use_original_as_scaffold, fragment_file, substitution_jobs):
        """
        Collect ALL inputs at the beginning in one interactive session.
        Fragment library and output directory are hardcoded - user not prompted.
        Returns dict with all configuration.
        """
        print("\n" + "="*70)
        print("ENZYME PROCESSOR PIPELINE - INPUT COLLECTION")
        print("="*70)
        print("\nPlease provide all inputs. Press Enter for defaults where shown.\n")
        
        # 1. Enzyme PDB file
        if enzyme_file is None:
            while True:
                enzyme_file = input("1. Enzyme PDB file path: ").strip()
                if os.path.exists(enzyme_file) and enzyme_file.endswith('.pdb'):
                    break
                print("   ❌ File not found or not a .pdb file. Try again.")
        
        # 2. Ligand structure file
        if molecule_file is None:
            while True:
                molecule_file = input("2. Ligand structure file (.mol/.sdf/.pdb): ").strip()
                ext = Path(molecule_file).suffix.lower()
                if os.path.exists(molecule_file) and ext in ['.mol', '.sdf', '.pdb']:
                    break
                print("   ❌ File not found or unsupported format. Try again.")
        
        # 3. Output directory - HARDCODED to script directory, no prompt
        if output_dir is None:
            # Use script's directory / enzyme_processor_results
            script_dir = os.path.dirname(os.path.abspath(__file__))
            output_dir = os.path.join(script_dir, 'enzyme_processor_results')
            os.makedirs(output_dir, exist_ok=True)
            print(f"\n3. Output directory: {output_dir}")
            print("   (automatically set to script directory)")
        
        # 4. Scaffold option
        if use_original_as_scaffold is None:
            print("\n4. Scaffold selection:")
            print("   1) Use molecule as-is (no fragmentation)")
            print("   2) Perform BRICS fragmentation to create scaffold")
            while True:
                choice = input("   Enter choice (1 or 2): ").strip()
                if choice in ['1', '2']:
                    use_original_as_scaffold = (choice == '1')
                    break
                print("   ❌ Invalid choice. Enter 1 or 2.")
        
        # 5. Fragment library - HARDCODED, no prompt
        if fragment_file is None:
            _base = "/Users/deboleenaadhikary/Downloads/Fragment_lib/350_frag_sorted.csv"
            _fallback = "/Users/deboleenaadhikary/Downloads/Fragment_lib/350_frag_cleaned.txt"
            fragment_file = _base if os.path.exists(_base) else _fallback
            fragment_file = _base if os.path.exists(_base) else _fallback
            fragment_file = "/Users/deboleenaadhikary/Downloads/Fragment_lib/350_frag_cleaned.txt"
        
        # Set default output directory if not provided
        if output_dir is None:
            script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
            output_dir = os.path.join(script_dir, 'enzyme_processor_results')
            print(f"\n5. Fragment library: {fragment_file}")
            if not os.path.exists(fragment_file):
                print(f"   ⚠️  WARNING: Fragment library not found at {fragment_file}")
                print("   Pipeline will continue but fragment substitution may fail.")
        
        # 6. Substitution jobs (attachment points + MW ranges)
        if substitution_jobs is None:
            print("\n6. Fragment substitution jobs:")
            print("   Define attachment points (SMILES with [*]) and MW ranges.")
            print("   Example: c1ccccc1[*] means attach at benzene ring\n")
            
            jobs = []
            while True:
                print(f"\n   Job {len(jobs) + 1}:")
                
                # Attachment SMILES
                while True:
                    att_smiles = input("   Attachment point SMILES (with [*], or 'done' to finish): ").strip()
                    if att_smiles.lower() == 'done':
                        break
                    if '[*]' in att_smiles or '[16*]' in att_smiles:
                        break
                    print("   ❌ Invalid. Must contain [*] or [16*]. Try again.")
                
                if att_smiles.lower() == 'done':
                    break
                
                # MW range
                print("   Molecular weight range:")
                while True:
                    mw_input = input("   Enter as 'min-max' (e.g., '0-200' or '0-350'): ").strip()
                    try:
                        if '-' in mw_input:
                            parts = mw_input.split('-')
                            mw_min = float(parts[0]) if parts[0] else None
                            mw_max = float(parts[1]) if parts[1] else None
                            mw_range = (mw_min, mw_max)
                            break
                        else:
                            print("   ❌ Invalid format. Use 'min-max' or 'All'.")
                    except ValueError:
                        print("   ❌ Invalid numbers. Try again.")
                
                jobs.append({
                    'attachment_smiles': att_smiles,
                    'mw_range': mw_range
                })
                print(f"   ✓ Job {len(jobs)} added")
            
            if not jobs:
                # Default job if none provided
                print("\n   No jobs defined. Using default: [*] attachment, 0-350 Da")
                jobs = [{'attachment_smiles': '[*]', 'mw_range': (0, 350)}]
            
            substitution_jobs = jobs
        
        print("\n" + "="*70)
        print("✓ All inputs collected. Starting pipeline...")
        print("="*70)
        
        return {
            'enzyme_file': enzyme_file,
            'molecule_file': molecule_file,
            'output_dir': output_dir,
            'use_original_as_scaffold': use_original_as_scaffold,
            'fragment_file': fragment_file,
            'substitution_jobs': substitution_jobs
        }
    
    def _validate_inputs(self):
        """Validate input files."""
        if not os.path.exists(self.enzyme_file):
            raise FileNotFoundError(f"Enzyme file not found: {self.enzyme_file}")
        
        if not os.path.exists(self.molecule_file):
            raise FileNotFoundError(f"Molecule file not found: {self.molecule_file}")
        
        if not self.enzyme_file.endswith('.pdb'):
            raise ValueError("Enzyme file must be a .pdb file")
        
        mol_ext = Path(self.molecule_file).suffix.lower()
        if mol_ext not in ['.mol', '.pdb', '.sdf']:
            raise ValueError(f"Molecule file must be .mol, .pdb, or .sdf (got {mol_ext})")
    
    def convert_molecule_to_smiles(self):
        """Convert molecule to SMILES using OpenBabel."""
        print("\n" + "="*60)
        print("STEP 1: Converting Molecule to SMILES")
        print("="*60)
        
        file_ext = Path(self.molecule_file).suffix.lower().replace('.', '')
        
        obConversion = ob.OBConversion()
        obConversion.SetInFormat(file_ext)
        obConversion.SetOutFormat("smi")
        
        mol = ob.OBMol()
        success = obConversion.ReadFile(mol, self.molecule_file)
        
        if not success:
            raise ValueError(f"Failed to read molecule file: {self.molecule_file}")
        
        raw_smiles = obConversion.WriteString(mol).strip()
        self.original_smiles = raw_smiles.split()[0] if raw_smiles else ''
        
        formula = mol.GetFormula()
        mol_weight = mol.GetMolWt()
        num_atoms = mol.NumAtoms()
        
        print(f"✓ Molecule file: {self.molecule_file}")
        print(f"✓ Format: {file_ext.upper()}")
        print(f"✓ SMILES (stored in memory): {self.original_smiles}")
        print(f"  Formula: {formula}")
        print(f"  Molecular Weight: {mol_weight:.2f}")
        print(f"  Number of Atoms: {num_atoms}")
        
        return self.original_smiles
    
    def process_scaffold_selection(self, use_original=None):
        """Ask user if molecule is scaffold or perform BRICS fragmentation."""
        print("\n" + "="*60)
        print("STEP 2: Scaffold Selection")
        print("="*60)
        
        if use_original is None and self.interactive:
            print("\nIs the molecule entered the scaffold itself, or would you like to make a scaffold?")
            print("  1) This should be my scaffold (use molecule as-is)")
            print("  2) Make a scaffold out of this (BRICS fragmentation)")
            
            while True:
                choice = input("\nEnter choice (1 or 2): ").strip()
                if choice == '1':
                    use_original = True
                    break
                elif choice == '2':
                    use_original = False
                    break
                else:
                    print("Invalid choice. Please enter 1 or 2.")
        elif use_original is None:
            use_original = True
        
        if use_original:
            print("\n✓ Using original molecule as scaffold")
            self.scaffold_smiles = self.original_smiles
            self.is_scaffold_original = True
            print(f"  Scaffold SMILES: {self.scaffold_smiles}")
        else:
            print("\n✓ Performing BRICS fragmentation...")
            self.scaffold_smiles = self._perform_brics_fragmentation()
            self.is_scaffold_original = False
        
        # Create RDKit mol object for scaffold
        self.scaffold_mol = Chem.MolFromSmiles(self.scaffold_smiles)
        if self.scaffold_mol is None:
            raise ValueError(f"Failed to parse scaffold SMILES: {self.scaffold_smiles}")
        
        return self.scaffold_smiles
    
    def _perform_brics_fragmentation(self):
        """Perform BRICS fragmentation and select highest MW fragment."""
        mol = Chem.MolFromSmiles(self.original_smiles)
        
        if mol is None:
            raise ValueError(f"Failed to parse SMILES: {self.original_smiles}")
        
        print("  Fragmenting molecule using BRICS...")
        fragments = BRICS.BRICSDecompose(mol)
        fragments_list = list(fragments)
        
        if not fragments_list:
            print("  ⚠ Warning: No fragments generated. Using original molecule as scaffold.")
            return self.original_smiles
        
        print(f"  Found {len(fragments_list)} fragments")
        
        fragment_data = []
        for i, frag_smiles in enumerate(fragments_list):
            clean_smiles = frag_smiles.replace('[*]', '')
            frag_mol = Chem.MolFromSmiles(clean_smiles)
            if frag_mol:
                mw = Descriptors.MolWt(frag_mol)
                fragment_data.append({
                    'index': i + 1,
                    'smiles': clean_smiles,
                    'original_smiles': frag_smiles,
                    'mol': frag_mol,
                    'mw': mw,
                    'num_atoms': frag_mol.GetNumAtoms()
                })
        
        fragment_data.sort(key=lambda x: x['mw'], reverse=True)
        self.brics_fragments = fragment_data
        
        print("\n  Fragments sorted by molecular weight (highest to lowest):")
        print("  " + "="*70)
        print(f"  {'#':<4} {'MW':<10} {'Atoms':<8} {'SMILES':<50}")
        print("  " + "-"*70)
        for frag in fragment_data[:10]:
            print(f"  {frag['index']:<4} {frag['mw']:<10.2f} {frag['num_atoms']:<8} {frag['smiles']:<50}")
        if len(fragment_data) > 10:
            print(f"  ... and {len(fragment_data) - 10} more fragments")
        
        scaffold = fragment_data[0]
        print("\n✓ Selected highest molecular weight fragment as scaffold:")
        print(f"  Scaffold SMILES: {scaffold['smiles']}")
        print(f"  Molecular Weight: {scaffold['mw']:.2f}")
        print(f"  Number of Atoms: {scaffold['num_atoms']}")
        
        return scaffold['smiles']
    
    def display_scaffold_structure(self):
        """Display scaffold structure with atom indices."""
        print("\n" + "="*60)
        print("STEP 3: Scaffold Structure Visualization")
        print("="*60)
        
        print(f"\nScaffold SMILES: {self.scaffold_smiles}")
        print(f"Number of atoms: {self.scaffold_mol.GetNumAtoms()}")
        
        # Display atom information
        print("\nAtom indices and types:")
        print("  " + "="*50)
        print(f"  {'Index':<8} {'Symbol':<10} {'Type':<15} {'Neighbors':<15}")
        print("  " + "-"*50)
        
        for atom in self.scaffold_mol.GetAtoms():
            idx = atom.GetIdx()
            symbol = atom.GetSymbol()
            hybrid = str(atom.GetHybridization())
            neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
            print(f"  {idx:<8} {symbol:<10} {hybrid:<15} {str(neighbors):<15}")
        
        # Try to save structure image
        try:
            img_path = os.path.join(self.output_dir, 'scaffold_structure.png')
            
            # Add atom indices to the molecule for visualization
            for atom in self.scaffold_mol.GetAtoms():
                atom.SetProp('atomNote', str(atom.GetIdx()))
            
            img = Draw.MolToImage(self.scaffold_mol, size=(600, 600))
            img.save(img_path)
            print(f"\n✓ Scaffold structure image saved to: {img_path}")
            print("  (Open this image to see atom positions)")
        except Exception as e:
            print(f"\n⚠ Could not save structure image: {e}")
    
    def mark_active_sites(self):
        """Let user mark active site reaction atoms."""
        print("\n" + "="*60, flush=True)
        print("STEP 4: Mark Active Site Reaction Atoms", flush=True)
        print("="*60, flush=True)
        
        if self.interactive:
            print("\nEnter active site atom indices separated by commas (e.g., 2,5,8)")
            print("Or press Enter if there are no specific active site atoms:")
            
            while True:
                user_input = input("\nActive site atoms: ").strip()
                
                if not user_input:
                    self.active_site_atoms = []
                    print("✓ No active site atoms marked")
                    break
                
                try:
                    indices = [int(x.strip()) for x in user_input.split(',')]
                    max_idx = self.scaffold_mol.GetNumAtoms() - 1
                    invalid = [i for i in indices if i < 0 or i > max_idx]
                    
                    if invalid:
                        print(f"❌ Invalid indices: {invalid}. Must be between 0 and {max_idx}")
                        continue
                    
                    self.active_site_atoms = indices
                    print(f"✓ Active site atoms marked: {self.active_site_atoms}")
                    break
                    
                except ValueError as e:
                    print(f"❌ Invalid input: {e}. Please enter comma-separated numbers.")
        else:
            # FIX: Keep the active sites passed from the web interface!
            print(f"✓ Active site atoms set from web interface: {self.active_site_atoms}", flush=True)
        
        return self.active_site_atoms
    
    def select_fragment_mw_range(self, job):
        """
        Let user select molecular weight range for fragments to be used.
        
        Parameters:
        -----------
        job : SubstitutionJob
            The job for which to select MW range
        """
        print(f"\n" + "-"*60)
        print(f"Fragment Molecular Weight Selection for Job {job.job_id}")
        print("-"*60)
        
        print("\nSelect the molecular weight range for fragments to attach:")
        print("  1) 0-200 Da")
        print("  2) 0-350 Da (optimum)")
        print("  3) Custom range (enter min-max)")

        mw_ranges = {
            '1': (0, 200),
            '2': (0, 350),
        }

        while True:
            choice = input("\nEnter choice (1-3, default=2): ").strip() or '2'

            if choice in mw_ranges:
                min_mw, max_mw = mw_ranges[choice]
                print(f"Selected: {min_mw}-{max_mw} Da")
                job.set_mw_range(min_mw, max_mw)
                break

            elif choice == '3':
                # Custom range
                print("\nEnter custom molecular weight range:")
                try:
                    min_input = input("  Minimum MW (Da): ").strip()
                    max_input = input("  Maximum MW (Da): ").strip()
                    
                    min_mw = float(min_input)
                    max_mw = float(max_input)
                    
                    if min_mw < 0:
                        print("❌ Minimum MW cannot be negative")
                        continue
                    
                    if max_mw <= min_mw:
                        print("❌ Maximum MW must be greater than minimum MW")
                        continue
                    
                    print(f"Selected: Custom range {min_mw}-{max_mw} Da")
                    job.set_mw_range(min_mw, max_mw)
                    break
                    
                except ValueError:
                    print("❌ Invalid input. Please enter numeric values.")
            
            else:
                print("❌ Invalid choice. Please enter a number from 1 to 6.")
    
    def create_substitution_jobs(self):
        """Let user create up to 3 substitution jobs with attachment points."""
        print("\n" + "="*60, flush=True)
        print("STEP 5: Create Substitution Jobs", flush=True)
        print("="*60, flush=True)
        
        if not self.interactive:
            # FIX: Actually process the jobs sent from the frontend!
            print("\nNon-interactive mode: Loading jobs from web interface...", flush=True)
            jobs = []
            for i, job_cfg in enumerate(self.substitution_jobs_config, 1):
                try:
                    job = SubstitutionJob(
                        job_id=i,
                        smiles_with_attachment=job_cfg['attachment_smiles'],
                        attachment_atoms=[], 
                        active_site_atoms=self.active_site_atoms
                    )
                    job.set_mw_range(job_cfg['mw_range'][0], job_cfg['mw_range'][1])
                    jobs.append(job)
                    print(f"✓ Job {i} initialized: {job_cfg['attachment_smiles']}", flush=True)
                except Exception as e:
                    print(f"❌ Error setting up Job {i}: {e}", flush=True)
            
            self.substitution_jobs = jobs
            return jobs
            
        print("\nYou can create up to 3 different substitution patterns (jobs).")
        print("Each job specifies where fragments can be attached using [*] markers.")
        
        print(f"\nYour scaffold SMILES (without attachment points):")
        print(f"  {self.scaffold_smiles}")
        
        print("\nHow to add attachment points:")
        print("  - Place [*] in parentheses next to the atom where you want substitution")
        print("  - Use explicit notation [CH], [C], [N] when needed for correct parsing")
        
        if 'c1ccccc1' in self.scaffold_smiles or 'c1cccc1' in self.scaffold_smiles:
            example_smiles = self.scaffold_smiles.replace('c1ccccc1', 'c1ccc([*])cc1', 1)
            print(f"  Original:  {self.scaffold_smiles}")
            print(f"  Modified:  {example_smiles}")
            print(f"                  ^^^ [*] added at para position")
        else:
            print(f"  Original:  {self.scaffold_smiles}")
            print(f"  Modified:  [Add [*] next to desired atom, e.g., C([*]) or c([*])]")
        
        print("\nIMPORTANT:")
        print("  - Attachment points CANNOT overlap with active sites")
        print(f"  - Active site atoms: {self.active_site_atoms if self.active_site_atoms else 'None'}")
        
        jobs = []
        for job_num in range(1, 4):
            print(f"\n{'='*60}")
            print(f"Job {job_num}")
            print(f"{'='*60}")
            
            create = input(f"\nCreate Job {job_num}? (y/n, default=n): ").strip().lower()
            
            if create != 'y':
                print(f"Skipping Job {job_num}")
                break
            
            while True:
                print(f"\nEnter SMILES with [*] attachment point(s) for Job {job_num}:")
                print(f"Your scaffold:  {self.scaffold_smiles}")
                print(f"Hint: Add ([*]) next to the atom where you want substitution")
                smiles_input = input("Modified SMILES: ").strip()
                
                if not smiles_input:
                    print("Empty SMILES, skipping this job")
                    break
                
                try:
                    job = SubstitutionJob(
                        job_id=job_num,
                        smiles_with_attachment=smiles_input,
                        attachment_atoms=[], 
                        active_site_atoms=self.active_site_atoms
                    )
                    
                    print(f"\n✅ Job {job_num} created successfully!")
                    
                    try:
                        img_path = os.path.join(self.output_dir, f'job{job_num}_structure.png')
                        for atom in job.mol.GetAtoms():
                            atom.SetProp('atomNote', str(atom.GetIdx()))
                        img = Draw.MolToImage(job.mol, size=(600, 600))
                        img.save(img_path)
                        print(f"  Structure image saved: {img_path}")
                    except:
                        pass
                    
                    self.select_fragment_mw_range(job)
                    jobs.append(job)
                    break
                    
                except Exception as e:
                    print(f"\n❌ Error: {e}")
                    print("Please try again with corrected SMILES")
                    retry = input(f"Retry Job {job_num}? (y/n): ").strip().lower()
                    if retry != 'y':
                        break
        
        self.substitution_jobs = jobs
        
        print(f"\n{'='*60}")
        print(f"Summary: {len(jobs)} substitution job(s) created")
        print(f"{'='*60}")
        
        for job in jobs:
            mw_str = f"{job.mw_range[0]}-{job.mw_range[1]}Da" if job.mw_range else "Not set"
            print(f"\nJob {job.job_id}:")
            print(f"  SMILES: {job.smiles_with_attachment}")
            print(f"  Attachment atoms: {job.attachment_atoms}")
            print(f"  Fragment MW range: {mw_str}")
        
        return jobs
    
    # ------------------------------------------------------------------
    # Volume-based fragment filtering helpers
    # ------------------------------------------------------------------

    def calculate_mol_volume(self, smiles):
        """
        Calculate the van der Waals volume (Å³) of a molecule from SMILES.
        Uses RDKit's AllChem.ComputeMolVolume with an embedded 3D conformer.

        Returns float volume in Å³, or None if calculation fails.
        """
        try:
            smiles = self._strip_dummy_atoms(smiles)
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            mol = Chem.AddHs(mol)
            result = AllChem.EmbedMolecule(mol, randomSeed=42, maxAttempts=200)
            if result != 0:
                # Fallback: use ETKDGv3
                params = AllChem.ETKDGv3()
                params.randomSeed = 42
                result = AllChem.EmbedMolecule(mol, params)
            if result != 0:
                return None
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=100)
            except Exception:
                pass
            vol = AllChem.ComputeMolVolume(mol, gridSpacing=0.2)
            return round(vol, 2)
        except Exception as e:
            print(f"  ⚠️  Volume calculation failed: {e}")
            return None

    def estimate_pocket_volume(self, pocket):
        """
        Estimate the binding pocket volume in Å³ from pocket data.

        Priority:
          1. fpocket 'Real volume' stored in pocket['volume'] (exact)
          2. Pocket bounding box:  size_x * size_y * size_z * 0.4
             (40% packing factor — a sphere inscribed in the box occupies ~52%
              of the box; biological pockets are typically 30-50% filled)

        Returns float volume in Å³.
        """
        # fpocket provides a direct volume estimate
        vol = pocket.get('volume', 0)
        if vol and vol > 10:
            return float(vol)

        # Fallback: bounding box heuristic
        sx = pocket.get('size_x', 20.0)
        sy = pocket.get('size_y', 20.0)
        sz = pocket.get('size_z', 20.0)
        return sx * sy * sz * 0.4

    # ── Sorted-CSV fragment library loader ───────────────────────────────────
    @staticmethod
    def _load_sorted_csv_fragments(csv_path, min_mw, max_mw):
        """
        Load only the fragments whose Avg_MW falls within [min_mw, max_mw] from
        a pre-sorted CSV (sorted ascending by Avg_MW column).

        Uses binary search (bisect) to locate the exact row window so rows
        outside the MW range are never read into memory.  This turns an O(n)
        full-library scan into an O(log n + k) operation where k is the number
        of qualifying fragments.

        Expected CSV columns: Fragment_SMILES, Avg_MW, Exact_MW

        Returns list of dicts: {index, smiles, mw, exact_mw}
        """
        import bisect, csv as _csv

        # ── Pass 1: read all Avg_MW values into a lightweight list for bisect ──
        avg_mws = []
        all_rows = []
        with open(csv_path, newline='', encoding='utf-8') as f:
            reader = _csv.DictReader(f)
            for row in reader:
                try:
                    avg_mws.append(float(row['Avg_MW']))
                    all_rows.append(row)
                except (ValueError, KeyError):
                    continue   # skip malformed rows

        # ── Binary search: find the slice [lo, hi) that satisfies the MW window ─
        lo = bisect.bisect_left(avg_mws, min_mw)
        hi = bisect.bisect_right(avg_mws, max_mw)

        result = []
        for global_idx, i in enumerate(range(lo, hi)):
            row = all_rows[i]
            smiles = row.get('Fragment_SMILES', '').strip()
            if not smiles:
                continue
            try:
                exact_mw = float(row.get('Exact_MW', avg_mws[i]))
            except (ValueError, TypeError):
                exact_mw = avg_mws[i]
            result.append({
                'index':    lo + global_idx,
                'smiles':   smiles,
                'mw':       avg_mws[i],
                'exact_mw': exact_mw,
            })
        return result, len(avg_mws)

    @staticmethod
    def _load_txt_fragments(txt_path):
        """
        Legacy loader: read all SMILES from a plain .txt file (one per line).
        Calculates MW on-the-fly via RDKit.  Kept for backward compatibility.
        Returns list of dicts: {index, smiles, mw}
        """
        fragments = []
        invalid_count = 0
        with open(txt_path, 'r') as f:
            raw_lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]

        for i, frag_smiles in enumerate(raw_lines):
            frag_mol = Chem.MolFromSmiles(frag_smiles)
            if frag_mol is None:
                invalid_count += 1
                continue
            try:
                canonical_smiles = Chem.MolToSmiles(frag_mol, canonical=True)
                verify_mol = Chem.MolFromSmiles(canonical_smiles)
                if verify_mol is None:
                    invalid_count += 1
                    continue
                frag_mol = verify_mol
            except Exception:
                invalid_count += 1
                continue
            try:
                rw = Chem.RWMol(frag_mol)
                dummy_idxs = sorted(
                    [a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 0],
                    reverse=True
                )
                for idx in dummy_idxs:
                    rw.RemoveAtom(idx)
                clean_mol = rw.GetMol()
                Chem.SanitizeMol(clean_mol)
                mw = Descriptors.MolWt(clean_mol)
            except Exception:
                mw = Descriptors.MolWt(frag_mol)
            fragments.append({'index': i, 'smiles': canonical_smiles, 'mw': mw})

        return fragments, invalid_count

    def substitute_fragments(self, fragment_file, max_fragment_volume=None):
        """
        Substitute fragments from the library into each substitution job.

        Fast path (sorted CSV):
            If `fragment_file` ends with '.csv', the file is treated as a
            pre-sorted fragment library (sorted ascending by Avg_MW).  Binary
            search locates the exact MW window for each job so rows outside
            the range are never loaded — O(log n + k) instead of O(n).
            MW values come directly from the CSV; no per-fragment RDKit
            calculation is needed.

        Legacy path (.txt):
            Plain SMILES file scanned fully with on-the-fly MW calculation
            (original behaviour, kept for backward compatibility).

        Parameters:
        -----------
        fragment_file : str
            Path to sorted CSV ('.csv') or legacy plain SMILES ('.txt').
        max_fragment_volume : float or None
            Maximum allowed VDW volume (Å³) per fragment (pocket-size guard).
        """
        print("\n" + "="*60)
        print("STEP 7: Fragment Substitution")
        print("="*60)

        if not self.substitution_jobs:
            print("\n⚠ No substitution jobs found. Skipping fragment substitution.")
            return {}

        if not os.path.exists(fragment_file):
            raise FileNotFoundError(f"Fragment file not found: {fragment_file}")

        use_csv = fragment_file.lower().endswith('.csv')
        print(f"\nFragment library : {fragment_file}")
        print(f"Mode             : {'sorted CSV (binary-search)' if use_csv else 'legacy .txt (full scan)'}")

        # ── If using CSV, we do per-job loading (each job may have a different
        #    MW window, so we binary-search independently for each job).
        # ── If using .txt, load everything once and filter in memory.
        txt_fragment_data = None   # populated only in legacy mode
        txt_total         = 0

        if not use_csv:
            print("\nLoading full fragment library…")
            txt_fragment_data, invalid_count = self._load_txt_fragments(fragment_file)
            txt_total = len(txt_fragment_data)
            print(f"✓ Loaded {txt_total} valid fragments")
            if invalid_count:
                print(f"  ⚠ Skipped {invalid_count} invalid SMILES")

            # Volume pre-filter (legacy mode only — CSV mode does it per-job below)
            if max_fragment_volume is not None and max_fragment_volume > 0:
                print(f"\n🔬 Volume filter: max VDW = {max_fragment_volume:.1f} Å³")
                vf, vs = [], 0
                for frag in txt_fragment_data:
                    fvol = self.calculate_mol_volume(frag['smiles'].replace('[*]', ''))
                    if fvol is None or fvol <= max_fragment_volume:
                        vf.append(frag)
                    else:
                        vs += 1
                print(f"   ✓ Kept {len(vf)} / {txt_total} ({vs} removed by volume filter)")
                txt_fragment_data = vf

        # ── Per-job processing ────────────────────────────────────────────────
        job_stats = {}

        for job in self.substitution_jobs:
            print(f"\n{'-'*60}")
            print(f"Processing Job {job.job_id}")
            print(f"{'-'*60}")
            print(f"  Scaffold : {job.smiles_with_attachment}")
            min_mw, max_mw = job.mw_range
            print(f"  MW range : {min_mw:.0f} – {max_mw:.0f} Da")

            # ── Retrieve fragments for this MW window ─────────────────
            if use_csv:
                raw_frags, total_in_lib = self._load_sorted_csv_fragments(
                    fragment_file, min_mw, max_mw
                )
                print(f"  Binary search on {total_in_lib:,} entries → "
                      f"{len(raw_frags):,} fragments in MW window "
                      f"(skipped {total_in_lib - len(raw_frags):,} outside range)")

                # Canonicalise SMILES and validate (only the window subset)
                fragment_data = []
                invalid_count = 0
                for frag in raw_frags:
                    frag_mol = Chem.MolFromSmiles(frag['smiles'])
                    if frag_mol is None:
                        invalid_count += 1
                        continue
                    try:
                        canonical = Chem.MolToSmiles(frag_mol, canonical=True)
                        if Chem.MolFromSmiles(canonical) is None:
                            invalid_count += 1
                            continue
                        frag['smiles'] = canonical
                    except Exception:
                        invalid_count += 1
                        continue
                    fragment_data.append(frag)

                if invalid_count:
                    print(f"  ⚠ Skipped {invalid_count} unparseable SMILES in window")

                # Optional volume filter on the (small) window subset
                if max_fragment_volume is not None and max_fragment_volume > 0:
                    vf, vs = [], 0
                    for frag in fragment_data:
                        fvol = self.calculate_mol_volume(frag['smiles'].replace('[*]', ''))
                        if fvol is None or fvol <= max_fragment_volume:
                            vf.append(frag)
                        else:
                            vs += 1
                    if vs:
                        print(f"  🔬 Volume filter removed {vs} more fragments")
                    fragment_data = vf

            else:
                # Legacy .txt mode — MW filter already applied above
                fragment_data = [
                    frag for frag in txt_fragment_data
                    if min_mw <= frag['mw'] <= max_mw
                ]
                total_in_lib = txt_total
                print(f"  MW filter → {len(fragment_data):,} / {total_in_lib:,} pass")

            if not fragment_data:
                print(f"  ⚠ No fragments in MW range, skipping job")
                job_stats[job.job_id] = {
                    'total_fragments': total_in_lib,
                    'filtered_fragments': 0,
                    'generated_ligands': 0,
                    'failed': 0,
                }
                continue

            # ── Fragment substitution ─────────────────────────────────
            ligands_generated = []
            failed_count = 0

            for frag in fragment_data:
                try:
                    combined_smiles = self._combine_scaffold_fragment(
                        job.smiles_with_attachment,
                        frag['smiles']
                    )
                    if combined_smiles:
                        ligands_generated.append({
                            'smiles':          combined_smiles,
                            'fragment_index':  frag['index'],
                            'fragment_smiles': frag['smiles'],
                            'fragment_mw':     frag['mw'],
                        })
                    else:
                        failed_count += 1
                except Exception as e:
                    failed_count += 1
                    print(f"  ❌ Exception with fragment {frag['index']}: {e}")

            print(f"  ✓ Generated {len(ligands_generated):,} ligands"
                  + (f"  (⚠ {failed_count} failed)" if failed_count else ""))

            if ligands_generated:
                self._save_ligands_to_files(job.job_id, ligands_generated)

            job_stats[job.job_id] = {
                'total_fragments':    total_in_lib,
                'filtered_fragments': len(fragment_data),
                'generated_ligands':  len(ligands_generated),
                'failed':             failed_count
            }
        
        # Print summary
        print(f"\n{'='*60}")
        print("Fragment Substitution Summary")
        print(f"{'='*60}")
        for job_id, stats in job_stats.items():
            print(f"\nJob {job_id}:")
            print(f"  Total fragments: {stats['total_fragments']}")
            print(f"  Filtered by MW: {stats['filtered_fragments']}")
            print(f"  Ligands generated: {stats['generated_ligands']}")
            if stats['failed'] > 0:
                print(f"  Failed: {stats['failed']}")
        
        return job_stats
    
    def _combine_scaffold_fragment(self, scaffold_smiles, fragment_smiles):
        """
        Combine scaffold and fragment properly handling multiple dummy atoms.
        Ensures ALL dummy atoms are removed and output is clean for OpenBabel.
        
        Parameters:
        -----------
        scaffold_smiles : str
            Scaffold SMILES with [*] attachment point (may have other dummies like [16*])
        fragment_smiles : str
            Fragment SMILES with [*] attachment point
        
        Returns:
        --------
        str : Combined SMILES with ALL dummies removed, clean for OpenBabel
        """
        try:
            from rdkit.Chem import RWMol, BondType
            
            # Parse molecules
            scaffold_mol = Chem.MolFromSmiles(scaffold_smiles)
            fragment_mol = Chem.MolFromSmiles(fragment_smiles)
            
            if scaffold_mol is None or fragment_mol is None:
                print(f"    ❌ Could not parse SMILES")
                return None
            
            def _find_attachment(mol, isotope_zero_only):
                for atom in mol.GetAtoms():
                    if atom.GetAtomicNum() != 0:
                        continue
                    if isotope_zero_only and atom.GetIsotope() != 0:
                        continue

                    neighbors = list(atom.GetNeighbors())
                    if not neighbors:
                        continue

                    attach_idx = neighbors[0].GetIdx()
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), attach_idx)
                    bond_type = bond.GetBondType() if bond else BondType.SINGLE
                    return atom.GetIdx(), attach_idx, bond_type
                return None, None, None

            # Find the substitution dummy [*] on the scaffold and the first dummy
            # on the fragment. Remaining dummies such as [16*] are removed later.
            scaf_subst_dummy, scaf_attach_atom, scaf_bond_type = _find_attachment(
                scaffold_mol, isotope_zero_only=True
            )
            frag_dummy, frag_attach_atom, frag_bond_type = _find_attachment(
                fragment_mol, isotope_zero_only=False
            )

            if None in [scaf_subst_dummy, scaf_attach_atom, frag_dummy, frag_attach_atom]:
                print(f"    ❌ Could not find attachment points")
                return None
            
            # Preserve full atom/bond metadata by combining the intact molecules
            # and only then removing the two attachment dummies.
            combo = Chem.CombineMols(scaffold_mol, fragment_mol)
            result = RWMol(combo)
            frag_offset = scaffold_mol.GetNumAtoms()
            frag_dummy_idx = frag_dummy + frag_offset
            frag_attach_idx = frag_attach_atom + frag_offset

            connect_bond_type = scaf_bond_type or BondType.SINGLE
            if frag_bond_type is not None and frag_bond_type != BondType.SINGLE:
                connect_bond_type = frag_bond_type

            result.AddBond(scaf_attach_atom, frag_attach_idx, connect_bond_type)

            for dummy_idx in sorted([scaf_subst_dummy, frag_dummy_idx], reverse=True):
                result.RemoveAtom(dummy_idx)

            # Now remove ALL remaining dummy atoms (like [16*])
            mol = result.GetMol()
            
            # Keep removing dummies until none left
            max_iterations = 10
            for iteration in range(max_iterations):
                dummies_found = []
                for atom in mol.GetAtoms():
                    if atom.GetAtomicNum() == 0:  # Any dummy atom
                        dummies_found.append(atom.GetIdx())
                
                if not dummies_found:
                    break  # No more dummies
                
                # Remove dummies (highest index first to preserve lower indices)
                dummies_found.sort(reverse=True)
                mol = RWMol(mol)
                for dummy_idx in dummies_found:
                    mol.RemoveAtom(dummy_idx)
                mol = mol.GetMol()

            # Fragment libraries can carry isotope-tagged atoms like [3OH] as
            # attachment-encoding artefacts. They should not survive into the
            # final dockable ligand.
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() != 0 and atom.GetIsotope() != 0:
                    atom.SetIsotope(0)
            
            # Sanitize with multiple attempts
            try:
                Chem.SanitizeMol(mol)
            except Chem.KekulizeException as e:
                print(f"    ⚠️ Kekulization failed: {e}")
                # Try alternative sanitization without kekulization
                try:
                    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE)
                    print(f"    ✓ Using alternative sanitization (no kekulize)")
                except Exception as e2:
                    print(f"    ❌ Sanitization failed completely: {e2}")
                    return None
            except Exception as e:
                print(f"    ❌ Sanitization error: {e}")
                return None
            
            # Try to assign stereochemistry
            try:
                Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            except Exception as e:
                print(f"    ⚠️ Could not assign stereochemistry: {e}")
                # Continue anyway
            
            # Get clean SMILES (no dummies!)
            try:
                output_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            except Exception as e:
                print(f"    ❌ Could not generate SMILES: {e}")
                return None
            
            # CRITICAL CHECK: Make sure no dummies remain
            if '[*]' in output_smiles or '(*)' in output_smiles:
                print(f"    ❌ ERROR: Dummy atoms still present in output!")
                return None
            
            # Check for disconnected fragments
            if '.' in output_smiles:
                print(f"    ❌ ERROR: Disconnected fragments")
                return None
            
            # Verify atom count
            expected = scaffold_mol.GetNumAtoms() + fragment_mol.GetNumAtoms()
            # Subtract all dummies
            for atom in scaffold_mol.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    expected -= 1
            for atom in fragment_mol.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    expected -= 1
            
            actual = mol.GetNumAtoms()
            
            if expected != actual:
                print(f"    ⚠️ WARNING: Atom count mismatch (expected {expected}, got {actual})")
            else:
                print(f"    ✓ Atom count correct: {actual} atoms")
            
            print(f"    ✓ Clean SMILES: {output_smiles[:80]}...")
            
            return output_smiles
            
        except Exception as e:
            print(f"    ❌ Exception: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def _save_ligands_to_files(self, job_id, ligands):
        """
        Save generated ligands to individual .smi files.
        
        Parameters:
        -----------
        job_id : int
            Job ID
        ligands : list
            List of ligand dictionaries
        """
        import shutil
        job_dir = os.path.join(self.output_dir, f'job{job_id}_ligands')

        # Wipe the whole directory so no stale .smi/.sdf/.pdbqt files
        # from a previous run survive into the current filtered set
        if os.path.exists(job_dir):
            shutil.rmtree(job_dir)
            print(f"  Cleared previous job{job_id} directory (stale ligand prevention)")
        os.makedirs(job_dir, exist_ok=True)

        print(f"\n  Saving ligands to: {job_dir}")

        # ── Deduplicate by canonical SMILES ──────────────────────────────────
        # The fragment library can contain chemically identical entries
        # (e.g. the same small group listed at different MW rows, or
        # stereoisomers that collapse to the same canonical form after
        # combining with the scaffold).  Keeping duplicates wastes docking
        # time and causes the same compound to appear multiple times in the
        # results CSV with different ligand_N names.
        seen_smiles = set()
        unique_ligands = []
        dup_count = 0
        for lig in ligands:
            try:
                mol = Chem.MolFromSmiles(lig['smiles'])
                canonical = Chem.MolToSmiles(mol, canonical=True) if mol else lig['smiles']
            except Exception:
                canonical = lig['smiles']
            if canonical in seen_smiles:
                dup_count += 1
                continue
            seen_smiles.add(canonical)
            lig = dict(lig)          # don't mutate caller's dict
            lig['smiles'] = canonical
            unique_ligands.append(lig)

        if dup_count:
            print(f"  ⚠ Removed {dup_count} duplicate SMILES "
                  f"({len(unique_ligands)} unique ligands remain)")
        ligands = unique_ligands
        # ─────────────────────────────────────────────────────────────────────

        for i, ligand in enumerate(ligands, 1):
            filename = f'ligand_{i}.smi'
            filepath = os.path.join(job_dir, filename)
            
            with open(filepath, 'w') as f:
                # Write SMILES with metadata as comment
                f.write(f"{ligand['smiles']}\n")
                f.write(f"# Fragment: {ligand['fragment_smiles']}\n")
                f.write(f"# Fragment MW: {ligand['fragment_mw']:.2f}\n")
                f.write(f"# Fragment Index: {ligand['fragment_index']}\n")
        
        print(f"  ✓ Saved {len(ligands)} .smi files")
    
    def convert_enzyme_to_pdbqt(self):
        """
        Convert enzyme PDB to PDBQT format for docking.
        Uses OpenBabel with -xr (rigid) and -xh (add hydrogens).
        
        Returns:
        --------
        str : Path to PDBQT file
        """
        print("\n" + "="*60)
        print("STEP 8: Converting Enzyme to PDBQT")
        print("="*60)
        
        if not self.enzyme_apo_file:
            print("\n⚠️ No apo enzyme PDB found, skipping conversion")
            return None
        
        # Output path
        pdbqt_path = self.enzyme_apo_file.replace('.pdb', '.pdbqt')
        
        print(f"\nConverting enzyme to PDBQT format...")
        print(f"  Input:  {self.enzyme_apo_file}")
        print(f"  Output: {pdbqt_path}")
        
        try:
            # OpenBabel command: -xr (rigid), -xh (add hydrogens)
            cmd = [
                'obabel',
                self.enzyme_apo_file,
                '-O', pdbqt_path,
                '-xr',  # Keep structure rigid
                '-xh'   # Add hydrogens
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            print(f"✓ Enzyme converted successfully")
            print(f"  Saved to: {pdbqt_path}")
            
            return pdbqt_path
            
        except subprocess.CalledProcessError as e:
            print(f"❌ OpenBabel conversion failed: {e}")
            if e.stderr:
                print(f"   Error: {e.stderr}")
            return None
        except FileNotFoundError:
            print(f"❌ OpenBabel not found. Please install: conda install -c conda-forge openbabel")
            return None
    
    def convert_ligands_to_pdbqt(self, job_id):
        """
        Convert ligands for a job to PDBQT format.
        Multiple fallback strategies for problematic molecules.
        
        Parameters:
        -----------
        job_id : int
            Job ID to process
        
        Returns:
        --------
        dict : Conversion statistics
        """
        job_ligands_dir = os.path.join(self.output_dir, f'job{job_id}_ligands')
        
        if not os.path.exists(job_ligands_dir):
            print(f"\n⚠️ No ligands directory for job {job_id}")
            return {
                'converted': 0,
                'failed': 0,
                'skipped': 0,
                'total': 0,
                'failed_names': [],
                'skipped_names': []
            }
        
        # Create subdirectories
        sdf_dir = os.path.join(job_ligands_dir, 'sdf_ligands')
        pdbqt_dir = os.path.join(job_ligands_dir, 'pdbqt_ligands')
        os.makedirs(sdf_dir, exist_ok=True)
        os.makedirs(pdbqt_dir, exist_ok=True)
        
        print(f"\n{'-'*60}")
        print(f"Converting Job {job_id} Ligands")
        print(f"{'-'*60}")
        
        # Find all .smi files
        smi_files = [f for f in os.listdir(job_ligands_dir) if f.endswith('.smi')]
        
        if not smi_files:
            print(f"  No .smi files found")
            return {
                'converted': 0,
                'failed': 0,
                'skipped': 0,
                'total': 0,
                'failed_names': [],
                'skipped_names': []
            }
        print(f"  Found {len(smi_files)} ligand files")
        
        converted = 0
        failed = 0
        skipped = 0
        failed_ligands = []
        skipped_ligands = []
        
        for smi_file in sorted(smi_files):
            smi_path = os.path.join(job_ligands_dir, smi_file)
            base_name = smi_file.replace('.smi', '')
            
            # Try multiple strategies in order
            success = False
            
            try:
                # Read SMILES
                with open(smi_path, 'r') as f:
                    smiles = f.readline().strip()
                
                if not smiles:
                    print(f"  ❌ {base_name}: Empty SMILES")
                    failed += 1
                    failed_ligands.append(base_name)
                    continue
                
                # Quick check: if molecule is too large, skip expensive strategies
                num_atoms = 0  # Initialize
                try:
                    mol = Chem.MolFromSmiles(smiles, sanitize=False)
                    if mol:
                        num_atoms = mol.GetNumAtoms()
                        
                        # Skip ultra-large molecules (> 300 atoms)
                        if num_atoms > 300:
                            print(f"  ⏭️  {base_name}: Too large ({num_atoms} atoms), skipping")
                            print(f"      Use external tools for molecules > 300 atoms")
                            skipped += 1
                            skipped_ligands.append(base_name)
                            self._save_diagnostic_info(smiles, base_name, job_ligands_dir, 
                                                      reason="Too large for automatic 3D generation")
                            continue
                except:
                    num_atoms = 0  # Reset on error
                
                # Strategy 1: RDKit 3D + MMFF optimization
                success = self._try_rdkit_3d_generation(smiles, base_name, sdf_dir, pdbqt_dir)
                
                # Strategy 2: RDKit 3D with ETKDGv3 (better for complex molecules)
                if not success:
                    success = self._try_etkdg_generation(smiles, base_name, sdf_dir, pdbqt_dir)
                
                # Strategy 3: RDKit basic 3D (no optimization)
                if not success:
                    success = self._try_basic_3d_generation(smiles, base_name, sdf_dir, pdbqt_dir)
                
                # Strategy 4: Direct OpenBabel conversion (last resort)
                if not success:
                    success = self._try_openbabel_direct(smi_path, base_name, pdbqt_dir)
                
                # Strategy 5: Unsanitized parsing (only for moderate-sized molecules or unknown size)
                if not success and (num_atoms == 0 or num_atoms <= 200):
                    success = self._try_unsanitized_parsing(smiles, base_name, sdf_dir, pdbqt_dir)
                
                if success:
                    converted += 1
                else:
                    failed += 1
                    failed_ligands.append(base_name)
                    print(f"  ❌ {base_name}: All conversion strategies failed")
                    # Save diagnostic info
                    self._save_diagnostic_info(smiles, base_name, job_ligands_dir,
                                              reason="All 3D generation strategies failed")
                    
            except Exception as e:
                print(f"  ❌ {base_name}: Unexpected error: {e}")
                failed += 1
                failed_ligands.append(base_name)
        
        print(f"\n  ✓ Converted: {converted}/{len(smi_files)}")
        if failed > 0:
            print(f"  ❌ Failed: {failed}")
            print(f"     Failed ligands: {', '.join(failed_ligands)}")
        if skipped > 0:
            print(f"  ⏭️  Skipped: {skipped} (too large for automatic conversion)")
            print(f"     Skipped ligands: {', '.join(skipped_ligands)}")
            print(f"     Check *_FAILED.txt files for external tool suggestions")
        
        print(f"  SDF files: {sdf_dir}")
        print(f"  PDBQT files: {pdbqt_dir}")
        
        return {
            'converted': converted, 
            'failed': failed, 
            'skipped': skipped,
            'total': len(smi_files), 
            'failed_names': failed_ligands,
            'skipped_names': skipped_ligands
        }
    
    @staticmethod
    def _strip_dummy_atoms(smiles):
        """
        Remove any wildcard/dummy [*] atoms from SMILES before 3D generation.
        These cause UFFTYPER warnings in OpenBabel and are meaningless for
        geometry optimisation (they should have been replaced by actual fragments).
        Replaces [*] with a hydrogen placeholder then round-trips through RDKit
        to produce a clean, sanitised SMILES.
        """
        if '[*]' not in smiles and '*' not in smiles:
            return smiles
        try:
            # Replace bare/bracketed wildcard with explicit H
            clean = smiles.replace('[*]', '[H]').replace('*', '[H]')
            mol = Chem.MolFromSmiles(clean)
            if mol is None:
                # Fall back: just strip them
                clean = smiles.replace('[*]', '').replace('*', '')
                mol = Chem.MolFromSmiles(clean)
            return Chem.MolToSmiles(mol) if mol is not None else smiles
        except Exception:
            return smiles

    def _prepare_ligand_3d_mol(self, smiles, allow_random_coords=False, optimize=True):
        """
        Build a fresh 3D conformer for a generated ligand.
        Uses a clean-graph → ETKDGv3 → force-field workflow, which is much more
        stable for puckered rings and heterocycles than plain EmbedMolecule.
        """
        smiles = self._strip_dummy_atoms(smiles)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        canonical = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        mol = Chem.MolFromSmiles(canonical)
        if mol is None:
            return None

        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        mol = Chem.AddHs(mol, addCoords=False)

        embed_params = []

        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.enforceChirality = True
        params.useSmallRingTorsions = True
        params.useMacrocycleTorsions = True
        params.useMacrocycle14config = True
        embed_params.append(params)

        params2 = AllChem.ETKDG()
        params2.randomSeed = 42
        params2.enforceChirality = True
        embed_params.append(params2)

        if allow_random_coords:
            params3 = AllChem.ETKDGv3()
            params3.randomSeed = 42
            params3.enforceChirality = True
            params3.useSmallRingTorsions = True
            params3.useMacrocycleTorsions = True
            params3.useMacrocycle14config = True
            params3.useRandomCoords = True
            embed_params.append(params3)

        embedded = False
        for params in embed_params:
            mol.RemoveAllConformers()
            if AllChem.EmbedMolecule(mol, params) == 0:
                embedded = True
                break

        if not embedded:
            return None

        if optimize:
            optimized = False
            try:
                if AllChem.MMFFHasAllMoleculeParams(mol):
                    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
                    optimized = True
            except Exception:
                optimized = False

            if not optimized:
                try:
                    AllChem.UFFOptimizeMolecule(mol, maxIters=500)
                except Exception:
                    pass

        try:
            Chem.AssignStereochemistryFrom3D(mol, replaceExistingTags=True)
        except Exception:
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        return mol

    @staticmethod
    def _write_ligand_3d_files(mol, base_name, sdf_dir, pdbqt_dir):
        """Write a prepared 3D ligand to SDF, then convert it to PDBQT."""
        sdf_path = os.path.join(sdf_dir, f'{base_name}.sdf')
        writer = Chem.SDWriter(sdf_path)
        writer.write(mol)
        writer.close()

        pdbqt_path = os.path.join(pdbqt_dir, f'{base_name}.pdbqt')
        cmd = ['obabel', sdf_path, '-O', pdbqt_path, '-xh', '--partialcharge', 'gasteiger']
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        return True

    def _try_rdkit_3d_generation(self, smiles, base_name, sdf_dir, pdbqt_dir):
        """Strategy 1: ETKDGv3 + force-field cleanup for generated ligands."""
        try:
            mol = self._prepare_ligand_3d_mol(
                smiles,
                allow_random_coords=False,
                optimize=True
            )
            if mol is None:
                return False
            return self._write_ligand_3d_files(mol, base_name, sdf_dir, pdbqt_dir)
        except Exception:
            return False
    
    def _try_etkdg_generation(self, smiles, base_name, sdf_dir, pdbqt_dir):
        """Strategy 2: ETKDGv3 with random-coordinate fallback."""
        try:
            mol = self._prepare_ligand_3d_mol(
                smiles,
                allow_random_coords=True,
                optimize=True
            )
            if mol is None:
                return False
            return self._write_ligand_3d_files(mol, base_name, sdf_dir, pdbqt_dir)
        except Exception:
            return False
    
    def _try_basic_3d_generation(self, smiles, base_name, sdf_dir, pdbqt_dir):
        """Strategy 3: ETKDGv3 embedding without force-field optimization."""
        try:
            mol = self._prepare_ligand_3d_mol(
                smiles,
                allow_random_coords=True,
                optimize=False
            )
            if mol is None:
                return False
            return self._write_ligand_3d_files(mol, base_name, sdf_dir, pdbqt_dir)
        except Exception:
            return False
    
    def _try_openbabel_direct(self, smi_path, base_name, pdbqt_dir):
        """Strategy 4: Direct OpenBabel conversion (last resort)."""
        try:
            pdbqt_path = os.path.join(pdbqt_dir, f'{base_name}.pdbqt')
            
            # Direct SMI → PDBQT with OpenBabel's 3D generation
            cmd = [
                'obabel',
                smi_path,
                '-O', pdbqt_path,
                '--gen3d',  # OpenBabel's 3D generation
                '-xh',
                '--partialcharge', 'gasteiger'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Verify the output file exists and has content
            if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 0:
                return True
            
            return False
            
        except Exception:
            return False
    
    def _try_unsanitized_parsing(self, smiles, base_name, sdf_dir, pdbqt_dir):
        """Strategy 5: Parse without sanitization for problematic SMILES."""
        try:
            smiles = self._strip_dummy_atoms(smiles)
            # Try parsing without sanitization
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is None:
                return False
            
            # Check size - skip if too large
            if mol.GetNumAtoms() > 200:
                return False
            
            # Try to fix common valence issues
            try:
                # Update explicit valence
                for atom in mol.GetAtoms():
                    atom.UpdatePropertyCache(strict=False)
                
                # Attempt sanitization in steps
                Chem.SanitizeMol(mol, 
                    sanitizeOps=Chem.SanitizeFlags.SANITIZE_FINDRADICALS |
                                Chem.SanitizeFlags.SANITIZE_SETAROMATICITY |
                                Chem.SanitizeFlags.SANITIZE_SETCONJUGATION |
                                Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION |
                                Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                    catchErrors=True)
            except:
                pass  # Continue with unsanitized molecule
            
            # Add hydrogens
            try:
                mol = Chem.AddHs(mol, addCoords=False)
            except:
                return False
            
            # Try embedding with very permissive parameters (but limited attempts)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.useRandomCoords = True
            params.useBasicKnowledge = False
            params.enforceChirality = False  # Don't enforce chirality
            params.useExpTorsionAnglePrefs = False  # Don't use torsion prefs
            
            result = AllChem.EmbedMolecule(mol, params)
            if result != 0:
                return False
            
            # Save and convert
            sdf_path = os.path.join(sdf_dir, f'{base_name}.sdf')
            writer = Chem.SDWriter(sdf_path)
            writer.write(mol)
            writer.close()
            
            pdbqt_path = os.path.join(pdbqt_dir, f'{base_name}.pdbqt')
            cmd = ['obabel', sdf_path, '-O', pdbqt_path, '-xh', '--partialcharge', 'gasteiger']
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            return True
            
        except Exception:
            return False
    
    
    def _save_diagnostic_info(self, smiles, base_name, job_dir, reason="Unknown"):
        """Save diagnostic information for failed conversions."""
        try:
            diag_file = os.path.join(job_dir, f'{base_name}_FAILED.txt')
            
            with open(diag_file, 'w') as f:
                f.write(f"Failed to convert: {base_name}\n")
                f.write(f"Reason: {reason}\n")
                f.write(f"SMILES: {smiles}\n\n")
                
                # Try to get molecular properties
                try:
                    mol = Chem.MolFromSmiles(smiles, sanitize=False)
                    if mol:
                        f.write(f"Number of atoms: {mol.GetNumAtoms()}\n")
                        f.write(f"Number of bonds: {mol.GetNumBonds()}\n")
                        
                        # Check for problematic atoms
                        problematic = []
                        for atom in mol.GetAtoms():
                            if atom.GetFormalCharge() != 0:
                                problematic.append(f"Atom {atom.GetIdx()}: charge {atom.GetFormalCharge()}")
                        
                        if problematic:
                            f.write(f"\nProblematic atoms:\n")
                            for p in problematic:
                                f.write(f"  {p}\n")
                except Exception as e:
                    f.write(f"\nCould not parse SMILES: {e}\n")
                
                f.write(f"\nSuggestions:\n")
                f.write(f"1. This molecule may be too large/complex for automatic 3D generation\n")
                f.write(f"2. Try using external tools (Schrödinger LigPrep, OpenEye)\n")
                f.write(f"3. Consider breaking into smaller fragments\n")
                f.write(f"4. Use 2D-to-3D conversion in docking software directly\n")
                
        except:
            pass  # Don't fail if diagnostics can't be saved
    
    # ==================================================================
    # STEP 9: Binding Pocket Detection — P2Rank (primary) / fpocket (fallback)
    # ==================================================================

    def detect_binding_pockets(self):
        """
        Detect binding pockets using P2Rank.
        P2Rank is run as:  <install_dir>/prank predict -f enzyme.pdb -o output/
        The script auto-locates 'prank' in common directories or asks the user.
        """
        print("\n" + "="*60)
        print("STEP 9: Binding Pocket Detection (P2Rank)")
        print("="*60)

        if not self.enzyme_apo_file or not os.path.exists(self.enzyme_apo_file):
            print("\n⚠️  No enzyme PDB found — skipping pocket detection")
            return []

        pockets = self._run_p2rank()

        if pockets:
            print(f"\n✓ P2Rank found {len(pockets)} pocket(s)")
            self._display_pockets(pockets, tool="P2Rank")
            return pockets

        # P2Rank found but returned nothing, OR not found at all
        print("\n⚠️  P2Rank returned no pockets — proceeding with BLIND docking")
        print("   (entire protein bounding box will be used as search space)")
        return []

    # ------------------------------------------------------------------
    # P2Rank
    # ------------------------------------------------------------------

    def _find_prank_exe(self):
        """
        Locate the 'prank' executable.
        P2Rank is typically downloaded as p2rank_X.X.X/ and run as ./prank from inside it.
        Searches in order:
          1. self._p2rank_dir if already set
          2. Common locations: output_dir/../p2rank_*, ~/Downloads/p2rank_*
          3. System PATH (unlikely but possible)
          4. Interactively asks the user
        Returns the full path to prank, or None.
        """
        import glob

        # Already found previously in this session
        if hasattr(self, '_p2rank_dir') and self._p2rank_dir:
            exe = os.path.join(self._p2rank_dir, 'prank')
            if os.path.isfile(exe):
                return exe

        # Search common locations for p2rank_* directories
        search_roots = [
            self.output_dir,
            os.path.dirname(self.output_dir),
            os.path.expanduser('~/Downloads'),
            os.path.expanduser('~'),
        ]
        for root in search_roots:
            for d in glob.glob(os.path.join(root, 'p2rank*')):
                exe = os.path.join(d, 'prank')
                if os.path.isfile(exe):
                    self._p2rank_dir = d
                    return exe

        # System PATH
        for name in ['prank', 'p2rank', 'p2rank.sh']:
            try:
                r = subprocess.run([name, '--version'],
                                   capture_output=True, text=True, timeout=10)
                if r.returncode == 0 or 'p2rank' in (r.stdout+r.stderr).lower():
                    self._p2rank_dir = None
                    return name
            except Exception:
                pass

        # Ask interactively
        if self.interactive:
            print("\n  P2Rank not found automatically.")
            print("  Enter the full path to your p2rank directory")
            print("  (e.g. /Users/you/Downloads/p2rank_2.5.1)")
            path = input("  Path: ").strip()
            exe  = os.path.join(path, 'prank')
            if os.path.isfile(exe):
                self._p2rank_dir = path
                return exe
            print(f"  ❌ 'prank' not found at {exe}")

        return None

    def _check_p2rank(self):
        """Return True if prank executable can be located."""
        return self._find_prank_exe() is not None

    def _run_p2rank(self):
        """
        Run P2Rank exactly as the user does manually:
            cd <p2rank_dir> && ./prank predict -f ./enzyme_apo.pdb -o <output_dir>

        Key fixes:
        - Explicitly injects Java 17 into the subprocess environment
          (shell exports like 'export JAVA_HOME=...' don't reach Python subprocesses)
        - Runs prank from its own directory (required for prank to find its lib/)
        - Copies the PDB next to the prank script so relative paths work cleanly
        """
        import shutil

        exe = self._find_prank_exe()
        if not exe:
            return []

        p2rank_dir = getattr(self, '_p2rank_dir', None)
        if not p2rank_dir:
            p2rank_dir = os.path.dirname(exe)

        # ── Build a subprocess environment with Java 17 injected ──────
        env = os.environ.copy()

        # Homebrew arm64 (Apple Silicon) Java 17 path
        java17_candidates = [
            '/opt/homebrew/opt/openjdk@17/bin',       # Apple Silicon brew
            '/usr/local/opt/openjdk@17/bin',           # Intel brew
            '/opt/homebrew/opt/openjdk/bin',           # latest brew openjdk
            '/usr/local/opt/openjdk/bin',
        ]
        java17_bin = None
        for candidate in java17_candidates:
            if os.path.isfile(os.path.join(candidate, 'java')):
                java17_bin = candidate
                break

        if java17_bin:
            env['JAVA_HOME'] = os.path.dirname(java17_bin)  # parent of bin/
            env['PATH'] = java17_bin + os.pathsep + env.get('PATH', '')
            print(f"  Java 17: {java17_bin}/java")
        else:
            print("  ⚠️  Java 17 not found in Homebrew paths — using system Java")

        # ── Copy PDB into p2rank directory (mirrors manual workflow) ──
        pdb_in_p2rank = os.path.join(p2rank_dir, 'enzyme_apo.pdb')
        shutil.copy(self.enzyme_apo_file, pdb_in_p2rank)

        # Output goes into Downloads/p2rank_results/
        out_dir = os.path.join(self.output_dir, 'p2rank_results')
        os.makedirs(out_dir, exist_ok=True)

        # ── Run: ./prank predict -f ./enzyme_apo.pdb -o <out_dir> ────
        # Use absolute exe path, cwd=p2rank_dir so lib/ and distlib/ are found
        cmd = [exe, 'predict', '-f', pdb_in_p2rank, '-o', out_dir]
        print(f"  Command: {' '.join(cmd)}")
        print(f"  Working dir: {p2rank_dir}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,
                cwd=p2rank_dir,
                env=env
            )
        except subprocess.TimeoutExpired:
            print("  ❌ P2Rank timed out after 300 s")
            return []
        except FileNotFoundError as e:
            print(f"  ❌ Cannot run prank: {e}")
            return []

        combined = result.stdout + result.stderr

        # ── Detect Java version errors ─────────────────────────────────
        if any(x in combined for x in
               ['ClassLoader', 'NoClassDefFoundError',
                'UnsupportedClassVersionError', 'Error: Unable to initialize']):
            print("  ❌ P2Rank crashed — Java version problem.")
            print(f"  Java used: {java17_bin or 'system'}")
            print("\n  Permanent fix — add to ~/.zshrc:")
            print('    export PATH="/opt/homebrew/opt/openjdk@17/bin:$PATH"')
            print('    export JAVA_HOME="/opt/homebrew/opt/openjdk@17"')
            print("\n  Then: source ~/.zshrc && re-run pipeline")
            return []

        # ── Find predictions CSV ───────────────────────────────────────
        pred_csv = None
        for root, dirs, files in os.walk(out_dir):
            for fn in files:
                if fn.endswith('_predictions.csv'):
                    pred_csv = os.path.join(root, fn)
                    break
            if pred_csv:
                break

        if not pred_csv:
            print(f"  ❌ Predictions CSV not found under {out_dir}")
            if combined.strip():
                print(f"  P2Rank output:\n  {combined.strip()[-500:]}")
            return []

        print(f"  ✓ Predictions CSV: {pred_csv}")
        return self._parse_p2rank_csv(pred_csv, os.path.dirname(pred_csv))

    def _parse_p2rank_csv(self, pred_csv, p2rank_out_dir):
        """
        Parse P2Rank *_predictions.csv.

        Actual P2Rank 2.x header (comma-separated, headers have leading spaces):
          name , rank , score , probability , sas_points , surf_atoms ,
          center_x , center_y , center_z , residue_ids , surf_atom_ids
        """
        import math

        pockets = []
        try:
            with open(pred_csv, newline='') as f:
                raw = f.read()

            lines = [l for l in raw.splitlines() if l.strip()]
            if not lines:
                print("  ⚠️  Predictions CSV is empty")
                return []

            # Parse header — strip all whitespace from each column name
            header = [h.strip() for h in lines[0].split(',')]

            for line in lines[1:]:
                vals = [v.strip() for v in line.split(',')]
                if len(vals) < len(header):
                    continue
                row = dict(zip(header, vals))

                try:
                    rank  = int(row.get('rank', 0))
                    score = float(row.get('score', 0))
                    prob  = float(row.get('probability', 0))
                    cx    = float(row.get('center_x', 0))
                    cy    = float(row.get('center_y', 0))
                    cz    = float(row.get('center_z', 0))
                    n_surf = int(float(row.get('surf_atoms',
                                   row.get('sas_points', 10))))

                    # ── Cα-distance box: use P2Rank residue_ids to find pocket
                    # residues, extract their Cα coordinates from the enzyme PDB,
                    # build a full pairwise distance matrix, and set the box
                    # dimension = the longest Cα–Cα distance found (no hard cap).
                    residue_ids_str = row.get('residue_ids', '')
                    box_dim = self._calpha_box_from_residues(
                        residue_ids_str, self.enzyme_apo_file
                    )
                    sx = sy = sz = box_dim

                    pockets.append({
                        'id':                 rank,
                        'name':               f'P2Rank Pocket {rank}',
                        'center_x':           cx,
                        'center_y':           cy,
                        'center_z':           cz,
                        'size_x':             round(sx, 1),
                        'size_y':             round(sy, 1),
                        'size_z':             round(sz, 1),
                        'volume':             score,
                        'druggability_score': prob,
                        'n_atoms':            n_surf,
                        'residue_ids':        residue_ids_str,
                    })
                except (ValueError, KeyError) as e:
                    continue

        except Exception as e:
            print(f"  ❌ CSV parse error: {e}")
            return []

        # Sort best (highest probability) first
        pockets.sort(key=lambda p: p['druggability_score'], reverse=True)

        if pockets:
            print(f"  Parsed {len(pockets)} pockets from CSV")
            print(f"  Best pocket: rank={pockets[0]['id']}, "
                  f"prob={pockets[0]['druggability_score']:.3f}, "
                  f"center=({pockets[0]['center_x']:.2f}, "
                  f"{pockets[0]['center_y']:.2f}, "
                  f"{pockets[0]['center_z']:.2f})")
        return pockets

    # ------------------------------------------------------------------
    # fpocket (fallback)
    # ------------------------------------------------------------------

    def _check_fpocket(self):
        """Return True if fpocket is on PATH."""
        try:
            r = subprocess.run(['fpocket', '--help'],
                               capture_output=True, text=True, timeout=5)
            return True
        except Exception:
            return False

    def _run_fpocket(self):
        """Run fpocket and parse its _info.txt output."""
        import shutil

        fp_dir = os.path.join(self.output_dir, 'fpocket_results')
        os.makedirs(fp_dir, exist_ok=True)

        pdb_copy = os.path.join(fp_dir, 'enzyme.pdb')
        shutil.copy(self.enzyme_apo_file, pdb_copy)

        try:
            subprocess.run(['fpocket', '-f', pdb_copy],
                           capture_output=True, text=True,
                           timeout=300, cwd=fp_dir)
        except subprocess.TimeoutExpired:
            print("  ❌ fpocket timed out")
            return []
        except FileNotFoundError:
            return []

        return self._parse_fpocket_results(fp_dir)

    def _parse_fpocket_results(self, fp_dir):
        """
        Parse fpocket _info.txt.

        Each pocket block:
            Pocket N :
                Score :                        0.39
                Druggability Score :           0.39
                Real volume (approximation) :  282.59
                x :   3.26
                y :   11.72
                z :   -4.90
        """
        import re

        pockets = []
        out_dirs = [d for d in os.listdir(fp_dir) if d.endswith('_out')]
        if not out_dirs:
            return pockets

        out_dir = os.path.join(fp_dir, out_dirs[0])
        info_files = [f for f in os.listdir(out_dir) if f.endswith('_info.txt')]
        if not info_files:
            return pockets

        info_path = os.path.join(out_dir, info_files[0])
        try:
            content = open(info_path).read()
            blocks  = re.split(r'Pocket\s+(\d+)\s+:', content)

            i = 1
            while i + 1 < len(blocks):
                pid  = int(blocks[i])
                body = blocks[i + 1]
                i   += 2

                p = {'id': pid, 'name': f'fpocket Pocket {pid}',
                     'size_x': 20.0, 'size_y': 20.0, 'size_z': 20.0}

                for line in body.splitlines():
                    line = line.strip()
                    if ':' not in line:
                        continue
                    k, _, v = line.partition(':')
                    k = k.strip().lower(); v = v.strip()
                    try:
                        if k == 'x':
                            p['center_x'] = float(v)
                        elif k == 'y':
                            p['center_y'] = float(v)
                        elif k == 'z':
                            p['center_z'] = float(v)
                        elif 'real volume' in k:
                            p['volume'] = float(v.split()[0])
                        elif 'druggability score' in k:
                            p['druggability_score'] = float(v)
                        elif k == 'score' and 'druggability_score' not in p:
                            p['druggability_score'] = float(v)
                    except (ValueError, IndexError):
                        pass

                if 'center_x' not in p:
                    # Derive from pocket PDB atoms
                    ppdb = os.path.join(out_dir, 'pockets',
                                        f'pocket{pid}_atm.pdb')
                    if os.path.exists(ppdb):
                        cx, cy, cz = self._centroid_from_pdb(ppdb)
                        p['center_x'] = cx; p['center_y'] = cy; p['center_z'] = cz
                        sx, sy, sz = self._box_dims_from_pdb(ppdb, padding=6.0)
                        p['size_x'] = sx
                        p['size_y'] = sy
                        p['size_z'] = sz

                if 'center_x' in p:
                    p.setdefault('volume', 0.0)
                    p.setdefault('druggability_score', 0.0)
                    p.setdefault('n_atoms', 0)

                    # Derive proper box from pocket PDB if available
                    ppdb = os.path.join(out_dir, 'pockets',
                                        f'pocket{pid}_atm.pdb')
                    if os.path.exists(ppdb) and p['size_x'] == 20.0:
                        sx, sy, sz = self._box_dims_from_pdb(ppdb, padding=6.0)
                        p['size_x'] = sx
                        p['size_y'] = sy
                        p['size_z'] = sz

                    pockets.append(p)

        except Exception as e:
            print(f"  ⚠️  fpocket parse error: {e}")

        pockets.sort(key=lambda p: p.get('druggability_score', 0), reverse=True)
        return pockets

    def _calpha_box_from_residues(self, residue_ids_str, enzyme_pdb,
                                   padding=4.0, fallback=20.0):
        """
        Compute docking-box dimension from P2Rank pocket residues.

        Strategy
        --------
        1. Parse the residue_ids string from P2Rank CSV
           (format: "A_45 A_47 B_102 ..." — chain_resnum tokens).
        2. Walk the enzyme PDB and collect the (x, y, z) of every CA
           (C-alpha) atom whose (chain, resnum) matches a pocket residue.
        3. Build the full pairwise Euclidean distance matrix for those CAs.
        4. The longest CA–CA distance found becomes the box side length,
           with a small padding added on each side (total += 2*padding/2).
           There is intentionally NO upper cap — the box adapts to
           whatever the pocket geometry demands.
        5. Falls back to `fallback` (Å) when residue_ids is absent,
           unparseable, or no matching CA atoms exist in the PDB.

        Parameters
        ----------
        residue_ids_str : str
            Raw value of the 'residue_ids' column, e.g. "A_45 A_47 B_102".
        enzyme_pdb : str
            Path to the apo-enzyme PDB file.
        padding : float
            Angstroms added to the longest distance (accounts for atom
            radius + a small solvent shell). Default 4 Å.
        fallback : float
            Box size returned when the calculation cannot be performed.

        Returns
        -------
        float
            Box side length in Å (longest CA–CA distance + padding).
        """
        import math

        # ── 1. Parse residue tokens: "A_45" → ('A', 45) ──────────────
        pocket_res = set()
        for token in residue_ids_str.split():
            token = token.strip()
            if not token:
                continue
            if '_' in token:
                parts = token.split('_')
                chain  = parts[0].upper()
                resnum_str = parts[1].lstrip('0') or '0'
                try:
                    resnum = int(resnum_str)
                    pocket_res.add((chain, resnum))
                except ValueError:
                    pass
            else:
                # No chain prefix — match any chain
                try:
                    pocket_res.add(('*', int(token)))
                except ValueError:
                    pass

        if not pocket_res:
            print(f"  ⚠️  No residue_ids parsed — using fallback box {fallback} Å")
            return fallback

        # ── 2. Extract Cα coords from PDB ─────────────────────────────
        ca_coords = []
        try:
            for line in open(enzyme_pdb):
                if not line.startswith('ATOM'):
                    continue
                atom_name = line[12:16].strip()
                if atom_name != 'CA':
                    continue
                chain  = line[21].strip().upper()
                try:
                    resnum = int(line[22:26].strip())
                    x      = float(line[30:38])
                    y      = float(line[38:46])
                    z      = float(line[46:54])
                except ValueError:
                    continue
                # Match against pocket residue set
                if (chain, resnum) in pocket_res or ('*', resnum) in pocket_res:
                    ca_coords.append((x, y, z))
        except Exception as e:
            print(f"  ⚠️  PDB read error during Cα extraction: {e}")
            return fallback

        if len(ca_coords) < 2:
            print(f"  ⚠️  Fewer than 2 Cα atoms found for pocket "
                  f"({len(ca_coords)} matched) — using fallback box {fallback} Å")
            return fallback

        # ── 3. Pairwise distance matrix → find maximum ────────────────
        max_dist = 0.0
        n = len(ca_coords)
        for i in range(n):
            xi, yi, zi = ca_coords[i]
            for j in range(i + 1, n):
                xj, yj, zj = ca_coords[j]
                d = math.sqrt((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)
                if d > max_dist:
                    max_dist = d

        box_dim = max_dist + padding
        print(f"  ✓ Cα-distance box: {n} pocket residues, "
              f"max Cα–Cα = {max_dist:.2f} Å → box = {box_dim:.1f} Å "
              f"(+ {padding} Å padding)")
        return round(box_dim, 1)

    def _centroid_from_pdb(self, pdb_file):
        """Return (cx, cy, cz) centroid of ATOM/HETATM coords in a PDB."""
        xs, ys, zs = [], [], []
        try:
            for line in open(pdb_file):
                if line.startswith(('ATOM', 'HETATM')):
                    try:
                        xs.append(float(line[30:38]))
                        ys.append(float(line[38:46]))
                        zs.append(float(line[46:54]))
                    except ValueError:
                        pass
        except Exception:
            pass
        if xs:
            return sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)
        return 0.0, 0.0, 0.0

    def _box_dims_from_pdb(self, pdb_file, padding=6.0):
        """Return (sx, sy, sz) bounding-box dimensions from a pocket PDB."""
        xs, ys, zs = [], [], []
        try:
            for line in open(pdb_file):
                if line.startswith(('ATOM', 'HETATM')):
                    try:
                        xs.append(float(line[30:38]))
                        ys.append(float(line[38:46]))
                        zs.append(float(line[46:54]))
                    except ValueError:
                        pass
        except Exception:
            pass
        if xs:
            return (max(xs)-min(xs)+padding,
                    max(ys)-min(ys)+padding,
                    max(zs)-min(zs)+padding)
        return 20.0, 20.0, 20.0

    # ------------------------------------------------------------------
    # Display & selection helpers
    # ------------------------------------------------------------------

    def _display_pockets(self, pockets, tool=""):
        """Print a formatted table of detected pockets."""
        label = f" ({tool})" if tool else ""
        print("\n" + "-"*68)
        print(f"Detected Pockets{label} — sorted best first:")
        print("-"*68)
        print(f"{'Rank':<6} {'Center (x, y, z)':<34} {'Score':<10} {'Box (Å)'}")
        print("-"*68)
        for p in pockets[:10]:   # show top 10
            coord = (f"({p['center_x']:.2f}, "
                     f"{p['center_y']:.2f}, "
                     f"{p['center_z']:.2f})")
            score = (f"{p.get('druggability_score', p.get('volume', 0)):.3f}")
            box   = (f"{p['size_x']:.0f}×{p['size_y']:.0f}×{p['size_z']:.0f}")
            print(f"{p['id']:<6} {coord:<34} {score:<10} {box}")
        if len(pockets) > 10:
            print(f"  ... and {len(pockets)-10} more pockets")
        print("-"*68)

    def select_best_pocket(self, pockets):
        """Return the top-ranked pocket (index 0 after sort)."""
        if not pockets:
            return None
        best = pockets[0]
        print(f"\n✓ Best pocket: Pocket {best['id']}")
        print(f"  Center : ({best['center_x']:.3f}, "
              f"{best['center_y']:.3f}, {best['center_z']:.3f})")
        print(f"  Box    : {best['size_x']:.0f} × "
              f"{best['size_y']:.0f} × {best['size_z']:.0f} Å")
        score_key = 'druggability_score' if best.get('druggability_score') else 'volume'
        print(f"  Score  : {best.get(score_key, 0):.3f}")
        return best


    def run_autodock_vina_docking(self, pockets, pdbqt_results):
        """Run AutoDock Vina for all ligands. Auto-falls back to blind docking."""
        print("\n" + "="*60)
        print("STEP 10: AutoDock Vina Molecular Docking")
        print("="*60)

        if not pdbqt_results or not pdbqt_results.get('enzyme_pdbqt'):
            print("\n⚠️ No enzyme PDBQT file, skipping docking")
            return None

        if not self._check_vina_installation():
            print("\n❌ AutoDock Vina not found!")
            print("   Install: conda install -c conda-forge vina")
            print("\n💡 Saving shell scripts for manual execution...")
            return self._generate_docking_scripts(pockets or [], pdbqt_results)

        # Decide: P2Rank/fpocket-guided or blind?
        pocket_ok, reason = self._validate_pockets(pockets or [])
        if not pocket_ok:
            print(f"\n⚠️ Pocket validation failed: {reason}")
            print("   → Switching to BLIND docking (covers entire protein)")
            blind = True
            # Create a single synthetic pocket entry from protein bbox
            cx, cy, cz, sx, sy, sz = self._get_protein_bbox(
                pdbqt_results['enzyme_pdbqt'])
            pocket = {'id': 1, 'name': 'Blind Docking',
                      'center_x': cx, 'center_y': cy, 'center_z': cz,
                      'size_x': sx, 'size_y': sy, 'size_z': sz,
                      'druggability_score': None, 'volume': None}
            selected_pockets = [pocket]
        else:
            blind = False
            selected_pockets = self._select_pockets_for_docking(pockets)
            if not selected_pockets:
                return None

        docking_results = {}
        for pocket in selected_pockets:
            print(f"\n{'='*60}")
            label = 'BLIND' if blind else f"Pocket {pocket['id']}"
            print(f"Docking at {label}: {pocket.get('name','')}")
            print(f"  Center: ({pocket['center_x']:.3f}, "
                  f"{pocket['center_y']:.3f}, {pocket['center_z']:.3f})")
            print(f"  Box:    {int(pocket['size_x'])} × "
                  f"{int(pocket['size_y'])} × {int(pocket['size_z'])} Å")
            print(f"{'='*60}")

            pocket_results = self._dock_at_pocket(pocket, pdbqt_results, blind=blind)
            docking_results[pocket['id']] = pocket_results

        self._display_docking_summary(docking_results)
        return docking_results

    def _check_vina_installation(self):
        """Check if AutoDock Vina is on PATH."""
        try:
            subprocess.run(['vina', '--version'],
                           capture_output=True, text=True, timeout=5)
            return True
        except Exception:
            return False

    def _select_pockets_for_docking(self, pockets):
        """
        Automatically dock at ALL pockets to find the best binding site.
        No user prompts - fully automatic.
        """
        if not pockets:
            return []
        
        # Dock at ALL pockets — will select best based on actual affinities later
        print(f"\n  → Docking at all {len(pockets)} pocket(s) to find best binding site")
        return pockets

    def _dock_at_pocket(self, pocket, pdbqt_results, blind=False):
        """Dock all ligands at one pocket, writing a config.txt per ligand."""
        enzyme_pdbqt = pdbqt_results['enzyme_pdbqt']
        pocket_dir   = os.path.join(self.output_dir,
                                    f'docking_pocket_{pocket["id"]}')
        # Remove stale files from previous runs to avoid mixing old logs/results.
        import shutil
        if os.path.exists(pocket_dir):
            shutil.rmtree(pocket_dir)
        os.makedirs(pocket_dir, exist_ok=True)
        print(f"\n  Output directory: {pocket_dir}")

        results = {'pocket_info': pocket, 'jobs': {},
                   'total_docked': 0, 'total_failed': 0}

        vina_broken    = False
        global_first   = True  # print config for very first ligand overall

        for job_id, job_stats in pdbqt_results.get('jobs', {}).items():
            if job_stats.get('converted', 0) == 0:
                continue

            pdbqt_dir = os.path.join(self.output_dir,
                                     f'job{job_id}_ligands', 'pdbqt_ligands')
            if not os.path.exists(pdbqt_dir):
                continue

            ligand_files = sorted(
                f for f in os.listdir(pdbqt_dir) if f.endswith('.pdbqt'))
            if not ligand_files:
                continue

            job_dir = os.path.join(pocket_dir, f'job{job_id}')
            os.makedirs(job_dir, exist_ok=True)

            print(f"\n{'-'*60}")
            print(f"Docking Job {job_id} Ligands ({len(ligand_files)} total)")
            print(f"{'-'*60}")

            # num_modes = number of substitution positions the user selected.
            # Each job corresponds to one clicked atom position in the viewer.
            num_modes = max(1, len(self.substitution_jobs))
            print(f"  Binding poses per ligand (num_modes): {num_modes}")

            docked = 0
            failed = 0
            first_err_shown = False

            for ligand_file in ligand_files:
                ligand_path = os.path.join(pdbqt_dir, ligand_file)
                ligand_name = ligand_file.replace('.pdbqt', '')
                print(f"  {ligand_name}...", end=' ', flush=True)

                try:
                    ok, err = self._run_vina_docking(
                        enzyme_pdbqt, ligand_path, pocket,
                        job_dir, ligand_name,
                        blind=blind, first_ligand=global_first,
                        num_modes=num_modes, energy_range=3)
                    global_first = False
                except FileNotFoundError:
                    vina_broken = True
                    break

                if ok:
                    print("✓")
                    docked += 1
                else:
                    print("❌")
                    failed += 1
                    if not first_err_shown and err:
                        # Show meaningful lines from the error (skip help text)
                        meaningful = [ln for ln in err.splitlines()
                                      if ln.strip() and
                                      not ln.strip().startswith('--') and
                                      'optional' not in ln.lower() and
                                      'display' not in ln.lower()]
                        if meaningful:
                            print(f"\n  ── Vina error ──")
                            print('\n'.join(f"  {l}" for l in meaningful[:10]))
                            print(f"  (full details in {ligand_name}_error.txt)")
                            print(f"  ────────────────\n")
                        first_err_shown = True

            if vina_broken:
                break

            results['jobs'][job_id] = {
                'docked': docked, 'failed': failed,
                'total': len(ligand_files), 'output_dir': job_dir}
            results['total_docked'] += docked
            results['total_failed'] += failed
            print(f"\n  ✓ Docked: {docked}/{len(ligand_files)}"
                  + (f"  ❌ Failed: {failed}" if failed else ""))

        return results
    
    def _validate_pockets(self, pockets):
        """
        Return (ok, reason).
        A pocket set is valid if ANY pocket has volume > 50 Å³ OR atom_count > 3.
        Some tools give atom counts but no volume — use atom count as fallback.
        """
        if not pockets:
            return False, "No pockets found"
        # Check volume first
        max_vol = max(p.get('volume', 0) or 0 for p in pockets)
        if max_vol >= 50:
            return True, "OK"
        # Fallback: accept if any pocket has >= 4 atoms
        max_atoms = max((p.get('atom_count', 0) or 0) + (p.get('n_atoms', 0) or 0) for p in pockets)
        if max_atoms >= 4:
            return True, "OK (using atom count)"
        return False, (f"No usable pocket data (max vol={max_vol:.1f} Å³, "
                       f"max atoms={max_atoms}). Switching to blind docking.")

    def _write_vina_config(self, config_path, receptor, ligand, out_file,
                           pocket, exhaustiveness=6, num_modes=3, energy_range=3):
        """Write vina config.txt. Only 'out' is needed — no 'log' (unsupported in v1.2.5).

        energy_range: max kcal/mol difference between best and worst reported pose.
        Vina only outputs a conformation if it lies within this window of the best pose
        AND is geometrically distinct (RMSD > ~1 Å).  Setting this to 3 instead of the
        Vina default of 2 gives small/rigid ligands enough room to produce a second or
        third distinct pose — otherwise those conformation columns would always be empty.
        """
        sx = int(round(pocket.get('size_x', 20)))
        sy = int(round(pocket.get('size_y', 20)))
        sz = int(round(pocket.get('size_z', 20)))
        with open(config_path, 'w') as f:
            f.write(f"receptor = {receptor}\n")
            f.write(f"ligand = {ligand}\n\n")
            f.write(f"center_x = {pocket['center_x']:.3f}\n")
            f.write(f"center_y = {pocket['center_y']:.3f}\n")
            f.write(f"center_z = {pocket['center_z']:.3f}\n\n")
            f.write(f"size_x = {sx}\n")
            f.write(f"size_y = {sy}\n")
            f.write(f"size_z = {sz}\n\n")
            f.write(f"exhaustiveness = {exhaustiveness}\n")
            f.write(f"num_modes = {num_modes}\n")
            f.write(f"energy_range = {energy_range}\n\n")
            f.write(f"out = {out_file}\n")

    def _get_protein_bbox(self, receptor_path):
        """Return (cx,cy,cz,sx,sy,sz) bounding box of the receptor."""
        xs, ys, zs = [], [], []
        for src in [receptor_path.replace('.pdbqt', '.pdb'), receptor_path]:
            if not os.path.exists(src):
                continue
            with open(src, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        try:
                            xs.append(float(line[30:38]))
                            ys.append(float(line[38:46]))
                            zs.append(float(line[46:54]))
                        except:
                            pass
            if xs:
                break
        if xs:
            cx = (max(xs) + min(xs)) / 2
            cy = (max(ys) + min(ys)) / 2
            cz = (max(zs) + min(zs)) / 2
            sx = max(xs) - min(xs) + 10
            sy = max(ys) - min(ys) + 10
            sz = max(zs) - min(zs) + 10
        else:
            cx, cy, cz, sx, sy, sz = 0.0, 0.0, 0.0, 60.0, 60.0, 60.0
        return cx, cy, cz, sx, sy, sz

    def _write_blind_vina_config(self, config_path, receptor, ligand, out_file,
                                 exhaustiveness=6, num_modes=3, energy_range=3):
        """Write blind-docking config.txt — box covers entire protein. No --log.

        energy_range=3: same reasoning as _write_vina_config — allows Vina to
        report additional distinct poses for rigid/small ligands that would
        otherwise produce empty Conformation_N columns.
        """
        cx, cy, cz, sx, sy, sz = self._get_protein_bbox(receptor)
        with open(config_path, 'w') as f:
            f.write(f"receptor = {receptor}\n")
            f.write(f"ligand = {ligand}\n\n")
            f.write(f"# Blind docking — box covers entire protein\n")
            f.write(f"center_x = {cx:.3f}\n")
            f.write(f"center_y = {cy:.3f}\n")
            f.write(f"center_z = {cz:.3f}\n\n")
            f.write(f"size_x = {int(sx)}\n")
            f.write(f"size_y = {int(sy)}\n")
            f.write(f"size_z = {int(sz)}\n\n")
            f.write(f"exhaustiveness = {exhaustiveness}\n")
            f.write(f"num_modes = {num_modes}\n")
            f.write(f"energy_range = {energy_range}\n\n")
            f.write(f"out = {out_file}\n")
        return {'center_x': cx, 'center_y': cy, 'center_z': cz,
                'size_x': sx, 'size_y': sy, 'size_z': sz}

    def _run_vina_docking(self, receptor, ligand, pocket, output_dir,
                          ligand_name, blind=False, first_ligand=False,
                          num_modes=1, energy_range=3):
        """
        Run AutoDock Vina v1.2.5.
        - ALL params in config.txt (except --log which v1.2.5 does NOT support)
        - vina stdout/stderr captured by Python and saved as the log file
        - num_modes: number of binding poses to generate — matches the number
          of substitution positions the user selected in the frontend.
        - energy_range: kcal/mol window for reporting poses (default 3 — wider
          than Vina's built-in default of 2 so small/rigid ligands consistently
          produce all requested poses instead of empty conformation columns).
        - Returns (True, None) on success or (False, error_str) on failure.
        """
        config_path = os.path.join(output_dir, f'{ligand_name}_config.txt')
        output_file = os.path.join(output_dir, f'{ligand_name}_docked.pdbqt')
        log_file    = os.path.join(output_dir, f'{ligand_name}_log.txt')

        try:
            if blind:
                self._write_blind_vina_config(
                    config_path, receptor, ligand, output_file,
                    num_modes=num_modes, energy_range=energy_range)
            else:
                self._write_vina_config(
                    config_path, receptor, ligand, output_file, pocket,
                    num_modes=num_modes, energy_range=energy_range)

            if first_ligand:
                print(f"\n  📄 Config for {ligand_name}:")
                with open(config_path) as cf:
                    for ln in cf:
                        print(f"     {ln}", end='')
                print()

            # vina --config only — no --log (unsupported in v1.2.5)
            cmd = ['vina', '--config', config_path]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

            # Save stdout+stderr as our log file (replaces --log)
            with open(log_file, 'w') as lf:
                lf.write(result.stdout)
                if result.stderr:
                    lf.write('\n' + result.stderr)

            ok = (result.returncode == 0
                  and os.path.exists(output_file)
                  and os.path.getsize(output_file) > 0)

            if not ok:
                err_msg = (result.stdout + result.stderr).strip()
                err_log = os.path.join(output_dir, f'{ligand_name}_error.txt')
                with open(err_log, 'w') as el:
                    el.write(f"returncode: {result.returncode}\n\n")
                    el.write(f"STDOUT:\n{result.stdout}\n")
                    el.write(f"STDERR:\n{result.stderr}\n\n")
                    el.write("CONFIG:\n")
                    with open(config_path) as cf:
                        el.write(cf.read())
                return False, err_msg[:600]

            return True, None

        except FileNotFoundError:
            print("\n  ❌ 'vina' not found in PATH. Install: conda install -c conda-forge vina")
            raise
        except subprocess.TimeoutExpired:
            return False, "Timeout after 600s"
        except Exception as e:
            return False, str(e)




    def _parse_conformations_from_pdbqt(self, pdbqt_path, max_confs):
        """
        Extract per-mode binding affinities directly from a docked PDBQT file.

        Every MODEL block written by Vina contains a header line of the form:
            REMARK VINA RESULT:   -4.557  0.000  0.000

        The first number is the binding affinity (kcal/mol) for that pose.
        Reading from the PDBQT is strictly more reliable than parsing stdout
        because:
          - Vina always writes this line for every mode it outputs
          - Stdout/log format has changed across Vina versions
          - The PDBQT file is the actual docking output, not a display artefact

        Returns a list of {'mode': N, 'affinity': float}, one entry per MODEL
        found (up to max_confs).  Never returns blanks — the list length equals
        however many modes Vina actually computed (which may be less than
        max_confs for rigid ligands, but every entry in the returned list is
        guaranteed to have a real numeric affinity).
        """
        conformations = []
        try:
            with open(pdbqt_path, 'r') as f:
                mode = 0
                for line in f:
                    if line.startswith('MODEL'):
                        mode += 1
                    elif line.startswith('REMARK VINA RESULT:'):
                        # Format: REMARK VINA RESULT:   -4.557  0.000  0.000
                        parts = line.split()
                        if len(parts) >= 4:
                            try:
                                affinity = float(parts[3])
                                conformations.append({
                                    'mode': mode if mode > 0 else len(conformations) + 1,
                                    'affinity': affinity
                                })
                                if len(conformations) >= max_confs:
                                    break
                            except (ValueError, IndexError):
                                continue
        except Exception:
            pass
        return conformations

    def _parse_vina_affinity(self, log_file):
        """Return mode-1 binding affinity from Vina log, or None."""
        try:
            with open(log_file, 'r') as f:
                for line in f:
                    parts = line.split()
                    if len(parts) >= 2 and parts[0] == '1':
                        try:
                            return float(parts[1])
                        except ValueError:
                            pass
        except Exception:
            pass
        return None


    def _parse_all_conformations(self, log_file, max_confs=3):
        """
        Parse top N conformations from Vina log file.
        
        Vina stdout format (what we save to log_file):
            mode |   affinity | dist from best mode
                 | (kcal/mol) | rmsd l.b.| rmsd u.b.
           -----+------------+----------+----------
              1       -5.912      0.000      0.000
              2       -5.719      1.234      2.567
              3       -5.631      1.789      3.012
        
        We look for lines starting with a digit followed by negative number.
        """
        conformations = []
        
        try:
            with open(log_file, 'r') as f:
                content = f.read()
            
            if not content.strip():
                return []
            
            lines = content.splitlines()
            
            # Find the results table by looking for the header
            # Accept variations: "mode", "affinity", case-insensitive
            start_idx = None
            for i, line in enumerate(lines):
                lower = line.lower()
                if 'mode' in lower and 'affinity' in lower:
                    # Found header — data starts after the separator line (----+----)
                    # Skip the header line and the next line (usually units or separator)
                    start_idx = i + 2
                    break
            
            if start_idx is None:
                # Fallback: look for lines that look like results (digit followed by negative number)
                # Example: "   1       -5.912"
                import re
                for line in lines:
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            mode = int(parts[0])
                            affinity = float(parts[1])
                            if affinity < 0:  # Valid affinity (always negative)
                                conformations.append({
                                    'mode': mode,
                                    'affinity': affinity
                                })
                                if len(conformations) >= max_confs:
                                    break
                        except (ValueError, IndexError):
                            continue
                return conformations
            
            # Parse from the identified start position
            for i in range(start_idx, min(start_idx + 30, len(lines))):
                if i >= len(lines):
                    break
                
                line = lines[i].strip()
                if not line:
                    continue
                
                parts = line.split()
                if len(parts) < 2:
                    continue
                
                try:
                    mode = int(parts[0])
                    affinity = float(parts[1])
                    
                    # Valid conformation: mode is 1-9, affinity is negative
                    if 1 <= mode <= 20 and affinity < 0:
                        conformations.append({
                            'mode': mode,
                            'affinity': affinity
                        })
                        if len(conformations) >= max_confs:
                            break
                except (ValueError, IndexError):
                    # Hit a non-data line — if we already have conformations, stop
                    if conformations:
                        break
                    continue
        
        except Exception as e:
            # Silently fail — caller will handle empty results
            pass
        
        return conformations
    
    def analyze_docking_results(self, docking_results, pdbqt_results):
        """
        Analyze docking results and generate ranked CSV report.
        
        Parameters:
        -----------
        docking_results : dict
            Docking results from run_autodock_vina_docking
        pdbqt_results : dict
            PDBQT conversion results
        
        Returns:
        --------
        str : Path to CSV file
        """
        print("\n" + "="*60)
        print("STEP 11: Docking Results Analysis")
        print("="*60)
        
        if not docking_results:
            print("\n⚠️ No docking results to analyze")
            return None
        
        # Collect all ligand data
        all_ligands = []
        
        for pocket_id, pocket_results in docking_results.items():
            for job_id, job_results in pocket_results.get('jobs', {}).items():
                output_dir = job_results.get('output_dir')
                
                if not output_dir or not os.path.exists(output_dir):
                    continue
                
                # Find all log files
                log_files = [f for f in os.listdir(output_dir) if f.endswith('_log.txt')]
                
                print(f"\n  Job {job_id}: Found {len(log_files)} log files in {output_dir}")
                
                for log_file in log_files:
                    ligand_name = log_file.replace('_log.txt', '')
                    log_path = os.path.join(output_dir, log_file)
                    pdbqt_path = os.path.join(output_dir, f'{ligand_name}_docked.pdbqt')
                    
                    # Check if PDBQT exists (confirms docking succeeded)
                    if not os.path.exists(pdbqt_path):
                        print(f"    ⚠️  {ligand_name}: docked PDBQT missing, skipping")
                        continue
                    
                    # Get SMILES and MW from original .smi file
                    smi_path = os.path.join(self.output_dir, f'job{job_id}_ligands', f'{ligand_name}.smi')
                    smiles, mw = self._get_smiles_and_mw(smi_path, pdbqt_path=pdbqt_path)

                    # ── Parse conformations ───────────────────────────────
                    # num_confs = number of substitution positions chosen
                    num_confs = max(1, len(self.substitution_jobs))

                    # Primary source: REMARK VINA RESULT lines in docked PDBQT.
                    # These are written by Vina for every pose it actually
                    # computed — they are the ground truth.
                    conformations = self._parse_conformations_from_pdbqt(
                        pdbqt_path, max_confs=num_confs)

                    # Fallback: stdout log (older Vina versions, blind docking)
                    if not conformations:
                        conformations = self._parse_all_conformations(
                            log_path, max_confs=num_confs)

                    if not conformations:
                        try:
                            with open(log_path) as f:
                                log_content = f.read()
                            if len(log_content) < 100:
                                print(f"    ⚠️  {ligand_name}: log too short ({len(log_content)} bytes)")
                            else:
                                print(f"    ⚠️  {ligand_name}: no conformations found in PDBQT or log")
                        except Exception as e:
                            print(f"    ⚠️  {ligand_name}: can't read log - {e}")
                        continue

                    # ── Strict padding: no blank conformation columns ─────
                    # If Vina produced fewer poses than requested (rigid/small
                    # ligands with limited conformational freedom), pad the
                    # missing slots with the weakest real affinity found.
                    # This is physically correct: it means "the weakest binding
                    # pose Vina could find for this ligand was X kcal/mol" and
                    # is far more informative than a blank cell.
                    actual_found = len(conformations)
                    if actual_found < num_confs:
                        worst_affinity = conformations[-1]['affinity']
                        reason = (
                            "rigid/few rotatable bonds — only 1 distinct "
                            "binding geometry exists"
                            if actual_found == 1
                            else f"only {actual_found} distinct poses within "
                                 f"energy_range found"
                        )
                        print(f"    ℹ️  {ligand_name}: {actual_found}/{num_confs} poses "
                              f"from Vina ({reason}). "
                              f"Padding missing pose(s) with worst found "
                              f"affinity ({worst_affinity:.3f} kcal/mol).")
                        for i in range(actual_found, num_confs):
                            conformations.append({
                                'mode':     i + 1,
                                'affinity': worst_affinity,
                                'padded':   True   # flag so CSV can mark it
                            })
                    
                    all_ligands.append({
                        'substrate_name': ligand_name,
                        'job_id': job_id,
                        'pocket_id': pocket_id,
                        'smiles': smiles,
                        'molecular_weight': mw,
                        'conformations': conformations,
                        'pdbqt_path': pdbqt_path,
                        'best_affinity': conformations[0]['affinity']
                    })
        
        if not all_ligands:
            print("\n⚠️ No valid docking results found")
            return None
        
        print(f"\n✓ Found {len(all_ligands)} docked ligands")
        
        # Sort ligands: first by best affinity (lower is better), then by MW (lower is better)
        all_ligands.sort(key=lambda x: (x['best_affinity'], x['molecular_weight']))
        
        # Generate CSV
        csv_path = self._generate_results_csv(all_ligands)
        
        # Auto-create ZIP files (both top 10 and all ligands) - no prompts
        zip_top10 = None
        zip_all = None
        
        if len(all_ligands) >= 10:
            zip_top10 = self._create_top_ligands_zip(all_ligands[:10], "top_10_ligands.zip")
        
        # Always create top_ligands.zip (all ligands) — used by the web UI download
        zip_all = self._create_top_ligands_zip(all_ligands, "top_ligands.zip")
        
        # Display summary
        self._display_analysis_summary(all_ligands, csv_path, zip_top10, zip_all)
        
        return csv_path
    
    def _get_smiles_and_mw(self, smi_path, pdbqt_path=None):
        """Get SMILES and molecular weight from .smi file, with PDBQT fallback."""
        # Primary source: generated .smi file from the same run
        try:
            with open(smi_path, 'r') as f:
                smiles = f.readline().strip()

            if smiles and not smiles.startswith('#'):
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    return smiles, Descriptors.MolWt(mol)
                return smiles, 0.0
        except Exception:
            pass

        # Fallback: recover ligand SMILES from docked PDBQT if .smi is missing
        if pdbqt_path and os.path.exists(pdbqt_path):
            try:
                ob_conv = ob.OBConversion()
                ob_conv.SetInFormat('pdbqt')
                ob_conv.SetOutFormat('smi')

                ob_mol = ob.OBMol()
                if ob_conv.ReadFile(ob_mol, pdbqt_path):
                    smi_line = ob_conv.WriteString(ob_mol).strip()
                    recovered = smi_line.split()[0] if smi_line else ''
                    if recovered:
                        mol = Chem.MolFromSmiles(recovered)
                        if mol:
                            return recovered, Descriptors.MolWt(mol)
                        return recovered, 0.0
            except Exception:
                pass

        return "Unknown", 0.0
    
    def _generate_results_csv(self, all_ligands):
        """
        Generate CSV with ranked docking results.

        Rules:
        - Only the top 3 ranked ligands are written (best binding affinity first).
        - The number of Conformation_N_Affinity columns equals the number of
          substitution positions the user selected (len(self.substitution_jobs)).
          If the user chose 1 position → 1 column.  3 positions → 3 columns.
        """
        csv_path = os.path.join(self.output_dir, 'docking_results_ranked.csv')

        # Number of conformation columns = number of substitution jobs
        num_confs  = max(1, len(self.substitution_jobs))
        # Only keep top 3 ligands
        top_ligands = all_ligands[:3]

        print(f"\n📊 Generating results CSV: {csv_path}")
        print(f"   Top ligands reported : {len(top_ligands)}")
        print(f"   Conformation columns : {num_confs}  "
              f"(= {num_confs} substitution position(s) selected)")

        with open(csv_path, 'w', newline='') as f:
            import csv
            writer = csv.writer(f)

            # ── Dynamic header ────────────────────────────────────────
            conf_headers = [f'Conformation_{i+1}_Affinity' for i in range(num_confs)]
            writer.writerow([
                'Rank',
                'Substrate_Name',
                'SMILES',
                'Molecular_Weight',
                *conf_headers,
                'PDBQT_File_Path',
                'Job_ID',
                'Pocket_ID',
            ])

            # ── Data rows (top 3 only) ────────────────────────────────
            for rank, ligand in enumerate(top_ligands, 1):
                confs = ligand['conformations']
                # Format each conformation value.
                # Padded slots (rigid ligand — fewer poses than requested) are
                # written as "VALUE*" so the reader knows this is the worst real
                # pose repeated, not a genuinely distinct binding conformation.
                conf_values = []
                for i in range(num_confs):
                    if i < len(confs):
                        val   = confs[i]['affinity']
                        label = f"{val}*" if confs[i].get('padded') else str(val)
                    else:
                        # Should never reach here after padding, but safe fallback
                        label = ''
                    conf_values.append(label)
                writer.writerow([
                    rank,
                    ligand['substrate_name'],
                    ligand['smiles'],
                    f"{ligand['molecular_weight']:.2f}",
                    *conf_values,
                    ligand['pdbqt_path'],
                    ligand['job_id'],
                    ligand['pocket_id'],
                ])

        print(f"✓ CSV written: {len(top_ligands)} ligand(s), "
              f"{num_confs} conformation column(s)")
        return csv_path
    
    def _create_top_ligands_zip(self, top_ligands, filename="top_ligands.zip"):
        """Create ZIP file containing top ligand PDBQT files."""
        import zipfile
        
        zip_path = os.path.join(self.output_dir, filename)
        
        print(f"\n📦 Creating ZIP file: {zip_path}")
        
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            # ── Include receptor (enzyme apo PDBQT) first ──────────────────
            enzyme_pdbqt = getattr(self, 'enzyme_apo_file', None)
            if enzyme_pdbqt:
                receptor_pdbqt = enzyme_pdbqt.replace('.pdb', '.pdbqt')
                if os.path.exists(receptor_pdbqt):
                    zipf.write(receptor_pdbqt, 'receptor/enzyme_apo.pdbqt')
                    print(f"  ✓ Added receptor: enzyme_apo.pdbqt")

            # ── Include ranked docked ligands ───────────────────────────────
            for i, ligand in enumerate(top_ligands, 1):
                if os.path.exists(ligand['pdbqt_path']):
                    arcname = f"ligands/rank_{i:03d}_{ligand['substrate_name']}_docked.pdbqt"
                    zipf.write(ligand['pdbqt_path'], arcname)

                    smi_path = os.path.join(
                        self.output_dir,
                        f"job{ligand['job_id']}_ligands",
                        f"{ligand['substrate_name']}.smi"
                    )
                    if os.path.exists(smi_path):
                        arcname_smi = f"ligands/rank_{i:03d}_{ligand['substrate_name']}.smi"
                        zipf.write(smi_path, arcname_smi)

        print(f"✓ ZIP created with receptor + {len(top_ligands)} ligand files")
        
        return zip_path
    
    def _display_analysis_summary(self, all_ligands, csv_path=None, zip_top10=None, zip_all=None):
        """Display summary of docking analysis."""
        print("\n" + "="*60)
        print("Docking Analysis Summary")
        print("="*60)
        
        print(f"\nTotal ligands analyzed: {len(all_ligands)}")
        
        # Output files
        if csv_path or zip_top10 or zip_all:
            print(f"\n📄 Output Files:")
        if csv_path:
            print(f"  CSV Report: {csv_path}")
        if zip_top10:
            print(f"  Top 10 ZIP: {zip_top10}")
        if zip_all:
            print(f"  All Ligands ZIP: {zip_all}")
        
        # Affinity statistics
        affinities = [lig['best_affinity'] for lig in all_ligands]
        
        print(f"\n📊 Binding Affinity Statistics:")
        print(f"  Best (most negative): {min(affinities):.2f} kcal/mol")
        print(f"  Worst (least negative): {max(affinities):.2f} kcal/mol")
        print(f"  Average: {sum(affinities)/len(affinities):.2f} kcal/mol")
        
        # Top 10 ligands
        print(f"\n🏆 Top 10 Best Binders:")
        print("-"*60)
        print(f"{'Rank':<6} {'Substrate':<20} {'Affinity':<12} {'MW':<10}")
        print("-"*60)
        
        for rank, ligand in enumerate(all_ligands[:10], 1):
            print(f"{rank:<6} {ligand['substrate_name']:<20} "
                  f"{ligand['best_affinity']:<12.2f} {ligand['molecular_weight']:<10.2f}")
        
        # Affinity distribution
        excellent = sum(1 for a in affinities if a <= -10)
        good = sum(1 for a in affinities if -10 < a <= -8)
        moderate = sum(1 for a in affinities if -8 < a <= -6)
        weak = sum(1 for a in affinities if a > -6)
        
        print(f"\n📈 Binding Affinity Distribution:")
        print(f"  Excellent (≤ -10 kcal/mol): {excellent} ligands")
        print(f"  Good (-10 to -8 kcal/mol): {good} ligands")
        print(f"  Moderate (-8 to -6 kcal/mol): {moderate} ligands")
        print(f"  Weak (> -6 kcal/mol): {weak} ligands")
        
        print("\n" + "="*60)
    
    def _display_docking_summary(self, docking_results):
        """Display docking results summary."""
        print("\n" + "="*60)
        print("Docking Summary")
        print("="*60)
        
        for pocket_id, results in docking_results.items():
            pocket_info = results['pocket_info']
            print(f"\nPocket {pocket_id}: {pocket_info.get('name', 'Unknown')}")
            print(f"  Center: ({pocket_info['center_x']:.2f}, {pocket_info['center_y']:.2f}, {pocket_info['center_z']:.2f})")
            
            for job_id, job_results in results.get('jobs', {}).items():
                print(f"\n  Job {job_id}:")
                print(f"    Docked: {job_results['docked']}/{job_results['total']} ligands")
                print(f"    Results: {job_results['output_dir']}")
            
            print(f"\n  Total: {results['total_docked']} ligands docked")
            if results['total_failed'] > 0:
                print(f"  Failed: {results['total_failed']}")
    
    def _generate_docking_scripts(self, pockets, pdbqt_results):
        """Generate bash scripts for manual docking execution."""
        print("\n📝 Generating docking scripts...")
        
        script_dir = os.path.join(self.output_dir, 'docking_scripts')
        os.makedirs(script_dir, exist_ok=True)
        
        for pocket in pockets:
            script_file = os.path.join(script_dir, f'dock_pocket_{pocket["id"]}.sh')
            
            with open(script_file, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write(f"# Docking script for Pocket {pocket['id']}\n")
                f.write(f"# {pocket.get('name', 'Unknown')}\n\n")
                
                f.write(f"RECEPTOR=\"{pdbqt_results['enzyme_pdbqt']}\"\n")
                f.write(f"CENTER_X={pocket['center_x']}\n")
                f.write(f"CENTER_Y={pocket['center_y']}\n")
                f.write(f"CENTER_Z={pocket['center_z']}\n")
                f.write(f"SIZE_X={pocket.get('size_x', 20)}\n")
                f.write(f"SIZE_Y={pocket.get('size_y', 20)}\n")
                f.write(f"SIZE_Z={pocket.get('size_z', 20)}\n\n")
                
                # Add docking commands for all jobs
                for job_id in pdbqt_results.get('jobs', {}).keys():
                    pdbqt_dir = os.path.join(self.output_dir, f'job{job_id}_ligands', 'pdbqt_ligands')
                    output_dir = os.path.join(self.output_dir, f'docking_pocket_{pocket["id"]}', f'job{job_id}')
                    
                    f.write(f"# Job {job_id}\n")
                    f.write(f"mkdir -p {output_dir}\n\n")
                    f.write(f"for LIGAND in {pdbqt_dir}/*.pdbqt; do\n")
                    f.write(f"    NAME=$(basename \"$LIGAND\" .pdbqt)\n")
                    f.write(f"    vina \\\n")
                    f.write(f"        --receptor $RECEPTOR \\\n")
                    f.write(f"        --ligand \"$LIGAND\" \\\n")
                    f.write(f"        --center_x $CENTER_X \\\n")
                    f.write(f"        --center_y $CENTER_Y \\\n")
                    f.write(f"        --center_z $CENTER_Z \\\n")
                    f.write(f"        --size_x $SIZE_X \\\n")
                    f.write(f"        --size_y $SIZE_Y \\\n")
                    f.write(f"        --size_z $SIZE_Z \\\n")
                    f.write(f"        --exhaustiveness 8 \\\n")
                    f.write(f"        --out {output_dir}/${{NAME}}_docked.pdbqt \\\n")
                    f.write(f"        --log {output_dir}/${{NAME}}_log.txt\n")
                    f.write(f"done\n\n")
            
            # Make script executable
            os.chmod(script_file, 0o755)
            print(f"  ✓ Created: {script_file}")
        
        print(f"\n💡 To run docking manually:")
        print(f"   cd {script_dir}")
        print(f"   ./dock_pocket_1.sh")
        
        return {'scripts_generated': True, 'script_dir': script_dir}
    
    def convert_all_to_pdbqt(self):
        """
        Convert enzyme and all ligands to PDBQT format.
        
        Returns:
        --------
        dict : Conversion results
        """
        print("\n" + "="*60)
        print("PDBQT CONVERSION FOR DOCKING")
        print("="*60)
        
        results = {
            'enzyme_pdbqt': None,
            'jobs': {}
        }
        
        # Convert enzyme
        enzyme_pdbqt = self.convert_enzyme_to_pdbqt()
        results['enzyme_pdbqt'] = enzyme_pdbqt
        
        # Convert ligands for each job
        if self.substitution_jobs:
            print("\n" + "="*60)
            print("Converting Ligands to PDBQT")
            print("="*60)
            
            for job in self.substitution_jobs:
                job_results = self.convert_ligands_to_pdbqt(job.job_id)
                results['jobs'][job.job_id] = job_results
        
        # Summary
        print("\n" + "="*60)
        print("PDBQT Conversion Summary")
        print("="*60)
        
        if enzyme_pdbqt:
            print(f"\n✓ Enzyme: {enzyme_pdbqt}")
        else:
            print(f"\n❌ Enzyme conversion failed")
        
        if results['jobs']:
            print(f"\nLigands:")
            total_converted = sum(r.get('converted', 0) for r in results['jobs'].values())
            total_failed = sum(r.get('failed', 0) for r in results['jobs'].values())
            total_skipped = sum(r.get('skipped', 0) for r in results['jobs'].values())
            total_ligands = sum(r.get('total', 0) for r in results['jobs'].values())
            
            for job_id, stats in results['jobs'].items():
                converted = stats.get('converted', 0)
                total = stats.get('total', 0)
                status_parts = [f"{converted}/{total} ligands"]
                if stats.get('skipped', 0) > 0:
                    status_parts.append(f"{stats['skipped']} skipped")
                print(f"  Job {job_id}: {', '.join(status_parts)}")
                
                failed_names = stats.get('failed_names') or []
                skipped_names = stats.get('skipped_names') or []
                if failed_names:
                    print(f"    Failed: {', '.join(failed_names)}")
                if skipped_names:
                    print(f"    Skipped (too large): {', '.join(skipped_names)}")
            
            print(f"\n  Total: {total_converted}/{total_ligands} ligands converted")
            if total_skipped > 0:
                print(f"  Skipped: {total_skipped} (molecules > 300 atoms)")
            if total_failed > 0:
                print(f"  Failed: {total_failed}")
                print(f"\n  💡 For failed/skipped ligands:")
                print(f"     - Check *_FAILED.txt files for diagnostic information")
                print(f"     - Ultra-large molecules (>300 atoms) are skipped automatically")
                print(f"     - Use external 3D generation tools (Schrödinger, OpenEye)")
                print(f"     - Or proceed with {total_converted} ligands for docking")
        
        return results
    
    def dock_ligands_with_vina(self, pockets, enzyme_pdbqt=None):
        """
        Dock all ligands into detected pockets using AutoDock Vina.
        
        Parameters:
        -----------
        pockets : list
            List of pocket dictionaries from detect_binding_pockets()
        enzyme_pdbqt : str, optional
            Path to enzyme PDBQT file
        
        Returns:
        --------
        dict : Docking results organized by pocket and ligand
        """
        print("\n" + "="*60)
        print("STEP 10: Molecular Docking (AutoDock Vina)")
        print("="*60)
        
        if not pockets:
            print("\n⚠️ No pockets provided, skipping docking")
            return {}
        
        if not enzyme_pdbqt:
            enzyme_pdbqt = self.enzyme_apo_file.replace('.pdb', '.pdbqt')
        
        if not os.path.exists(enzyme_pdbqt):
            print(f"❌ Enzyme PDBQT not found: {enzyme_pdbqt}")
            return {}
        
        # Check if vina is available
        try:
            subprocess.run(['vina', '--help'], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("❌ AutoDock Vina not found. Install with:")
            print("   conda install -c conda-forge autodock-vina")
            return {}
        
        # Collect all ligand PDBQT files
        ligand_files = []
        for job_dir in ['job1_ligands', 'job2_ligands', 'job3_ligands']:
            pdbqt_dir = os.path.join(self.output_dir, job_dir, 'pdbqt_ligands')
            if os.path.exists(pdbqt_dir):
                ligands = [os.path.join(pdbqt_dir, f) for f in os.listdir(pdbqt_dir) 
                          if f.endswith('.pdbqt')]
                ligand_files.extend(ligands)
        
        if not ligand_files:
            print("❌ No ligand PDBQT files found")
            return {}
        
        print(f"\nDocking Setup:")
        print(f"  Enzyme: {enzyme_pdbqt}")
        print(f"  Ligands: {len(ligand_files)} files")
        print(f"  Pockets: {len(pockets)} binding sites")
        print(f"  Total docking runs: {len(ligand_files) * len(pockets)}")
        
        # Create docking directory
        docking_dir = os.path.join(self.output_dir, 'docking_results')
        os.makedirs(docking_dir, exist_ok=True)
        
        # Run docking for each pocket
        all_results = {}
        
        for pocket_idx, pocket in enumerate(pockets, 1):
            print(f"\n{'-'*60}")
            print(f"Docking to Pocket {pocket_idx}")
            print(f"{'-'*60}")
            print(f"  Center: ({pocket['center_x']:.2f}, {pocket['center_y']:.2f}, {pocket['center_z']:.2f})")
            print(f"  Box size: {pocket['size_x']}×{pocket['size_y']}×{pocket['size_z']} Ų")
            
            pocket_dir = os.path.join(docking_dir, f'pocket_{pocket_idx}')
            os.makedirs(pocket_dir, exist_ok=True)
            
            pocket_results = self._dock_to_pocket(
                enzyme_pdbqt, 
                ligand_files, 
                pocket, 
                pocket_dir,
                pocket_idx
            )
            
            all_results[f'pocket_{pocket_idx}'] = pocket_results
        
        # Generate summary
        self._generate_docking_summary(all_results, docking_dir)
        
        return all_results
    
    def _dock_to_pocket(self, enzyme_pdbqt, ligand_files, pocket, output_dir, pocket_idx):
        """Dock all ligands to a single pocket."""
        results = []
        
        total = len(ligand_files)
        for i, ligand_path in enumerate(ligand_files, 1):
            ligand_name = os.path.basename(ligand_path).replace('.pdbqt', '')
            
            # Output files
            out_pdbqt = os.path.join(output_dir, f'{ligand_name}_docked.pdbqt')
            log_file = os.path.join(output_dir, f'{ligand_name}_log.txt')
            
            print(f"  [{i}/{total}] Docking {ligand_name}...", end=' ', flush=True)
            
            try:
                # Run AutoDock Vina
                cmd = [
                    'vina',
                    '--receptor', enzyme_pdbqt,
                    '--ligand', ligand_path,
                    '--center_x', str(pocket['center_x']),
                    '--center_y', str(pocket['center_y']),
                    '--center_z', str(pocket['center_z']),
                    '--size_x', str(pocket['size_x']),
                    '--size_y', str(pocket['size_y']),
                    '--size_z', str(pocket['size_z']),
                    '--exhaustiveness', '8',
                    '--out', out_pdbqt,
                    '--log', log_file
                ]
                
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=300  # 5 minute timeout per ligand
                )
                
                # Parse binding affinity from log
                affinity = self._parse_vina_log(log_file)
                
                if affinity is not None:
                    print(f"✓ {affinity:.1f} kcal/mol")
                    results.append({
                        'ligand': ligand_name,
                        'affinity': affinity,
                        'output': out_pdbqt,
                        'log': log_file,
                        'success': True
                    })
                else:
                    print(f"⚠️ Failed to parse results")
                    results.append({
                        'ligand': ligand_name,
                        'affinity': None,
                        'success': False
                    })
                    
            except subprocess.TimeoutExpired:
                print(f"⏱️ Timeout")
                results.append({
                    'ligand': ligand_name,
                    'affinity': None,
                    'success': False,
                    'error': 'timeout'
                })
            except Exception as e:
                print(f"❌ Error: {e}")
                results.append({
                    'ligand': ligand_name,
                    'affinity': None,
                    'success': False,
                    'error': str(e)
                })
        
        return results
    
    def _parse_vina_log(self, log_file):
        """
        Parse AutoDock Vina log file to extract top 3 binding affinities.
        
        Returns:
        --------
        list : List of affinities for top 3 conformations, or empty list if failed
        """
        try:
            affinities = []
            with open(log_file, 'r') as f:
                in_results = False
                for line in f:
                    # Look for results table
                    if 'mode |   affinity | dist from best mode' in line:
                        in_results = True
                        continue
                    
                    if in_results:
                        # Parse result lines: "   1       -7.5      0.000      0.000"
                        line = line.strip()
                        if not line or line.startswith('-'):
                            break
                        
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                mode = int(parts[0])
                                affinity = float(parts[1])
                                affinities.append(affinity)
                                
                                # Get top 3 conformations
                                if len(affinities) >= 3:
                                    break
                            except ValueError:
                                continue
            
            return affinities if affinities else []
        except:
            return []
    
    def _generate_docking_summary(self, all_results, docking_dir):
        """Generate summary of docking results."""
        print(f"\n{'='*60}")
        print("Docking Summary")
        print(f"{'='*60}")
        
        summary_file = os.path.join(docking_dir, 'docking_summary.txt')
        csv_file = os.path.join(docking_dir, 'docking_summary.csv')
        
        # Prepare data
        all_data = []
        for pocket_name, results in all_results.items():
            for result in results:
                if result['success'] and result['affinity'] is not None:
                    all_data.append({
                        'pocket': pocket_name,
                        'ligand': result['ligand'],
                        'affinity': result['affinity'],
                        'output': result.get('output', '')
                    })
        
        if not all_data:
            print("\n⚠️ No successful docking results")
            return
        
        # Sort by affinity (best first)
        all_data.sort(key=lambda x: x['affinity'])
        
        # Write summary file
        with open(summary_file, 'w') as f:
            f.write("AutoDock Vina Docking Results Summary\n")
            f.write("="*60 + "\n\n")
            
            f.write("Top 10 Best Binding Affinities:\n")
            f.write("-"*60 + "\n")
            for i, data in enumerate(all_data[:10], 1):
                f.write(f"{i:2d}. {data['ligand']:30s} ({data['pocket']}): {data['affinity']:6.1f} kcal/mol\n")
            
            f.write("\n" + "="*60 + "\n\n")
            f.write("Results by Pocket:\n")
            f.write("-"*60 + "\n")
            
            for pocket_name in sorted(all_results.keys()):
                pocket_data = [d for d in all_data if d['pocket'] == pocket_name]
                if pocket_data:
                    f.write(f"\n{pocket_name.upper()}:\n")
                    for data in sorted(pocket_data, key=lambda x: x['affinity'])[:5]:
                        f.write(f"  {data['ligand']:30s}: {data['affinity']:6.1f} kcal/mol\n")
        
        # Write CSV file
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['pocket', 'ligand', 'affinity', 'output'])
            writer.writeheader()
            writer.writerows(all_data)
        
        # Print summary to console
        print(f"\n✓ Docked: {len(all_data)} ligand-pocket combinations")
        print(f"\nTop 5 Best Binding Affinities:")
        for i, data in enumerate(all_data[:5], 1):
            print(f"  {i}. {data['ligand']:25s} → {data['pocket']:10s}: {data['affinity']:6.1f} kcal/mol")
        
        print(f"\n📄 Full results:")
        print(f"   Summary: {summary_file}")
        print(f"   CSV: {csv_file}")
        print(f"   Docked structures: {docking_dir}/pocket_*/")
    
    def prepare_apo_protein(self):
        """Clean enzyme PDB using BioPython."""
        print("\n" + "="*60)
        print("STEP 6: Preparing Apo Protein")
        print("="*60)
        
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('enzyme', self.enzyme_file)
        
        print(f"\nLoading enzyme: {self.enzyme_file}")
        
        all_chains = []
        for model in structure:
            for chain in model:
                all_chains.append(chain.id)
        self.chains_info = list(set(all_chains))
        print(f"Available chains: {self.chains_info}")
        
        selector = ChainAProteinSelect()
        self.enzyme_apo_file = os.path.join(self.output_dir, 'enzyme_apo.pdb')
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(self.enzyme_apo_file, selector)
        
        self.removed_hetero = selector.removed_hetero
        self.removed_chains = list(set(selector.removed_chains))
        
        print(f"\n✓ Saved clean structure: {self.enzyme_apo_file}")
        print(f"  Chains removed: {self.removed_chains if self.removed_chains else 'None'}")
        print(f"  Heteroatoms removed: {len(self.removed_hetero)}")
        
        return self.enzyme_apo_file
    
    def get_results(self):
        """Get all processing results."""
        return {
            'original_smiles': self.original_smiles,
            'scaffold_smiles': self.scaffold_smiles,
            'scaffold_mol': self.scaffold_mol,
            'is_scaffold_original': self.is_scaffold_original,
            'brics_fragments': self.brics_fragments,
            'active_site_atoms': self.active_site_atoms,
            'substitution_jobs': self.substitution_jobs,
            'fragment_library': self.fragment_library,
            'enzyme_original': self.enzyme_file,
            'enzyme_apo': self.enzyme_apo_file,
            'molecule_file': self.molecule_file,
            'chains_original': self.chains_info,
            'chains_removed': self.removed_chains,
            'removed_heteroatoms': self.removed_hetero,
            'output_directory': self.output_dir
        }
    
    def process_all(self, use_original_as_scaffold=None, fragment_file=None, progress_cb=None):
        """
        Run complete processing pipeline.
        progress_cb(pct: int, stage: str) is called at each milestone so the
        web UI can show real percentage progress instead of a fake timer.
        """
        def _rpt(pct, stage=''):
            msg = f"[PROGRESS] {pct}% — {stage}"
            try:
                sys.__stdout__.write(msg + '\n')
                sys.__stdout__.flush()
            except Exception:
                pass
            if callable(progress_cb):
                progress_cb(pct, stage)


        # Step 1-5: Scaffold and job creation
        _rpt(2, 'Converting molecule to SMILES…')
        self.convert_molecule_to_smiles()
        _rpt(5, 'Selecting scaffold…')
        self.process_scaffold_selection(use_original=use_original_as_scaffold)
        self.display_scaffold_structure()
        _rpt(8, 'Marking active sites…')
        self.mark_active_sites()
        _rpt(10, 'Creating substitution jobs…')
        self.create_substitution_jobs()

        # Step 6: Prepare apo protein
        _rpt(14, 'Preparing apo protein…')
        self.prepare_apo_protein()

        # Step 6.5: Early pocket detection for volume-based fragment filtering
        pockets = []
        max_fragment_volume = None
        if self.enzyme_apo_file and os.path.exists(self.enzyme_apo_file):
            _rpt(20, 'Detecting binding pockets…')
            pockets = self.detect_binding_pockets()

            if pockets:
                best_pocket = pockets[0]
                pocket_volume = self.estimate_pocket_volume(best_pocket)
                scaffold_volume = self.calculate_mol_volume(self.scaffold_smiles) or 0.0
                base_remaining_volume = max(0.0, pocket_volume - scaffold_volume)
                relaxed_remaining_volume = max(0.0, (pocket_volume * 1.25) - scaffold_volume)

                print(f"\n📐 Volume analysis:")
                print(f"   Pocket volume  (best pocket): {pocket_volume:.1f} Å³")
                print(f"   Scaffold volume (VDW):         {scaffold_volume:.1f} Å³")
                print(f"   Raw remaining space (Y-X):      {base_remaining_volume:.1f} Å³")
                print(f"   Relaxed limit ((Y*1.25)-X):     {relaxed_remaining_volume:.1f} Å³")

                if relaxed_remaining_volume > 0:
                    max_fragment_volume = relaxed_remaining_volume
                    print(f"   ✓ Fragments > {relaxed_remaining_volume:.1f} Å³ will be excluded")
                else:
                    print(f"   ⚠️  Scaffold already fills/exceeds relaxed limit — volume filter disabled")
            else:
                print("\n  No pockets detected — volume filter disabled; using MW filter only")

        # Step 7: Fragment substitution
        fragment_stats = {}
        if self.substitution_jobs:
            fragment_file = getattr(self, 'fragment_file', None)

            if fragment_file is None and self.interactive:
                fragment_file = input("\nEnter path to fragments .txt file: ").strip()

            if fragment_file:
                _rpt(28, 'Running fragment substitution…')
                fragment_stats = self.substitute_fragments(
                    fragment_file,
                    max_fragment_volume=max_fragment_volume
                )

        # Step 8: Convert to PDBQT
        pdbqt_results = {}
        if self.enzyme_apo_file or fragment_stats:
            _rpt(45, 'Converting to PDBQT…')
            pdbqt_results = self.convert_all_to_pdbqt()

        # Step 9: Molecular docking
        docking_results = {}
        if pdbqt_results.get('enzyme_pdbqt'):
            _rpt(58, 'Running molecular docking…')

            if not pockets:
                pockets = self.detect_binding_pockets()

            docking_results = self.run_autodock_vina_docking(pockets, pdbqt_results)

        # Step 10: CSV + ZIP generation
        results_csv = None
        if docking_results:
            _rpt(92, 'Generating results CSV & ZIP…')
            results_csv = self.analyze_docking_results(docking_results, pdbqt_results)

        _rpt(100, 'Complete')

        results = self.get_results()
        results['fragment_stats']  = fragment_stats
        results['pdbqt_results']   = pdbqt_results
        results['pockets']         = pockets
        results['docking_results'] = docking_results
        results['results_csv']     = results_csv

        return results


def main():
    """Main function."""
    print("="*60)
    print("Enzyme-Ligand Processing Pipeline")
    print("With Active Site Marking and Fragment Substitution")
    print("="*60)
    
    # Set output directory to script directory automatically
    script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
    output_dir = os.path.join(script_dir, 'enzyme_processor_results')
    
    if len(sys.argv) >= 3:
        enzyme_file = sys.argv[1]
        molecule_file = sys.argv[2]
        # Command-line arg 3 can override output_dir if provided
        if len(sys.argv) > 3:
            output_dir = sys.argv[3]
    else:
        print("\nEnter input files:")
        enzyme_file = input("Enzyme PDB file path: ").strip()
        molecule_file = input("Molecule file path (.mol/.pdb/.sdf): ").strip()
        # Output directory is automatic - no prompt
        print(f"\nOutput directory: {output_dir}")
        print("(automatically set to script directory)")
    
    try:
        processor = EnzymeLigandProcessor(enzyme_file, molecule_file, output_dir, interactive=True)
        results = processor.process_all()
        
        print("\n" + "="*60)
        print("✅ PROCESSING COMPLETE!")
        print("="*60)
        print(f"\n📋 Results:")
        print(f"   Scaffold SMILES: {results['scaffold_smiles']}")
        print(f"   Active site atoms: {results['active_site_atoms']}")
        print(f"   Substitution jobs: {len(results['substitution_jobs'])}")
        
        for job in results['substitution_jobs']:
            mw_str = f"{job.mw_range[0]}-{job.mw_range[1]}Da" if job.mw_range else "Not set"
            print(f"     Job {job.job_id}: {len(job.attachment_atoms)} attachment(s), MW range: {mw_str}")
        
        print(f"   Apo Protein: {results['enzyme_apo']}")
        
        # Show fragment substitution results if available
        if 'fragment_stats' in results and results['fragment_stats']:
            print(f"\n💡 Fragment Substitution Results:")
            for job_id, stats in results['fragment_stats'].items():
                if stats['generated_ligands'] > 0:
                    output_dir_path = os.path.join(output_dir, f'job{job_id}_ligands')
                    print(f"   Job {job_id}: {stats['generated_ligands']} ligands → {output_dir_path}/")
                    print(f"            (ligand_1.smi, ligand_2.smi, ... ligand_{stats['generated_ligands']}.smi)")
        
        # Show PDBQT conversion results if available
        if 'pdbqt_results' in results and results['pdbqt_results']:
            pdbqt = results['pdbqt_results']
            print(f"\n🎯 PDBQT Files for Docking:")
            
            if pdbqt.get('enzyme_pdbqt'):
                print(f"   Enzyme: {pdbqt['enzyme_pdbqt']}")
            
            if pdbqt.get('jobs'):
                print(f"   Ligands:")
                for job_id, stats in pdbqt['jobs'].items():
                    if stats['converted'] > 0:
                        pdbqt_dir = os.path.join(output_dir, f'job{job_id}_ligands', 'pdbqt_ligands')
                        print(f"     Job {job_id}: {stats['converted']} ligands → {pdbqt_dir}/")
                        print(f"              (ligand_1.pdbqt, ligand_2.pdbqt, ... ligand_{stats['converted']}.pdbqt)")
        
        # Show pocket detection results if available
        if 'pockets' in results and results['pockets']:
            print(f"\n🔍 Detected Binding Pockets:")
            for pocket in results['pockets']:
                print(f"   Pocket {pocket['id']}: {pocket.get('name', 'Unknown')}")
                print(f"      Center: ({pocket['center_x']:.2f}, {pocket['center_y']:.2f}, {pocket['center_z']:.2f})")
                if pocket.get('druggability_score') is not None:
                    print(f"      Druggability: {pocket['druggability_score']:.2f}")
        
        # Show docking results if available
        if 'docking_results' in results and results['docking_results']:
            print(f"\n⚗️  Docking Results:")
            
            for pocket_id, pocket_results in results['docking_results'].items():
                pocket_info = pocket_results.get('pocket_info', {})
                print(f"\n   Pocket {pocket_id}: {pocket_info.get('name', 'Unknown')}")
                
                for job_id, job_results in pocket_results.get('jobs', {}).items():
                    if job_results['docked'] > 0:
                        print(f"     Job {job_id}: {job_results['docked']} ligands docked → {job_results['output_dir']}/")
                
                print(f"   Total: {pocket_results['total_docked']} ligands docked")
        
        # Show analysis results if available
        if 'results_csv' in results and results['results_csv']:
            print(f"\n📊 Docking Analysis Results:")
            print(f"   Ranked CSV: {results['results_csv']}")
            print(f"   Contains: Top 3 conformations per ligand")
            print(f"   Sorted by: Binding affinity, then molecular weight")
            
            # Check if ZIP was created
            zip_path = os.path.join(output_dir, 'top_ligands.zip')
            if os.path.exists(zip_path):
                print(f"   Top Ligands ZIP: {zip_path}")
        
        return results
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)




app = Flask(__name__, template_folder='.')
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, 'enzyme_processor_results')
UPLOADS_DIR = os.path.join(OUTPUT_DIR, 'uploads')
RUNTIME_LOG = os.path.join(OUTPUT_DIR, 'backend_runtime.log')

class JobLogStream(io.TextIOBase):
    """Capture print() output from a pipeline thread and forward it to _job_log()."""

    def __init__(self, job_id):
        super().__init__()
        self.job_id = job_id
        self._buffer = ""

    def write(self, s):
        if s is None:
            return 0
        s = str(s)
        if not s:
            return 0

        s = s.replace('\r\n', '\n').replace('\r', '\n')
        self._buffer += s

        while '\n' in self._buffer:
            line, self._buffer = self._buffer.split('\n', 1)
            line = line.strip()
            if line:
                _job_log(self.job_id, line)

        return len(s)

    def flush(self):
        if self._buffer.strip():
            _job_log(self.job_id, self._buffer.strip())
            self._buffer = ""

# ── Thread-safe job store for pipeline progress tracking ────────────────────
_jobs = {}          # job_id -> {pct, stage, done, error, results}
_jobs_lock = threading.Lock()

def _job_set(job_id, **kwargs):
    with _jobs_lock:
        _jobs.setdefault(job_id, {}).update(kwargs)

def _job_get(job_id):
    with _jobs_lock:
        return dict(_jobs.get(job_id, {}))

def _append_runtime_log(line):
    try:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        with open(RUNTIME_LOG, 'a', encoding='utf-8') as fh:
            fh.write(line + '\n')
    except Exception:
        pass

def _job_log(job_id, message):
    ts = time.strftime('%H:%M:%S')
    line = f"[{ts}] {message}"
    full_line = f"[{ts}] [JOB {job_id}] {message}"

    # Write directly to the real console/stdout.
    try:
        sys.__stdout__.write(full_line + '\n')
        sys.__stdout__.flush()
    except Exception:
        pass

    # Store per-job log lines for /pipeline-log
    with _jobs_lock:
        job = _jobs.setdefault(job_id, {})
        job.setdefault('logs', []).append(line)

    _append_runtime_log(full_line)


@app.after_request
def disable_http_cache(resp):
    # Prevent stale HTML/JS so frontend always reflects latest backend logic.
    resp.headers['Cache-Control'] = 'no-store, no-cache, must-revalidate, max-age=0'
    resp.headers['Pragma'] = 'no-cache'
    resp.headers['Expires'] = '0'
    return resp


@app.route('/')
def index():
    template_name = (
        'enzyme_processor_patched.html'
        if os.path.exists(os.path.join(SCRIPT_DIR, 'enzyme_processor_patched.html'))
        else 'enzyme_processor.html'
    )
    return render_template(template_name)


@app.route('/generate-scaffold', methods=['POST'])
def generate_scaffold():
    """
    Accepts a ligand file + scaffold choice.
    Returns: { success, sdf_block, scaffold_smiles }

    Pipeline (mirrors the user's reference code pattern):
      1. Parse the uploaded file with RDKit (SDMolSupplier / MolFromPDBFile / OBabel fallback)
      2. Get a canonical, artefact-free SMILES:  Chem.MolToSmiles(mol, isomericSmiles=True)
      3. Re-parse from that SMILES:  Chem.MolFromSmiles(smiles)
         This guarantees a clean mol object with no file-format dummies, no stale
         conformers, and no leftover BRICS attachment atoms — regardless of the input file.
      4. If BRICS mode: decompose the clean mol, pick largest fragment, strip
         ALL [n*] dummies with re.sub before re-parsing as a fresh clean mol.
      5. AddHs → ETKDG(seed=42) → UFFOptimizeMolecule → RemoveHs
         (UFF is more tolerant of unusual atom types than MMFF)
      6. Serialize to SDF and return.
    """
    def _embed_and_optimise(mol):
        """Add Hs, embed with ETKDG(seed=42), UFF-optimise, remove Hs.
        Falls back through progressively more permissive strategies."""
        mol = Chem.AddHs(mol)

        # Strategy 1: ETKDGv3 with seed
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        ok = AllChem.EmbedMolecule(mol, params)

        # Strategy 2: basic ETKDG with seed (user's reference code approach)
        if ok != 0:
            params2 = AllChem.ETKDG()
            params2.randomSeed = 42
            ok = AllChem.EmbedMolecule(mol, params2)

        # Strategy 3: random coords as last resort
        if ok != 0:
            params3 = AllChem.ETKDGv3()
            params3.randomSeed = 42
            params3.useRandomCoords = True
            ok = AllChem.EmbedMolecule(mol, params3)

        if ok == 0:
            # UFF is more robust than MMFF for unusual substituents / small fragments
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            except Exception:
                pass

        return Chem.RemoveHs(mol)

    def _clean_smiles_to_mol(smiles_str):
        """Parse a SMILES (possibly with BRICS [n*] dummies) and return a
        sanitised mol with all dummy atoms removed.

        Capping strategy — WHY we use C (methyl) not [H]:
          BRICS marks each cut bond with a numbered dummy [n*].  The dummy sits
          at the position where the adjacent fragment connected.  Replacing it
          with [H] removes that bonding position entirely, making the scaffold
          one carbon shorter than the original chain — this is the "missing
          carbon" the user observes.

          Replacing [n*] with C (a plain methyl cap) restores the carbon count:
          the attachment site becomes a CH₃ group, which is chemically neutral,
          does not affect the scaffold's interaction with the active site, and
          preserves the chain length seen in the original ligand.

          [*]   → C  (plain attachment dummy → methyl cap)
          [1*]  → C  (BRICS cut on bond type 1)
          [4*]  → C  etc.
        """
        # Replace every isotope-numbered BRICS dummy with plain carbon
        clean = re.sub(r'\[\d+\*\]', 'C', smiles_str)
        # Replace bare [*] or bare * with carbon too
        clean = re.sub(r'\[\*\]', 'C', clean)
        clean = re.sub(r'(?<!\[)\*(?!\])', 'C', clean)
        mol = Chem.MolFromSmiles(clean)
        if mol is not None:
            return mol
        # Fallback: try with [H] cap in case the methyl creates valence conflict
        clean2 = re.sub(r'\[\d+\*\]', '[H]', smiles_str)
        clean2 = re.sub(r'\[\*\]', '[H]', clean2)
        clean2 = re.sub(r'(?<!\[)\*(?!\])', '[H]', clean2)
        return Chem.MolFromSmiles(clean2)  # None if still invalid

    try:
        ligand_file = request.files.get('ligandFile')
        if not ligand_file:
            return jsonify({'success': False, 'error': 'No ligand file provided'}), 400

        scaffold_choice = request.form.get('scaffoldChoice', '1')

        suffix = Path(secure_filename(ligand_file.filename)).suffix.lower()
        with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
            ligand_file.save(tmp.name)
            tmp_path = tmp.name

        # ── Step 1: Parse uploaded file ─────────────────────────────────────
        rdmol = None
        if suffix in ('.sdf', '.mol'):
            # removeHs=True gives us a clean heavy-atom mol straight away
            supplier = Chem.SDMolSupplier(tmp_path, removeHs=True, sanitize=True)
            rdmol = next((m for m in supplier if m is not None), None)
        elif suffix == '.pdb':
            rdmol = Chem.MolFromPDBFile(tmp_path, removeHs=True, sanitize=True)

        if rdmol is None:
            # OpenBabel fallback for unsupported formats
            obConv = ob.OBConversion()
            obConv.SetInFormat(suffix.replace('.', ''))
            obConv.SetOutFormat('smi')
            obmol = ob.OBMol()
            obConv.ReadFile(obmol, tmp_path)
            raw_smi = obConv.WriteString(obmol).strip().split()[0]
            rdmol = Chem.MolFromSmiles(raw_smi)
            if rdmol is None:
                os.unlink(tmp_path)
                return jsonify({'success': False, 'error': 'Could not parse ligand file'}), 400

        os.unlink(tmp_path)

        # ── Step 2: Canonical clean SMILES from parsed mol ───────────────────
        # Chem.MolToSmiles strips file-format artefacts and gives the true
        # chemical graph — exactly what the reference code does.
        canonical_smiles = Chem.MolToSmiles(rdmol, isomericSmiles=True)

        # ── Step 3: Re-parse from SMILES → clean mol (no conformers, no dummies)
        # This mirrors: mol_clean = Chem.MolFromSmiles(smiles) in the reference code.
        base_mol = Chem.MolFromSmiles(canonical_smiles)
        if base_mol is None:
            return jsonify({'success': False, 'error': 'Canonical SMILES round-trip failed'}), 400

        # ── Step 4: Scaffold selection ───────────────────────────────────────
        if scaffold_choice == '2':
            # BRICS decompose → pick largest fragment → STRIP dummies
            frags = BRICS.BRICSDecompose(base_mol)
            if frags:
                # Rank by heavy-atom count (after stripping dummies)
                best_frag_smiles = None
                best_ha_count = -1
                for frag_smi in frags:
                    clean_mol = _clean_smiles_to_mol(frag_smi)
                    if clean_mol is None:
                        continue
                    ha = clean_mol.GetNumHeavyAtoms()
                    if ha > best_ha_count:
                        best_ha_count = ha
                        best_frag_smiles = frag_smi

                if best_frag_smiles:
                    scaffold_mol = _clean_smiles_to_mol(best_frag_smiles)
                    if scaffold_mol is None:
                        scaffold_mol = base_mol
                else:
                    scaffold_mol = base_mol
            else:
                scaffold_mol = base_mol
        else:
            scaffold_mol = base_mol

        # ── Step 5: Generate clean 3D coordinates ────────────────────────────
        # Always regenerate from scratch (ETKDG seed=42 + UFF) so the displayed
        # structure is consistent regardless of what coords (if any) the file had.
        # This is the exact workflow from the reference code.
        scaffold_mol = _embed_and_optimise(scaffold_mol)

        # ── Step 6: Scaffold SMILES for display (clean — no [n*] dummies) ────
        scaffold_smiles = Chem.MolToSmiles(scaffold_mol, isomericSmiles=True, canonical=True)

        # ── Step 7: Serialize to SDF block for the 3D viewer ─────────────────
        sdf_block = Chem.MolToMolBlock(scaffold_mol)

        return jsonify({
            'success': True,
            'sdf_block': sdf_block,
            'scaffold_smiles': scaffold_smiles
        })

    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/generate-job-smiles', methods=['POST'])
def generate_job_smiles():
    """
    Receives the scaffold SDF block (stored client-side from /generate-scaffold)
    and the 0-based atom_idx that the user clicked in the 3Dmol viewer.

    Approach:
      1. Parse the scaffold SDF exactly as shown in the viewer, so atom_idx stays in
         the same atom space as 3Dmol.
      2. Validate that the clicked atom is a carbon with at least one available H.
      3. Expand hydrogens explicitly, replace one H on the clicked carbon with [*],
         then collapse the remaining hydrogens again. This preserves the graph and
         handles both implicit-H and explicit-H-count carbons correctly.
      4. Serialize from the scaffold's canonical start atom, after moving the dummy
         atom to the front of the atom order so RDKit keeps [*] as a branch instead
         of making it the main-chain continuation.
      5. Normalise bare '*' → '[*]' and validate round-trip.
    """
    def _get_canonical_root_atom_idx(mol):
        """
        Ask RDKit which atom it used first when writing the canonical scaffold SMILES.
        This avoids fragile isomorphism remapping and preserves the scaffold's natural
        left-to-right reading direction in the returned job SMILES.
        """
        Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        props = mol.GetPropsAsDict(includePrivate=True, includeComputed=True)
        output_order = list(props.get('_smilesAtomOutputOrder', []))
        if not output_order:
            return 0
        return int(output_order[0])

    def _add_attachment_point(mol, atom_idx):
        """
        Replace one hydrogen on atom_idx with a [*] dummy.
        This is graph-safe for both:
          - carbons with implicit Hs (CH3/CH2/aromatic CH), and
          - carbons whose heavy-atom MolBlock stores H count as explicitHs=1
            (common for chiral centres in SDF/MolBlock round-trips).
        """
        atom = mol.GetAtomWithIdx(atom_idx)

        # Safety 1: must be carbon
        if atom.GetAtomicNum() != 6:
            raise ValueError(
                f'Atom {atom_idx} is {atom.GetSymbol()} — only carbon atoms '
                'can be substitution sites.'
            )

        # Safety 2: must have at least one implicit H available.
        #
        # WHY NOT degree >= 4:
        #   GetDegree() counts bond *partners*, not bond order.
        #   A carbonyl carbon C(=O) has degree=3 (three neighbours: e.g. N, O, C)
        #   but its valence is already full: 1(N) + 2(=O) + 1(C) = 4.
        #   Adding [*] would make it 5-valent — chemically invalid.
        #   RDKit skips valence enforcement on atoms bonded to wildcards during
        #   SanitizeMol, so no error is raised and a wrong SMILES is silently emitted.
        #
        # GetTotalNumHs() correctly returns 0 for:
        #   - Carbonyl / ester / amide C (degree 3, double bond fills valence)
        #   - Quaternary C (degree 4)
        # And returns > 0 for valid sites:
        #   - Terminal CH3 (3 Hs)  •  CH2 (2 Hs)  •  CH chiral centre (1 H)
        if atom.GetTotalNumHs() == 0:
            raise ValueError(
                f'Atom {atom_idx} (C, degree={atom.GetDegree()}) has no available H — '
                'it is likely a carbonyl, ester, or fully-substituted carbon. '
                'Please select a CH, CH₂, or CH₃ carbon instead.'
            )

        mol_with_h = Chem.AddHs(Chem.Mol(mol), addCoords=True)
        rw = Chem.RWMol(mol_with_h)
        expanded_atom = rw.GetAtomWithIdx(atom_idx)
        h_neighbors = [
            nbr.GetIdx() for nbr in expanded_atom.GetNeighbors()
            if nbr.GetAtomicNum() == 1
        ]
        if not h_neighbors:
            raise ValueError(
                f'Atom {atom_idx} has no removable H after explicit-H expansion. '
                'Please select a CH, CH₂, or CH₃ carbon instead.'
            )

        # Reuse one attached hydrogen atom as the substitution dummy. This consumes
        # a real hydrogen slot instead of adding a fifth bond to the carbon.
        dummy_idx = h_neighbors[0]
        dummy_atom = rw.GetAtomWithIdx(dummy_idx)
        dummy_atom.SetAtomicNum(0)
        dummy_atom.SetIsotope(0)
        dummy_atom.SetFormalCharge(0)
        dummy_atom.SetNoImplicit(True)
        dummy_atom.SetNumExplicitHs(0)
        dummy_atom.SetNumRadicalElectrons(0)
        dummy_atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

        mol_with_dummy = rw.GetMol()
        Chem.SanitizeMol(mol_with_dummy)

        result_mol = Chem.RemoveHs(mol_with_dummy)
        Chem.SanitizeMol(result_mol)
        Chem.AssignStereochemistry(result_mol, cleanIt=True, force=True)

        subst_dummies = [
            a.GetIdx() for a in result_mol.GetAtoms()
            if a.GetAtomicNum() == 0 and a.GetIsotope() == 0
        ]
        if len(subst_dummies) != 1:
            raise ValueError(
                f'Expected exactly one [*] after H replacement, found {len(subst_dummies)}.'
            )

        return result_mol, subst_dummies[0]

    def _renumber_dummy_first(mol):
        """
        With RDKit's non-canonical rooted writer, moving the dummy atom to the front
        of the atom order keeps it as a branch on the selected carbon instead of the
        main-chain continuation.
        """
        dummy_atoms = [
            atom.GetIdx() for atom in mol.GetAtoms()
            if atom.GetAtomicNum() == 0 and atom.GetIsotope() == 0
        ]
        if len(dummy_atoms) != 1:
            raise ValueError(
                f'Expected exactly one substitution dummy, found {len(dummy_atoms)}.'
            )

        dummy_idx = dummy_atoms[0]
        new_order = [dummy_idx] + [idx for idx in range(mol.GetNumAtoms()) if idx != dummy_idx]
        renumbered = Chem.RenumberAtoms(mol, new_order)
        renumber_map = {old_idx: new_idx for new_idx, old_idx in enumerate(new_order)}
        return renumbered, renumber_map

    try:
        data     = request.get_json(silent=True) or {}
        sdf_block = data.get('sdf_block', '')
        atom_idx  = int(data.get('atom_idx', -1))

        if not sdf_block:
            return jsonify({'success': False, 'error': 'No sdf_block provided'}), 400

        # Parse the SDF that was sent to the 3Dmol viewer.
        # The SDF from /generate-scaffold is already heavy-atom only, so
        # atom indices here match exactly what atom.index returns in 3Dmol.
        scaffold_mol = Chem.MolFromMolBlock(sdf_block, removeHs=True, sanitize=True)
        if scaffold_mol is None:
            return jsonify({'success': False, 'error': 'Could not parse sdf_block'}), 400

        num_atoms = scaffold_mol.GetNumAtoms()
        if atom_idx < 0 or atom_idx >= num_atoms:
            return jsonify({
                'success': False,
                'error': f'atom_idx {atom_idx} is out of range (mol has {num_atoms} heavy atoms)'
            }), 400

        root_atom_idx = _get_canonical_root_atom_idx(Chem.Mol(scaffold_mol))

        # ── Core: replace one H with [*] at the exact clicked atom ───────────
        try:
            result_mol, _ = _add_attachment_point(scaffold_mol, atom_idx)
        except ValueError as ve:
            return jsonify({'success': False, 'error': str(ve)}), 400

        # Keep the substitution dummy as a branch during rooted SMILES writing.
        result_mol, renumber_map = _renumber_dummy_first(result_mol)
        rooted_atom_idx = renumber_map.get(root_atom_idx, 0)

        job_smiles = Chem.MolToSmiles(
            result_mol,
            rootedAtAtom=rooted_atom_idx,
            canonical=False,
            isomericSmiles=True
        )

        # Normalise bare '*' to '[*]' (RDKit sometimes emits plain * for isotope-0 dummy)
        job_smiles = re.sub(r'(?<!\[)\*(?!\])', '[*]', job_smiles)

        # ── Validate round-trip and semantics ────────────────────────────────
        check = Chem.MolFromSmiles(job_smiles)
        if check is None:
            return jsonify({
                'success': False,
                'error': f'Generated SMILES failed round-trip validation: {job_smiles}'
            }), 500

        subst_dummies = [
            a for a in check.GetAtoms()
            if a.GetAtomicNum() == 0 and a.GetIsotope() == 0
        ]
        if len(subst_dummies) != 1:
            return jsonify({
                'success': False,
                'error': f'Expected exactly one [*], found {len(subst_dummies)} in: {job_smiles}'
            }), 500

        if not any(n.GetAtomicNum() == 6 for n in subst_dummies[0].GetNeighbors()):
            return jsonify({
                'success': False,
                'error': '[*] is not bonded to a carbon atom in the generated SMILES.'
            }), 500

        return jsonify({'success': True, 'job_smiles': job_smiles})

    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'success': False, 'error': str(e)}), 500



@app.route('/run-pipeline', methods=['POST'])
def run_pipeline():
    """Legacy blocking route — kept for backward compat. Prefer /start-pipeline."""
    return start_pipeline()


@app.route('/start-pipeline', methods=['POST'])
def start_pipeline():
    """
    Non-blocking: save files, spawn pipeline thread, return job_id immediately.
    Frontend polls /pipeline-status?job_id=xxx for real-time progress.
    """
    try:
        output_dir  = OUTPUT_DIR
        uploads_dir = UPLOADS_DIR
        os.makedirs(uploads_dir, exist_ok=True)
        _job_log('SYSTEM', 'Received /start-pipeline request')

        enzyme_file = request.files['enzymeFile']
        ligand_file = request.files['ligandFile']
        enzyme_path = os.path.join(uploads_dir, secure_filename(enzyme_file.filename))
        ligand_path = os.path.join(uploads_dir, secure_filename(ligand_file.filename))
        enzyme_file.save(enzyme_path)
        ligand_file.save(ligand_path)

        scaffold_choice = request.form.get('scaffoldChoice', '1')
        use_original    = (scaffold_choice == '1')

        raw_jobs = json.loads(request.form.get('jobs', '[]'))
        substitution_jobs = []
        for job in raw_jobs:
            min_mw = float(job['min']) if job.get('min') else 0.0
            max_mw = float(job['max']) if job.get('max') else 350.0
            substitution_jobs.append({
                'attachment_smiles': job['smiles'],
                'mw_range': (min_mw, max_mw)
            })

        job_id = str(uuid.uuid4())
        _job_set(job_id, pct=1, stage='Starting…', done=False, error=None, results=None)
        _job_log(
            job_id,
            f"Queued pipeline request (enzyme={os.path.basename(enzyme_path)}, "
            f"ligand={os.path.basename(ligand_path)}, jobs={len(substitution_jobs)})"
        )

        def _run():
            started_at = time.time()
            stream = JobLogStream(job_id)
            try:
                _job_log(job_id, "Worker thread started")
                with redirect_stdout(stream), redirect_stderr(stream):
                    last_progress = {'pct': None, 'stage': None}

                    def _cb(pct, stage=''):
                        _job_set(job_id, pct=pct, stage=stage)
                        if pct != last_progress['pct'] or stage != last_progress['stage']:
                            stage_msg = f" — {stage}" if stage else ""
                            _job_log(job_id, f"Progress {pct}%{stage_msg}")
                            last_progress['pct'] = pct
                            last_progress['stage'] = stage

                    _cb(1, 'Initialising processor…')

                    processor = EnzymeLigandProcessor(
                        enzyme_file=enzyme_path,
                        molecule_file=ligand_path,
                        output_dir=output_dir,
                        use_original_as_scaffold=use_original,
                        substitution_jobs=substitution_jobs,
                        interactive=False
                    )

                    processor.active_site_atoms = []
                    _job_log(job_id, "Processor ready; starting process_all()")
                    results = processor.process_all(
                        use_original_as_scaffold=use_original,
                        progress_cb=_cb
                    )

                    _job_log(job_id, f"Pipeline core completed in {time.time() - started_at:.1f}s")

                    ligands_generated = 0
                    if results.get('docking_results'):
                        for v in results['docking_results'].values():
                            ligands_generated += v.get('total_docked', 0)

                    best_affinity = 'N/A'
                    results_csv = results.get('results_csv', '')
                    if results_csv and os.path.exists(results_csv):
                        try:
                            with open(results_csv) as f:
                                reader = csv.DictReader(f)
                                vals = []
                                for row in reader:
                                    for col in ('Conformation_1_Affinity', 'affinity',
                                                'Affinity', 'binding_affinity', 'score', 'Score'):
                                        if row.get(col):
                                            try:
                                                vals.append(float(row[col]))
                                            except ValueError:
                                                pass
                                            break
                                if vals:
                                    best_affinity = f"{min(vals):.2f}"
                        except Exception as e:
                            _job_log(job_id, f"Could not parse affinity from CSV: {e}")

                    csv_rel = zip_rel = ''
                    if results_csv:
                        try:
                            csv_rel = os.path.relpath(results_csv, output_dir)
                        except ValueError:
                            csv_rel = os.path.basename(results_csv)

                    zip_path = os.path.join(output_dir, 'top_ligands.zip')
                    if os.path.exists(zip_path):
                        zip_rel = 'top_ligands.zip'

                    _job_set(
                        job_id,
                        pct=100,
                        stage='Complete',
                        done=True,
                        error=None,
                        results={
                            'success': True,
                            'scaffold': results.get('scaffold_smiles', 'N/A'),
                            'ligands': ligands_generated,
                            'affinity': f"{best_affinity} kcal/mol",
                            'csv': csv_rel,
                            'zip': zip_rel,
                        }
                    )
                    _job_log(job_id, "Job completed successfully")

            except Exception as exc:
                import traceback
                tb = traceback.format_exc()
                _job_log(job_id, tb)
                _job_log(job_id, f"Pipeline failed: {exc}")
                _job_set(job_id, pct=0, stage='Error', done=True, error=str(exc), results=None)
            
            finally:
                stream.flush()

        threading.Thread(target=_run, daemon=True).start()
        return jsonify({'success': True, 'job_id': job_id})

    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'success': False, 'error': str(e)}), 500


@app.route('/pipeline-status')
def pipeline_status():
    """Return current progress for a running pipeline job."""
    job_id = request.args.get('job_id', '')
    job = _job_get(job_id)
    if not job:
        return jsonify({'success': False, 'error': 'Unknown job_id'}), 404
    return jsonify({
        'success': True,
        'pct':     job.get('pct', 0),
        'stage':   job.get('stage', ''),
        'done':    job.get('done', False),
        'error':   job.get('error'),
        'results': job.get('results'),
    })


@app.route('/pipeline-log')
def pipeline_log():
    """
    Return log lines accumulated for a running/completed job.
    Supports cursor-based incremental fetch so the frontend only
    receives new lines each poll.

    Query params:
        job_id  – required
        offset  – optional int; return lines from this index onward (default 0)
    """
    job_id = request.args.get('job_id', '')
    try:
        offset = int(request.args.get('offset', 0))
    except (ValueError, TypeError):
        offset = 0

    job = _job_get(job_id)
    if not job:
        return jsonify({'success': False, 'error': 'Unknown job_id'}), 404

    all_logs = job.get('logs', [])
    new_lines = all_logs[offset:]
    return jsonify({
        'success': True,
        'lines':   new_lines,
        'total':   len(all_logs),   # next offset the client should send
        'done':    job.get('done', False),
    })


@app.route('/download')
def download_file():
    filename = request.args.get('file', '')
    if not filename:
        abort(400)
    output_dir = OUTPUT_DIR
    full_path  = os.path.realpath(os.path.join(output_dir, filename))
    if not full_path.startswith(os.path.realpath(output_dir)):
        abort(403)
    return send_from_directory(output_dir, filename, as_attachment=True)




if __name__ == "__main__":
    import socket
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(line_buffering=True)

    def _find_free_port(start=5001, max_tries=20):
        # Skip 5000 — macOS AirPlay Receiver squats on it.
        # Plain bind (no SO_REUSEADDR) gives a truthful answer.
        for port in range(start, start + max_tries):
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                try:
                    s.bind(("127.0.0.1", port))
                    return port
                except OSError:
                    continue
        raise RuntimeError(
            f"No free port found in range {start}-{start + max_tries - 1}."
        )

    PORT = _find_free_port(5001)

    _actual_template = (
        'enzyme_processor_patched.html'
        if os.path.exists(os.path.join(SCRIPT_DIR, 'enzyme_processor_patched.html'))
        else 'enzyme_processor.html'
    )
    print("=" * 60)
    print("Enzyme Drug Discovery Pipeline - Server Starting")
    print("=" * 60)
    print(f"\nServer URL: http://127.0.0.1:{PORT}")
    print(f"Template:   {os.path.join(SCRIPT_DIR, _actual_template)}")
    print(f"Runtime log: {RUNTIME_LOG}")
    print("\nNote: Debug mode enabled, but reloader disabled for Spyder compatibility")
    print("Tip: to reclaim port 5000, disable AirPlay Receiver in")
    print("     System Settings > General > AirDrop & Handoff")
    print("Press CTRL+C to stop the server")
    print("=" * 60)

    app.run(debug=True, port=PORT, use_reloader=False)
