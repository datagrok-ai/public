#language: python
#name: molToHelmConverterPy
#description: Converts molecules to HELM notation based on monomer library
#input: dataframe moleculesDataframe
#input: column moleculesColumn {semType: Molecule}
#input: string libraryJSON
#output: dataframe result_helm {action:join(moleculesDataframe)} [Sequences, in HELM format]
molListToProcess = moleculesDataframe[moleculesColumn].tolist()
import pandas as pd
import numpy as np
"""
Aggregated file combining all modules from the logics folder.
Generated automatically - do not edit manually.
"""

# External library imports
from collections import defaultdict
from enum import Enum
from itertools import combinations
from rdkit import Chem
from rdkit import RDLogger
from typing import Dict
from typing import List
from typing import Optional
from typing import Tuple
import json
import os

# ============================================================================
# Content from: fragment_graph.py
# ============================================================================

from rdkit import Chem
from typing import Optional, List, Dict, Tuple
from enum import Enum


class LinkageType(Enum):
    """Types of linkages between fragments"""
    PEPTIDE = "peptide"
    DISULFIDE = "disulfide"
    ESTER = "ester"
    ETHER = "ether"
    THIOETHER = "thioether"
    UNKNOWN = "unknown"


class FragmentNode:
    """Represents a single molecular fragment (amino acid/monomer)"""
    
    def __init__(self, fragment_id: int, mol: Chem.Mol):
        self.id = fragment_id
        self.mol = mol  # RDKit molecule object
        self.smiles = Chem.MolToSmiles(mol, canonical=True) if mol else ""
        self.monomer = None  # Will be filled by matcher - MonomerData object
        self.is_c_terminal = False
        self.is_n_terminal = False
    
    def __repr__(self):
        monomer_name = self.monomer.symbol if self.monomer else "Unknown"
        return f"FragmentNode(id={self.id}, monomer={monomer_name}, smiles={self.smiles[:20]}...)"


class FragmentLink:
    """Represents a connection between two fragments"""
    
    def __init__(
        self, 
        from_node_id: int, 
        to_node_id: int,
        linkage_type: LinkageType,
        from_atom_idx: Optional[int] = None,
        to_atom_idx: Optional[int] = None
    ):
        self.from_node_id = from_node_id
        self.to_node_id = to_node_id
        self.linkage_type = linkage_type
        self.from_atom_idx = from_atom_idx  # Atom index in from_node's molecule
        self.to_atom_idx = to_atom_idx      # Atom index in to_node's molecule
    
    def __repr__(self):
        return f"FragmentLink({self.from_node_id} --{self.linkage_type.value}--> {self.to_node_id})"


class FragmentGraph:
    """
    Graph structure representing a molecule as fragments and their connections.
    
    Supports:
    - Linear peptides (chain of peptide bonds)
    - Cyclic peptides (peptide bond from last to first)
    - Disulfide bridges (additional S-S links)
    - Branched structures (multiple links per fragment)
    """
    
    def __init__(self):
        self.nodes: Dict[int, FragmentNode] = {}  # node_id -> FragmentNode
        self.links: List[FragmentLink] = []
    
    def add_node(self, node: FragmentNode) -> int:
        """Add a fragment node to the graph"""
        self.nodes[node.id] = node
        return node.id
    
    def add_link(self, link: FragmentLink):
        """Add a linkage between two nodes"""
        if link.from_node_id not in self.nodes or link.to_node_id not in self.nodes:
            raise ValueError(f"Cannot add link: nodes {link.from_node_id} or {link.to_node_id} not in graph")
        self.links.append(link)
    
    def get_node(self, node_id: int) -> Optional[FragmentNode]:
        """Get a node by ID"""
        return self.nodes.get(node_id)
    
    def get_neighbors(self, node_id: int) -> List[Tuple[int, LinkageType]]:
        """Get all neighbors of a node with their linkage types"""
        neighbors = []
        for link in self.links:
            if link.from_node_id == node_id:
                neighbors.append((link.to_node_id, link.linkage_type))
            elif link.to_node_id == node_id:
                neighbors.append((link.from_node_id, link.linkage_type))
        return neighbors
    
    def get_ordered_nodes(self) -> List[FragmentNode]:
        """
        Get nodes in sequential order (for linear/cyclic peptides).
        For branched structures, returns a depth-first traversal.
        """
        if not self.nodes:
            return []
        
        # Find starting node (N-terminal for peptides)
        start_node_id = None
        for node_id, node in self.nodes.items():
            if node.is_n_terminal:
                start_node_id = node_id
                break
        
        # If no N-terminal found, use first node
        if start_node_id is None:
            start_node_id = min(self.nodes.keys())
        
        # Traverse the graph
        ordered = []
        visited = set()
        self._traverse_from_node(start_node_id, visited, ordered)
        
        return ordered
    
    def _traverse_from_node(self, node_id: int, visited: set, ordered: list):
        """Helper for depth-first traversal"""
        if node_id in visited:
            return
        
        visited.add(node_id)
        ordered.append(self.nodes[node_id])
        
        # Get peptide bond neighbors first (to maintain chain order)
        peptide_neighbors = []
        other_neighbors = []
        
        for link in self.links:
            if link.from_node_id == node_id and link.to_node_id not in visited:
                if link.linkage_type == LinkageType.PEPTIDE:
                    peptide_neighbors.append(link.to_node_id)
                else:
                    other_neighbors.append(link.to_node_id)
        
        # Visit peptide bonds first, then others
        for neighbor_id in peptide_neighbors + other_neighbors:
            self._traverse_from_node(neighbor_id, visited, ordered)
    
    def get_fragment_sequence(self) -> List[str]:
        """Get sequence of monomer symbols (for matched fragments)"""
        ordered_nodes = self.get_ordered_nodes()
        return [
            node.monomer.symbol if node.monomer else f"X{node.id}"
            for node in ordered_nodes
        ]
    
    def is_cyclic(self) -> bool:
        """
        Detect if the peptide is cyclic.
        A cyclic peptide has a peptide bond connecting the last residue back to near the beginning.
        Handles cases where N-terminal caps (like 'ac' from Lys_Ac) create an extra fragment at position 0.
        """
        if len(self.nodes) < 3:
            return False
        
        # Get ordered nodes
        ordered = self.get_ordered_nodes()
        if len(ordered) < 3:
            return False
        
        # Get the last node ID
        last_id = ordered[-1].id
        
        # For a cyclic peptide, the last residue should connect back to one of the first few residues
        # (usually first, but could be second if there's an N-terminal cap like 'ac')
        # Check if last node has a peptide bond to any of the first 3 nodes
        first_few_ids = [ordered[i].id for i in range(min(3, len(ordered)))]
        
        for link in self.links:
            if link.linkage_type == LinkageType.PEPTIDE:
                # Check if link connects last node to one of the first few nodes
                if (link.from_node_id == last_id and link.to_node_id in first_few_ids) or \
                   (link.to_node_id == last_id and link.from_node_id in first_few_ids):
                    return True
        
        return False
    
    def __len__(self):
        return len(self.nodes)
    
    def __repr__(self):
        return f"FragmentGraph(nodes={len(self.nodes)}, links={len(self.links)})"
    
    def to_dict(self) -> dict:
        """Convert graph to dictionary for serialization"""
        return {
            "nodes": [
                {
                    "id": node.id,
                    "smiles": node.smiles,
                    "monomer": node.monomer.symbol if node.monomer else None,
                    "is_n_terminal": node.is_n_terminal,
                    "is_c_terminal": node.is_c_terminal
                }
                for node in self.nodes.values()
            ],
            "links": [
                {
                    "from": link.from_node_id,
                    "to": link.to_node_id,
                    "type": link.linkage_type.value,
                    "from_atom": link.from_atom_idx,
                    "to_atom": link.to_atom_idx
                }
                for link in self.links
            ]
        }

# ============================================================================
# Content from: fragment_processor.py
# ============================================================================

from rdkit import Chem


class BondDetector:
    #GENERALIZATION ITEM: BOND PATTERNS SHOULD BE DERIVED FROM LIBRARY
    def __init__(self):
        # True peptide bond: C and N both in backbone (each bonded to carbons)
        # First carbon can be aliphatic or aromatic (for amino acids like NMe2Abz)
        # Carbonyl carbon is sp2 (X3)
        # Exclude if carbonyl is in a small ring (r5 or r6) to avoid cleaving lactams like Pyr
        # !r5 = not in 5-membered ring, !r6 = not in 6-membered ring
        # This preserves lactams but allows large macrocycles and proline (C=O outside ring)
        # Nitrogen can be X2 (proline, imino) or X3 (standard amino, N-methyl)
        # N-C bond can be single (-) or double (=) for imine bonds in dehydro amino acids
        # Alpha carbon after N can be sp3 (X4) or sp2 (X3) for dehydroamino acids
        self.peptide_bond = Chem.MolFromSmarts('[#6]-[C;X3;!r5;!r6](=[O;X1])-[N;X2,X3]~[C;X3,X4]')
        # True disulfide bond: S-S where each S is bonded to carbon (cysteine residues)
        self.disulfide_bond = Chem.MolFromSmarts('[C;X4]-[S;X2]-[S;X2]-[C;X4]')
        # Primary amine at N-terminus (can be NH2 or NH3+), alpha-C can be sp3 or sp2
        self.primary_amine = Chem.MolFromSmarts('[N;H2,H3;X3,X4]-[C;X3,X4]')

    def find_cleavable_bonds(self, mol: Chem.Mol):
        """
        Find all cleavable bonds in the molecule.
        
        Returns:
            List of tuples: (atom1_idx, atom2_idx, LinkageType)
        """
        try:
            all_bonds = []

            # Find peptide bonds
            peptide_bonds = self._find_peptide_bonds(mol)
            all_bonds.extend([(bond[0], bond[1], LinkageType.PEPTIDE) for bond in peptide_bonds])
            
            # Find disulfide bonds
            disulfide_bonds = self._find_disulfide_bonds(mol)
            all_bonds.extend([(bond[0], bond[1], LinkageType.DISULFIDE) for bond in disulfide_bonds])

            # Order peptide bonds from N to C (keep disulfide bonds unordered)
            peptide_only = [(b[0], b[1]) for b in all_bonds if b[2] == LinkageType.PEPTIDE]
            ordered_peptide = self._order_bonds_from_n_to_c(mol, peptide_only)
            
            # Rebuild with types
            ordered_bonds = [(b[0], b[1], LinkageType.PEPTIDE) for b in ordered_peptide]
            ordered_bonds.extend([b for b in all_bonds if b[2] != LinkageType.PEPTIDE])

            return ordered_bonds

        except Exception:
            return []

    def _find_peptide_bonds(self, mol: Chem.Mol):
        bonds = []
        try:
            matches = mol.GetSubstructMatches(self.peptide_bond)
            for match in matches:
                if len(match) >= 5:
                    # Pattern: [C;X3,X4]-[C;X3](=[O;X1])-[N;X2,X3]~[C;X3,X4]
                    # match[0]=alpha-C (sp2 or sp3), match[1]=carbonyl-C, match[2]=O, match[3]=N, match[4]=next-alpha-C (sp2 or sp3)
                    c_atom = match[1]  # Carbonyl carbon
                    n_atom = match[3]  # Nitrogen
                    bonds.append((c_atom, n_atom))
        except Exception:
            pass
        return bonds
    
    def _find_disulfide_bonds(self, mol: Chem.Mol):
        """Find disulfide bonds (S-S linkages)"""
        bonds = []
        try:
            matches = mol.GetSubstructMatches(self.disulfide_bond)
            for match in matches:
                if len(match) >= 4:
                    # Pattern: [C;X4]-[S;X2]-[S;X2]-[C;X4]
                    # match[0]=C, match[1]=S, match[2]=S, match[3]=C
                    s1_atom = match[1]  # First sulfur
                    s2_atom = match[2]  # Second sulfur
                    bonds.append((s1_atom, s2_atom))
        except Exception:
            pass
        return bonds

    def _order_bonds_from_n_to_c(self, mol: Chem.Mol, bonds):
        if not bonds:
            return bonds

        n_terminal = self._find_n_terminal(mol)
        if n_terminal is None:
            return bonds

        ordered = []
        visited = set()
        current = n_terminal

        while current is not None and len(ordered) < len(bonds):
            next_bond = None
            for bond in bonds:
                if bond not in visited and bond[1] == current:
                    next_bond = bond
                    break

            if next_bond is None:
                break

            ordered.append(next_bond)
            visited.add(next_bond)
            current = next_bond[0]

        for bond in bonds:
            if bond not in visited:
                ordered.append(bond)

        return ordered

    def _find_n_terminal(self, mol: Chem.Mol):
        try:
            matches = mol.GetSubstructMatches(self.primary_amine)
            if matches:
                return matches[0][0]

            max_h = -1
            n_term = None
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 7:
                    h_count = atom.GetTotalNumHs()
                    if h_count > max_h:
                        max_h = h_count
                        n_term = atom.GetIdx()
            return n_term

        except Exception:
            return None


class FragmentProcessor:
    def __init__(self, monomer_library):
        self.monomer_library = monomer_library
        self.bond_detector = BondDetector()

    def process_molecule(self, mol: Chem.Mol) -> FragmentGraph:
        """
        Process a molecule into a fragment graph.
        
        Args:
            mol: RDKit molecule object
        
        Returns:
            FragmentGraph object containing fragments and their connections
        """
        graph = FragmentGraph()
        # Store original molecule for fragment recovery
        graph.original_mol = mol
        
        try:
            bonds_to_cleave = self.bond_detector.find_cleavable_bonds(mol)

            if not bonds_to_cleave:
                # Single fragment (no cleavable bonds)
                node = FragmentNode(0, mol)
                node.is_n_terminal = True
                node.is_c_terminal = True
                graph.add_node(node)
                return graph

            # Extract bond info for fragmentation
            bond_indices = []
            bond_info = []  # (bond_idx, atom1, atom2, linkage_type)
            seen_bonds = set()  # Track which bonds we've already added
            
            for atom1, atom2, linkage_type in bonds_to_cleave:
                bond = mol.GetBondBetweenAtoms(atom1, atom2)
                if bond:
                    bond_idx = bond.GetIdx()
                    if bond_idx not in seen_bonds:
                        bond_indices.append(bond_idx)
                        bond_info.append((bond_idx, atom1, atom2, linkage_type))
                        seen_bonds.add(bond_idx)
                    # Skip duplicate bonds silently
                # Skip invalid bonds silently

            if not bond_indices:
                # No valid bonds found
                node = FragmentNode(0, mol)
                node.is_n_terminal = True
                node.is_c_terminal = True
                graph.add_node(node)
                return graph

            # Fragment the molecule
            fragmented_mol = Chem.FragmentOnBonds(mol, bond_indices, addDummies=True)
            
            # Get fragments as molecules
            fragments = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)
            
            # Get atom mappings separately (which original atoms are in which fragment)
            atom_mappings = Chem.GetMolFrags(fragmented_mol, asMols=False, fragsMolAtomMapping=True)
            
            # Store bond cleavage info for recovery - we'll use this to selectively re-fragment
            graph.cleaved_bond_indices = bond_indices
            graph.bond_info = bond_info
            graph.atom_mappings = atom_mappings
            print(f"DEBUG: Created {len(fragments)} fragments, cleaved {len(bond_indices)} bonds")

            # Create nodes for each fragment
            fragment_nodes = []
            for i, frag in enumerate(fragments):
                clean_frag = self._clean_fragment(frag)
                if clean_frag and clean_frag.GetNumAtoms() >= 3:
                    is_c_terminal = (i == len(fragments) - 1)
                    is_n_terminal = (i == 0)
                    # No normalization! Use fragment as-is
                    node = FragmentNode(i, clean_frag)
                    node.is_c_terminal = is_c_terminal
                    node.is_n_terminal = is_n_terminal
                    graph.add_node(node)
                    fragment_nodes.append((i, node))

            # Create links between fragments based on the actual cleaved bonds
            # Build mapping: original atom index → (fragment_idx, new_atom_idx_in_fragment)
            atom_to_fragment_and_idx = {}
            for frag_idx, original_atom_indices in enumerate(atom_mappings):
                for new_idx_in_frag, original_atom_idx in enumerate(original_atom_indices):
                    atom_to_fragment_and_idx[original_atom_idx] = (frag_idx, new_idx_in_frag)
            
            print(f"DEBUG: Processing {len(bond_info)} cleaved bonds to create links")
            print(f"DEBUG: atom_to_fragment_and_idx has {len(atom_to_fragment_and_idx)} entries")
            
            # For each cleaved bond, determine which fragments it connects
            link_count = 0
            for bond_idx, atom1_orig, atom2_orig, linkage_type in bond_info:
                # Find which fragments contain these atoms and their new indices
                frag1_info = atom_to_fragment_and_idx.get(atom1_orig)
                frag2_info = atom_to_fragment_and_idx.get(atom2_orig)
                
                if frag1_info is None or frag2_info is None:
                    print(f"DEBUG: Skipping bond atoms {atom1_orig}-{atom2_orig}: not found in fragments")
                    continue
                
                frag1, atom1_new = frag1_info
                frag2, atom2_new = frag2_info
                    
                # Create link even if both atoms are in same fragment (internal bond like in Phe_4Sdihydroorotamido)
                # This creates a "self-link" that will be used during recovery to reconstruct the monomer
                link = FragmentLink(frag1, frag2, linkage_type, 
                                   from_atom_idx=atom1_new, to_atom_idx=atom2_new)
                graph.add_link(link)
                link_count += 1
                
                if frag1 == frag2:
                    print(f"DEBUG: Link {link_count}: {linkage_type.value.upper()} SELF-LINK frag{frag1} "
                          f"orig_atoms({atom1_orig}<->{atom2_orig}) frag_atoms({atom1_new}<->{atom2_new})")
                else:
                    print(f"DEBUG: Link {link_count}: {linkage_type.value.upper()} frag{frag1}<->frag{frag2} "
                          f"orig_atoms({atom1_orig}<->{atom2_orig}) frag_atoms({atom1_new}<->{atom2_new})")
            
            print(f"DEBUG: Created {link_count} links total")

            return graph

        except Exception as e:
            # Fallback: single node with original molecule
            node = FragmentNode(0, mol)
            node.is_n_terminal = True
            node.is_c_terminal = True
            graph.add_node(node)
            return graph

    def _clean_fragment(self, mol: Chem.Mol):
        try:
            mol_copy = Chem.Mol(mol)
            atoms_to_remove = []

            for atom in mol_copy.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    atoms_to_remove.append(atom.GetIdx())

            atoms_to_remove.sort(reverse=True)
            if atoms_to_remove:
                emol = Chem.EditableMol(mol_copy)
                for atom_idx in atoms_to_remove:
                    emol.RemoveAtom(atom_idx)
                return emol.GetMol()

            return mol_copy

        except Exception:
            return None

    def _reconstruct_fragment_with_links(self, node_ids: list, graph: FragmentGraph, 
                                         links_to_exclude: list) -> Chem.Mol:
        """
        Reconstruct a molecule by combining multiple fragment nodes, using link information.
        
        Args:
            node_ids: List of node IDs to merge
            graph: The fragment graph
            links_to_exclude: List of FragmentLink objects connecting the nodes to merge
        
        Returns:
            Combined RDKit molecule, or None if reconstruction fails
        """
        if not node_ids or not hasattr(graph, 'original_mol'):
            return None
        
        if not hasattr(graph, 'cleaved_bond_indices') or not hasattr(graph, 'bond_info'):
            return None
        
        try:
            # Find which bond indices correspond to the links we want to exclude
            bonds_to_exclude_indices = []
            
            for link in links_to_exclude:
                # Find the bond_info entry that matches this link's original atoms
                # We need to find which bond connected these fragments
                for bond_list_idx, (bond_idx, atom1, atom2, linkage_type) in enumerate(graph.bond_info):
                    # Check if this bond connects the fragments in this link
                    if hasattr(graph, 'atom_mappings'):
                        # Find which fragments contain these atoms
                        frag1 = None
                        frag2 = None
                        for frag_idx, atom_indices in enumerate(graph.atom_mappings):
                            if atom1 in atom_indices:
                                frag1 = frag_idx
                            if atom2 in atom_indices:
                                frag2 = frag_idx
                        
                        # If this bond connects the two fragments in the link, exclude it
                        if (frag1 == link.from_node_id and frag2 == link.to_node_id) or \
                           (frag1 == link.to_node_id and frag2 == link.from_node_id):
                            bonds_to_exclude_indices.append(bond_list_idx)
                            print(f"DEBUG: Excluding {linkage_type.value} bond at index {bond_list_idx} (atoms {atom1}<->{atom2})")
                            break
            
            # Create new bond list excluding the bonds we want to keep
            new_bond_indices = [
                bond_idx for i, bond_idx in enumerate(graph.cleaved_bond_indices)
                if i not in bonds_to_exclude_indices
            ]
            
            print(f"DEBUG reconstruct: Original had {len(graph.cleaved_bond_indices)} cleaved bonds, "
                  f"excluding {len(bonds_to_exclude_indices)} bonds, new list has {len(new_bond_indices)} bonds")
            
            # Re-fragment with the modified bond list
            if not new_bond_indices:
                # No bonds to cleave - return whole molecule
                return graph.original_mol
            
            fragmented_mol = Chem.FragmentOnBonds(graph.original_mol, new_bond_indices, addDummies=True)
            fragments = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)
            new_atom_mappings = Chem.GetMolFrags(fragmented_mol, asMols=False, fragsMolAtomMapping=True)
            
            # Find which new fragment contains atoms from our target nodes
            # Look for the fragment that contains atoms from the first node we want to merge
            sorted_nodes = sorted(node_ids)
            first_node_atoms = set(graph.atom_mappings[sorted_nodes[0]])
            
            target_fragment_idx = None
            for new_frag_idx, new_atoms in enumerate(new_atom_mappings):
                # Check if this new fragment contains any atoms from our first target node
                if first_node_atoms & set(new_atoms):
                    target_fragment_idx = new_frag_idx
                    break
            
            print(f"DEBUG reconstruct: Got {len(fragments)} fragments after re-fragmentation, "
                  f"target_fragment_idx={target_fragment_idx}")
            
            if target_fragment_idx is not None and target_fragment_idx < len(fragments):
                clean_frag = self._clean_fragment(fragments[target_fragment_idx])
                return clean_frag if clean_frag else fragments[target_fragment_idx]
            
            return None
        
        except Exception as e:
            print(f"DEBUG reconstruct: Exception: {e}")
            return None
    
    def _reconstruct_fragment(self, node_ids: list, graph: FragmentGraph) -> Chem.Mol:
        """
        Reconstruct a molecule by combining multiple fragment nodes.
        Re-fragments the original molecule, excluding bonds between the nodes to merge.
        """
        if not node_ids or not hasattr(graph, 'original_mol') or not hasattr(graph, 'cleaved_bond_indices'):
            return None
        
        try:
            # Sort node IDs to ensure consistent ordering
            sorted_nodes = sorted(node_ids)
            
            # Identify which bonds to exclude (bonds between consecutive merged nodes)
            bonds_to_exclude = set()
            for i in range(len(sorted_nodes) - 1):
                # We want to keep the bond between node i and node i+1
                # This bond would be at position sorted_nodes[i] in the cleaved_bond_indices
                if sorted_nodes[i] + 1 == sorted_nodes[i + 1]:
                    # Consecutive nodes - exclude the bond between them
                    if sorted_nodes[i] < len(graph.cleaved_bond_indices):
                        bonds_to_exclude.add(sorted_nodes[i])
            
            # Create new bond list excluding the bonds we want to keep
            new_bond_indices = [
                bond_idx for i, bond_idx in enumerate(graph.cleaved_bond_indices)
                if i not in bonds_to_exclude
            ]
            
            print(f"DEBUG reconstruct: Original had {len(graph.cleaved_bond_indices)} cleaved bonds, "
                  f"excluding {len(bonds_to_exclude)} bonds, new list has {len(new_bond_indices)} bonds")
            
            # Re-fragment with the modified bond list
            if not new_bond_indices:
                # No bonds to cleave - return whole molecule
                return graph.original_mol
            
            fragmented_mol = Chem.FragmentOnBonds(graph.original_mol, new_bond_indices, addDummies=True)
            fragments_tuple = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)
            fragments = list(fragments_tuple)
            
            # Find which fragment corresponds to our merged nodes
            # The merged nodes should be at the position of the first node ID in sorted order
            target_idx = sorted_nodes[0]
            
            # Account for excluded bonds shifting fragment indices
            adjusted_idx = target_idx - sum(1 for excluded_idx in bonds_to_exclude if excluded_idx < target_idx)
            
            print(f"DEBUG reconstruct: Got {len(fragments)} fragments after re-fragmentation, "
                  f"target_idx={target_idx}, adjusted_idx={adjusted_idx}")
            
            if adjusted_idx < len(fragments):
                clean_frag = self._clean_fragment(fragments[adjusted_idx])
                return clean_frag if clean_frag else fragments[adjusted_idx]
            
            return None
        
        except Exception as e:
            print(f"DEBUG reconstruct: Exception: {e}")
            return None

    def _merge_nodes_in_graph(self, graph: FragmentGraph, nodes_to_merge: list, 
                              new_node: FragmentNode) -> None:
        """
        Remove old nodes, add new merged node, update all links.
        Preserves terminal flags from edge nodes.
        """
        if not nodes_to_merge:
            return
        
        # Sort node IDs to identify edge nodes
        sorted_nodes = sorted(nodes_to_merge)
        leftmost = sorted_nodes[0]
        rightmost = sorted_nodes[-1]
        
        # Preserve terminal flags
        if leftmost in graph.nodes:
            new_node.is_n_terminal = graph.nodes[leftmost].is_n_terminal
        if rightmost in graph.nodes:
            new_node.is_c_terminal = graph.nodes[rightmost].is_c_terminal
        
        # Update links: replace references to merged nodes
        updated_links = []
        nodes_to_merge_set = set(nodes_to_merge)
        
        for link in graph.links:
            from_in = link.from_node_id in nodes_to_merge_set
            to_in = link.to_node_id in nodes_to_merge_set
            
            # Skip internal links between merged nodes
            if from_in and to_in:
                continue
            
            # Update link if one end is being merged
            new_from = new_node.id if from_in else link.from_node_id
            new_to = new_node.id if to_in else link.to_node_id
            
            updated_links.append(FragmentLink(new_from, new_to, link.linkage_type))
        
        # Remove old nodes
        for node_id in nodes_to_merge:
            if node_id in graph.nodes:
                del graph.nodes[node_id]
        
        # Add new node and update links
        graph.add_node(new_node)
        graph.links = updated_links

    def recover_unmatched_fragments(self, graph: FragmentGraph, matcher) -> bool:
        """
        Try to recover unmatched fragments by merging with neighbors based on graph links.
        Returns True if any merges were successful.
        """
        # Identify unmatched nodes
        unmatched_nodes = []
        for node_id, node in graph.nodes.items():
            if node.monomer and node.monomer.symbol.startswith("X"):
                unmatched_nodes.append(node_id)
        
        if not unmatched_nodes:
            return False
        
        print(f"DEBUG: Found {len(unmatched_nodes)} unmatched nodes: {unmatched_nodes}")
        
        had_changes = False
        
        # Try to recover each unmatched node
        for node_id in unmatched_nodes:
            # Check if node still exists (might have been merged already)
            if node_id not in graph.nodes:
                continue
            
            # Get neighbors from graph links (returns list of (neighbor_id, linkage_type))
            neighbors = graph.get_neighbors(node_id)
            
            if not neighbors:
                print(f"DEBUG: Node {node_id} has no neighbors")
                continue
            
            print(f"DEBUG: Node {node_id} neighbors: {[(n[0], n[1].value) for n in neighbors]}")
            
            # Try merging with each individual neighbor first
            for neighbor_id, linkage_type in neighbors:
                if neighbor_id not in graph.nodes:
                    continue
                    
                nodes_to_merge = sorted([node_id, neighbor_id])
                print(f"DEBUG: Trying to merge nodes {nodes_to_merge} (via {linkage_type.value} bond)")
                
                # Find the links between nodes we're merging
                links_to_exclude = []
                for link in graph.links:
                    from_in = link.from_node_id in nodes_to_merge
                    to_in = link.to_node_id in nodes_to_merge
                    if from_in and to_in:
                        links_to_exclude.append(link)
                
                # Reconstruct combined molecule
                combined_mol = self._reconstruct_fragment_with_links(nodes_to_merge, graph, links_to_exclude)
                if not combined_mol:
                    print(f"DEBUG: Failed to reconstruct molecule for {nodes_to_merge}")
                    continue
                
                print(f"DEBUG: Reconstructed mol with {combined_mol.GetNumAtoms()} atoms")
                
                # Count expected connections for this merged fragment
                # Get all unique neighbors of the merged set (excluding internal connections)
                all_neighbors = set()
                for nid in nodes_to_merge:
                    if nid in graph.nodes:
                        node_neighbors = graph.get_neighbors(nid)
                        for neighbor_id, _ in node_neighbors:
                            if neighbor_id not in nodes_to_merge:
                                all_neighbors.add(neighbor_id)
                
                num_connections = len(all_neighbors)
                print(f"DEBUG: Expecting {num_connections} connections")
                
                # Try to match the combined fragment
                monomer = matcher.find_exact_match(combined_mol, num_connections)
                
                if monomer:
                    print(f"DEBUG: SUCCESS! Matched to {monomer.symbol}")
                    # Success! Create new merged node
                    new_node_id = min(nodes_to_merge)
                    new_node = FragmentNode(new_node_id, combined_mol)
                    new_node.monomer = monomer
                    
                    # Merge nodes in graph
                    self._merge_nodes_in_graph(graph, nodes_to_merge, new_node)
                    
                    had_changes = True
                    break  # Stop trying other neighbors for this node
                else:
                    print(f"DEBUG: No match found for merge {nodes_to_merge}")
            
            if had_changes:
                break  # Restart from beginning after a successful merge
        
        return had_changes

# ============================================================================
# Content from: helm_generator.py
# ============================================================================

class HELMGenerator:
    """
    Generates HELM notation from fragment graphs or monomer lists.
    
    Supports:
    - Linear peptides
    - Cyclic peptides
    - Disulfide bridges
    - Custom linkages
    """
    
    def __init__(self):
        #GENERALIZATION ITEM: POLYMER TYPES SHOULD BE DERIVED FROM LIBRARY
        self.polymer_types = {
            "peptide": "PEPTIDE",
            "rna": "RNA",
            "dna": "DNA",
            "chemical": "CHEM"
        }

    def generate_helm_from_graph(self, graph: FragmentGraph) -> str:
        """
        Generate HELM notation from a FragmentGraph.
        
        Args:
            graph: FragmentGraph containing matched monomers and their connections
        
        Returns:
            HELM notation string
        """
        if len(graph) == 0:
            return ""
        
        # Get ordered sequence of monomers
        ordered_nodes = graph.get_ordered_nodes()
        sequence_symbols = [node.monomer.symbol if node.monomer else "X" for node in ordered_nodes]
        
        # Check if cyclic
        is_cyclic = graph.is_cyclic()
        
        # Generate sequence notation
        if is_cyclic:
            # Cyclic: wrap multi-letter monomers in brackets, single-letter ones stay as-is
            formatted_symbols = [f"[{symbol}]" if len(symbol) > 1 else symbol for symbol in sequence_symbols]
            sequence = ".".join(formatted_symbols)
        else:
            # Linear: no brackets
            sequence = ".".join(sequence_symbols)
        
        # Collect non-sequential connections (disulfide bridges, cyclic bonds, etc.)
        connections = []
        
        if is_cyclic:
            # Find the actual cyclic peptide bond (last residue connects back to beginning)
            # This handles cases where N-terminal caps (like 'ac') are at position 1
            last_id = ordered_nodes[-1].id
            first_few_ids = [ordered_nodes[i].id for i in range(min(3, len(ordered_nodes)))]
            
            for link in graph.links:
                if link.linkage_type == LinkageType.PEPTIDE:
                    # Check if this is the cyclic bond (last to one of first few)
                    is_cyclic_bond = False
                    from_id, to_id = None, None
                    
                    if link.from_node_id == last_id and link.to_node_id in first_few_ids:
                        from_id, to_id = link.from_node_id, link.to_node_id
                        is_cyclic_bond = True
                    elif link.to_node_id == last_id and link.from_node_id in first_few_ids:
                        from_id, to_id = link.to_node_id, link.from_node_id
                        is_cyclic_bond = True
                    
                    if is_cyclic_bond:
                        # Find positions (1-indexed)
                        from_pos = next((i + 1 for i, n in enumerate(ordered_nodes) if n.id == from_id), None)
                        to_pos = next((i + 1 for i, n in enumerate(ordered_nodes) if n.id == to_id), None)
                        
                        if from_pos and to_pos:
                            connections.append(f"PEPTIDE1,PEPTIDE1,{from_pos}:R2-{to_pos}:R1")
                            break
        
        # Add disulfide bridges
        for link in graph.links:
            if link.linkage_type == LinkageType.DISULFIDE:
                # Get positions in ordered sequence (1-indexed)
                from_pos = None
                to_pos = None
                for i, node in enumerate(ordered_nodes):
                    if node.id == link.from_node_id:
                        from_pos = i + 1
                    if node.id == link.to_node_id:
                        to_pos = i + 1
                
                if from_pos and to_pos:
                    # Format: PEPTIDE1,PEPTIDE1,from_pos:R3-to_pos:R3
                    connections.append(f"PEPTIDE1,PEPTIDE1,{from_pos}:R3-{to_pos}:R3")
        
        # Generate final HELM notation
        if connections:
            connection_str = "|".join(connections)
            helm = f"PEPTIDE1{{{sequence}}}${connection_str}$$$V2.0"
        else:
            helm = f"PEPTIDE1{{{sequence}}}$$$$"
        
        return helm

    def generate_helm_notation(self, monomers) -> str:
        """
        Legacy method: Generate HELM notation from a list of monomers.
        Kept for backward compatibility.
        
        Args:
            monomers: List of MonomerData objects
        
        Returns:
            HELM notation string
        """
        if not monomers:
            return ""

        sequence = ".".join([monomer.symbol for monomer in monomers])
        helm = f"PEPTIDE1{{{sequence}}}$$$$"

        return helm

# ============================================================================
# Content from: monomer_library.py
# ============================================================================

from rdkit import Chem
from rdkit import RDLogger
from collections import defaultdict
from itertools import combinations
import json
import os

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.warning')

class MonomerData:
    def __init__(self):
        self.symbol = ""
        self.name = ""
        self.mol = None
        self.smiles = ""  # Original SMILES with R-groups
        self.r_groups = {}  # R-group label -> cap SMILES
        self.r_group_count = 0
        self.capped_smiles_cache = {}  # Cache: frozenset of removed R-groups -> canonical SMILES

    def __repr__(self):
        return f"Monomer({self.symbol}: {self.name}, R-groups: {self.r_group_count})"
    
    def get_capped_smiles_for_removed_rgroups(self, removed_rgroups: frozenset) -> str:
        """
        Get canonical SMILES with specific R-groups removed (lazy generation with caching).
        
        Args:
            removed_rgroups: frozenset of R-group labels that were removed (e.g., {'R1', 'R2'})
        
        Returns:
            Canonical SMILES with those R-groups removed, or empty string on error
            
        Example:
            For monomer with R1, R2:
            - get_capped_smiles_for_removed_rgroups({'R1'}) → SMILES with R1 removed, R2 kept
            - get_capped_smiles_for_removed_rgroups({'R2'}) → SMILES with R2 removed, R1 kept
            - get_capped_smiles_for_removed_rgroups({'R1', 'R2'}) → SMILES with both removed
        """
        # Check cache first
        if removed_rgroups in self.capped_smiles_cache:
            return self.capped_smiles_cache[removed_rgroups]
        
        # Generate on demand
        smiles = self._get_smiles_with_rgroups_removed(removed_rgroups)
        
        # Cache for future use
        self.capped_smiles_cache[removed_rgroups] = smiles
        
        return smiles
    
    def _get_smiles_with_rgroups_removed(self, removed_rgroups: frozenset) -> str:
        """
        Generate canonical SMILES with specific R-groups removed and others capped.
        
        Args:
            removed_rgroups: Set of R-group labels where bonds were broken (e.g., {'R1', 'R2'})
        
        Returns:
            Canonical SMILES string matching fragment structure
            
        Logic:
            - R-groups in removed_rgroups: Remove dummy atom (bond was broken)
            - R-groups NOT in removed_rgroups: Cap according to library (e.g., OH, H)
            - Final SMILES has NO [*:X] markers to match fragment SMILES
        """
        try:
            mol_copy = Chem.Mol(self.mol)
            
            # Identify which R-groups to cap vs remove
            kept_rgroups = set(self.r_groups.keys()) - removed_rgroups
            
            # Process each R-group
            # IMPORTANT: SMILES [*:1] uses atom map numbers, not isotopes!
            dummy_atoms_to_process = []
            for atom in mol_copy.GetAtoms():
                if atom.GetAtomicNum() == 0:  # Dummy atom (R-group)
                    map_num = atom.GetAtomMapNum()
                    if map_num > 0:
                        r_label = f"R{map_num}"
                        if r_label in removed_rgroups:
                            # Just remove this dummy atom
                            dummy_atoms_to_process.append((atom.GetIdx(), 'remove', r_label))
                        elif r_label in kept_rgroups:
                            # Need to cap this R-group
                            cap_smiles = self.r_groups.get(r_label, '')
                            dummy_atoms_to_process.append((atom.GetIdx(), 'cap', cap_smiles))
            
            # Apply caps to kept R-groups, remove others
            # Process in two passes: first cap, then remove
            # Cap R-groups: Replace [*:X] with the cap group (e.g., H or OH)
            for atom_idx, action, data in sorted(dummy_atoms_to_process, reverse=True):
                if action == 'cap':
                    cap_smiles = data
                    # For R1 cap '[*:1][H]', we just remove [*:1] (implicit H added)
                    # For R2 cap 'O[*:2]', we need to add O when removing [*:2]
                    # Simplified: check if cap has O
                    if 'O' in cap_smiles and '[*:' in cap_smiles:
                        # R2-like cap: need to add OH group
                        # Get the neighbor atom of the dummy
                        atom = mol_copy.GetAtomWithIdx(atom_idx)
                        neighbors = atom.GetNeighbors()
                        if neighbors:
                            neighbor = neighbors[0]
                            # Add OH to the neighbor before removing dummy
                            emol = Chem.EditableMol(mol_copy)
                            new_o_idx = emol.AddAtom(Chem.Atom(8))  # Oxygen
                            emol.AddBond(neighbor.GetIdx(), new_o_idx, Chem.BondType.SINGLE)
                            emol.RemoveAtom(atom_idx)
                            mol_copy = emol.GetMol()
                    else:
                        # R1-like cap: just remove dummy (implicit H)
                        emol = Chem.EditableMol(mol_copy)
                        emol.RemoveAtom(atom_idx)
                        mol_copy = emol.GetMol()
                elif action == 'remove':
                    # Just remove the dummy atom
                    emol = Chem.EditableMol(mol_copy)
                    emol.RemoveAtom(atom_idx)
                    mol_copy = emol.GetMol()
            
            if mol_copy:
                # Sanitize to add implicit hydrogens where needed
                Chem.SanitizeMol(mol_copy)
                # Generate canonical SMILES without any R-group markers
                return Chem.MolToSmiles(mol_copy, canonical=True)
            return ""
        except Exception as e:
            return ""


class MonomerLibrary:
    def __init__(self):
        self.monomers = {}
        self.smiles_to_monomer = {}
        self.name_to_monomer = {}
        self.symbol_to_monomer = {}

    def load_from_helm_json(self, json_path: str) -> None:
        if not os.path.exists(json_path):
            return

        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
        except Exception:
            return

        successful = 0
        for monomer_dict in data:
            try:
                monomer = self._parse_monomer(monomer_dict)
                if monomer and monomer.mol is not None:
                    self.monomers[monomer.symbol] = monomer
                    self.symbol_to_monomer[monomer.symbol] = monomer

                    clean_name = monomer.name.lower().replace(" ", "").replace("-", "").replace("_", "")
                    self.name_to_monomer[clean_name] = monomer

                    successful += 1
            except Exception:
                continue

    def _parse_monomer(self, monomer_dict: dict):
        # IMPORTANT: Only load PEPTIDE monomers (amino acids)
        # The library contains RNA, CHEM, etc. with overlapping symbols (A, C, G, T, U)
        polymer_type = monomer_dict.get('polymerType', '')
        if polymer_type != 'PEPTIDE':
            return None
        
        monomer = MonomerData()
        monomer.symbol = monomer_dict.get('symbol', '')
        monomer.name = monomer_dict.get('name', '')

        if not monomer.symbol:
            return None

        smiles = monomer_dict.get('smiles', '')
        molfile = monomer_dict.get('molfile', '')

        if smiles:
            try:
                monomer.mol = Chem.MolFromSmiles(smiles)
                monomer.smiles = smiles
            except Exception:
                monomer.mol = None

        if monomer.mol is None and molfile:
            try:
                monomer.mol = Chem.MolFromMolBlock(molfile)
                if monomer.mol:
                    monomer.smiles = Chem.MolToSmiles(monomer.mol)
            except Exception:
                monomer.mol = None

        if monomer.mol is None:
            return None
        
        # Parse R-groups
        rgroups_list = monomer_dict.get('rgroups', [])
        for rgroup in rgroups_list:
            label = rgroup.get('label', '')
            cap_smiles = rgroup.get('capGroupSMILES', '')
            if label and cap_smiles:
                monomer.r_groups[label] = cap_smiles
        
        monomer.r_group_count = len(monomer.r_groups)

        return monomer

    def find_monomer_by_fragment_smiles(self, fragment_smiles: str, num_connections: int):
        """
        Find monomer by matching fragment SMILES with on-demand R-group removal.
        
        Args:
            fragment_smiles: Canonical SMILES of the fragment  
            num_connections: Number of connections this fragment has in the graph
        
        Returns:
            MonomerData if match found, None otherwise
            
        Logic:
            - Fragment with N connections → N R-groups were removed during fragmentation
            - For monomer with M R-groups, try all C(M,N) combinations of which N R-groups were removed
            - Generate SMILES for each combination on-demand (with caching)
            
        Example:
            Fragment has 1 connection, monomer has R1, R2:
            - Try removing R1 → check if SMILES matches
            - Try removing R2 → check if SMILES matches
        """
        # Search through all monomers
        for symbol, monomer in self.monomers.items():
            # Skip if monomer doesn't have enough R-groups
            if monomer.r_group_count < num_connections:
                continue
            
            # Generate all combinations of num_connections R-groups that could have been removed
            r_group_labels = list(monomer.r_groups.keys())
            
            # For each combination of R-groups that could have been removed
            for removed_combo in combinations(r_group_labels, num_connections):
                removed_set = frozenset(removed_combo)
                
                # Generate SMILES with these R-groups removed (lazy, cached)
                candidate_smiles = monomer.get_capped_smiles_for_removed_rgroups(removed_set)
                
                # Check if it matches the fragment
                if candidate_smiles == fragment_smiles:
                    return monomer
        
        return None

    def find_monomer_by_symbol(self, symbol: str):
        return self.symbol_to_monomer.get(symbol)

# ============================================================================
# Content from: monomer_matcher.py
# ============================================================================

from rdkit import Chem


class MonomerMatcher:
    """
    Matches molecular fragments to monomers using graph-aware R-group analysis.
    
    Revolutionary approach:
    - No hardcoded mappings
    - No complex normalization
    - Direct string comparison of canonical SMILES
    - Graph topology determines which R-groups are capped
    """
    
    def __init__(self, monomer_library: MonomerLibrary):
        self.monomer_library = monomer_library

    def find_exact_match(self, fragment: Chem.Mol, num_connections: int = 0):
        """
        Find exact match for a fragment based on graph topology.
        
        Args:
            fragment: RDKit molecule object representing a fragment
            num_connections: Number of connections this fragment has in the graph
        
        Returns:
            MonomerData object if match found, None otherwise
        """
        try:
            # Get canonical SMILES of the fragment
            frag_smiles = Chem.MolToSmiles(fragment, canonical=True)
            if not frag_smiles:
                return None

            # Use the library's new graph-aware matching
            match = self.monomer_library.find_monomer_by_fragment_smiles(
                frag_smiles, num_connections
            )
            
            return match

        except Exception:
            return None
    
    def match_graph(self, graph: FragmentGraph):
        """
        Match all fragments in a graph to monomers.
        
        Args:
            graph: FragmentGraph with unmatched nodes
        
        Returns:
            Number of successfully matched nodes
        """
        matched_count = 0
        
        for node_id, node in graph.nodes.items():
            # Count connections for this node
            neighbors = graph.get_neighbors(node_id)
            num_connections = len(neighbors)
            
            # Find matching monomer
            monomer = self.find_exact_match(node.mol, num_connections)
            
            if monomer:
                node.monomer = monomer
                matched_count += 1
        
        return matched_count

# ============================================================================
# Content from: pipeline.py
# ============================================================================

from rdkit import Chem
import os
import json

# Global variables for caching
_MONOMER_LIBRARY = None
_PROCESSOR = None
_MATCHER = None
_HELM_GENERATOR = None


def _load_monomer_library():
    global _MONOMER_LIBRARY
    if _MONOMER_LIBRARY is None:
        # Define path to library relative to current directory
        current_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(current_dir)
        library_path = os.path.join(project_root, "libraries", "HELMCoreLibrary.json")

        if not os.path.exists(library_path):
            return None

        print("Loading monomer library...")
        _MONOMER_LIBRARY = MonomerLibrary()
        _MONOMER_LIBRARY.load_from_helm_json(library_path)

        if not _MONOMER_LIBRARY.monomers:
            return None
        
        print(f"Monomer library loaded: {len(_MONOMER_LIBRARY.monomers)} monomers")

    return _MONOMER_LIBRARY


def _get_processors():
    """
    Get or create singleton instances of processors.
    Returns tuple: (processor, matcher, helm_generator)
    """
    global _PROCESSOR, _MATCHER, _HELM_GENERATOR
    
    if _PROCESSOR is None or _MATCHER is None or _HELM_GENERATOR is None:
        library = _load_monomer_library()
        if not library:
            return None, None, None
        
        _PROCESSOR = FragmentProcessor(library)
        _MATCHER = MonomerMatcher(library)
        _HELM_GENERATOR = HELMGenerator()
    
    return _PROCESSOR, _MATCHER, _HELM_GENERATOR


def preload_library():
    """
    Preload the monomer library and initialize processors once at the start.
    Returns True if successful, False otherwise.
    """
    library = _load_monomer_library()
    if library is None:
        return False
    
    # Initialize processors
    processor, matcher, generator = _get_processors()
    return processor is not None


def convert_molecules_batch(molfiles: list, library_json: str = None) -> list:
    """
    Convert a batch of molecules from molfile format to HELM notation.
    
    Args:
        molfiles: List of molfile strings
        library_json: Optional monomer library as JSON string.
                     If None, uses default cached library from HELMCoreLibrary.json
    
    Returns:
        List of tuples: (success: bool, helm_notation: str)
        success is True if molecule was successfully converted, False otherwise
    """
    # Determine which library to use
    if library_json is None:
        # Use cached global library
        global _PROCESSOR
        if _PROCESSOR is None:
            print("Initializing monomer library and processors...")
            if not preload_library():
                print("ERROR: Failed to load monomer library")
                return [(False, "Library initialization failed") for _ in molfiles]
            print()
        
        # Use shared processor instances
        processor, matcher, helm_generator = _get_processors()
        if not processor:
            return [(False, "") for _ in molfiles]
    else:
        # Load custom library from provided JSON string (no caching)
        try:
            library_data = json.loads(library_json)
        except Exception as e:
            print(f"ERROR: Failed to parse library JSON: {str(e)}")
            return [(False, "Invalid JSON") for _ in molfiles]
        
        print(f"Loading custom library from JSON string...")
        library = MonomerLibrary()
        
        # Parse the library data
        successful = 0
        for monomer_dict in library_data:
            try:
                monomer = library._parse_monomer(monomer_dict)
                if monomer and monomer.mol is not None:
                    library.monomers[monomer.symbol] = monomer
                    library.symbol_to_monomer[monomer.symbol] = monomer
                    clean_name = monomer.name.lower().replace(" ", "").replace("-", "").replace("_", "")
                    library.name_to_monomer[clean_name] = monomer
                    successful += 1
            except Exception:
                continue
        
        if not library.monomers:
            print("ERROR: No monomers loaded from custom library")
            return [(False, "Library loading failed") for _ in molfiles]
        
        print(f"Custom library loaded: {len(library.monomers)} monomers")
        
        # Create processor instances for this library
        processor = FragmentProcessor(library)
        matcher = MonomerMatcher(library)
        helm_generator = HELMGenerator()
    
    results = []
    
    for i in range(len(molfiles)):
        molfile = molfiles[i]
        mol = Chem.MolFromMolBlock(molfile)
        if not mol:
            results.append((False, ""))
            continue
        
        try:
            # Process molecule into fragment graph
            graph = processor.process_molecule(mol)
            
            # Match each fragment to a monomer using graph topology
            unknown_count = 0
            for node_id, node in graph.nodes.items():
                # Count connections for this node
                neighbors = graph.get_neighbors(node_id)
                num_connections = len(neighbors)
                
                # Find matching monomer
                monomer = matcher.find_exact_match(node.mol, num_connections)
                if monomer:
                    node.monomer = monomer
                else:
                    unknown_count += 1
                    mock_monomer = MonomerData()
                    mock_monomer.symbol = f"X{unknown_count}"
                    mock_monomer.name = f"Unknown_{unknown_count}"
                    node.monomer = mock_monomer
            
            # Try to recover unmatched fragments by merging with neighbors
            max_recovery_attempts = 3  # Prevent infinite loops
            for attempt in range(max_recovery_attempts):
                had_changes = processor.recover_unmatched_fragments(graph, matcher)
                if not had_changes:
                    break
            
            if len(graph.nodes) > 0:
                helm_notation = helm_generator.generate_helm_from_graph(graph)
                results.append((True, helm_notation))
            else:
                results.append((False, ""))
        except Exception as e:
            results.append((False, f"Error: {str(e)}"))
    
    return results

res_helm_list = convert_molecules_batch(molListToProcess, library_json=libraryJSON)
result_helm = pd.DataFrame(map(lambda x: x[1], res_helm_list), columns=["regenerated sequences"])