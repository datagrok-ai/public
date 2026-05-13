#language: python
#name: molToHelmConverterPy
#description: Converts molecules to HELM notation based on monomer library
#input: dataframe moleculesDataframe
#input: column moleculesColumn {semType: Molecule}
#input: file libraryFile
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
import re

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
        """Helper for depth-first traversal with bidirectional link support"""
        if node_id in visited:
            return

        visited.add(node_id)
        ordered.append(self.nodes[node_id])

        # Follow links in BOTH directions but prefer the canonical (from→to)
        # direction. Link direction depends on bond detection order and is not
        # guaranteed to match backbone direction (e.g. FC01 stapled peptides).
        peptide_fwd = []
        peptide_bwd = []
        other_fwd = []
        other_bwd = []

        for link in self.links:
            if link.from_node_id == node_id and link.to_node_id not in visited:
                if link.linkage_type == LinkageType.PEPTIDE:
                    peptide_fwd.append(link.to_node_id)
                else:
                    other_fwd.append(link.to_node_id)
            elif link.to_node_id == node_id and link.from_node_id not in visited:
                if link.linkage_type == LinkageType.PEPTIDE:
                    peptide_bwd.append(link.from_node_id)
                else:
                    other_bwd.append(link.from_node_id)

        # Forward first, backward as fallback
        for neighbor_id in peptide_fwd + peptide_bwd + other_fwd + other_bwd:
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
        
        # Check if any of the last few residues connect back to any of the first few.
        # Checking multiple positions on each end handles branch nodes (like 'ac')
        # that the bidirectional traversal may place at the edges.
        first_few_ids = set(ordered[i].id for i in range(min(3, len(ordered))))
        last_few_ids = set(ordered[-i - 1].id for i in range(min(3, len(ordered))))

        for link in self.links:
            if link.linkage_type == LinkageType.PEPTIDE:
                if (link.from_node_id in last_few_ids and link.to_node_id in first_few_ids) or \
                   (link.to_node_id in last_few_ids and link.from_node_id in first_few_ids):
                    return True

        return False
    
    def find_all_cycles(self) -> List[List[int]]:
        """
        Find all cycles in the graph using DFS.
        Returns list of cycles, where each cycle is a list of node IDs.
        """
        cycles = []
        visited = set()
        rec_stack = set()
        parent = {}
        
        def dfs(node_id: int, path: List[int]):
            visited.add(node_id)
            rec_stack.add(node_id)
            path.append(node_id)
            
            # Get peptide bond neighbors
            neighbors = [n_id for n_id, link_type in self.get_neighbors(node_id) 
                        if link_type == LinkageType.PEPTIDE]
            
            for neighbor_id in neighbors:
                if neighbor_id not in visited:
                    parent[neighbor_id] = node_id
                    dfs(neighbor_id, path[:])
                elif neighbor_id in rec_stack and neighbor_id != parent.get(node_id):
                    # Found a cycle - extract it from path
                    cycle_start_idx = path.index(neighbor_id)
                    cycle = path[cycle_start_idx:] + [neighbor_id]
                    # Normalize cycle (start from smallest ID)
                    min_idx = cycle.index(min(cycle[:-1]))  # Don't include duplicate last element
                    normalized = cycle[min_idx:-1] + cycle[:min_idx]
                    if normalized not in cycles:
                        cycles.append(normalized)
            
            rec_stack.remove(node_id)
        
        # Try starting DFS from each unvisited node
        for node_id in self.nodes.keys():
            if node_id not in visited:
                parent[node_id] = None
                dfs(node_id, [])
        
        return cycles
    
    def get_connected_components(self) -> List[List[int]]:
        """
        Find all connected components in the graph.
        Returns list of components, where each component is a list of node IDs.
        """
        visited = set()
        components = []
        
        def dfs_component(node_id: int, component: List[int]):
            visited.add(node_id)
            component.append(node_id)
            neighbors = self.get_neighbors(node_id)
            for neighbor_id, _ in neighbors:
                if neighbor_id not in visited:
                    dfs_component(neighbor_id, component)
        
        for node_id in self.nodes.keys():
            if node_id not in visited:
                component = []
                dfs_component(node_id, component)
                components.append(sorted(component))
        
        return components
    
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
        # Alpha carbon after N can be sp3 (X4) or sp2 (X3) for dehydroamino acids, or aromatic (#6 includes both)
        self.peptide_bond = Chem.MolFromSmarts('[#6]-[C;X3;!r5;!r6](=[O;X1])-[N;X2,X3]~[#6;X3,X4]')
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

            # Filter out internal amide bonds in CHEM linkers like FC01.
            # FC01 pattern: C(=O)-N-ArRing-N-C(=O) — two amide bonds connect to the
            # same aromatic ring via the alpha-C position (match[4]).
            # Real aromatic amino acids (3Abz) only have ONE such bond per ring.
            skip_indices = set()
            ring_info = mol.GetRingInfo()
            rings = ring_info.AtomRings()

            # Map: ring_frozenset -> list of match indices where match[4] is aromatic on that ring
            # Only consider small rings (5-6 atoms) — large macrocycles should not be filtered
            ring_to_matches = {}
            for i, match in enumerate(matches):
                if len(match) < 5:
                    continue
                alpha_c_atom = mol.GetAtomWithIdx(match[4])
                if alpha_c_atom.GetIsAromatic():
                    for ring in rings:
                        if match[4] in ring and len(ring) <= 6:
                            ring_key = frozenset(ring)
                            if ring_key not in ring_to_matches:
                                ring_to_matches[ring_key] = []
                            ring_to_matches[ring_key].append(i)
                            break

            # If 2+ matches share an aromatic ring at their alpha-C position, skip them
            for ring_key, match_indices in ring_to_matches.items():
                if len(match_indices) >= 2:
                    skip_indices.update(match_indices)

            for i, match in enumerate(matches):
                if len(match) >= 5 and i not in skip_indices:
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


    def _find_staple_sidechain_bonds(self, mol, existing_bonds):
        """
        Find non-backbone bonds to cleave in macrocycles.

        Handles three types of macrocyclic cross-links:
        1. RCMtrans/RCMcis (stapled peptides): C=C double bond in the linker.
           Cleaves one hop away on each side to keep the correct R3 chain length.
        2. FC01-type (thioether staples): C-S bonds in the linker.
        3. Alkyl cross-links (bi-cyclic peptides): pure C-C chains connecting
           two amino acid side chains (R3-R3). Detected by finding non-backbone
           segments in large macrocycles and cleaving at their midpoint.
        """
        ring_info = mol.GetRingInfo()
        large_rings = [set(ring) for ring in ring_info.AtomRings() if len(ring) > 10]
        if not large_rings:
            return []

        existing_atom_pairs = set()
        for a1, a2, _ in existing_bonds:
            existing_atom_pairs.add((min(a1, a2), max(a1, a2)))
        # Also track existing bond atoms for backbone detection
        existing_bond_atoms = set()
        for a1, a2, _ in existing_bonds:
            existing_bond_atoms.add(a1)
            existing_bond_atoms.add(a2)

        additional_bonds = []
        seen = set()

        for ring in large_rings:
            ring_list = list(ring)

            # --- Type 1: C=C double bonds (RCMtrans/RCMcis) ---
            # Only if molecule has quaternary alpha-methyl C (staple monomer signature)
            quat_alpha = Chem.MolFromSmarts('[N;X2,X3]-[C;X4;H0](-[C;X3](=[O;X1]))-[CH3]')
            has_quat = quat_alpha and mol.HasSubstructMatch(quat_alpha)
            if has_quat:
                for bond in mol.GetBonds():
                    a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                    if a1 not in ring or a2 not in ring:
                        continue
                    if (bond.GetBondTypeAsDouble() >= 2 and
                            mol.GetAtomWithIdx(a1).GetAtomicNum() == 6 and
                            mol.GetAtomWithIdx(a2).GetAtomicNum() == 6):
                        for cc_atom_idx in (a1, a2):
                            other_cc = a2 if cc_atom_idx == a1 else a1
                            cc_atom = mol.GetAtomWithIdx(cc_atom_idx)
                            for nbr in cc_atom.GetNeighbors():
                                nbr_idx = nbr.GetIdx()
                                if nbr_idx == other_cc or nbr_idx not in ring:
                                    continue
                                for nbr2 in nbr.GetNeighbors():
                                    nbr2_idx = nbr2.GetIdx()
                                    if nbr2_idx == cc_atom_idx or nbr2_idx not in ring:
                                        continue
                                    pair = (min(nbr_idx, nbr2_idx), max(nbr_idx, nbr2_idx))
                                    if pair not in existing_atom_pairs and pair not in seen:
                                        seen.add(pair)
                                        additional_bonds.append((nbr_idx, nbr2_idx, LinkageType.UNKNOWN))

            # --- Type 2: C-S thioether bonds (FC01) ---
            # Only true thioethers (S bonded to C on both sides), NOT disulfide-adjacent
            for bond in mol.GetBonds():
                a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if a1 not in ring or a2 not in ring:
                    continue
                at1, at2 = mol.GetAtomWithIdx(a1), mol.GetAtomWithIdx(a2)
                if ((at1.GetAtomicNum() == 6 and at2.GetAtomicNum() == 16) or
                        (at1.GetAtomicNum() == 16 and at2.GetAtomicNum() == 6)):
                    s_atom = at2 if at2.GetAtomicNum() == 16 else at1
                    # Skip if S is bonded to another S (disulfide bridge path)
                    if any(n.GetAtomicNum() == 16 for n in s_atom.GetNeighbors()):
                        continue
                    pair = (min(a1, a2), max(a1, a2))
                    if pair not in existing_atom_pairs and pair not in seen:
                        seen.add(pair)
                        additional_bonds.append((a1, a2, LinkageType.UNKNOWN))

        # --- Type 3: Alkyl cross-link paths (bi-cyclic R3-R3) ---
        # Find pairs of alpha-C atoms connected by pure carbon chains (no N/O/S
        # in the path). These are R3-R3 cross-links between different cycles.
        # Cleave at the midpoint of each such chain.
        alpha_c_pat = Chem.MolFromSmarts('[N]-[C;X4]-[C;X3](=[O])')
        if alpha_c_pat:
            ac_matches = mol.GetSubstructMatches(alpha_c_pat)
            alpha_c_set = list(set(m[1] for m in ac_matches))
            for i, ac1 in enumerate(alpha_c_set):
                for ac2 in alpha_c_set[i + 1:]:
                    path = Chem.GetShortestPath(mol, ac1, ac2)
                    if not path or len(path) < 4 or len(path) > 12:
                        continue
                    # All middle atoms must be C with no N/O/S neighbors
                    middle_ok = True
                    for mid_idx in path[1:-1]:
                        atom = mol.GetAtomWithIdx(mid_idx)
                        if atom.GetAtomicNum() != 6:
                            middle_ok = False
                            break
                        if any(n.GetAtomicNum() in (7, 8, 16) for n in atom.GetNeighbors()):
                            middle_ok = False
                            break
                    if middle_ok:
                        mid = len(path) // 2
                        pair = (min(path[mid - 1], path[mid]), max(path[mid - 1], path[mid]))
                        if pair not in existing_atom_pairs and pair not in seen:
                            seen.add(pair)
                            additional_bonds.append((path[mid - 1], path[mid], LinkageType.UNKNOWN))

        return additional_bonds

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

            # Detect R3 side-chain bonds for staple monomers (R8, S5, etc.)
            r3_bonds = self._find_staple_sidechain_bonds(mol, bonds_to_cleave)
            if r3_bonds:
                bonds_to_cleave.extend(r3_bonds)

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
            graph.uncleaned_fragments = fragments  # Keep fragments with dummy atoms for R-group SMILES

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
            if node.monomer and node.monomer.is_unknown:
                unmatched_nodes.append(node_id)
        
        if not unmatched_nodes:
            return False
        
        had_changes = False
        
        # Try to recover each unmatched node
        for node_id in unmatched_nodes:
            # Check if node still exists (might have been merged already)
            if node_id not in graph.nodes:
                continue
            
            # Get neighbors from graph links (returns list of (neighbor_id, linkage_type))
            neighbors = graph.get_neighbors(node_id)
            
            if not neighbors:
                continue
            
            # Try merging with each individual neighbor first
            for neighbor_id, linkage_type in neighbors:
                if neighbor_id not in graph.nodes:
                    continue
                    
                nodes_to_merge = sorted([node_id, neighbor_id])
                
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
                
                # Try to match the combined fragment (exact match only)
                monomer = matcher.find_exact_match(combined_mol, num_connections)

                if monomer:
                    # Success! Create new merged node
                    new_node_id = min(nodes_to_merge)
                    new_node = FragmentNode(new_node_id, combined_mol)
                    new_node.monomer = monomer
                    
                    # Merge nodes in graph
                    self._merge_nodes_in_graph(graph, nodes_to_merge, new_node)
                    
                    had_changes = True
                    break  # Stop trying other neighbors for this node
            
            if had_changes:
                break  # Restart from beginning after a successful merge
        
        return had_changes
    
    def recover_unmatched_with_stereo_agnostic(self, graph: FragmentGraph, matcher) -> int:
        """
        Separate recovery procedure: Try to match remaining unmatched fragments 
        using stereochemistry-agnostic comparison.
        
        This handles poor quality input data where stereochemistry is not assigned.
        Only called after regular recovery attempts have finished.
        
        Args:
            graph: FragmentGraph with some unmatched nodes
            matcher: MonomerMatcher instance
        
        Returns:
            Number of fragments that were successfully matched
        """
        from rdkit import Chem
        
        # Find all unmatched nodes (nodes with mock/unknown monomers)
        unmatched_nodes = []
        for node_id, node in graph.nodes.items():
            if node.monomer and node.monomer.is_unknown:
                unmatched_nodes.append(node_id)
        
        if not unmatched_nodes:
            return 0
        
        print(f"DEBUG: Attempting stereo-agnostic recovery for {len(unmatched_nodes)} unmatched nodes")
        
        matched_count = 0
        
        for node_id in unmatched_nodes:
            if node_id not in graph.nodes:
                continue
            
            node = graph.nodes[node_id]
            
            # Get fragment SMILES
            fragment_smiles = Chem.MolToSmiles(node.mol, canonical=True)
            
            # Count connections
            neighbors = graph.get_neighbors(node_id)
            num_connections = len(neighbors)
            
            # Try stereo-agnostic matching
            monomer = matcher.monomer_library.find_monomer_by_fragment_smiles_no_stereo(
                fragment_smiles, num_connections
            )
            
            if monomer:
                print(f"DEBUG: Stereo-agnostic match for node {node_id}: {monomer.symbol}")
                node.monomer = monomer
                matched_count += 1
            else:
                print(f"DEBUG: No stereo-agnostic match for node {node_id}")
        
        return matched_count

    def recover_unmatched_by_merging_stereo_agnostic(self, graph: FragmentGraph, matcher) -> bool:
        """
        Final recovery pass: merge pairs of BOTH-unmatched neighbor fragments and
        try stereo-agnostic matching on the combined result.

        This handles monomers like Phe_4Sdihydroorotamido that have internal amide
        bonds which get incorrectly cleaved, producing two unmatched fragments.

        Only merges when BOTH fragments in a pair are unmatched — never touches
        already-matched nodes to avoid regressions.

        Returns True if any merges were successful.
        """
        def _is_unmatched(node):
            return (node.monomer and
                    node.monomer.is_unknown)

        unmatched_ids = [nid for nid, node in graph.nodes.items() if _is_unmatched(node)]
        if not unmatched_ids:
            return False

        had_changes = False

        for node_id in unmatched_ids:
            if node_id not in graph.nodes:
                continue
            if not _is_unmatched(graph.nodes[node_id]):
                continue

            neighbors = graph.get_neighbors(node_id)
            for neighbor_id, linkage_type in neighbors:
                if neighbor_id not in graph.nodes:
                    continue
                # Only merge with another unmatched neighbor
                if not _is_unmatched(graph.nodes[neighbor_id]):
                    continue

                nodes_to_merge = sorted([node_id, neighbor_id])

                # Find internal links between the merge candidates
                links_to_exclude = []
                for link in graph.links:
                    if (link.from_node_id in nodes_to_merge and
                            link.to_node_id in nodes_to_merge):
                        links_to_exclude.append(link)

                combined_mol = self._reconstruct_fragment_with_links(
                    nodes_to_merge, graph, links_to_exclude)
                if not combined_mol:
                    continue

                # Count external connections for merged fragment
                all_neighbors = set()
                for nid in nodes_to_merge:
                    if nid in graph.nodes:
                        for nbr_id, _ in graph.get_neighbors(nid):
                            if nbr_id not in nodes_to_merge:
                                all_neighbors.add(nbr_id)
                num_connections = len(all_neighbors)

                # Try exact match first, then stereo-agnostic
                monomer = matcher.find_exact_match(combined_mol, num_connections)
                if not monomer:
                    combined_smiles = Chem.MolToSmiles(combined_mol, canonical=True)
                    monomer = matcher.monomer_library.find_monomer_by_fragment_smiles_no_stereo(
                        combined_smiles, num_connections)

                if monomer:
                    new_node_id = min(nodes_to_merge)
                    new_node = FragmentNode(new_node_id, combined_mol)
                    new_node.monomer = monomer
                    self._merge_nodes_in_graph(graph, nodes_to_merge, new_node)
                    had_changes = True
                    break  # Restart from outer loop

            if had_changes:
                break

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
    - Multi-chain structures (BILN peptides)
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
        
        Supports multi-chain structures (BILN peptides):
        - Detects all cycles (rings) using SSSR-like algorithm
        - Each cycle becomes a separate PEPTIDE chain
        - R1-R2 connections define backbone within each chain
        - R3 connections link chains together
        
        Args:
            graph: FragmentGraph containing matched monomers and their connections
        
        Returns:
            HELM notation string
        """
        if len(graph) == 0:
            return ""
        
        # Find all cycles in the graph (each cycle will be a separate PEPTIDE chain)
        cycles = graph.find_all_cycles()
        
        # Decision: Use multi-chain HELM only if:
        # 1. Multiple cycles exist (BILN-style structure), OR
        # 2. There are standalone nodes not in any cycle (attached fragments)
        
        if not cycles:
            # No cycles - simple linear peptide
            return self._generate_simple_helm(graph)
        
        if len(cycles) == 1:
            # Single cycle (with or without standalone branch nodes like 'ac')
            # _generate_simple_helm handles branches correctly with proper R-group detection
            return self._generate_simple_helm(graph)

        # Multi-chain structure detected (multiple cycles)
        return self._generate_multi_chain_helm(graph, cycles)
    
    def _generate_simple_helm(self, graph: FragmentGraph) -> str:
        """
        Generate HELM for simple linear or single-cycle peptides.
        This is the original implementation for backward compatibility.
        """
        # Get ordered sequence of monomers (backbone)
        ordered_nodes_raw = graph.get_ordered_nodes()
        
        # Check if cyclic
        is_cyclic = graph.is_cyclic()
        
        # Filter backbone: nodes that are part of R1-R2 chain are backbone
        # Nodes lacking R1 (like 'ac' acetyl cap) are branches regardless of position
        backbone_nodes = []
        for node in ordered_nodes_raw:
            is_branch = False
            if node.monomer and len(ordered_nodes_raw) > 1:
                has_r1 = 'R1' in node.monomer.r_groups
                if not has_r1 and not node.monomer.is_unknown:
                    is_branch = True
            if not is_branch:
                backbone_nodes.append(node)
        
        ordered_nodes = backbone_nodes
        sequence_symbols = [node.monomer.symbol if node.monomer else "X" for node in ordered_nodes]
        
        # Detect branch nodes (nodes not in backbone)
        ordered_node_ids = {node.id for node in ordered_nodes}
        branch_nodes = [(node_id, node) for node_id, node in graph.nodes.items() 
                       if node_id not in ordered_node_ids]
        
        # Generate sequence notation — always bracket multi-char symbols (HELM spec requirement,
        # also needed for inline SMILES like [*:1]NC(CC(=O)O)C(=O)[*:2])
        formatted_symbols = [f"[{symbol}]" if len(symbol) > 1 else symbol for symbol in sequence_symbols]
        sequence = ".".join(formatted_symbols)
        
        # Collect non-sequential connections (disulfide bridges, cyclic bonds, etc.)
        connections = []
        
        if is_cyclic:
            # Find the actual cyclic peptide bond (last residue connects back to beginning)
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
        
        # Handle branch nodes (side chain modifications)
        # Create separate PEPTIDE chains for each branch
        branch_chains = []
        if branch_nodes:
            for branch_idx, (branch_node_id, branch_node) in enumerate(branch_nodes, start=2):
                branch_chain_name = f"PEPTIDE{branch_idx}"
                branch_symbol = branch_node.monomer.symbol if branch_node.monomer else f"X{branch_node_id}"
                
                # Format branch chain (single monomer)
                # In cyclic peptides, always use brackets for consistency with reference HELM
                if is_cyclic:
                    branch_chains.append(f"{branch_chain_name}{{[{branch_symbol}]}}")
                else:
                    branch_chains.append(f"{branch_chain_name}{{{branch_symbol}}}")
                
                # Find which backbone node this branch connects to
                for link in graph.links:
                    backbone_node_id = None
                    if link.from_node_id == branch_node_id and link.to_node_id in ordered_node_ids:
                        backbone_node_id = link.to_node_id
                    elif link.to_node_id == branch_node_id and link.from_node_id in ordered_node_ids:
                        backbone_node_id = link.from_node_id
                    
                    if backbone_node_id is not None:
                        # Find position of backbone node (1-indexed)
                        backbone_pos = next((i + 1 for i, n in enumerate(ordered_nodes) if n.id == backbone_node_id), None)
                        if backbone_pos:
                            # Determine which R-group the branch uses
                            branch_r_group = "R1"
                            if branch_node.monomer:
                                if 'R1' in branch_node.monomer.r_groups:
                                    branch_r_group = "R1"
                                elif 'R2' in branch_node.monomer.r_groups:
                                    branch_r_group = "R2"
                            
                            # Connection: backbone position R3 (side chain) to branch position 1 R-group
                            connections.append(f"PEPTIDE1,{branch_chain_name},{backbone_pos}:R3-1:{branch_r_group}")
                            break
        
        # Generate final HELM notation
        all_chains = [f"PEPTIDE1{{{sequence}}}"] + branch_chains
        helm_chains = "|".join(all_chains)
        
        if connections:
            connection_str = "|".join(connections)
            helm = f"{helm_chains}${connection_str}$$$V2.0"
        else:
            helm = f"{helm_chains}$$$$V2.0"
        
        return helm
    
    def _generate_multi_chain_helm(self, graph: FragmentGraph, cycles: list) -> str:
        """
        Generate HELM for multi-chain structures (BILN peptides).
        
        Strategy:
        1. Each cycle becomes a separate PEPTIDE chain
        2. Nodes not in cycles become additional chains
        3. R3 connections between chains are added as cross-links
        """
        # Identify which nodes belong to which cycles
        nodes_in_cycles = set()
        for cycle in cycles:
            nodes_in_cycles.update(cycle)
        
        # Find standalone nodes (not in any cycle)
        standalone_nodes = [nid for nid in graph.nodes.keys() if nid not in nodes_in_cycles]
        
        # Build PEPTIDE chains
        chains = []
        chain_node_map = {}  # Maps node_id -> (chain_idx, position_in_chain)
        
        # Add cyclic chains
        for cycle_idx, cycle in enumerate(cycles, start=1):
            chain_name = f"PEPTIDE{cycle_idx}"
            # Create sequence from cycle nodes
            sequence_symbols = []
            for pos, node_id in enumerate(cycle):
                node = graph.nodes[node_id]
                symbol = node.monomer.symbol if node.monomer else f"X{node_id}"
                sequence_symbols.append(symbol)
                chain_node_map[node_id] = (cycle_idx, pos + 1)  # 1-indexed position
            
            # Format with brackets for multi-letter symbols
            formatted = [f"[{s}]" if len(s) > 1 else s for s in sequence_symbols]
            sequence = ".".join(formatted)
            chains.append(f"{chain_name}{{{sequence}}}")
        
        # Add standalone chains (linear fragments not in cycles)
        next_chain_idx = len(cycles) + 1
        for node_id in standalone_nodes:
            chain_name = f"PEPTIDE{next_chain_idx}"
            node = graph.nodes[node_id]
            symbol = node.monomer.symbol if node.monomer else f"X{node_id}"
            chains.append(f"{chain_name}{{{symbol}}}")
            chain_node_map[node_id] = (next_chain_idx, 1)
            next_chain_idx += 1
        
        # Build connections
        connections = []
        
        # Add cyclic connections (R1-R2 within each cycle)
        for cycle_idx, cycle in enumerate(cycles, start=1):
            if len(cycle) >= 3:
                # Connect last to first
                chain_name = f"PEPTIDE{cycle_idx}"
                last_pos = len(cycle)
                connections.append(f"{chain_name},{chain_name},{last_pos}:R2-1:R1")
        
        # Add inter-chain connections (R3 links) and disulfide bridges
        processed_links = set()
        for link in graph.links:
            link_key = tuple(sorted([link.from_node_id, link.to_node_id]))
            if link_key in processed_links:
                continue
            
            from_chain_info = chain_node_map.get(link.from_node_id)
            to_chain_info = chain_node_map.get(link.to_node_id)
            
            if not from_chain_info or not to_chain_info:
                continue
            
            from_chain, from_pos = from_chain_info
            to_chain, to_pos = to_chain_info
            
            # Skip intra-cycle backbone peptide bonds (already handled by R1-R2 connection)
            if from_chain == to_chain and link.linkage_type == LinkageType.PEPTIDE:
                # Check if this is a sequential bond within the cycle
                cycle = cycles[from_chain - 1] if from_chain <= len(cycles) else []
                # Sequential bonds: adjacent positions or last-to-first
                if abs(from_pos - to_pos) == 1 or (from_pos == 1 and to_pos == len(cycle)) or (to_pos == 1 and from_pos == len(cycle)):
                    processed_links.add(link_key)
                    continue
            
            # Add cross-chain connections or intra-chain disulfide bridges
            if link.linkage_type == LinkageType.DISULFIDE:
                # Disulfide uses R3 (side chain cysteine)
                r_group = "R3"
            elif link.linkage_type == LinkageType.PEPTIDE:
                # Cross-chain peptide bond (side chain R3 connection)
                r_group = "R3"
            else:
                r_group = "R3"
            
            from_chain_name = f"PEPTIDE{from_chain}"
            to_chain_name = f"PEPTIDE{to_chain}"
            connections.append(f"{from_chain_name},{to_chain_name},{from_pos}:{r_group}-{to_pos}:{r_group}")
            processed_links.add(link_key)
        
        # Generate final HELM
        helm_chains = "|".join(chains)
        if connections:
            connection_str = "|".join(connections)
            helm = f"{helm_chains}${connection_str}$$$V2.0"
        else:
            helm = f"{helm_chains}$$$$V2.0"
        
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
import re

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.warning')

def remove_stereochemistry_from_smiles(smiles: str) -> str:
    """
    Remove stereochemistry markers from SMILES string.
    Converts [C@@H], [C@H] to C, etc.

    This is used for matching when input molecules don't have stereochemistry defined.
    Only strips brackets from SMILES organic subset atoms (B,C,N,O,P,S,F,Cl,Br,I).
    Atoms like Se, Te, etc. must keep their brackets to remain valid SMILES.
    """
    if not smiles:
        return smiles

    # SMILES organic subset: atoms that can appear without brackets
    organic_subset = {'B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I'}

    # Remove @ symbols (stereochemistry markers)
    smiles_no_stereo = re.sub(r'(@+)', '', smiles)

    # Remove explicit H and brackets only for organic subset atoms
    # [C@@H] -> [CH] -> C, but [SeH] must stay as [SeH]
    def _simplify_bracket(match):
        atom = match.group(1)  # e.g. 'C', 'Se', 'N'
        has_h = match.group(2)  # 'H' or ''
        if atom in organic_subset:
            return atom  # Strip brackets (and H) for organic subset
        elif has_h:
            return f'[{atom}H]'  # Keep brackets and H for non-organic atoms
        else:
            return f'[{atom}]'  # Keep brackets for non-organic atoms

    smiles_no_stereo = re.sub(r'\[([A-Z][a-z]?)(H?)\]', _simplify_bracket, smiles_no_stereo)

    return smiles_no_stereo

class MonomerData:
    def __init__(self):
        self.symbol = ""
        self.name = ""
        self.mol = None
        self.smiles = ""  # Original SMILES with R-groups
        self.r_groups = {}  # R-group label -> cap SMILES
        self.r_group_count = 0
        self.capped_smiles_cache = {}  # Cache: frozenset of removed R-groups -> canonical SMILES
        self.is_unknown = False  # True for unmatched fragments with inline SMILES

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


def _canonicalize_no_stereo(smiles: str) -> str:
    """
    Remove stereochemistry and re-canonicalize through RDKit.
    This ensures consistent canonical SMILES regardless of how the molecule was constructed.
    String-only stereo removal can produce non-canonical SMILES.
    """
    no_stereo = remove_stereochemistry_from_smiles(smiles)
    mol = Chem.MolFromSmiles(no_stereo)
    if mol:
        return Chem.MolToSmiles(mol, canonical=True)
    return no_stereo  # Fallback to string version if parse fails


class MonomerLibrary:
    def __init__(self):
        self.monomers = {}
        self.smiles_to_monomer = {}
        self.name_to_monomer = {}
        self.symbol_to_monomer = {}
        # Hash indices for O(1) matching (built after loading)
        self._smiles_index = {}         # canonical_smiles -> MonomerData
        self._smiles_no_stereo_index = {}  # stereo-free_smiles -> MonomerData

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

        # Build hash indices for O(1) matching
        self._build_smiles_indices()

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

    def _build_smiles_indices(self):
        """
        Pre-compute all possible capped SMILES for every monomer and build
        hash indices for O(1) lookup. Called once after loading all monomers.

        For each monomer with M R-groups, generates capped SMILES for all
        possible R-group removal combinations (up to 2^M - 1 entries, typically 1-7).

        Deduplicates monomers with identical SMILES+R-groups to avoid redundant
        capping computations (important for large libraries with variants).
        """
        self._smiles_index = {}
        self._smiles_no_stereo_index = {}

        # Dedup: group monomers by (smiles, r_group_keys) to avoid recomputing
        # identical capped forms for monomers with the same structure
        seen_structures = {}  # (smiles, r_group_frozenset) -> list of capped entries

        for symbol, monomer in self.monomers.items():
            if monomer.r_group_count == 0:
                continue

            r_group_labels = list(monomer.r_groups.keys())
            struct_key = (monomer.smiles, frozenset(monomer.r_groups.items()))

            if struct_key in seen_structures:
                # Reuse cached capped SMILES from an identical monomer
                for capped_smiles, n_removed in seen_structures[struct_key]:
                    key = (capped_smiles, n_removed)
                    if key not in self._smiles_index:
                        self._smiles_index[key] = monomer
                    ns_canonical = _canonicalize_no_stereo(capped_smiles)
                    if ns_canonical:
                        ns_key = (ns_canonical, n_removed)
                        if ns_key not in self._smiles_no_stereo_index:
                            self._smiles_no_stereo_index[ns_key] = monomer
                continue

            # First time seeing this structure — compute capped SMILES
            cached_entries = []

            for n_removed in range(1, monomer.r_group_count + 1):
                for removed_combo in combinations(r_group_labels, n_removed):
                    removed_set = frozenset(removed_combo)
                    capped_smiles = monomer.get_capped_smiles_for_removed_rgroups(removed_set)

                    if not capped_smiles:
                        continue

                    cached_entries.append((capped_smiles, n_removed))

                    key = (capped_smiles, n_removed)
                    if key not in self._smiles_index:
                        self._smiles_index[key] = monomer

                    ns_canonical = _canonicalize_no_stereo(capped_smiles)
                    if ns_canonical:
                        ns_key = (ns_canonical, n_removed)
                        if ns_key not in self._smiles_no_stereo_index:
                            self._smiles_no_stereo_index[ns_key] = monomer

            seen_structures[struct_key] = cached_entries

    def find_monomer_by_fragment_smiles(self, fragment_smiles: str, num_connections: int):
        """
        Find monomer by matching fragment SMILES. O(1) hash lookup.
        """
        return self._smiles_index.get((fragment_smiles, num_connections))

    def find_monomer_by_fragment_smiles_no_stereo(self, fragment_smiles: str, num_connections: int):
        """
        Find monomer by matching fragment SMILES WITHOUT stereochemistry.
        Used in recovery for handling poor quality input data. O(1) hash lookup.
        """
        ns_canonical = _canonicalize_no_stereo(fragment_smiles)
        if not ns_canonical:
            return None

        return self._smiles_no_stereo_index.get((ns_canonical, num_connections))

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

def _generate_rgroup_smiles(graph, node_id):
    """
    Generate SMILES with R-group markers ([*:1], [*:2], ...) for an unmatched fragment.
    Uses the uncleaned fragment (with dummy atoms from FragmentOnBonds) stored in the graph.
    Falls back to plain SMILES from the cleaned mol if uncleaned data isn't available.
    """
    # Try to use uncleaned fragment with dummy atoms
    if hasattr(graph, 'uncleaned_fragments') and node_id < len(graph.uncleaned_fragments):
        uncleaned = graph.uncleaned_fragments[node_id]
        try:
            mol = Chem.RWMol(Chem.Mol(uncleaned))
            r_num = 1
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    atom.SetIsotope(0)
                    atom.SetAtomMapNum(r_num)
                    r_num += 1
            return Chem.MolToSmiles(mol)
        except Exception:
            pass

    # Fallback: plain SMILES from cleaned mol (no R-groups)
    node = graph.nodes.get(node_id)
    if node and node.mol:
        return Chem.MolToSmiles(node.mol, canonical=True)
    return "?"


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


def convert_molecules_batch(molecules: list, library_json: str = None, input_type: str = "auto") -> list:
    """
    Convert a batch of molecules to HELM notation.
    
    Args:
        molecules: List of molecule strings (molfiles or SMILES)
        library_json: Optional monomer library as JSON string.
                     If None, uses default cached library from HELMCoreLibrary.json
        input_type: Type of input molecules - "molfile", "smiles", or "auto" (default).
                   "auto" will attempt to detect the format automatically.
    
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
                return [(False, "Library initialization failed") for _ in molecules]
            print()
        
        # Use shared processor instances
        processor, matcher, helm_generator = _get_processors()
        if not processor:
            return [(False, "") for _ in molecules]
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
            return [(False, "Library loading failed") for _ in molecules]
        
        print(f"Custom library loaded: {len(library.monomers)} monomers")

        # Build hash indices for O(1) matching
        library._build_smiles_indices()

        # Create processor instances for this library
        processor = FragmentProcessor(library)
        matcher = MonomerMatcher(library)
        helm_generator = HELMGenerator()
    
    # Helper function to detect molecule format
    def _is_molfile(mol_string: str) -> bool:
        """Check if string is a molfile (starts with RDKit molfile markers or has multiple lines)"""
        if not mol_string:
            return False
        lines = mol_string.strip().split('\n')
        # Molfiles typically have multiple lines and specific format
        if len(lines) > 3:
            # Check for V2000 or V3000 molfile markers
            if 'V2000' in mol_string or 'V3000' in mol_string:
                return True
            # Check for typical molfile structure (counts line format)
            if len(lines) > 3:
                counts_line = lines[3] if len(lines) > 3 else ""
                # Molfile counts line has specific format with atom/bond counts
                if len(counts_line) >= 6 and counts_line[:6].replace(' ', '').isdigit():
                    return True
        return False
    
    results = []
    
    for i in range(len(molecules)):
        mol_string = molecules[i]
        
        # Determine input type and parse molecule
        if input_type == "auto":
            # Auto-detect format
            if _is_molfile(mol_string):
                mol = Chem.MolFromMolBlock(mol_string)
            else:
                # Assume SMILES if not molfile
                mol = Chem.MolFromSmiles(mol_string)
        elif input_type == "molfile":
            mol = Chem.MolFromMolBlock(mol_string)
        elif input_type == "smiles":
            mol = Chem.MolFromSmiles(mol_string)
        else:
            results.append((False, f"Invalid input_type: {input_type}"))
            continue
        
        if not mol:
            results.append((False, ""))
            continue
        
        try:
            # Process molecule into fragment graph
            graph = processor.process_molecule(mol)
            
            # Match each fragment to a monomer using graph topology
            for node_id, node in graph.nodes.items():
                # Count connections for this node
                neighbors = graph.get_neighbors(node_id)
                num_connections = len(neighbors)

                # Find matching monomer
                monomer = matcher.find_exact_match(node.mol, num_connections)
                if monomer:
                    node.monomer = monomer
                else:
                    # Generate inline SMILES with R-group markers for unmatched fragments
                    mock_monomer = MonomerData()
                    mock_monomer.is_unknown = True
                    mock_monomer.symbol = _generate_rgroup_smiles(graph, node_id)
                    mock_monomer.name = "Unknown"
                    mock_monomer.r_groups = {f'R{j+1}': '' for j in range(num_connections)}
                    mock_monomer.r_group_count = num_connections
                    node.monomer = mock_monomer
            
            # Try to recover unmatched fragments by merging with neighbors
            max_recovery_attempts = 3  # Prevent infinite loops
            for attempt in range(max_recovery_attempts):
                had_changes = processor.recover_unmatched_fragments(graph, matcher)
                if not had_changes:
                    break
            
            # After regular recovery, try stereo-agnostic matching for remaining unmatched fragments
            # This handles poor quality data with missing stereochemistry
            stereo_matched = processor.recover_unmatched_with_stereo_agnostic(graph, matcher)
            if stereo_matched > 0:
                print(f"DEBUG: Stereo-agnostic recovery matched {stereo_matched} additional fragments")

            # Final pass: merge pairs of both-unmatched neighbor fragments
            # with stereo-agnostic matching (handles split monomers like Phe_4Sdihydroorotamido)
            for attempt in range(max_recovery_attempts):
                had_changes = processor.recover_unmatched_by_merging_stereo_agnostic(graph, matcher)
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


def convert_molfiles_to_helm(molfiles: list, library_json: str = None) -> list:
    """
    Convert a batch of molfiles to HELM notation.
    Convenience wrapper for convert_molecules_batch with input_type="molfile".
    
    Args:
        molfiles: List of molfile strings
        library_json: Optional monomer library as JSON string
    
    Returns:
        List of tuples: (success: bool, helm_notation: str)
    """
    return convert_molecules_batch(molfiles, library_json=library_json, input_type="molfile")


def convert_smiles_to_helm(smiles_list: list, library_json: str = None) -> list:
    """
    Convert a batch of SMILES to HELM notation.
    Convenience wrapper for convert_molecules_batch with input_type="smiles".
    
    Args:
        smiles_list: List of SMILES strings
        library_json: Optional monomer library as JSON string
    
    Returns:
        List of tuples: (success: bool, helm_notation: str)
    """
    return convert_molecules_batch(smiles_list, library_json=library_json, input_type="smiles")

global libraryJSON
with open(libraryFile) as f:
    libraryJSON = f.read()

res_helm_list = convert_molecules_batch(molListToProcess, library_json=libraryJSON)
result_helm = pd.DataFrame(map(lambda x: x[1], res_helm_list), columns=["regenerated sequences"])