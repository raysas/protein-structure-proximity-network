# Protein Structure Proximity Network

This project aims to create a network out of a protein structure stored in the standard PDB format.

The proximity network (in this context) of a protein is defined as a graph $G=(V, E)$ where $V$ is the set of amino acid residues (indexed) within the protein sequence and $E$ is the set of edges connecting these residues. An edge is established between two residues when they are close to each other, i.e., distance is small and will be taken $<5Å$ by default.

## Dependencies:

- Biopython: _PDB module to extract and parse PDB files_
- Nglview: _for protein structure visualization_
- NetworkX: _for graph construction and analysis_
- Pyvis: _for netwokk visualization_

<!-- <iframe src="test/karate_graph.html" width="100%" height="600"></iframe> -->

## References

Hamelryck, T., Manderick, B. (2003) PDB parser and structure class implemented in Python. Bioinformatics 19: 2308–2310

