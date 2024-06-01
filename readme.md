## Protein Structure Proximity Network

This project aims to create a network out of a protein structure stored in the standard PDB format.

The proximity network (in this context) of a protein is defined as a graph $G=(V, E)$ where $V$ is the set of amino acid residues (indexed) within the protein sequence and $E$ is the set of edges connecting these residues. An edge is established between two residues when they are close to each other, i.e., distance is small and will be taken $=<8Å$ by default.

### Dependencies:

- Biopython: _PDB module to extract and parse PDB files_
- NetworkX: _for graph construction and analysis_
- Pyvis: _for netwokk visualization_
- Pandas
- Numpy
- Seaborn
- Matplotlib

<!-- <iframe src="test/karate_graph.html" width="100%" height="600"></iframe> -->

### Installation

First clone this repository locally:

```bash
git clone https://github.com/raysas/protein-structure-proximity-network.git
```

### Running the code

The workflow takes from a use a pdb id to generate the network. All that has to be done is to run the following command on the terminal:

```bash
python code/generate_network.py <pdb_id>
```

Or if you wish to set the distance threshold to a different value, you can run:

```bash
python code/generate_network.py <pdb_id> <distance_threshold>
```

_Make sure you input a valid PDB id. Check [here](https://www.rcsb.org/docs/general-help/identifiers-in-pdb) for more info._

If you wish you can run each function from extracting the pdb file to visualzing tee network by importing the python module and calling the functions.

```python
import sys
sys.path.append('<path-to-this-repo>/code')

from generate_network import *
```

***p.s.** this is a work in progress and the code is not yet optimized. More features to be added soon listed in the [to-do](./TODO.md) list*

### Output

In the data folder you will find a new directory of named by pdb_id after running


data
 └──  pdb_id
        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── log.txt
        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├──pdb_id.pdb
        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├──pdb_id_contact_map.png
        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├──pdb_id_t_network.graphml
        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; └──pdb_id_t_network_viz.html

After retreiving teh PDB structure, it will extract residues and coordinates information to generate a contact map. The contact map is a heatmap of the distance between residues. Note that the coordinates for each residues are defined by the coordinates of their alpha carbon.

<p align='center'>
<img src='./images/6xdc_contact_map.png' alt='contact map' width=70%>
</p>

The network is then generated (by default threshold distance=8) and saved in graphml format. The network is also visualized in an interactive html file. The layout is set based on the x and y coordinates for each residues.

<p align='center'>
<img src='./images/6xdc_net.png' alt='protein proximity network' width=70%>
</p>


## References

Hamelryck, T., Manderick, B. (2003) PDB parser and structure class implemented in Python. Bioinformatics 19: 2308–2310

