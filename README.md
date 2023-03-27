# parKVFinder

![GitHub release (latest by date)](https://img.shields.io/github/v/release/LBC-LNBio/parKVFinder)
![GitHub all releases](https://img.shields.io/github/downloads/LBC-LNBio/parKVFinder/total)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/LBC-LNBio/parKVFinder/testing.yml)

Parallel KVFinder (parKVFinder) is a powerful open-source software designed to detect and comprehensively characterize biomolecular cavities of any type. It provides detailed information on the spatial, depth, constitutional, and hydropathy characteristics of each cavity.

The spatial characterization of each cavity includes information on its shape, volume, and surface area. The depth characterization determines the depth of each cavity point, shown in the B-factors, and calculates the average and maximum depth per cavity. The constitutional characterization identifies the amino acids that form the identified cavities. Lastly, the hydropathy characterization maps the Eisenberg & Weiss hydrophobicity scale onto the surface points of each cavity, shown in the Q-factor, and estimates the average hydropathy for each cavity.

With its extensive range of features, parKVFinder is an invaluable tool for various biomolecular applications, such as drug discovery, durg design, and protein-protein interaction studies.

## Wiki

We provide a [wiki](https://github.com/LBC-LNBio/parKVFinder/wiki), this page was built to help you get started with our cavity detection software.

1. [Download & Installation](https://github.com/LBC-LNBio/parKVFinder/wiki/parKVFinder-Installation)
    - [GCC Installation](https://github.com/LBC-LNBio/parKVFinder/wiki/GCC-Installation)
2. [PyMOL Plug-in Installation](https://github.com/LBC-LNBio/parKVFinder/wiki/PyMOL-Plugin-Installation)
    - [PyMOL Installation](https://github.com/LBC-LNBio/parKVFinder/wiki/PyMOL-Installation)
3. [Tutorial](https://github.com/LBC-LNBio/parKVFinder/wiki/parKVFinder-Tutorial)
4. [Manual](https://github.com/LBC-LNBio/parKVFinder/wiki/parKVFinder-Manual)
5. [About](https://github.com/LBC-LNBio/parKVFinder/wiki/About)

If you are planning on using parKVFinder on Windows, refer to this parKVFinder-win [repository](https://github.com/LBC-LNBio/parKVFinder-win).


## Citation

If you use <ins>parKVFinder</ins>, please cite:

João Victor da Silva Guerra, Helder Veras Ribeiro Filho, Leandro Oliveira Bortot, Rodrigo Vargas Honorato, José Geraldo de Carvalho Pereira, Paulo Sergio Lopes de Oliveira. ParKVFinder: A thread-level parallel approach in biomolecular cavity detection. SoftwareX (2020). https://doi.org/10.1016/j.softx.2020.100606.

If you use <ins>depth and hydropathy characterization</ins>, please also cite:

Guerra, J.V.d., Ribeiro-Filho, H.V., Jara, G.E. et al. pyKVFinder: an efficient and integrable Python package for biomolecular cavity detection and characterization in data science. BMC Bioinformatics 22, 607 (2021). https://doi.org/10.1186/s12859-021-04519-4.

## License

The software is licensed under the terms of the GNU General Public License version 3 (GPL3) and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

---
