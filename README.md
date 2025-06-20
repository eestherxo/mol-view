# 🔬MolView
A web-based molecular visualization tool that allows users to input SMILES strings, molecule names or upload molecular files to interactively view 3D molecular structures.

---

## 💡Features
- Parse and display molecular structures from:
    - SMILES Strings
    - IUPAC/Common names
    - Molecular files 
- Supported file formats:
    - SDF (Structure Data File)
    - MOL (MDL Molfile)
    - XYZ (XYZ coordinate file)

---

## 🌐Live Demo 
👉 [Try it on Streamlit](https://mol-view.streamlit.app) 

[![Live Demo](https://img.shields.io/badge/Live-Demo-blue?style=for-the-badge&logo=streamlit)](https://mol-view.streamlit.app)

--- 

## 🧰Tech Stack 
- RDKit: Molecular parsing, SMILES handling and 3D coordinate generation 
- Cirpy: Molecule name conversion to SMILES
- Open Babel: File format conversion 
- Streamlit: UI and deployment
- 3Dmol.js: Interactive 3D visualization of molecular structures

--- 

## ▶️Run Locally
1. Clone the repository:
```bash
git clone https://github.com/eestherxo/mol-view.git
cd mol-view
```

2. Install the required dependencies:
```bash
conda env create -f environment.yml
conda activate mol-view
```

3. Run the application:
```bash
streamlit run main.py
```
