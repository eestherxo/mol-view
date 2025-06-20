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
👉 [Try it on Streamlit](https://your-app-name.streamlit.app) 

[![Live Demo](https://img.shields.io/badge/Live-Demo-blue?style=for-the-badge&logo=streamlit)](https://your-app-name.streamlit.app)

--- 

## 🧰Tech Stack 
- RDKit: Molecular parsing, SMILES handling and 3D coordinate generation 
- Cirpy: Molecule name conversion to SMILES
- Open Babel: File format conversion 
- 3Dmol.js: Interactive 3D visualization of molecular structures

--- 

## ▶️Run Locally
```bash
streamlit run main.py
```