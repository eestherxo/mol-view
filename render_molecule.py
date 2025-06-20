from typing import BinaryIO
import tempfile
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import pybel
import cirpy


def mol_from_file(molecule_file: BinaryIO) -> Chem.Mol:
    """
    Convert an uploaded molecular file (MOL, SDF, XYZ) to an RDKit molecule object.

    Args:
        molecule_file: The uploaded file object

    Returns:
        RDKit molecule object.
    """
    file_extension = molecule_file.name.split(".")[-1].lower()
    with tempfile.NamedTemporaryFile(delete=False, suffix=f".{file_extension}") as temp_file:
        temp_file.write(molecule_file.read())
        temp_path = temp_file.name

    if file_extension == "mol":
        mol = Chem.MolFromMolFile(temp_path)
        if mol is None:
            raise ValueError("Could not parse MOL file")

    elif file_extension == "sdf":
        # load the first molecule from an SDF file
        supplier = Chem.SDMolSupplier(temp_path)
        if len(supplier) == 0 or supplier[0] is None:
            raise ValueError("Could not parse SDF file")
        return supplier[0]

    elif file_extension == "xyz":
        pybel_mol = next(pybel.readfile("xyz", temp_path))
        mol_block = pybel_mol.write("mol")
        mol = Chem.MolFromMolBlock(mol_block)
        if mol is None:
            raise ValueError("Could not convert XYZ to RDKit molecule")

    else:
        raise ValueError(f"Unsupported file format: {file_extension}")

    return mol


def mol_from_smiles(molecule: str) -> Chem.Mol:
    """
    Convert a SMILES string or Molecule name to a RDKit molecule object
    Args:
        molecule (str): The SMILES representation or name of the molecule.

    Returns:
        RDKit molecule object.
    """
    # Convert SMILES to RDKit molecule object
    mol = Chem.MolFromSmiles(molecule.strip())
    if mol is None:
        try:
            cirpy_smiles = cirpy.resolve(molecule.strip(), 'smiles')
            if cirpy_smiles:
                mol = Chem.MolFromSmiles(cirpy_smiles)
        except Exception as e:
            raise ValueError(f"Could not convert SMILES or name to RDKit molecule: {str(e)}")
    return mol


def render_molecule(molecule: Chem.Mol) -> str:
    """
    Render a RDKit molecule object  using 3Dmol.js.

    Args:
        molecule: A RDKit molecule object

    Returns:
        str: HTML string containing the rendered 3D model of the molecule.
    """
    if molecule is None:
        raise ValueError("No molecule provided")

    mol = Chem.AddHs(molecule)

    # Generate 3D coordinates for the molecule
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)

    # Convert RDKit molecule to 3Dmol.js format
    block = Chem.MolToMolBlock(mol)

    # Get molecule name safely
    try:
        # Try to get SMILES representation if possible, but handle gracefully if not
        smiles = Chem.MolToSmiles(molecule) if molecule else "Unknown"
        molecule_name = smiles
        # Only try to resolve name if we have a valid SMILES
        if smiles and smiles != "Unknown":
            try:
                resolved_name = cirpy.resolve(smiles, "iupac_name")
                if resolved_name:
                    molecule_name = resolved_name
            except:
                pass  # If cirpy fails, just use the SMILES
    except:
        molecule_name = "Unknown Molecule"

    # Create a 3Dmol.js view
    html_template = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>3D Molecular Viewer</title>
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
        <style>
            body {{
                margin: 0;
                padding: 0;
                font-family: Arial, sans-serif;
            }}
            #viewer {{
                width: 100%;
                height: 100vh;
                position: relative;
            }}
            .controls {{
                position: absolute;
                top: 10px;
                left: 10px;
                z-index: 1000;
                background: rgba(255, 255, 255, 0.9);
                padding: 10px;
                border-radius: 5px;
                box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
            }}
            .control-button {{
                margin: 2px;
                padding: 5px 10px;
                border: none;
                border-radius: 3px;
                background: #3498db;
                color: white;
                cursor: pointer;
                font-size: 12px;
            }}
            .control-button:hover {{
                background-color: #2980b9;
            }}
            .info {{
                position: absolute;
                bottom: 10px;
                left: 10px;
                z-index: 1000;
                background: rgba(255, 255, 255, 0.9);
                padding: 5px 10px;
                border-radius: 3px;
                font-size: 12px;
                color: #333;
            }}
        </style>
    </head>
    <body>
    <div id="viewer"></div>
    <div class="controls">
        <div style="margin-bottom: 5px; font-weight: bold; color: #333;">Molecule: {molecule_name}</div>
        <button class="control-button" onclick="resetView()">Reset View</button>
        <button class="control-button" onclick="toggleSpin()">Toggle Spin</button>
    </div>
    <div class="info">
        Use mouse to rotate, scroll to zoom
    </div>
    <script>
        let viewer = $3Dmol.createViewer("viewer", {{
            defaultcolors: $3Dmol.rasmolElementColors
        }});
        
        let spinning = false;
        
        const molData = `{block}`;

        function loadMolecule() {{
            viewer.clear();
            viewer.addModel(molData, 'mol');
            viewer.setStyle({{'stick': {{'radius':0.15}}, 'sphere': {{'scale': 0.25}}}});
            viewer.setBackgroundColor('white');
            viewer.zoomTo();
            viewer.render();
        }}
        
        function resetView() {{
            viewer.zoomTo();
            viewer.render();
        }}
        
        function toggleSpin() {{
            if (spinning) {{
                viewer.spin(false);
                spinning = false;
            }} else {{
                viewer.spin("y", 1);
                spinning = true;
            }}
        }}

        loadMolecule();
    </script>
    </body>
    </html>
    """

    return html_template
