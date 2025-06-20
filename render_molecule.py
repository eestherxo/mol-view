from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import cirpy

def render_molecule(molecule: str) -> str:
    """
    Render a molecule from its SMILES representation or name using RDKit and 3Dmol.js.

    Args:
        molecule (str): The SMILES representation or name of the molecule.

    Returns:
        str: HTML string containing the rendered 3D model of the molecule.
    """
    # Convert SMILES to RDKit molecule object
    mol = Chem.MolFromSmiles(molecule)
    if mol is None:
        mol = cirpy.resolve(molecule, 'smiles')
        mol = Chem.MolFromSmiles(mol)
    mol = Chem.AddHs(mol)

    # Generate 3D coordinates for the molecule
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)

    # Convert RDKit molecule to 3Dmol.js format
    block = Chem.MolToMolBlock(mol)

    # Create a 3Dmol.js view
    # viewer = py3Dmol.view(width=800, height=600)
    # viewer.addModel(block, "mol")
    # viewer.setStyle({'stick': {'radius':0.15}, 'sphere': {'scale': 0.25}})
    # viewer.setBackgroundColor('white')
    # viewer.zoomTo()

    html_template = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>3D Molecular Viewer</title>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.0/3Dmol-min.js"></script>
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
        <div style="margin-bottom: 5px; font-weight: bold; color: #333;">Molecule: {cirpy.resolve(molecule, "smiles")}</div>
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
