import streamlit as st
import streamlit.components as components
from render_molecule import render_molecule, mol_from_file, mol_from_smiles

# Layout configuration
st.set_page_config("3D MolView", layout="centered")
st.title("3D MolView")
st.divider()
# col1, col2 = st.columns([1,2], gap="medium", border=True)

# Initialize session state for molecule
if "molecule" not in st.session_state:
    st.session_state.molecule = ""

def on_visualize():
    if molecule_file is not None:
        # Process file upload if available
        try:
            st.session_state.molecule = mol_from_file(molecule_file)
        except Exception as e:
            st.error(f"Could not process the uploaded file: {str(e)}")
    elif text_input:
        # Process text input if no file but text is available
        try:
            st.session_state.molecule = mol_from_smiles(text_input)
        except Exception as e:
            st.error(f"Could not process the SMILES input: {str(e)}")
    else:
        st.warning("Please enter a SMILES string or upload a file")


with st.container():
    st.header("Molecular Viewer")

    # Only render when a molecule is specified
    if st.session_state.molecule:
        try:
            st.components.v1.html(render_molecule(st.session_state.molecule), height=600, width=800)
        except Exception as e:
            st.error(f"Error rendering molecule: {str(e)}")
    else:
        st.info("Ready for Visualization", width=350)

with st.container():
    st.header("Molecule Input")
    text_input = st.text_input("Enter SMILES String or Molecule Name:", placeholder="e.g. CCO or Ethanol")
    molecule_file = st.file_uploader("Upload Molecular File", type=["MOL", "SDF", "XYZ"])
    st.button("Visualize", on_click=on_visualize, use_container_width=True)