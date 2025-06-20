import streamlit as st
import streamlit.components as components
from render_molecule import render_molecule, mol_from_file, mol_from_smiles

# Layout configuration
st.set_page_config("3D MolView", layout="wide")
st.title("3D MolView")

# Create columns without borders initially
col1, col2 = st.columns([1,2], gap="medium")

# Initialize session state for molecule
if "molecule" not in st.session_state:
    st.session_state.molecule = None

# The visualize function accesses the input values through the form
def on_visualize():
    text_input = st.session_state.get("text_input", "")
    molecule_file = st.session_state.get("molecule_file")

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

with col1:
    with st.form("molecule_input_form"):
        st.header("Molecule Input")
        st.text_input("Enter SMILES String or Molecule Name:",
                     placeholder="e.g. CCO or Ethanol",
                     key="text_input")
        st.file_uploader("Upload Molecular File",
                        type=["MOL", "SDF", "XYZ"],
                        key="molecule_file")
        st.form_submit_button("Visualize", on_click=on_visualize, use_container_width=True)

# Add a container with border to column 2
with col2:
    with st.container(border=True):
        st.header("Molecular Viewer")

        # Only render when a molecule is specified
        if st.session_state.molecule:
            try:
                st.components.v1.html(render_molecule(st.session_state.molecule), height=600, width=800)
            except Exception as e:
                st.error(f"Error rendering molecule: {str(e)}")
        else:
            st.info("Ready for Visualization", width=300)
