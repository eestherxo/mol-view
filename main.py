import streamlit as st
import streamlit.components as components
from render_molecule import render_molecule

# Layout configuration
st.set_page_config("3D MolView", layout="wide")
st.title("3D MolView")
col1, col2 = st.columns([1,2], gap="medium", border=True)

# Initialize session state for molecule
if "molecule" not in st.session_state:
    st.session_state.molecule = ""

def on_visualize():
    st.session_state.molecule = molecule_input

with col1:
    st.header("Molecule Input")
    molecule_input = st.text_input("Enter SMILES String or Molecule Name:", placeholder="e.g. CCO or Ethanol")
    st.button("Visualize", on_click=on_visualize, use_container_width=True)
    st.file_uploader("Upload Molecular File", type=["MOL", "SDF", "XYZ"])
    
with col2:
    st.header("Molecular Viewer")
    
    # Only render when a molecule is specified
    if st.session_state.molecule:
        try:
            st.components.v1.html(render_molecule(st.session_state.molecule), height=600, width=800)
        except Exception as e:
            st.error(f"Error rendering molecule: {str(e)}")
    else:
        st.info("Ready for Visualization", width=350)
