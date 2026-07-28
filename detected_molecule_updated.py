import streamlit as st
import requests
import random
import base64
from io import BytesIO
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Draw, AllChem, rdMolDescriptors
import py3Dmol
import stmol
from streamlit_ketcher import st_ketcher

# ------------------- Helper functions --------------------

def mol_to_image(smiles, size=(400, 200)):
    """Convert SMILES to 2D molecular image."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=size)
    return img

def get_molecule_name(smiles):
    """Query PubChem for the common name (Title) of a SMILES."""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/Title/JSON"
        resp = requests.get(url, timeout=5)
        if resp.status_code == 200:
            data = resp.json()
            return data['PropertyTable']['Properties'][0]['Title']
        else:
            return None
    except:
        return None

def get_molecule_info(smiles):
    """Get comprehensive molecule info from PubChem."""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/MolecularFormula,MolecularWeight,InChI,InChIKey,IUPACName/JSON"
        resp = requests.get(url, timeout=5)
        if resp.status_code == 200:
            data = resp.json()
            props = data['PropertyTable']['Properties'][0]
            return {
                'formula': props.get('MolecularFormula', 'N/A'),
                'weight': props.get('MolecularWeight', 'N/A'),
                'inchi': props.get('InChI', 'N/A'),
                'inchikey': props.get('InChIKey', 'N/A'),
                'iupac': props.get('IUPACName', 'N/A')
            }
        return None
    except:
        return None

def compute_properties(smiles):
    """Compute molecular properties from SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    props = {
        "Molecular Weight": round(Descriptors.MolWt(mol), 2),
        "LogP (cLogP)": round(Descriptors.MolLogP(mol), 2),
        "H‑bond Donors": Lipinski.NumHDonors(mol),
        "H‑bond Acceptors": Lipinski.NumHAcceptors(mol),
        "Rotatable Bonds": Lipinski.NumRotatableBonds(mol),
        "Heavy Atoms": Descriptors.HeavyAtomCount(mol),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "Fraction Csp3": round(Descriptors.FractionCSP3(mol), 2),
        "Aromatic Rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "Saturated Rings": rdMolDescriptors.CalcNumSaturatedRings(mol),
        "Num Rings": rdMolDescriptors.CalcNumRings(mol),
        "Num Atoms": mol.GetNumAtoms(),
        "Num Bonds": mol.GetNumBonds(),
    }
    return props

def lipinski_check(props):
    """Check Lipinski's Rule of 5 violations."""
    violations = []
    if props["Molecular Weight"] > 500:
        violations.append("MW > 500")
    if props["LogP (cLogP)"] > 5:
        violations.append("LogP > 5")
    if props["H‑bond Donors"] > 5:
        violations.append("H‑bond donors > 5")
    if props["H‑bond Acceptors"] > 10:
        violations.append("H‑bond acceptors > 10")
    passed = len(violations) == 0
    return passed, violations

def drug_likeness_summary(props):
    """Return a fun summary string."""
    passed, violations = lipinski_check(props)
    if passed:
        return "✅ Likely good oral bioavailability (Lipinski compliant)."
    else:
        return f"⚠️ {len(violations)} Lipinski violation(s). May have absorption issues."

def create_3d_viewer(smiles):
    """Create and return 3D viewer HTML."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    # Add hydrogens and generate 3D conformer
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_h, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol_h)
    
    # Convert to PDB block
    pdb = Chem.MolToPDBBlock(mol_h)
    
    # Create 3D viewer
    viewer = py3Dmol.view(width=500, height=400)
    viewer.addModel(pdb, "pdb")
    viewer.setStyle({"stick": {"color": "spectrum"}, "sphere": {"scale": 0.2}})
    viewer.zoomTo()
    
    # Add controls
    viewer.addButton('Reset View', 'view.reset()')
    viewer.addButton('Spin On/Off', 'spin = !spin; if(spin) view.spin(true); else view.spin(false);')
    
    return viewer

def random_molecule():
    """Return a random SMILES from a curated list of drugs."""
    drugs = {
        "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "Paracetamol": "CC(=O)NC1=CC=C(C=C1)O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Penicillin G": "CC1(C(=O)N2C(C(=O)NC(CC3=CC=CC=C3)C(=O)O)SC2(C)C)N1C(=O)CC4=CC=CC=C4",
        "Dopamine": "C1=CC(=C(C=C1CCN)O)O",
        "Serotonin": "C1=CC2=C(C=C1O)NC=C2CCN",
        "Glucose": "C(C1C(C(C(C(O1)O)O)O)O)O",
        "Ethanol": "CCO",
        "Cisplatin": "N.N.Cl[Pt]Cl",
        "Morphine": "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O",
        "Nicotine": "CN1CCCC1C2=CN=CC=C2",
    }
    name, smi = random.choice(list(drugs.items()))
    return name, smi

def get_property_explanation(prop_name):
    """Return detailed explanation for each molecular property."""
    explanations = {
        "Molecular Weight": """
        **Molecular Weight (MW)** - The total mass of all atoms in the molecule.
        - **Why it matters**: Smaller molecules (<500 Da) are generally better absorbed orally
        - **Lipinski's Rule**: Should be ≤ 500 Da for good oral bioavailability
        - **Drug-like range**: 150-500 Da
        """,
        "LogP (cLogP)": """
        **LogP (Partition Coefficient)** - Measures how lipophilic (fat-loving) a molecule is.
        - **What it means**: Higher LogP = more hydrophobic, better at crossing cell membranes
        - **Lipinski's Rule**: Should be ≤ 5 for good oral bioavailability
        - **Drug-like range**: 1-3 is optimal for oral drugs
        - **Too high (>5)**: May have poor solubility and toxicity issues
        """,
        "H‑bond Donors": """
        **Hydrogen Bond Donors** - Atoms that can donate hydrogen in hydrogen bonds.
        - **Examples**: -OH, -NH, -NH₂ groups
        - **Lipinski's Rule**: Should be ≤ 5
        - **Why it matters**: Too many donors reduce membrane permeability
        """,
        "H‑bond Acceptors": """
        **Hydrogen Bond Acceptors** - Atoms that can accept hydrogen in hydrogen bonds.
        - **Examples**: Oxygen and nitrogen atoms in functional groups
        - **Lipinski's Rule**: Should be ≤ 10
        - **Why it matters**: Too many acceptors can reduce oral absorption
        """,
        "Rotatable Bonds": """
        **Rotatable Bonds** - Single bonds that can rotate around their axis.
        - **What it means**: More rotatable bonds = more flexible molecule
        - **Drug-like range**: ≤ 10 for good bioavailability
        - **Why it matters**: Flexibility affects how well the drug can bind to targets
        """,
        "Heavy Atoms": """
        **Heavy Atoms** - All atoms except hydrogen.
        - **What it measures**: Molecular size and complexity
        - **Drug-like range**: 10-30 heavy atoms
        - **Why it matters**: Influences synthesis difficulty and drug-likeness
        """,
        "TPSA": """
        **Topological Polar Surface Area** - Surface area of polar atoms (oxygen, nitrogen).
        - **What it measures**: Ability to cross cell membranes
        - **Drug-like range**: ≤ 140 Å² for good brain permeability
        - **Optimal range**: 20-100 Å² for oral drugs
        - **Why it matters**: Low TPSA = better membrane crossing ability
        """,
        "Fraction Csp3": """
        **Fraction of sp³ Hybridized Carbons** - Carbons with tetrahedral geometry.
        - **What it measures**: Molecular 3D complexity
        - **Drug-like range**: > 0.3 for better drug-likeness
        - **Why it matters**: More sp³ carbons = more 3D structure, better binding
        """,
        "Aromatic Rings": """
        **Aromatic Rings** - Ring structures with delocalized electrons (benzene rings).
        - **Examples**: Benzene, pyridine, furan
        - **Drug-like range**: 1-3 aromatic rings
        - **Why it matters**: Aromatic rings influence drug binding and metabolism
        """,
        "Saturated Rings": """
        **Saturated Rings** - Rings with only single bonds (no double bonds).
        - **Examples**: Cyclohexane, piperidine
        - **Why it matters**: Adds 3D character and complexity to molecules
        """,
        "Num Rings": """
        **Total Number of Rings** - All ring structures combined.
        - **Drug-like range**: 1-4 rings
        - **Why it matters**: Rings add rigidity and 3D shape to drug molecules
        """
    }
    return explanations.get(prop_name, "Property description not available.")

def display_smiles_help():
    """Display help for writing valid SMILES strings."""
    st.sidebar.markdown("""
    ### ✏️ SMILES Writing Guide
    
    **Basic Rules:**
    - Atoms are represented by element symbols (C, N, O, etc.)
    - Single bonds are implicit (just write atoms together)
    - Double bonds use `=` (e.g., `C=C`)
    - Triple bonds use `#` (e.g., `C#C`)
    - Rings use numbers (e.g., `C1CC1`)
    - Branches use parentheses (e.g., `CC(C)C`)
    
    **Examples:**
    - Ethanol: `CCO`
    - Benzene: `c1ccccc1` (lowercase `c` for aromatic)
    - Aspirin: `CC(=O)OC1=CC=CC=C1C(=O)O`
    
    **Common Mistakes:**
    - ❌ Wrong: `C6H6` (use `c1ccccc1`)
    - ❌ Wrong: `CH3CH2OH` (use `CCO`)
    - ❌ Wrong: `Benzene` (must be SMILES string)
    
    **Tips:**
    - Start from a known SMILES string
    - Use PubChem to find SMILES
    - Check your string on SMILES validator sites
    """)

# ------------------- Streamlit App -------------------------

st.set_page_config(page_title="🧬 DetectED Molecule Explorer", layout="wide")

# Custom CSS for better UI and font colors
st.markdown("""
<style>
    /* Main background and text */
    .stApp {
        background-color: #0e1117;
    }
    
    /* Header styling */
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        background: linear-gradient(45deg, #4facfe, #00f2fe);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        padding: 1rem 0;
    }
    
    /* White text for all main content */
    .stMarkdown, .stText, .stInfo, .stSuccess, .stWarning, .stError {
        color: white !important;
    }
    
    /* Black text specifically for metric labels and values */
    .stMetric label, .stMetric div, .stMetric [data-testid="stMetricLabel"], 
    .stMetric [data-testid="stMetricValue"], .stMetric [data-testid="stMetricDelta"] {
        color: black !important;
    }
    
    /* Override metric background to be white for readability */
    .stMetric {
        background-color: white !important;
        border-radius: 8px;
        padding: 0.5rem 1rem !important;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        margin: 0.2rem 0;
    }
    
    /* Property cards with white background and black text */
    .property-card {
        background-color: white !important;
        border-radius: 10px;
        padding: 1rem;
        margin: 0.5rem 0;
        color: black !important;
    }
    
    .property-card * {
        color: black !important;
    }
    
    /* Make expanders have white text */
    .streamlit-expanderHeader, .streamlit-expanderContent {
        color: white !important;
    }
    
    /* Sidebar text */
    .css-1d391kg, .css-1vq4p4l {
        color: white !important;
    }
    
    /* Tab text */
    .stTabs [data-baseweb="tab"] {
        color: white !important;
    }
    
    /* Code blocks */
    .stCodeBlock {
        background-color: #1e1e1e !important;
        color: #d4d4d4 !important;
    }
    
    /* Buttons */
    .stButton button {
        background-color: #4facfe !important;
        color: white !important;
        border-radius: 8px !important;
    }
    
    /* Input fields */
    .stTextInput input {
        background-color: #1e1e1e !important;
        color: white !important;
        border: 1px solid #4facfe !important;
    }
    
    /* Select boxes */
    .stSelectbox div {
        color: white !important;
    }
    
    /* Make sure tooltips are readable */
    .stTooltipContent {
        color: black !important;
    }
    
    /* Info boxes keep white text */
    .stInfo, .stSuccess, .stWarning, .stError {
        color: white !important;
    }
    
    /* Specifically target metric values to be black */
    div[data-testid="stMetricValue"] {
        color: black !important;
        font-weight: bold;
    }
    
    div[data-testid="stMetricLabel"] {
        color: black !important;
    }
    
    /* Make sure the 3D viewer has proper spacing */
    .stExpander {
        background-color: transparent;
    }
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="main-header">🧬 DetectED Molecule Explorer</div>', unsafe_allow_html=True)
st.markdown("*For teaching cheminformatics – early disease detection & AI*")

# Sidebar
with st.sidebar:
    st.markdown("### 🔬 DetectED Nonprofit")
    st.markdown("*Advancing disease detection through education*")
    
    st.markdown("---")
    st.markdown("### 🎯 Quick Actions")
    
    # Quick example buttons
    col1, col2 = st.columns(2)
    with col1:
        if st.button("💊 Aspirin"):
            st.session_state.quick_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
            st.session_state.quick_name = "Aspirin"
    with col2:
        if st.button("☕ Caffeine"):
            st.session_state.quick_smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
            st.session_state.quick_name = "Caffeine"
    
    col3, col4 = st.columns(2)
    with col3:
        if st.button("🍷 Ethanol"):
            st.session_state.quick_smiles = "CCO"
            st.session_state.quick_name = "Ethanol"
    with col4:
        if st.button("🎲 Random"):
            name, smi = random_molecule()
            st.session_state.quick_smiles = smi
            st.session_state.quick_name = name
    
    st.markdown("---")
    display_smiles_help()
    
    st.markdown("---")
    st.caption("Built with ❤️ by DetectED – Early Disease Detection & AI")

# ------------------- Session state for history ------------
if 'history' not in st.session_state:
    st.session_state.history = []
if 'current_smiles' not in st.session_state:
    st.session_state.current_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
if 'quick_smiles' not in st.session_state:
    st.session_state.quick_smiles = None

# ------------------- Main input area ---------------------

# Check for quick action
if st.session_state.quick_smiles:
    smiles = st.session_state.quick_smiles
    st.session_state.current_smiles = smiles
    st.session_state.quick_smiles = None
else:
    smiles = st.session_state.current_smiles

st.markdown("---")
st.markdown("### ✏️ Draw or Input Your Molecule")

tab1, tab2, tab3 = st.tabs(["✏️ Draw Structure", "⌨️ Type SMILES", "📚 Examples"])

with tab1:
    st.markdown("**Draw your molecule using the editor below:**")
    st.markdown("*Click on atoms to add elements, drag to create bonds. Click 'Apply' when done.*")
    
    try:
        drawn_smiles = st_ketcher(smiles)
        if drawn_smiles and drawn_smiles != smiles:
            smiles = drawn_smiles
            st.session_state.current_smiles = smiles
            # Check if valid
            mol_check = Chem.MolFromSmiles(smiles)
            if mol_check:
                st.success(f"✅ Molecule drawn successfully! SMILES: `{smiles}`")
            else:
                st.error("⚠️ Invalid structure drawn. Please try again.")
    except Exception as e:
        st.error(f"Error loading molecular editor: {str(e)}")
        st.info("Please make sure 'streamlit-ketcher' is installed: `pip install streamlit-ketcher`")

with tab2:
    st.markdown("**Enter SMILES string manually:**")
    manual_smiles = st.text_input("SMILES:", value=smiles)
    if st.button("🔍 Analyze SMILES"):
        if manual_smiles.strip():
            mol_check = Chem.MolFromSmiles(manual_smiles)
            if mol_check:
                smiles = manual_smiles
                st.session_state.current_smiles = smiles
                st.success(f"✅ Valid SMILES: `{smiles}`")
            else:
                st.error("❌ Invalid SMILES string. Please check the format.")

with tab3:
    st.markdown("**Select from common molecules:**")
    examples = {
        "💊 Aspirin (Pain relief)": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "💊 Ibuprofen (Anti-inflammatory)": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "💊 Paracetamol (Pain relief)": "CC(=O)NC1=CC=C(C=C1)O",
        "☕ Caffeine (Stimulant)": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "🧪 Penicillin G (Antibiotic)": "CC1(C(=O)N2C(C(=O)NC(CC3=CC=CC=C3)C(=O)O)SC2(C)C)N1C(=O)CC4=CC=CC=C4",
        "🧬 Dopamine (Neurotransmitter)": "C1=CC(=C(C=C1CCN)O)O",
        "🧬 Serotonin (Hormone)": "C1=CC2=C(C=C1O)NC=C2CCN",
        "🍷 Ethanol (Alcohol)": "CCO",
        "🧪 Morphine (Analgesic)": "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O",
        "🚬 Nicotine (Stimulant)": "CN1CCCC1C2=CN=CC=C2",
    }
    
    cols = st.columns(3)
    for idx, (name, smi) in enumerate(examples.items()):
        with cols[idx % 3]:
            if st.button(name, key=f"ex_{idx}"):
                smiles = smi
                st.session_state.current_smiles = smiles
                st.success(f"Loaded: {name}")

# ------------------- Process the molecule -----------------
if smiles.strip():
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        st.error("❌ Invalid SMILES string. Please check the format and try again.")
        st.info("💡 Tip: Use the Draw Structure tab or examples dropdown for valid molecules.")
    else:
        # ---- Get name ----
        name = get_molecule_name(smiles)
        pubchem_info = get_molecule_info(smiles)
        
        if name:
            st.markdown(f"## 🧪 {name}")
        else:
            st.markdown("## 🧪 Unknown compound")

        # Add to history
        if not st.session_state.history or st.session_state.history[-1][0] != smiles:
            st.session_state.history.append((smiles, name))

        # ---- Display 2D structure and properties ----
        img = mol_to_image(smiles)
        if img:
            col_img, col_props = st.columns([1, 1])
            with col_img:
                st.image(img, caption="2D Structure", use_container_width=True)
                # Download button
                buffered = BytesIO()
                img.save(buffered, format="PNG")
                img_str = base64.b64encode(buffered.getvalue()).decode()
                href = f'<a href="data:image/png;base64,{img_str}" download="molecule.png" style="color: #4facfe;">⬇️ Download 2D Structure</a>'
                st.markdown(href, unsafe_allow_html=True)
                
                # Show SMILES
                st.code(f"SMILES: {smiles}", language="text")
                
                # PubChem info
                if pubchem_info:
                    with st.expander("📋 PubChem Information"):
                        st.write(f"**Formula:** {pubchem_info['formula']}")
                        st.write(f"**Molecular Weight:** {pubchem_info['weight']} Da")
                        st.write(f"**IUPAC Name:** {pubchem_info['iupac']}")
                        st.write(f"**InChI:** `{pubchem_info['inchi'][:100]}...`")
                        st.write(f"**InChIKey:** `{pubchem_info['inchikey']}`")
            
            with col_props:
                # ---- Properties ----
                props = compute_properties(smiles)
                if props:
                    st.markdown("### 📊 Molecular Properties")
                    st.markdown("*Hover over each property for detailed explanation*")
                    
                    # Property display with tooltips - using custom CSS for black text
                    cols = st.columns(2)
                    with cols[0]:
                        st.metric("MW", f"{props['Molecular Weight']} Da", 
                                 help=get_property_explanation("Molecular Weight"))
                        st.metric("LogP", props["LogP (cLogP)"], 
                                 help=get_property_explanation("LogP (cLogP)"))
                        st.metric("H‑bond Donors", props["H‑bond Donors"], 
                                 help=get_property_explanation("H‑bond Donors"))
                        st.metric("TPSA", f"{props['TPSA']} Å²", 
                                 help=get_property_explanation("TPSA"))
                    with cols[1]:
                        st.metric("H‑bond Acceptors", props["H‑bond Acceptors"], 
                                 help=get_property_explanation("H‑bond Acceptors"))
                        st.metric("Rotatable Bonds", props["Rotatable Bonds"], 
                                 help=get_property_explanation("Rotatable Bonds"))
                        st.metric("Heavy Atoms", props['Heavy Atoms'], 
                                 help=get_property_explanation("Heavy Atoms"))
                        st.metric("Total Rings", props['Num Rings'], 
                                 help=get_property_explanation("Num Rings"))
                    
                    with st.expander("🔍 Additional Descriptors"):
                        cols_add = st.columns(2)
                        with cols_add[0]:
                            st.metric("Fraction Csp³", f"{props['Fraction Csp3']:.3f}", 
                                     help=get_property_explanation("Fraction Csp3"))
                            st.metric("Aromatic Rings", props['Aromatic Rings'], 
                                     help=get_property_explanation("Aromatic Rings"))
                        with cols_add[1]:
                            st.metric("Saturated Rings", props['Saturated Rings'], 
                                     help=get_property_explanation("Saturated Rings"))
                            st.metric("Atoms/Bonds", f"{props['Num Atoms']} atoms, {props['Num Bonds']} bonds")

        # ---- Lipinski & fun summary ----
        st.markdown("### 💊 Drug‑likeness Assessment")
        
        # Lipinski explanation
        with st.expander("📚 What is Lipinski's Rule of 5?"):
            st.markdown("""
            **Lipinski's Rule of 5** is a rule of thumb to evaluate drug-likeness:
            
            1. **Molecular Weight** ≤ 500 Da
            2. **LogP** ≤ 5 (lipophilicity)
            3. **H-bond Donors** ≤ 5 (NH, OH groups)
            4. **H-bond Acceptors** ≤ 10 (N, O atoms)
            
            **Why it matters:**
            - Molecules that violate these rules often have poor oral bioavailability
            - Used in early drug discovery to filter candidates
            - Not absolute - many approved drugs violate one or more rules
            
            **Limitations:**
            - Doesn't consider synthetic accessibility
            - Doesn't predict toxicity or efficacy
            - Some drugs (e.g., antibiotics) intentionally violate rules
            """)
        
        # Display Lipinski results
        passed, violations = lipinski_check(props)
        col1, col2 = st.columns([2, 1])
        with col1:
            if passed:
                st.success("✅ **Passes Lipinski Rule of 5** (good oral bioavailability likely).")
            else:
                st.warning("⚠️ **Fails Lipinski** – violations:")
                for v in violations:
                    st.write(f"  • {v}")
            st.info(drug_likeness_summary(props))
        with col2:
            # Summary badge
            if passed:
                st.markdown("""
                <div style="background-color: #1e3c72; border-radius: 10px; padding: 1rem; text-align: center; color: white;">
                    <h3 style="color: #4facfe;">✅ Drug-like</h3>
                </div>
                """, unsafe_allow_html=True)
            else:
                st.markdown(f"""
                <div style="background-color: #3d1e1e; border-radius: 10px; padding: 1rem; text-align: center; color: white;">
                    <h3 style="color: #ff6b6b;">⚠️ {len(violations)} Violations</h3>
                </div>
                """, unsafe_allow_html=True)

        # ---- 3D viewer - ORIGINAL WORKING VERSION ----
        with st.expander("🧬 3D Structure Viewer (Interactive)"):
            st.markdown("*Drag to rotate, scroll to zoom. The molecule is optimized in 3D space.*")
            try:
                import streamlit.components.v1 as components
                
                # Use the original working approach from your code
                mol_h = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol_h, randomSeed=42)
                AllChem.MMFFOptimizeMolecule(mol_h)
                pdb = Chem.MolToPDBBlock(mol_h)
                
                # Create 3D viewer exactly as in your original code
                viewer = py3Dmol.view(width=400, height=400)
                viewer.addModel(pdb, "pdb")
                viewer.setStyle({"stick": {}, "sphere": {"scale": 0.2}})
                viewer.zoomTo()
                viewer.spin(False)
                
                # Display in Streamlit using the original method
                viewer_html = viewer._make_html()
                components.html(viewer_html, height=450, width=450)
                
                # Also try stmol as fallback
                try:
                    stmol.showmol(
                        st,
                        pdb,
                        width=400,
                        height=400,
                        style={"stick": {}, "sphere": {"scale": 0.2}},
                        zoom=1.2,
                        spin=False
                    )
                except:
                    pass  # Fallback already handled above
                
            except Exception as e:
                st.error(f"Error displaying 3D structure: {str(e)}")
                st.info("Make sure you have py3Dmol and stmol installed correctly.")

        # ---- Similarity search (simplified) ----
        with st.expander("🔍 Find Similar Molecules"):
            st.markdown("**Compare with known drugs using Tanimoto similarity:**")
            reference_drugs = {
                "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
                "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
                "Paracetamol": "CC(=O)NC1=CC=C(C=C1)O",
                "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                "Morphine": "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O",
            }
            
            if mol:
                fp_current = Chem.RDKFingerprint(mol)
                similarities = []
                for drug_name, drug_smi in reference_drugs.items():
                    drug_mol = Chem.MolFromSmiles(drug_smi)
                    if drug_mol:
                        fp_drug = Chem.RDKFingerprint(drug_mol)
                        sim = Chem.DataStructs.TanimotoSimilarity(fp_current, fp_drug)
                        similarities.append((drug_name, sim))
                
                similarities.sort(key=lambda x: x[1], reverse=True)
                
                st.write("**Top 3 most similar drugs in our database:**")
                for name, sim in similarities[:3]:
                    st.progress(sim, text=f"{name}: {sim:.2%} similarity")

        # ---- Fun fact ----
        with st.expander("📖 Did you know?"):
            facts = {
                "Aspirin": "Aspirin was first synthesized in 1897 and is one of the most widely used medications.",
                "Ibuprofen": "Ibuprofen was discovered in the 1960s while researching for a new rheumatoid arthritis treatment.",
                "Caffeine": "Caffeine is the world's most widely consumed psychoactive substance.",
                "Paracetamol": "Paracetamol is also known as acetaminophen and is used for pain relief and fever.",
                "Penicillin G": "Penicillin was the first antibiotic discovered by Alexander Fleming in 1928.",
                "Dopamine": "Dopamine is a neurotransmitter that plays a role in reward and motor control.",
                "Serotonin": "Serotonin is a neurotransmitter that contributes to well-being and happiness.",
                "Ethanol": "Ethanol has been consumed by humans for thousands of years and is produced by fermentation.",
                "Morphine": "Morphine is one of the oldest known painkillers, isolated from opium in 1804.",
                "Nicotine": "Nicotine is named after Jean Nicot, who introduced tobacco to France in the 16th century.",
            }
            if name and name in facts:
                st.write(f"🧪 **{name}:** {facts[name]}")
            else:
                st.write("🧪 This molecule is part of our teaching collection – explore its properties!")
                if pubchem_info and pubchem_info['iupac'] != 'N/A':
                    st.write(f"📌 IUPAC Name: {pubchem_info['iupac']}")

        # ---- Export options ----
        with st.expander("💾 Export Data"):
            st.markdown("**Export molecular data:**")
            col1, col2, col3 = st.columns(3)
            
            with col1:
                # SMILES download
                st.download_button(
                    label="📋 Download SMILES",
                    data=smiles,
                    file_name=f"{name if name else 'molecule'}.smi",
                    mime="text/plain"
                )
            
            with col2:
                # Properties as CSV
                if props:
                    import pandas as pd
                    df = pd.DataFrame([props])
                    csv = df.to_csv(index=False)
                    st.download_button(
                        label="📊 Download Properties (CSV)",
                        data=csv,
                        file_name=f"{name if name else 'molecule'}_properties.csv",
                        mime="text/csv"
                    )
            
            with col3:
                # SDF download
                if mol:
                    sdf = Chem.MolToMolBlock(mol)
                    st.download_button(
                        label="🧪 Download SDF",
                        data=sdf,
                        file_name=f"{name if name else 'molecule'}.sdf",
                        mime="chemical/x-mdl-sdfile"
                    )

        # ---- Interactive Quiz ----
        with st.expander("🧠 Interactive Quiz"):
            st.markdown("**Test your molecular knowledge:**")
            
            quiz_options = [
                "What does LogP measure?",
                "What is the Lipinski limit for H-bond donors?",
                "What does TPSA stand for?",
            ]
            
            quiz_answers = {
                "What does LogP measure?": "Lipophilicity (fat-liking ability)",
                "What is the Lipinski limit for H-bond donors?": "≤ 5 donors",
                "What does TPSA stand for?": "Topological Polar Surface Area",
            }
            
            selected_quiz = st.selectbox("Choose a quiz question:", quiz_options)
            if st.button("Show Answer"):
                st.info(f"**Answer:** {quiz_answers[selected_quiz]}")
            
            if st.button("🎲 Random Drug Fact"):
                facts_list = [
                    "Penicillin was discovered by accident in 1928 when Alexander Fleming noticed mold killing bacteria.",
                    "Aspirin's chemical name is acetylsalicylic acid.",
                    "Caffeine is found in coffee, tea, and chocolate.",
                    "The molecular formula of ethanol is C₂H₅OH.",
                    "Dopamine is involved in the brain's reward system.",
                    "Paracetamol was discovered in 1877.",
                    "Morphine is still used as a powerful painkiller in hospitals.",
                ]
                st.success(random.choice(facts_list))

# ---- History section ----
st.markdown("---")
st.markdown("### 📜 Recent Molecules")

if st.session_state.history:
    cols = st.columns(min(len(st.session_state.history), 5))
    for idx, (smi, name) in enumerate(st.session_state.history[-5:]):
        display_name = name if name else f"Mol {idx+1}"
        with cols[idx % len(cols)]:
            if st.button(f"{display_name}", key=f"hist_{idx}"):
                st.session_state.current_smiles = smi
                st.rerun()
    
    if st.button("🗑️ Clear History", key="clear_history"):
        st.session_state.history = []
        st.rerun()
else:
    st.write("No molecules viewed yet. Start exploring!")

# ---- Footer ----
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #888; padding: 1rem;">
    <p>🧬 DetectED Molecule Explorer v2.0 | Built with Streamlit, RDKit, and ❤️</p>
    <p style="font-size: 0.8rem;">For educational purposes in cheminformatics and AI</p>
</div>
""", unsafe_allow_html=True)
