import streamlit as st
import requests
import random
import base64
from io import BytesIO
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Draw, AllChem, rdMolDescriptors

# ------------------- Optional 3D viewer -------------------
try:
    import py3Dmol
    import stmol
    has_3d = True
except ImportError:
    has_3d = False
    st.info("For 3D view, install: pip install py3Dmol stmol")

# ------------------- Helper functions --------------------

def mol_to_image(smiles, size=(400, 200)):
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

def compute_properties(smiles):
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
    }
    return props

def lipinski_check(props):
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

def show_3d_structure(smiles):
    """Display interactive 3D using py3Dmol."""
    if not has_3d:
        st.warning("Install py3Dmol and stmol to view 3D structures.")
        return
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return
    # Add hydrogens and generate 3D conformer
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_h, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol_h)
    # Convert to PDB block
    pdb = Chem.MolToPDBBlock(mol_h)
    stmol.showmol(
        st,
        pdb,
        width=400,
        height=300,
        style={"stick": {}, "sphere": {"scale": 0.2}},
        zoom=1.2,
        spin=False
    )

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
        "Cisplatin": "N.N.Cl[Pt]Cl",  # works but no 2D
    }
    name, smi = random.choice(list(drugs.items()))
    return name, smi

# ------------------- Streamlit App -------------------------

st.set_page_config(page_title="DetectED Molecule Explorer", layout="wide")
st.title("🧬 DetectED Molecule Explorer")
st.markdown("*For teaching cheminformatics – early disease detection & AI*")
st.sidebar.image("https://via.placeholder.com/150?text=DetectED", use_container_width=True)  # Replace with your logo
st.sidebar.markdown("### 🎯 DetectED Nonprofit")

# ------------------- Session state for history ------------
if 'history' not in st.session_state:
    st.session_state.history = []

# ------------------- Main input area ---------------------
col1, col2 = st.columns([2, 1])

with col1:
    input_method = st.radio("Choose input method:", ("Select from examples", "Type SMILES", "Random molecule"))

    if input_method == "Select from examples":
        examples = {
            "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
            "Paracetamol": "CC(=O)NC1=CC=C(C=C1)O",
            "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "Penicillin G": "CC1(C(=O)N2C(C(=O)NC(CC3=CC=CC=C3)C(=O)O)SC2(C)C)N1C(=O)CC4=CC=CC=C4",
            "Dopamine": "C1=CC(=C(C=C1CCN)O)O",
            "Serotonin": "C1=CC2=C(C=C1O)NC=C2CCN",
            "Glucose": "C(C1C(C(C(C(O1)O)O)O)O)O",
            "Ethanol": "CCO",
        }
        chosen = st.selectbox("Pick a molecule", list(examples.keys()))
        smiles = examples[chosen]
        st.info(f"SMILES: `{smiles}`")
    elif input_method == "Type SMILES":
        smiles = st.text_input("Enter SMILES string:", value="CC(=O)OC1=CC=CC=C1C(=O)O")
        if st.button("Show molecule"):
            pass
    else:  # Random
        if st.button("🎲 Give me a random molecule"):
            name, smi = random_molecule()
            smiles = smi
            st.session_state.random_name = name
        else:
            smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # default aspirin
        if 'random_name' in st.session_state:
            st.success(f"Random: **{st.session_state.random_name}**")

with col2:
    st.markdown("### 📜 History")
    if st.session_state.history:
        for idx, (smi, name) in enumerate(reversed(st.session_state.history[-5:])):
            st.write(f"{idx+1}. {name or smi[:10]}...")
        if st.button("Clear history"):
            st.session_state.history = []
    else:
        st.write("No molecules viewed yet.")

# ------------------- Process the molecule -----------------
if smiles.strip():
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        st.error("Invalid SMILES string. Please check and try again.")
    else:
        # ---- Get name ----
        name = get_molecule_name(smiles)
        if name:
            st.subheader(f"🧪 {name}")
        else:
            st.subheader("🧪 Unknown compound")

        # Add to history
        st.session_state.history.append((smiles, name))

        # ---- Display 2D structure ----
        img = mol_to_image(smiles)
        if img:
            col_img, col_props = st.columns([1, 1])
            with col_img:
                st.image(img, caption="2D Structure", use_container_width=True)
                # Download button
                buffered = BytesIO()
                img.save(buffered, format="PNG")
                img_str = base64.b64encode(buffered.getvalue()).decode()
                href = f'<a href="data:image/png;base64,{img_str}" download="molecule.png">⬇️ Download image</a>'
                st.markdown(href, unsafe_allow_html=True)
            with col_props:
                # ---- Properties ----
                props = compute_properties(smiles)
                if props:
                    st.markdown("### 📊 Properties")
                    cols = st.columns(2)
                    with cols[0]:
                        st.metric("MW", f"{props['Molecular Weight']} Da")
                        st.metric("LogP", props["LogP (cLogP)"])
                        st.metric("H‑bond Donors", props["H‑bond Donors"])
                    with cols[1]:
                        st.metric("H‑bond Acceptors", props["H‑bond Acceptors"])
                        st.metric("Rotatable Bonds", props["Rotatable Bonds"])
                        st.metric("TPSA", f"{props['TPSA']} Å²")
                    with st.expander("🔍 More descriptors"):
                        st.write(f"**Fraction Csp³**: {props['Fraction Csp3']}")
                        st.write(f"**Aromatic Rings**: {props['Aromatic Rings']}")
                        st.write(f"**Saturated Rings**: {props['Saturated Rings']}")
                        st.write(f"**Total Rings**: {props['Num Rings']}")
                        st.write(f"**Heavy Atoms**: {props['Heavy Atoms']}")

        # ---- Lipinski & fun summary ----
        st.markdown("### 💊 Drug‑likeness")
        passed, violations = lipinski_check(props)
        if passed:
            st.success("✅ **Passes Lipinski Rule of 5** (good oral bioavailability likely).")
        else:
            st.warning("❌ **Fails Lipinski** – violations:")
            for v in violations:
                st.write(f"- {v}")
        st.info(drug_likeness_summary(props))

        # ---- 3D viewer (optional) ----
        if has_3d:
            with st.expander("🧬 View 3D structure (interactive)"):
                show_3d_structure(smiles)
        else:
            with st.expander("🧬 3D structure"):
                st.warning("Install py3Dmol and stmol for 3D view.")

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
            }
            if name and name in facts:
                st.write(facts[name])
            else:
                st.write("This molecule is part of our teaching collection – explore its properties!")

        # ---- Quiz button ----
        if st.button("🧠 Quiz Me!"):
            st.balloons()
            st.success("Try to guess the molecule's name from its SMILES! (or ask your students)")

# Sidebar footer
st.sidebar.markdown("---")
st.sidebar.caption("Built with ❤️ by DetectED – Early Disease Detection & AI")