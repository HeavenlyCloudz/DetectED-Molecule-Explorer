import gradio as gr
import requests
import random
import base64
from io import BytesIO
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Draw, AllChem, rdMolDescriptors
import py3Dmol

# ------------------- Helper functions --------------------

def mol_to_image(smiles, size=(400, 200)):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=size)
    return img

def get_molecule_name(smiles):
    """Query PubChem for the common name."""
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
    passed, violations = lipinski_check(props)
    if passed:
        return "✅ Likely good oral bioavailability (Lipinski compliant)."
    else:
        return f"⚠️ {len(violations)} Lipinski violation(s). May have absorption issues."

def get_3d_html(smiles):
    """Generate interactive 3D viewer HTML using py3Dmol."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "<p>Invalid SMILES</p>"
    # Generate 3D conformer
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_h, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol_h)
    pdb = Chem.MolToPDBBlock(mol_h)

    view = py3Dmol.view(width=400, height=300)
    view.addModel(pdb, "pdb")
    view.setStyle({"stick": {}, "sphere": {"scale": 0.2}})
    view.zoomTo()
    return view._repr_html_()  # returns HTML string for embedding

def random_molecule():
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
    }
    name, smi = random.choice(list(drugs.items()))
    return name, smi

# ------------------- Gradio App -------------------------

def process_smiles(smiles, history_state):
    """Main processing function called on any input change."""
    if not smiles.strip():
        return (None, None, "", "", "", "", "", history_state)
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, None, "Invalid SMILES", "", "", "", "", history_state)
    
    # Get name
    name = get_molecule_name(smiles)
    display_name = name if name else "Unknown compound"
    
    # Generate image
    img = mol_to_image(smiles)
    img_pil = img if img else None
    
    # Compute properties
    props = compute_properties(smiles)
    if not props:
        return (None, None, "Could not compute properties", "", "", "", "", history_state)
    
    # Lipinski
    passed, violations = lipinski_check(props)
    lipinski_text = "✅ Passes Lipinski Rule of 5" if passed else "❌ Fails Lipinski – violations: " + ", ".join(violations)
    summary = drug_likeness_summary(props)
    
    # Build property strings
    prop_text = (
        f"**MW:** {props['Molecular Weight']} Da\n"
        f"**LogP:** {props['LogP (cLogP)']}\n"
        f"**H‑bond Donors:** {props['H‑bond Donors']}\n"
        f"**H‑bond Acceptors:** {props['H‑bond Acceptors']}\n"
        f"**Rotatable Bonds:** {props['Rotatable Bonds']}\n"
        f"**TPSA:** {props['TPSA']} Å²\n"
        f"**Fraction Csp³:** {props['Fraction Csp3']}\n"
        f"**Aromatic Rings:** {props['Aromatic Rings']}\n"
        f"**Saturated Rings:** {props['Saturated Rings']}\n"
        f"**Total Rings:** {props['Num Rings']}\n"
        f"**Heavy Atoms:** {props['Heavy Atoms']}"
    )
    
    # 3D HTML
    try:
        html_3d = get_3d_html(smiles)
    except:
        html_3d = "<p>3D generation failed</p>"
    
    # Fun fact
    facts = {
        "Aspirin": "Aspirin was first synthesized in 1897 and is one of the most widely used medications.",
        "Ibuprofen": "Ibuprofen was discovered in the 1960s while researching for a new rheumatoid arthritis treatment.",
        "Caffeine": "Caffeine is the world's most widely consumed psychoactive substance.",
        "Paracetamol": "Paracetamol is also known as acetaminophen and is used for pain relief and fever.",
        "Penicillin G": "Penicillin was the first antibiotic discovered by Alexander Fleming in 1928.",
        "Dopamine": "Dopamine is a neurotransmitter that plays a role in reward and motor control.",
        "Serotonin": "Serotonin is a neurotransmitter that contributes to well-being and happiness.",
    }
    fun_fact = facts.get(name, "This molecule is part of our teaching collection – explore its properties!")
    
    # Update history
    history = history_state if history_state else []
    history.append((smiles, display_name))
    history_text = "\n".join([f"{i+1}. {name}" for i, (_, name) in enumerate(history[-5:])])
    
    return (img_pil, display_name, prop_text, lipinski_text, summary, html_3d, fun_fact, history)

def random_action(history_state):
    name, smi = random_molecule()
    return smi, history_state

# Gradio UI
with gr.Blocks(title="DetectED Molecule Explorer", theme=gr.themes.Soft()) as demo:
    gr.Markdown("# 🧬 DetectED Molecule Explorer")
    gr.Markdown("*For teaching cheminformatics – early disease detection & AI*")
    
    # State for history
    history_state = gr.State([])
    
    with gr.Row():
        with gr.Column(scale=2):
            input_method = gr.Radio(
                choices=["Select from examples", "Type SMILES", "Random molecule"],
                value="Select from examples",
                label="Choose input method"
            )
            example_dropdown = gr.Dropdown(
                choices=[
                    "Aspirin", "Ibuprofen", "Paracetamol", "Caffeine",
                    "Penicillin G", "Dopamine", "Serotonin", "Glucose", "Ethanol"
                ],
                value="Aspirin",
                label="Select a molecule",
                visible=True
            )
            smiles_input = gr.Textbox(
                label="Enter SMILES string",
                value="CC(=O)OC1=CC=CC=C1C(=O)O",
                visible=False
            )
            random_button = gr.Button("🎲 Give me a random molecule", visible=False)
            # Hidden textbox to hold the current SMILES (we'll update it)
            current_smiles = gr.Textbox(value="CC(=O)OC1=CC=CC=C1C(=O)O", label="SMILES", visible=False)
        
        with gr.Column(scale=1):
            gr.Markdown("### 📜 History")
            history_display = gr.Textbox(label="", lines=5, interactive=False, placeholder="No molecules viewed yet.")
            clear_history = gr.Button("Clear history")
    
    # Outputs
    with gr.Row():
        with gr.Column(scale=1):
            molecule_name = gr.Textbox(label="Molecule Name", interactive=False)
            mol_image = gr.Image(label="2D Structure", type="pil")
        with gr.Column(scale=1):
            properties = gr.Markdown("### 📊 Properties\n")
            lipinski = gr.Markdown("### 💊 Drug‑likeness\n")
            summary = gr.Markdown()
    
    with gr.Row():
        with gr.Column():
            gr.Markdown("### 🧬 3D Structure (interactive)")
            mol_3d = gr.HTML()
        with gr.Column():
            fun_fact = gr.Markdown("### 📖 Did you know?\n")
            quiz_button = gr.Button("🧠 Quiz Me!")
            quiz_output = gr.Markdown()
    
    # Update visibility based on input method
    def update_visibility(method):
        return {
            example_dropdown: gr.update(visible=(method == "Select from examples")),
            smiles_input: gr.update(visible=(method == "Type SMILES")),
            random_button: gr.update(visible=(method == "Random molecule"))
        }
    input_method.change(update_visibility, inputs=input_method, outputs=[example_dropdown, smiles_input, random_button])
    
    # When example changes
    def set_smiles_from_example(choice):
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
        return examples[choice]
    example_dropdown.change(set_smiles_from_example, inputs=example_dropdown, outputs=current_smiles)
    
    # When SMILES typed
    def set_smiles_from_text(text):
        return text
    smiles_input.change(set_smiles_from_text, inputs=smiles_input, outputs=current_smiles)
    
    # Random button
    def random_and_set(history):
        name, smi = random_molecule()
        return smi, history
    random_button.click(random_and_set, inputs=[history_state], outputs=[current_smiles, history_state])
    
    # Process when current_smiles changes
    current_smiles.change(
        process_smiles,
        inputs=[current_smiles, history_state],
        outputs=[mol_image, molecule_name, properties, lipinski, summary, mol_3d, fun_fact, history_display]
    )
    
    # Clear history
    def clear_history_fn():
        return [], ""
    clear_history.click(clear_history_fn, outputs=[history_state, history_display])
    
    # Quiz button
    def quiz_response():
        return "🎉 Try to guess the molecule's name from its SMILES! (or ask your students)"
    quiz_button.click(quiz_response, outputs=quiz_output)
    
    # Initial load
    demo.load(
        process_smiles,
        inputs=[current_smiles, history_state],
        outputs=[mol_image, molecule_name, properties, lipinski, summary, mol_3d, fun_fact, history_display]
    )

if __name__ == "__main__":
    demo.launch()
