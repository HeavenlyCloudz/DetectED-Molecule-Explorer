# 🧬 DetectED Molecule Explorer

An interactive cheminformatics teaching app built with **RDKit**, **Streamlit**, and **Ketcher** for **DetectED Labs**, an extension of DetectED. Our goal is to make early disease detection and AI education accessible, this time through hands-on molecular exploration!

---

## ✨ Features

- **Multiple Input Methods** – Choose from three ways to explore molecules:
  - ✏️ **Draw Structure** – Interactive Ketcher molecular editor for drawing molecules visually
  - ⌨️ **Type SMILES** – Enter any SMILES string manually
  - 📚 **Examples** – Pick from curated examples (Aspirin, Ibuprofen, Caffeine, Morphine, Nicotine, etc.)
  - 🎲 **Random Molecule** – Surprise button for exploration and classroom engagement
  - 💊 **Quick Actions** – Sidebar shortcuts for common molecules

- **Auto‑Name Lookup** – Fetches the common drug name and comprehensive data from PubChem (when available)

- **2D Structure Viewer** – Displays high‑quality 2D molecular images with download option

- **3D Interactive Viewer** – Rotate, zoom, and explore molecules in 3D with:
  - Full rotation and zoom controls
  - Spin on/off toggle
  - Reset view button
  - Spectrum color representation

- **Comprehensive Molecular Properties** – Computes key descriptors:
  - Molecular Weight, LogP (cLogP), TPSA
  - H‑bond Donors/Acceptors, Rotatable Bonds
  - Fraction Csp³, Aromatic/Saturated/Total Rings
  - Heavy Atoms, Total Atoms, Total Bonds
  - Detailed tooltips explaining each property

- **PubChem Integration** – Retrieves additional information:
  - Molecular Formula
  - IUPAC Name
  - InChI and InChIKey
  - Verified Molecular Weight

- **Lipinski Rule of 5** – Automatic pass/fail check with detailed violation breakdown and educational explanation

- **Drug‑Likeness Assessment** – Plain‑language interpretation of bioavailability potential with visual status badges

- **Similarity Search** – Compare molecules against a database of known drugs using Tanimoto similarity

- **Viewing History** – Tracks recently viewed molecules for easy revisiting with quick-click navigation

- **Fun Facts** – Historical and biological context for common drugs

- **Interactive Quiz Mode** – Test knowledge with:
  - Multiple-choice style questions
  - Random drug facts generator
  - Educational reinforcement

- **Export Options** – Download molecular data in multiple formats:
  - 📋 SMILES (.smi)
  - 📊 Properties as CSV
  - 🧪 SDF (Structure Data File)
  - 🖼️ 2D Structure as PNG

- **Comprehensive Documentation** – Built-in help system for:
  - SMILES writing guide with examples
  - Property explanations with tooltips
  - Lipinski's Rule of 5 educational content

- **Dark Theme Optimized** – Professional dark interface with:
  - White text for readability
  - Black text on white backgrounds for metrics
  - Color-coded status indicators

---

## ▶️ How to Use
Input Method
  - Draw Structure – Use the Ketcher editor to draw your molecule visually
  - Type SMILES – Enter a SMILES string manually
  - Examples – Select from 10+ common drug molecules
  - Quick Actions – Use sidebar buttons for instant access
  - Random – Get a surprise molecule

Exploration
  - View the 2D structure and automatically fetched name
  - Explore molecular properties with detailed tooltips
  - Check Lipinski compliance and drug-likeness
  - Toggle the 3D viewer for interactive exploration
  - Compare against similar drugs in our database

Export & Share
  - Download 2D structure as PNG
  - Export SMILES, SDF, or properties CSV
  - Share SMILES strings with colleagues

Educational Tools
  - Test yourself with the interactive quiz
  - Learn from fun facts and historical context
  - Explore property explanations with tooltips
  - Review viewing history for continued learning

---

## 📚 Teaching Applications
This app is designed for classroom and workshop settings to teach:
- Cheminformatics Basics – SMILES notation, molecular representations, structure drawing
- Drug Discovery Principles – Lipinski's Rule of 5, drug‑likeness assessment
- Computational Chemistry – Descriptors and property prediction
- Molecular Visualization – 2D and 3D structure exploration
- Data Integration – Working with PubChem API and chemical databases
- Critical Thinking – Interpreting property data for drug design

---

## 🎯 Classroom Activities
- Drug Discovery Workshop – Use Lipinski's Rule to evaluate hypothetical drug candidates
- Structure-Activity Relationships – Compare similar molecules and their properties
- 3D Visualization – Explore molecular geometry and conformation
- Quiz Challenges – Test student knowledge with built-in quiz mode
- Research Projects – Export data for further analysis

---

## 🙏 Acknowledgements
- RDKit – Open‑source cheminformatics toolkit
- Streamlit – Rapid web app development framework
- Ketcher – Professional molecular structure editor
- py3Dmol – 3D molecular visualization
- PubChem – Chemical information database API
- Pandas – Data handling and export

