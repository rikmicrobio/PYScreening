
# PYScreening ðŸŽ¯
## 1. **Introduction** ðŸ“˜
### Welcome to the repository for our latest project on virtual screening of small molecules using **Lipinski's Rule of Five**. This project aims to develop a comprehensive pipeline for filtering large sets of molecules based on **Lipinski's Rule**, a well-established guideline for predicting the oral bioavailability of drug candidates. The workflow includes:

- ðŸ“„ **Reading SMILES strings** from input files
- ðŸ§¬ **Converting SMILES to RDKit molecules**
- âš™ï¸ **Applying Lipinski's Rule of Five**
- ðŸ’¾ **Saving filtered molecules** in CSV and SDF formats

This repository includes the complete source code, documentation, and example datasets to reproduce the results. We believe this tool will be valuable for researchers in the fields of **drug discovery** and **computational chemistry**.

![pic_info](https://user-images.githubusercontent.com/45164491/213928955-f94c8260-fd60-457d-be9d-500662abe62e.png)

*Picture Courtesy: Chagas, C. M., Moss, S., & Alisaraie, L. (2018). Drug metabolites and their effects on the development of adverse reactions: Revisiting Lipinskiâ€™s Rule of Five. International Journal of Pharmaceutics, 549(1-2), 133-149.*

---

## 2. **Python Code Breakdown** ðŸ§‘â€ðŸ’» (Step 1: Screening)
### This Python script performs virtual screening of small molecules based on **Lipinskiâ€™s Rule of Five**. Below are the key steps involved:

### 1ï¸âƒ£ **Importing Required Modules** ðŸ“¦
We begin by importing necessary modules:
- `csv` for handling CSV files
- `rdkit` for working with chemical structures
- `Lipinski` and `Descriptors` from `rdkit.Chem` for computing molecular properties

### 2ï¸âƒ£ **Reading SMILES Strings** ðŸ”¬
- The script reads a list of **SMILES strings** from a CSV file using `csv.reader()`.
- These strings represent the molecular structures to be screened.

### 3ï¸âƒ£ **SMILES to RDKit Molecules** âš—ï¸
- Converts SMILES strings to **RDKit molecule objects** using `Chem.MolFromSmiles()`.
- Invalid SMILES strings are discarded.

### 4ï¸âƒ£ **Applying Lipinski's Rule of Five** ðŸ§®
- Filters molecules based on the following criteria:
  - â‰¤ 5 hydrogen bond donors
  - â‰¤ 10 hydrogen bond acceptors
  - Molecular weight â‰¤ 500 Daltons
  - LogP â‰¤ 5

- **Lipinski's Rule** helps in predicting the oral bioavailability of drug candidates. 

### 5ï¸âƒ£ **Saving Filtered Molecules to CSV** ðŸ’¾
- Writes the filtered molecules to a **CSV file** using `Chem.MolToSmiles()`.
- Provides a concise record of molecules that passed the screening.

### 6ï¸âƒ£ **Final Output** ðŸ§¾
- Prints the total number of molecules that passed the filter.

> ðŸ’¡ This script is ideal for quickly screening large datasets of molecules, filtering out candidates that donâ€™t meet **Lipinskiâ€™s Rule** for oral bioavailability.

---

## 3. **Python Code Breakdown** ðŸ§‘â€ðŸ’» (Step 2: Saving as SDF)
### This Python script generates **3D structures** of the filtered molecules and saves them in **SDF format** for further analysis. 

### 1ï¸âƒ£ **Importing OS Module** ðŸ—‚ï¸
- The `os` module is used to create directories and handle file operations.

### 2ï¸âƒ£ **Creating Output Folder** ðŸ“
- A folder named `Screened_Molecules_3D` is created to store the output **SDF files**.

### 3ï¸âƒ£ **Saving Molecules in SDF Format** ðŸ’½
- The filtered molecules are converted into **SDF format** using `Chem.MolToMolBlock()`.
- Files are named sequentially (`mol_1.sdf`, `mol_2.sdf`, etc.) and saved in the created folder.

### 4ï¸âƒ£ **Storing SDF Files** ðŸ“‚
- Ensures all the filtered molecules are stored properly for future use in 3D visualization, **3D-QSAR**, or other **molecular modeling** techniques.

> ðŸŽ¯ **Utility**: This script provides an easy way to convert filtered molecules into **3D SDF format**, which is crucial for many downstream applications in drug discovery.

---

## 4. **Contact Information** ðŸ“§

For more details, please reach out:

- **Email**: [rikgangulybioinfo@gmail.com](mailto:rikgangulybioinfo@gmail.com)
- **Lab**: Computational Biology Laboratory,  
  North-Eastern Hill University, Shillong, India

---

## 5. **License** ðŸ“œ
This project is licensed under the MIT License. The MIT License allows:

- **Permission to Use**: Free usage, modification, and distribution for any purpose, personal or commercial.
- **No Liability**: The project is provided "as is," without warranty of any kind, express or implied.
- **Attribution**: You must retain the original copyright notice in any copies or substantial portions of the software.

Please see the LICENSE file for more details.

