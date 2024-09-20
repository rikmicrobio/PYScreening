# Re-defining the readme_content variable without emojis to make it simpler and avoid encoding issues.

readme_content_simple = """
# PYScreening
## 1. Introduction
Welcome to the repository for our latest project on virtual screening of small molecules using **Lipinski's Rule of Five**. This project aims to develop a comprehensive pipeline for filtering large sets of molecules based on **Lipinski's Rule**, a well-established guideline for predicting the oral bioavailability of drug candidates. The workflow includes:

- Reading SMILES strings from input files
- Converting SMILES to RDKit molecules
- Applying Lipinski's Rule of Five
- Saving filtered molecules in CSV and SDF formats

This repository includes the complete source code, documentation, and example datasets to reproduce the results. We believe this tool will be valuable for researchers in the fields of **drug discovery** and **computational chemistry**.

![pic_info](https://user-images.githubusercontent.com/45164491/213928955-f94c8260-fd60-457d-be9d-500662abe62e.png)

*Picture Courtesy: Chagas, C. M., Moss, S., & Alisaraie, L. (2018). Drug metabolites and their effects on the development of adverse reactions: Revisiting Lipinski’s Rule of Five. International Journal of Pharmaceutics, 549(1-2), 133-149.*

---

## 2. Python Code Breakdown (Step 1: Screening)
This Python script performs virtual screening of small molecules based on **Lipinski’s Rule of Five**. Below are the key steps involved:

1. **Importing Required Modules**
   - The script imports necessary modules like `csv`, `rdkit`, and functions from `Lipinski` and `Descriptors`.

2. **Reading SMILES Strings**
   - The script reads a list of **SMILES strings** from a CSV file.

3. **Converting SMILES to RDKit Molecules**
   - Converts SMILES strings to **RDKit molecule objects** using `Chem.MolFromSmiles()`.

4. **Applying Lipinski's Rule of Five**
   - Filters molecules based on hydrogen bond donors, acceptors, molecular weight, and logP.

5. **Saving Filtered Molecules to CSV**
   - Writes the filtered molecules to a **CSV file** using `Chem.MolToSmiles()`.

6. **Final Output**
   - Prints the number of molecules that passed the filter.

---

## 3. Python Code Breakdown (Step 2: Saving as SDF)
This Python script generates **3D structures** of the filtered molecules and saves them in **SDF format**.

1. **Importing OS Module**
   - Uses the `os` module to handle file operations and create directories.

2. **Creating Output Folder**
   - Creates a folder named `Screened_Molecules_3D` to store **SDF files**.

3. **Saving Molecules in SDF Format**
   - Converts molecules to **SDF format** using `Chem.MolToMolBlock()`.

4. **Storing SDF Files**
   - Stores all filtered molecules in SDF format for future use.

---

## 4. Contact Information

For more details, please reach out:

- **Email**: rikgangulybioinfo@gmail.com
- **Lab**: Computational Biology Laboratory,  
  North-Eastern Hill University, Shillong, India

---

## 5. License
This project is licensed under the MIT License. The MIT License allows:

- Free usage, modification, and distribution for any purpose.
- The project is provided "as is," without warranty of any kind.
- You must retain the original copyright notice in any copies.

Please see the LICENSE file for more details.
"""

# Save the simplified version as a Readme.md file
file_path = "/mnt/data/Readme_simple.md"
with open(file_path, "w") as f:
    f.write(readme_content_simple)

file_path
