# Docking

Docking [package](https://datagrok.ai/help/develop/develop#packages) is a robust bioinformatics tool for the [Datagrok](https://datagrok.ai) platform. It is designed for conducting molecular docking simulations. Molecular docking is a crucial technique in computational biology and drug discovery. It allows researchers to predict the binding interactions between small molecules (ligands) and macromolecules (receptors). The primary purpose of the Docking package is to provide a user-friendly interface for performing these simulations efficiently and accurately.

The Docking package seamlessly integrates with [Autodock GPU](https://catalog.ngc.nvidia.com/orgs/hpc/containers/autodock), a widely-used software for molecular docking simulations. This integration ensures high-performance and reliable results in predicting the binding configurations of ligands to target macromolecules.

## How To

### Prepare data

Create or load a dataframe containig the ligands to be docked. For demonstration purposes, consider using the provided demo data located under **System:AppData/Docking/demo_files**.

### Prepare macromolecule

Utilize [AutoDock tools](https://www.chem.uwec.edu/chem491_w01/Chem491-Molecules%20and%20Medicine%202008/AutoDock%20Tutorial.pdf) to prepare the target protein and docking configurations files. Locate these files in a folder under **System:AppData/Docking/targets**.

### Run docking

1. Navigate to Chem > Autodock. A dialog appears.

2. Configure the parameters:

* `Ligands`: Specify the column within the provided dataframe that contains the small molecules to be docked.
* `Target`: Provide the folder containing the docking configurations and the macromolecule for simulation.
* `Conformations`: Define the number of conformations for the simulation.

When the results are ready, you can:

1. **View docking poses**: The added column containing docking poses allows users to visualize the spatial orientation of each ligand within the binding pocket of the receptor.
2. **Check free energy**: The column with free energy numbers provides users with quantitative information about the stability of each ligand-receptor complex.
3. **Zoom to pocket**: By clicking on a pose, a Molstar viewer appears in the context panel. It automatically zooms to the binding pocket of interest, providing a closer look at the ligand-receptor interaction.
4. **Explore additional properties**: Calculated properties for the selected pose allow users to explore binding affinities, interaction energies, or other relevant information for detailed analysis.