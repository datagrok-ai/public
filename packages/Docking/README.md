# Docking

Docking [package](https://datagrok.ai/help/develop/develop#packages) is a robust bioinformatics tool for the [Datagrok](https://datagrok.ai) platform. It is designed for conducting molecular docking simulations. Molecular docking is a crucial technique in computational biology and drug discovery. It allows researchers to predict the binding interactions between small molecules (ligands) and macromolecules (receptors). The primary purpose of the Docking package is to provide a user-friendly interface for performing these simulations efficiently and accurately.

## Features

1. **Autodock GPU integration**

The Docking package seamlessly integrates with [Autodock GPU](https://catalog.ngc.nvidia.com/orgs/hpc/containers/autodock), a widely-used software for molecular docking simulations. This integration ensures high-performance and reliable results in predicting the binding configurations of ligands to target macromolecules.

2. **Simplified workflow**

Users can easily provide input data, specify ligands, and define the target folder containing docking configurations and macromolecules. The package offers flexibility through adjustable parameters such as the number of docking confirmations, allowing researchers to customize simulations based on their specific experimental requirements.

3. **Analysis**

To run the analysis go to **Chem > Autodock**. A dialog opens.

* `Table`: A dataframe where each row represents a unique molecular configuration.
* `Ligands`: Column within the provided dataframe that contains small molecules to be docked.
* `Target`: Folder that contains docking configurations and the macromolecule for simulation. This folder is located under *System:AppData/Docking/targets*.
* `Conformations`: Number of conformations.