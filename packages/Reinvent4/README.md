# Reinvent4

The [Reinvent4](https://datagrok.ai/help/develop/develop#packages) is a plugin that seamlessly integrates the [Reinvent4](https://github.com/MolecularAI/REINVENT4) molecular design tool with the [Datagrok](https://datagrok.ai) platform.  

[Reinvent4](https://github.com/MolecularAI/REINVENT4) is a powerful tool for **de novo molecular design**, **scaffold hopping**, **R-group replacement**, **linker design**, **molecule optimization**, and other small molecule design tasks.  

## Getting Started  

### Configuring Reinvent4  

The plugin includes several pre-configured [configuration files](https://github.com/datagrok-ai/public/tree/master/packages/Reinvent4/files/reinvent). If you need a custom configuration, follow the official [Reinvent4 instructions](https://github.com/MolecularAI/REINVENT4/blob/main/configs/toml/scoring_components_example.toml).  

Once your configuration is ready, place the files in a folder under **System:AppData/Reinvent4/reinvent**.

The folder name will automatically appear as the **configuration name** in the Datagrok plugin UI.  

### Generating Molecules  

1. Navigate to **Chem > Generate molecules...**.  
2. In the dialog, configure the following parameters:  
   - **Ligand seed**: The starting point for ligand generation.  
   - **Optimize**: The optimization criteria (folder where all configuration files will be stored).  
3. Click **Run** to start the molecular generation process.  