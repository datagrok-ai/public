---
title: "Getting started"
sidebar_position: 0
format: mdx
---

```mdx-code-block
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';
import BrowserWindow from '@site/src/components/browser-window';
```

Setting up docking configurations can be challenging, but we've streamlined this process.
In this section, we'll guide you through each key concept, from preparing configurations
to storing them and beyond.

### Prepare macromolecule (target)

<div style = {{ display: 'flex' }}>

<div style = {{ width: '60%' }}>

Autodock plugin contains several pre-created [targets](https://github.com/datagrok-ai/public/tree/master/packages/Docking/files/targets).
To add your own to the list:

* Prepare the macromolecule and AutoDock grid parameter file using [AutoDock tools](https://ccsb.scripps.edu/mgltools/downloads/).
* Ensure that the macromolecule is in the `PDBQT` format.

For detailed instructions, review the Autodock tutorials
([first](https://www.chem.uwec.edu/chem491_w01/Chem491-Molecules%20and%20Medicine%202008/AutoDock%20Tutorial.pdf),
[second](https://omicstutorials.com/a-comprehensive-bioinformatics-tutorial-mastering-ligand-protein-docking-with-autodock/)).

Put these files in a folder under **System:AppData/Docking/targets**.
The name of the folder will appear as the target name in the Datagrok plugin UI.

</div>

<div style = {{ width: '40%', display: 'flex', gap: '20px', 'justify-content': 'end' }}>

```mdx-code-block
<BrowserWindow url="" bodyStyle={{'padding': '0px'}}>
```

<img
  src={require('./_pics/targets-folder.png').default}
  style= {{ height: '160px', 'border-radius': '5px'}}
/>

```mdx-code-block
</BrowserWindow>
```

</div>

</div>

:::tip Atomic maps

The AutoDock calculation fails if for a particular ligand corresponding atomic maps not available.
To run docking on a big ligand dataset, we suggest including all available atomic map types:
`A C HD N NA OA SA CL`.

:::

### Prepare data

Create or load a dataframe that contains ligands to be docked. 
For demonstration purposes, consider using the provided [demo data](https://github.com/datagrok-ai/public/blob/master/packages/Docking/files/demo_files/demo_dataset.csv).

### Run docking simulations

* Navigate to **Chem > Autodock**. A dialog appears.
* Configure the parameters:
  * **Ligands**: Specify the column within the provided dataframe that contains the small molecules to be docked.
  * **Target**: Choose the folder containing the docking configurations and the macromolecule for simulation.
  * **Conformations**: Define the number of conformations to search.
* Run the calculations

```mdx-code-block
<BrowserWindow url='' bodyStyle={{'padding': '0px'}}>
```

![docking-simulations](./_pics/docking-simulations.gif)

```mdx-code-block
</BrowserWindow>
```


:::tip Performance note

During the first run for a target Autodock calculates the macromolecule grids.
Grid calculation works on a single CPU core, so it can take about a minute. 
Datagrok caches calculated grids, so subsequent runs with the same target will be much faster.

:::