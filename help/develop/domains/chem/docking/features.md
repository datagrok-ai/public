---
title: "Features"
sidebar_position: 0
format: mdx
---

```mdx-code-block
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';
import BrowserWindow from '@site/src/components/browser-window';
```

Upon completing your docking simulations (please refer to the [getting started](./getting-started) section for details),
a range of functionalities becomes available. In this section, we will guide you through all of them.

### View docking results

<div style = {{ display: 'flex' }}>

<div style = {{ width: '50%' }}>

Upon completion of docking simulations, users can review the results comprising **docking poses** and **binding energy** values.

The column displaying binding energy provides quantitative insights into the stability of ligand-receptor complexes, where
lower values indicate stronger binding.

This column is color-coded, ranging from green for lower values to red for higher values,
facilitating straightforward interpretation.

</div>

<div style = {{ width: '50%', display: 'flex', gap: '20px', 'justify-content': 'end' }}>

```mdx-code-block
<BrowserWindow url="" bodyStyle={{'padding': '0px'}}>
```

<img
  src={require('./_pics/docking-results.png').default}
  style= {{ height: '200px', 'border-radius': '5px'}}
/>

```mdx-code-block
</BrowserWindow>
```

</div>

</div>

### Zoom into binding pocket

Clicking on a pose activates a Molstar viewer in the context panel, automatically zooming to the binding pocket.
This allows for a closer examination of the molecular interactions within the pocket.

```mdx-code-block
<BrowserWindow url='' bodyStyle={{'padding': '0px'}}>
```

![docking-simulations](./_pics/zoom-into-pocket.gif)

```mdx-code-block
</BrowserWindow>
```

### Explore additional properties

In the docking panel, we display additional properties returned by Autodock, such as interaction energies and binding affinities.
These properties can be added to the entire dataframe by simply clicking on the plus icon.

```mdx-code-block
<BrowserWindow url='' bodyStyle={{'padding': '0px'}}>
```

![docking-simulations](./_pics/additional-properties.gif)

```mdx-code-block
</BrowserWindow>
```

### Download pose in PDB and CIF formats

<div style = {{ display: 'flex' }}>

<div style = {{ width: '50%' }}>

To download the pose in the context of the protein, simply right click to open the context menu.
Within this menu, you'll find a **Download** option under which you can select either **as CIF** or **as PDB**. 

</div>

<div style = {{ width: '50%', display: 'flex', gap: '20px', 'justify-content': 'end' }}>

```mdx-code-block
<BrowserWindow url="" bodyStyle={{'padding': '0px'}}>
```

<img
  src={require('./_pics/molstar-download.gif').default}
  style= {{ height: '350px', 'width': '', 'border-radius': '5px'}}
/>

```mdx-code-block
</BrowserWindow>
```

</div>

</div>