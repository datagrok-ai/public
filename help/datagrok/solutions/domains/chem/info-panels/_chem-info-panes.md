---
title: Info panes for chemical data
---

[Overview: info panes PLACEHOLDER]

The info panes are displayed automatically in the **Context Panel** when you select a molecule in a dataset.

:::tip

To add info pane as a column, hover over the first descriptor shown, click the **more actions** icon > **Add as column**. This adds a column with

:::

:::note developers

You can [create custom info panes](../../../../../develop/how-to/add-info-panel.md).

:::

## Actions

The following actions are available from the **Actions** info pane:

* **Add to variables**: Adds a cell value as a variable which can be then be used in [scripts](../../../../../compute/scripting.md)
* **Copy value**: Copies a cell value
* **Add to favorites**: [To be deleted????]
* **Sketch**: Opens a [sketcher](../chem.md)

## Chemistry

The following info panes and actions are available form the **Chemistry** info pane:

* **Gasteiger Partial Charges**: [RDKit-based script](https://dev.datagrok.ai/script/7acf813d-4f65-51f2-bc3f-503cde26c460)
* **Descriptors**: The following descriptors are shown by default:
  * FractionCSP3
  * HeavyAtomCount
  * NHOHCount

  To calculate other descriptors, click the **Select** button and choose from the list. The info pane updates accordingly.
* **Properties**: Shows molecular properties, such as formula, molecular weight, LogP, and others.

## Biology

The following info panes and actions are available form the **Biology** info pane:

* **Drug Likeness**:
  
  <details>
  <summary> Drug likeness calculation details </summary>
  
  The approach for druglikeness calculation involves a list of substructure fragments with associated druglikeness scores, which are used to determine the druglikeness score of a molecule by summing up the scores of the fragments present in the molecule and dividing it by the square root of the number of fragments. The frequency of each fragment in traded drugs and non-drug-like compounds is used to eliminate redundant fragments, and then the druglikeness score is calculated by comparing the frequency of the fragment in traded drugs to non-drug-like compounds. 

  `d = sum(Vi)/sqrt(n)`

  You can find more information about this method on the [Open Molecules website](https://openmolecules.org/properties/properties.html#druglikeness).

  </details>

* **Structural Alerts**:
* **Toxicity**:

## Structure

## Databases
