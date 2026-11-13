# Compound summary

A compound's side panel opens with a short read-only summary of it: **ID**, **Name** and **Structure**, captions on
the left and values on the right, as in a form. The structure is the platform's drawing of the molecule,
`grok.chem.drawMolecule(compound.smiles, 160, 100)`, an element of 160×100 pixels.

The plugin calls:

```ts
export function compoundSummary(compound: {id: string, name: string, smiles: string}): HTMLElement | Control
```
