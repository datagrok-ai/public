# Sketcher

**Apps** > **Chem** > **Sketcher** draws one molecule or query in the Crux sketcher, which fills the view. It opens
in molecule mode. For the query tools (atom lists, query bonds, generic atoms and the query sections of Atom
Properties), turn on **Query mode** in the sketcher's gear menu; a link to a query opens in query mode.

The context panel opens beside the sketcher. Whatever you draw is the current object, so the panel shows the
molecule's panes (its structure, properties, identifiers and the rest) and follows every change. The help panel
closes to give the sketcher room; opened again, it shows this page.

The app always shows Crux. The sketcher you picked for other molecule editors stays as it is.

The app opens on an example, atorvastatin, unless its address names a molecule.

## The link

The page address is a link to the drawing: copy it from the browser's address bar and send it, and the app opens on
the same molecule. It is the app's address with one parameter, `smiles`:

```
/apps/Chem/Sketcher?smiles=CC(C)c1c(C(%3DO)Nc2ccccc2)...
```

- A molecule goes into the link as CXSMILES without coordinates: SMILES, plus what only CXSMILES holds, such as
  enhanced stereo groups, R labels and attachment points. The drawing's layout is not part of it: the app lays the
  molecule out again when it opens the link.
- A drawing with query features that SMILES cannot hold (an atom list, a query bond, a custom query) goes in as
  CXSMARTS: SMARTS, plus what only the CX part holds, such as enhanced stereo groups, attachment points and link
  nodes, so a shared query reopens with them. A query whose drawing writes another SMARTS through the CX route (one
  with R labels or explicit hydrogens drawn, for example) goes in as plain SMARTS. Atom aliases are not kept in a
  query's link.
- An empty drawing gives an empty parameter (`?smiles=`), which opens an empty sketcher. The address without the
  parameter opens the example.
- The address changes as you draw; it does not add a step to the browser's history, so **Back** leaves the app.
- A drawing too large for an address (over 6,000 characters in the link) is not written into it: the app says so
  once, and the address shows the app alone. Use **Export...** to share such a molecule.

You can also write a link yourself, with any SMILES, CXSMILES or SMARTS, URL-encoded. A text the sketcher cannot read
opens an empty sketcher with the molecule shown as malformed.

## Export

**Export...** in the ribbon opens the Export dialog. It previews what it saves, fitted to the dialog: resize the
dialog to see it larger.

- **Style**, for SVG, PNG and CDXML:
  - "ACS Document 1996": the molecule laid out afresh, as Crux lays out a molecule read from its text, and drawn in
    the ACS 1996 style: black, bonds 0.508 cm (14.4 pt), lines 0.6 pt, labels 10 pt. The SVG is sized in points, so a
    document places it at print size.
  - "As drawn": the molecule as you drew it on the canvas, in its coordinates, with the sketcher's look (its atom
    colours, its bond length).
- **Format**: SVG, PNG at 300 dpi, CDXML for ChemDraw (its ACS Document 1996 document style, in either layout),
  MOL V2000, MOL V3000, SDF, SMILES, CXSMILES, SMARTS, InChI and InChIKey. The MOL and SDF files keep the drawing's
  coordinates. A format the drawing has no text in says why: a query has no SMILES and no InChI, and InChI does not
  hold R labels. "MOL V2000" is written as V3000 when V2000 cannot hold the molecule, and the dialog says so.
- **Transparent background**: for SVG and PNG.

**Save** downloads the file. **Copy** puts it on the clipboard: the PNG as a picture, every other format as text.
