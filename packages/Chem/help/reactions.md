# Reactions

Chem draws, edits and exports chemical reactions and multi-step routes in three places: the Reaction Sketcher
app, the Reaction panel of a reaction cell, and a reaction cell's editor.

The Reaction panel needs Chem's **Chemical Engine** setting set to Crux. The Reaction Sketcher, a reaction cell's
editor and its **Export reaction...** work with either setting.

Two example tables come with Chem: open **Files** > **App Data** > **Chem**, then `test-reactions.csv` (reactions and
routes as reaction SMILES) or `reactions-cdxml.csv` (ChemDraw pages). A ChemDraw page is fitted to its cell. A page
with more than its reaction, such as a key of reagents beside it, shows its reaction in the cell, with the texts beside
it; the editor and the export keep the whole page. A cell too small to show it readably says so: double-click it to
open the page in the editor. A page's reaction is
made of its steps: a molecule on no step is left out of it, and the Reaction panel and **Export reaction...** say
how many were left out.

## The Reaction Sketcher

Open it from **Apps** > **Chem** > **Reactions** > **Reaction Sketcher**. It holds one reaction or one route. It
opens with the context panel and the help panel closed, so that the drawing has the room; open them again as usual.

On the left, the sketcher draws the reaction. Use its arrow tool for the arrow. **Text**, under the sketcher, holds
the same reaction as text, and the two follow each other:

- draw, and Text is written again;
- type in Text, and the drawing follows when you pause or leave the box. The sketcher's **Undo** takes back what the
  text changed in the drawing.

Between the sketcher and Text, **Reactant**, **Reagent** and **Product** add a molecule to the drawing in that place:
type or paste its SMILES or a molecule, then press Enter or **+**, or click **Draw...** to draw it (for **Next step**,
**Draw...** opens on the route's last product, ready to change into the next one). **Reagent** also
takes a reagent's name, such as HATU or DMF, and suggests names in a list as you type; the name is added as its
structure. A
reactant goes after the reactants, a reagent over the arrow, a product after the arrow. **Next step** builds a route:
it draws an arrow after the last product, which becomes the new step's reactant, and adds the molecule as that step's
product. Each addition is one edit, which Undo takes back.

A reaction or route opened from text, with no coordinates of its own, is drawn arranged: a long route wraps into rows.

Text reads:

- reaction SMILES or SMARTS (`reactants>agents>products`), with a CX part too;
- a route: `A>>B>>C`, or steps separated by `--**--`, Chem's route notation. The sketcher's **Copy** > **SMILES**
  writes a route this way, and reaction cells read it;
- an RXN block or an RD file.

A molecule with no arrow is not a reaction yet: draw an arrow, or type `reactants>>products`.

Under Text, a line says what is read ("One reaction.", "A route of 5 steps.") or what is wrong. An orange line names
what Text does not keep of the drawing. For example, a reagent drawn as a label such as HATU is written as `*` in
reaction SMILES. Draw its structure to keep it. Conditions written over the arrow are not kept either.

**Map numbers**, at the top of the drawing on the right, shows or hides the atom map numbers in its drawings. It is
off until you switch it on.

The ribbon has four commands:

- **Examples**: nine reactions and routes to start from. An example's title stays while you edit its reaction; it
  goes when you clear the drawing or put a different reaction in its place.
- **Reset to example**: the default route again; your reaction is lost.
- **Export...**: the Export panel (see [Export](#export)).
- **Add to table**: puts the reaction into a table, as text in a reaction column. With a table open that has a
  reaction column, it adds a row to that table. Otherwise it makes a new table, Reactions, with the reaction and its
  title, and opens it.

The app keeps your reaction in the browser and opens on it next time. When the drawing holds what a reaction text
cannot (a catalytic cycle, texts, arrow styles, a label such as HATU), the drawing is kept too, and the app opens on
it as you drew it; typing in Text makes the text the reaction again. If the kept drawing cannot be read within a few
seconds, the sketcher draws the text instead, and the line under Text says so. If the sketcher cannot open the kept text
(for example, a route whose step has no product), it says so and why over the empty drawing. The text stays the
reaction: it is drawn on the right and exported as it is. Correct it in Text, or draw a new reaction in its place.

While you drag in the sketcher, the drawing on the right keeps its picture, and it is drawn again when you let go.

## The Reaction panel

The panel shows for a reaction cell in the context panel, and on the right of the Reaction Sketcher:

- the reaction drawn large. **Full screen** draws it as large as the window allows, its legend under it, with a
  route's steps in a list of their own under it; close it with its × or Escape;
- a route's **Steps**, one by one;
- a template's **Template notes**: what no picture draws;
- **Mapping**: the reaction's atoms are mapped automatically, a step at a time. A reaction that already has map
  numbers keeps them, and a template is not mapped. The reaction centre is coloured from the map numbers: green for a
  bond made, red for a bond broken, blue for a bond whose order changed. A legend under the picture names the colours
  it uses.

## Editing a reaction cell

Double-click a reaction cell, or right-click it and choose **Edit reaction...**, to edit it in the sketcher. This
works with either Chemical Engine setting. **OK** writes the reaction back in the cell's own notation. The full-screen icon in the dialog's title bar, or a double-click
on the title, makes the dialog fill the window; click it again to go back.
The sketcher fills the dialog as you resize it. A reaction with no coordinates of its own opens arranged. OK without a change writes nothing, and **Cancel** never
writes. A cell that holds a ChemDraw page (CDXML) is written back as CDXML, with its texts.

Before you click OK, a line at the bottom of the dialog says what the cell's notation will not keep of the drawing,
such as a label written as `*` or the texts over the arrow. A red line means OK cannot write the reaction; the line
says why. A reaction the editor cannot open, such as a route that is not one chain of steps, is said so in place of the
sketcher, and OK stays off.

## Export

**Export...** opens the Export panel. It is in the Reaction Sketcher's ribbon and on the Reaction panel; a reaction cell's
menu has it as **Export reaction...**.

The panel previews the scheme, fitted to the panel: resize the panel to see it larger. What it saves is at print size:

- **Style**: "ACS Document 1996" (black, Arial 10 pt, bonds 0.508 cm) or "Screen".
- **Width**: a journal's single column (8.3 cm), double column (17.8 cm), or free. A route wraps to the width.
- **Fit each step on one line**: each step of a route on its own row, numbered, all at one size: the largest at which
  the widest step fits the width, so the bonds are one length across the figure. Under half of print size, a warning
  says how small the text is and which width to choose.
- **As drawn**: the reaction as it is drawn, not arranged: the Reaction Sketcher's drawing, or a cell's ChemDraw page
  or RXN block with coordinates. It is drawn at the style's print size, in its colours, with the map numbers as the
  panel says. The width, the arrangements and the reaction centre's colours do not apply.
- **Format**: SVG, PNG at 300 or 600 dpi (its background transparent if you like), or CDXML for ChemDraw.
  **Copy as picture** puts the PNG on the clipboard. In the Reaction Sketcher, CDXML is the drawing as you drew it, with
  its texts, such as the conditions over the arrow. SVG, PNG and the preview draw the texts over and under the
  arrows too.
- **Reaction centre colours**: off for ACS 1996. A picture with the colours carries their legend.
- **Map numbers**: off for ACS 1996; for Screen, on when the picture you export from shows them. It applies to CDXML
  too: on, a mapped reaction's page keeps its atom mapping for ChemDraw; off, the page leaves it out. In a picture
  drawn smaller than 70% of print size (46% for Screen) the map numbers are small: the panel says so and which width
  draws them larger.

**Save** downloads what the preview shows.
