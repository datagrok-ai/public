# Tile Viewer — manual checklist

Exploratory only. Everything that could be automated has moved into the section's
Playwright scenarios, so what is left here is what a person still has to judge by eye.
Automated elsewhere, do not repeat here: the Edit Form designer (field removal, label
deletion, RESET, CLOSE AND APPLY) and the layout round-trip in
`tile-viewer-selection-form-editor.md`; dragging a tile between lanes, the lanes ladder
and the persistence peak in `tile-viewer-lanes-persist.md`; the property surface, the
auto-generate states and the context-menu inventory in `tile-viewer.md`.

All scenarios start with: Close all → open demog → add a Tile Viewer.

## Colour coding on tiles

Colour coding is applied from the GRID's column menu — the Tile Viewer's own context menu
has no Colour Coding item. Whether a tile field reflects it, and how legibly, is a visual
judgement, which is why it stays here.

1. In the grid, open the AGE column's menu and switch colour coding to Linear
2. Look at the tiles: the AGE field's background should form a gradient across rows, and the
   value text should stay readable against it
3. Switch the RACE column to Categorical colour coding
4. Look at the tiles: each RACE value should get a visibly distinct background, and no two
   adjacent categories should be hard to tell apart
5. Turn colour coding off on both columns — the tile fields return to a plain background

## Edit Form with a different source table

The SketchView has no table field on the current build, so this flow cannot be driven as
written and is kept as a manual probe of the broken-binding behaviour.

1. Open demog and spgi-100
2. On the spgi-100 view add a Tile Viewer, open Edit Form..., design a card, CLOSE AND APPLY
3. Switch the viewer's Table to demog
4. Look at the tiles: report anything that renders as a broken or empty binding placeholder

## Picking the lanes column from the property panel

The automated scenarios set the lanes column through the viewer's property object, because
the picker is a canvas grid with no selector for an individual row. The UI pick itself is
therefore exercised nowhere else, and this is the step that covers it.

1. Open the viewer's Properties → Data and click the Lanes column picker
2. Choose RACE: the tiles regroup into one lane per race, each lane headed by its value
3. Choose a numeric column such as AGE: report whether the picker offers it at all, and if
   it does, whether the resulting lanes are usable or one lane per distinct number
4. Clear the picker: the tiles collapse back into a single lane

## Hamburger menu walk

Nothing on this surface is automated: no scenario opens the hamburger menu, and none
measures how long a menu takes to open. This walk is therefore the ONLY coverage the
surface has, and a defect missed here is missed everywhere.

1. Open the viewer's hamburger menu and walk every top-level entry and submenu
2. Report any entry that opens nothing, opens the wrong thing, or leaves the page unresponsive
3. On a table with MANY distinct lane values, time the hamburger OPENING itself — that is
   where the cost sits (GROK-19356), not in walking the submenus. Report a noticeable delay
   between the click and the menu appearing
4. Note that a group row (Properties, Lanes) does not expand on hover, and clicking it closes
   the whole menu — report it if either behaves differently for you, since that is what makes
   this surface awkward to automate and worth a human pass

---
{
  "order": 200,
  "datasets": ["System:DemoFiles/demog.csv"]
}
