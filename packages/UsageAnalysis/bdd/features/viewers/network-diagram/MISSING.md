# Network diagram — what the features cannot say yet

Written during the viewer-gaps round of 23 Oct 2026, against
`files/TestTrack/Viewers/NetworkDiagram/network-diagram.md` and the viewer's "Automation surface"
(`core/client/d4/lib/src/viewers/network_diagram/CLAUDE.md`). The status reports nodes (areas, names,
`rows of node`) but nothing about edges, the empty canvas or the view's scale.

## Still open

- Clicking an edge selects the rows behind it (md "Selection gestures", step 6; Select Edges On Click off): needs
  `edge "<from> -> <to>"` areas and a `rows of edge` reading from vis's edge geometry.
- A double-click on empty canvas clears the selection and fits the view (md step 7): needs an `empty space` area and a
  `view scale` reading (`Network.getScale()`).
- Show Arrows: needs an `arrow heads` reading (the edges vis draws with an arrow), to claim 4 on demog-1000 and 0 off.

## Won't do

- The core#1041 crashes (one column as both nodes, an unmerged node hovered, the node columns removed, a click then a
  close) and the edge widths after a filter: xamgle `d4 | Viewers | Network diagram | Filter, hover and close edge
  cases` and `d4 | Viewers | Network diagram | After a filter the edge width and color aggregate the rows left…`.
- No rows and deleted rows: xamgle `d4 | Viewers | Data Manipulation | Network diagram`.
- More than 10000 nodes ("Too many nodes to display"), a wrong `onNodeExpandFunction`, a negative edge width: low value.
