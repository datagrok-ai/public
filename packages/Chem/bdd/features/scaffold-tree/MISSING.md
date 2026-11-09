# Scaffold Tree — what the features do not claim

## Still open

- public#2550 (a node's counter did not follow the table's filter): `hits of node <i>` counts the whole table, not the
  label the node draws; needs a `shown of node <i>` reading in `ScaffoldTreeViewer.getWidgetStatus()` (Chem
  `src/widgets/scaffold-tree.ts`); claim: the reading drops after `user filters rows where "Stereo Category" is "S_UNKN"`.

## Won't do

- public#2504 (the filter panel's hover text for the Scaffold Tree filter): tooltip wording.
- public#4091 (a deep node's add/edit/remove icons with six Scaffold Trees in a view): an edge case; the icon areas
  would need clipping to the visible tree.
- public#2139, public#2626: Chem `src/tests/scaffold-tree-tests.ts` ('edit invalid structure', skipped under
  GROK-16714; 'parent node drawn with an explicit H opens its child with no warning').
