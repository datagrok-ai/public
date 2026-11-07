# Scaffold Tree — what the features cannot say yet

Written during the round of 26 Oct 2026 that turned closed issues into regression scenarios
(`scaffold-tree-regressions.feature`), with the existing steps only: no library, binding or code change.
Every scaffold there is sketched by hand, so none of them runs the Python tree generation.

## Still open

1. **The count a node shows** (public#2550: the node counters did not follow the table's filter). `hits of node
   <i>` is the node's bitset count over the whole table, not the label the node draws. Wanted in
   `ScaffoldTreeViewer.getWidgetStatus()` (Chem `src/widgets/scaffold-tree.ts`): a `shown of node <i>` reading
   of the counter label.

   ```gherkin
   Scenario: Node counts follow the table's filter (public#2550)
     # a pyridine node sketched by hand, as in scaffold-tree-regressions.feature
     Then the "hits of node 1" reading of Scaffold Tree viewer should be 17
     When user remembers the "shown of node 1" reading of Scaffold Tree viewer
     And user filters rows where "Stereo Category" is "S_UNKN"
     # NEW SIGNAL (Chem): the counter label the node draws
     Then the "shown of node 1" reading of Scaffold Tree viewer should be lower than remembered
   ```

   The exact count is to be measured once the reading exists.

2. **The filter panel's hover text** (public#2504: the panel's ? did not describe the Scaffold Tree filter as it
   does other filters). Needs a step that hovers the filter indicator of the filter panel and reads its tooltip
   text; the panel's `summary of <caption>` reading is what the counter's tooltip lists, but nothing yet claims the
   tooltip itself.

   ```gherkin
   Scenario: The panel's hover describes the Scaffold Tree filter (public#2504)
     # a Scaffold Tree filter card with a checked pyridine node
     # NEW STEP: When user hovers over the filter indicator of filter panel
     Then tooltip should contain the text "Structure"
   ```

3. **The icons of a deep node** (public#4091: with more than five Scaffold Trees in a view, a node at a deep level of
   the hierarchy, drawn toward the right edge, showed no add/edit/remove icons on hover). The `add icon of node <i>` /
   `edit icon of node <i>` areas are reported whenever the icon has a size, so an icon pushed past the viewer's right
   edge (clipped by the scrolling tree, not hidden) is still reported, and `should have an "edit icon of node 9" area`
   would pass with the bug back. Wanted in `ScaffoldTreeViewer.getWidgetStatus()` (Chem `src/widgets/scaffold-tree.ts`):
   a `depth of node <i>` reading and the icon areas clipped to the tree's visible box (dropped when outside it). Not
   run on the stand (2026-11-05): the scaffold-tree features are red on master today.

   ```gherkin
   Scenario: A deep node shows its icons on hover (public#4091)
     # six Scaffold Tree viewers on SPGI, a tree generated and expanded
     # NEW READING (Chem): depth of node <i>; icon areas clipped to the visible tree
     When user hovers over the "node 9" area of Scaffold Tree viewer
     Then the "depth of node 9" reading of Scaffold Tree viewer should be at least 4
     And Scaffold Tree viewer should have an "edit icon of node 9" area
   ```

## Package tests, not scenarios

- **public#2139** (a scaffold made invalid in the sketcher and fixed again did not filter): `src/tests/scaffold-tree-tests.ts`
  'edit invalid structure', skipped under GROK-16714.
- **public#2626** (a child of a scaffold drawn with an explicit H opened with the "not a superstructure" warning): the
  same file, 'parent node drawn with an explicit H opens its child with no warning'.
