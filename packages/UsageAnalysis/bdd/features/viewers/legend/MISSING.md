# Legend — what the features cannot say yet

Written during the Grid / Filter panel / Legend gap round of 23 Oct 2026, against the old
`playwright-public/legend` specs, `files/TestTrack/Viewers/Legend/*.md`, the legend's DOM contract
(`core/client/d4/lib/src/viewer_base/legend_mixin.dart`, `widgets/legend/`) and the library's legend
steps (`bindings/tiers/viewers/steps.ts`, `src/runtime/viewer-legend.ts`). Added in that round:
`legend-stability.feature` (zoom and pan never move the legend, a folded corner legend stays folded
through Position Auto and resizes, a 1000-category column docks and lists 872 items).

Pushed back: the survey's "collapsed corner legend stays collapsed after a layout round trip" — the
chevron's fold is not a viewer setting (`legend-placement.feature` says why), and the old spec's
claim was about Position Auto and resizes, which `legend-stability.feature` now makes. The mini icon
click is already driven (`user clicks on mini legend icon of <viewer>`).

## Still open

- **Per-section readings.** A scatter plot coloured by one column and marked by another has two
  legend sections (`main`, `extra`); the legend publishes one total (`data-legend-items`), so
  `legend-sections.test.ts`'s "changing the markers column leaves the colour section's items
  untouched" cannot be said. Needs `data-legend-items-<section>` (or a section's item keys) on the
  legend root:

  ```gherkin
  Scenario: The colour and the markers are independent sections
    When user sets "Markers" property of scatter plot viewer to "SEX"
    Then the "main" section of the legend of scatter plot viewer should list 4 items
    And the "extra" section of the legend of scatter plot viewer should list 2 items
    When user sets "Markers" property of scatter plot viewer to "DIS_POP"
    Then the "main" section of the legend of scatter plot viewer should list the same items as before
    And the "extra" section of the legend of scatter plot viewer should list 6 items
  ```
- **The inner splitter between two sections** (`LegendView.initInnerSplitter`): no name on it and
  no step drags it; `the legend of X should be taller than before` reads the whole legend. Needs a
  name (`legend-inner-splitter`) and per-section heights.
- **The markers column selector in the legend header** (`legend-selectors.test.ts`, "picking a
  column works, including after a tooltip cycle"): the element is named (`legend markers selector`),
  but no step picks a column in a selector that is not a viewer's own hit area
  (`user picks {string} in the column selector at the {string} area of {widget}` needs an area).
  A `legend markers selector` area in the viewer status, or a generic "picks X in {element}" for a
  column combobox, would let the scenario be written.
- **A label that does not fit** (`legend-sizing-and-chrome.test.ts`): whether an item's label is
  cut (`scrollWidth > clientWidth` on `.d4-legend-value`) and the docked legend's width cap. The
  geometry step `should not be clipped` reads the element box, not text overflow inside it. Needs a
  `data-legend-truncated` count or a per-item overflow flag.
- **Scrolling a long list** (`legend-stability.test.ts`, "hundreds of categories stay virtualised"):
  the list renders a window of its items, so `user scrolls {element} to the middle of its list`
  cannot reach an item that is not in the page yet. Needs the first rendered item key published
  (`data-legend-first-item`) and a step that scrolls the legend list by a page.
