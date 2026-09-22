/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-tooltip-labels.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/grid.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, filterPasses, filterPassesAll, noneSelected, selectedPassFilter, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, dragSelectionOverArea, hoverArea, noErrors, painted, pointerAway, propertiesShouldBe, readingAtLeast, readingHigher, readingIs, setProperties, setProperty, showsRows, tooltipColumns, tooltipSomeColumns, tooltipValue} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot tooltip and marker labels", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-tooltip-labels.feature", import.meta.url);
  test("Scatter plot tooltip and marker labels", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(22, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(23, "And properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["Show Tooltip","inherit from table"],["Label Columns",""]]), [["Show Tooltip","inherit from table"],["Label Columns",""]]);
    await run.scenario("Hovering a marker makes its row hovered and shows the table's tooltip", async () => {
      await session.step(28, "Then the \"hovered row\" reading of scatter plot viewer should be 0", () => readingIs(page, "hovered row", el("scatter plot viewer"), 0));
      await session.step(29, "When user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(30, "Then the \"hovered row\" reading of scatter plot viewer should be 11", () => readingIs(page, "hovered row", el("scatter plot viewer"), 11));
      await session.step(31, "And tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(32, "And the tooltip should show some columns", () => tooltipSomeColumns(page));
      await session.step(33, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
      await session.step(34, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A custom column list replaces what the tooltip shows", async () => {
      await session.step(38, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Show Tooltip","show custom tooltip"],["Row Tooltip","AGE\nSEX"],["Data Values","Do not add"],["showLabels","Always"]]), [["Show Tooltip","show custom tooltip"],["Row Tooltip","AGE\nSEX"],["Data Values","Do not add"],["showLabels","Always"]]);
      await session.step(43, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(44, "Then the \"hovered row\" reading of scatter plot viewer should be 11", () => readingIs(page, "hovered row", el("scatter plot viewer"), 11));
      await session.step(45, "And the tooltip should show columns \"AGE, SEX\"", () => tooltipColumns(page, "AGE, SEX"));
      await session.step(46, "And the tooltip should show \"AGE\" as \"46\"", () => tooltipValue(page, "AGE", "46"));
      await session.step(47, "And the tooltip should show \"SEX\" as \"F\"", () => tooltipValue(page, "SEX", "F"));
      await session.step(48, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
      await session.step(49, "And user sets \"Show Tooltip\" property of scatter plot viewer to \"do not show\"", () => setProperty(page, "Show Tooltip", el("scatter plot viewer"), "do not show"));
      await session.step(50, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(51, "Then the \"hovered row\" reading of scatter plot viewer should be 11", () => readingIs(page, "hovered row", el("scatter plot viewer"), 11));
      await session.step(52, "And tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(53, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
      await session.step(54, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Show Tooltip","inherit from table"],["Row Tooltip",""],["Data Values","Merge"],["showLabels","Auto"]]), [["Show Tooltip","inherit from table"],["Row Tooltip",""],["Data Values","Merge"],["showLabels","Auto"]]);
      await session.step(59, "And user hovers over the \"marker of row 11\" area of scatter plot viewer", () => hoverArea(page, "marker of row 11", el("scatter plot viewer")));
      await session.step(60, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(61, "And the tooltip should show some columns", () => tooltipSomeColumns(page));
      await session.step(62, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
      await session.step(63, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Labels are drawn for every row, then only for the selected ones", async () => {
      await session.step(66, "Then no rows should be selected", () => noneSelected(page));
      await session.step(67, "And the \"labels shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "labels shown", el("scatter plot viewer"), 0));
      await session.step(68, "When user sets \"Label Columns\" property of scatter plot viewer to \"AGE\"", () => setProperty(page, "Label Columns", el("scatter plot viewer"), "AGE"));
      await session.step(69, "Then the \"labels shown\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "labels shown", el("scatter plot viewer"), 1));
      await session.step(70, "When user sets \"Show Labels For\" property of scatter plot viewer to \"Selected\"", () => setProperty(page, "Show Labels For", el("scatter plot viewer"), "Selected"));
      await session.step(71, "Then the \"labels shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "labels shown", el("scatter plot viewer"), 0));
      await session.step(72, "When user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(73, "Then some rows should be selected", () => someSelected(page));
      await session.step(74, "And the \"labels shown\" reading of scatter plot viewer should be higher than before", () => readingHigher(page, "labels shown", el("scatter plot viewer")));
      await session.step(75, "When user clicks on the \"empty space\" area of scatter plot viewer", () => clickArea(page, "empty space", el("scatter plot viewer")));
      await session.step(76, "Then no rows should be selected", () => noneSelected(page));
      await session.step(77, "And the \"labels shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "labels shown", el("scatter plot viewer"), 0));
      await session.step(78, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Label Columns",""],["Show Labels For","All"]]), [["Label Columns",""],["Show Labels For","All"]]);
      await session.step(81, "Then the \"labels shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "labels shown", el("scatter plot viewer"), 0));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Labels on a datetime axis survive a filter that keeps one category and a selection", async () => {
      await session.step(85, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X","STARTED"],["Y","RACE"],["Label Columns","USUBJID"],["Show Labels For","Selected"]]), [["X","STARTED"],["Y","RACE"],["Label Columns","USUBJID"],["Show Labels For","Selected"]]);
      await session.step(90, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
      await session.step(91, "When user adds a categorical filter on \"RACE\" keeping \"Caucasian\"", () => addCategoricalFilter(page, "RACE", "Caucasian"));
      await session.step(92, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(93, "And scatter plot viewer should show 896 rows", () => showsRows(page, el("scatter plot viewer"), 896));
      await session.step(94, "When user drags a selection box over the \"view\" area of scatter plot viewer", () => dragSelectionOverArea(page, "view", el("scatter plot viewer")));
      await session.step(95, "Then some rows should be selected", () => someSelected(page));
      await session.step(96, "And every selected row should pass the filter", () => selectedPassFilter(page));
      await session.step(97, "And the \"labels shown\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "labels shown", el("scatter plot viewer"), 1));
      await session.step(98, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(99, "And no errors should have been logged", () => noErrors(page));
      await session.step(100, "When user adds a categorical filter on \"RACE\" keeping \"Asian, Black, Caucasian, Other\"", () => addCategoricalFilter(page, "RACE", "Asian, Black, Caucasian, Other"));
      await session.step(101, "And user clicks on the \"empty space\" area of scatter plot viewer", () => clickArea(page, "empty space", el("scatter plot viewer")));
      await session.step(102, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(103, "And no rows should be selected", () => noneSelected(page));
      await session.step(104, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X","WEIGHT"],["Y","HEIGHT"],["Label Columns",""],["Show Labels For","All"]]), [["X","WEIGHT"],["Y","HEIGHT"],["Label Columns",""],["Show Labels For","All"]]);
      await session.step(109, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(110, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
