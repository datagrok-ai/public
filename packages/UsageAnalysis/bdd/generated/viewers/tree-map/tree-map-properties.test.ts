/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tree-map/tree-map-properties.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tree-map]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, filterIsExactlyCategory, filterPasses, filterPassesAll, onlyOfSelected, selectWhereIs, selectWhereOneOf, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, areaRepainted, clickArea, hasArea, hasNoArea, noErrors, painted, readingIs, readingReads, resizeTo, restoreSize, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {clickEmptySpace} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Tree map selection band, Row Source Selected, a filtering bar chart, Auto Layout and the categorical gate", () => {
  const session = feature(test, "features/viewers/tree-map/tree-map-properties.feature", import.meta.url);
  test("Tree map selection band, Row Source Selected, a filtering bar chart, Auto Layout and the categorical gate", {tag: ["@journey", "@viewers", "@realizes:viewers.tree-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a tree map viewer with:", () => addViewerWith(page, "tree map", [["splitByColumnNames","RACE"]]), [["splitByColumnNames","RACE"]]);
    await session.step(20, "Then the \"split columns\" reading of tree map viewer should be \"RACE\"", () => readingReads(page, "split columns", el("tree map viewer"), "RACE"));
    await session.step(21, "And the \"leaves\" reading of tree map viewer should be 4", () => readingIs(page, "leaves", el("tree map viewer"), 4));
    await session.step(22, "And the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
    await session.step(23, "And the \"rows of Asian\" reading of tree map viewer should be 15", () => readingIs(page, "rows of Asian", el("tree map viewer"), 15));
    await session.step(24, "And tree map viewer should be painted", () => painted(page, el("tree map viewer")));
    await run.scenario("Show Selected Rows off takes the selection band away and keeps the rows selected", async () => {
      await session.step(27, "Given user clears the row selection", () => clearSelection(page));
      await session.step(28, "When user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
      await session.step(29, "Then 15 rows should be selected", () => selectedRowCount(page, 15));
      await session.step(30, "And the \"selected rows of Asian\" reading of tree map viewer should be 15", () => readingIs(page, "selected rows of Asian", el("tree map viewer"), 15));
      await session.step(31, "And tree map viewer should have a \"selection of Asian\" area", () => hasArea(page, el("tree map viewer"), "selection of Asian"));
      await session.step(32, "And tree map viewer should not have a \"selection of Caucasian\" area", () => hasNoArea(page, el("tree map viewer"), "selection of Caucasian"));
      await session.step(33, "When user sets \"showSelectedRows\" property of tree map viewer to \"false\"", () => setProperty(page, "showSelectedRows", el("tree map viewer"), "false"));
      await session.step(34, "Then tree map viewer should not have a \"selection of Asian\" area", () => hasNoArea(page, el("tree map viewer"), "selection of Asian"));
      await session.step(35, "And the \"leaf Asian\" area of tree map viewer should have repainted", () => areaRepainted(page, "leaf Asian", el("tree map viewer")));
      await session.step(36, "And 15 rows should be selected", () => selectedRowCount(page, 15));
      await session.step(37, "And only rows where \"RACE\" is \"Asian\" should be selected", () => onlyOfSelected(page, "RACE", "Asian"));
      await session.step(38, "And the \"selected rows of Asian\" reading of tree map viewer should be 15", () => readingIs(page, "selected rows of Asian", el("tree map viewer"), 15));
      await session.step(39, "And tree map viewer should have a \"leaf Asian\" area", () => hasArea(page, el("tree map viewer"), "leaf Asian"));
      await session.step(40, "When user sets \"showSelectedRows\" property of tree map viewer to \"true\"", () => setProperty(page, "showSelectedRows", el("tree map viewer"), "true"));
      await session.step(41, "Then tree map viewer should have a \"selection of Asian\" area", () => hasArea(page, el("tree map viewer"), "selection of Asian"));
      await session.step(42, "And the \"leaf Asian\" area of tree map viewer should have repainted", () => areaRepainted(page, "leaf Asian", el("tree map viewer")));
      await session.step(43, "When user clears the row selection", () => clearSelection(page));
      await session.step(44, "Then tree map viewer should not have a \"selection of Asian\" area", () => hasNoArea(page, el("tree map viewer"), "selection of Asian"));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source Selected builds the map from the selection and draws no band", async () => {
      await session.step(48, "Given user clears the row selection", () => clearSelection(page));
      await session.step(49, "When user selects rows where \"RACE\" is one of \"Asian, Black\"", () => selectWhereOneOf(page, "RACE", "Asian, Black"));
      await session.step(50, "Then 42 rows should be selected", () => selectedRowCount(page, 42));
      await session.step(51, "And tree map viewer should have a \"selection of Asian\" area", () => hasArea(page, el("tree map viewer"), "selection of Asian"));
      await session.step(52, "And tree map viewer should have a \"selection of Black\" area", () => hasArea(page, el("tree map viewer"), "selection of Black"));
      await session.step(53, "When user sets \"rowSource\" property of tree map viewer to \"Selected\"", () => setProperty(page, "rowSource", el("tree map viewer"), "Selected"));
      await session.step(54, "Then the \"row source\" reading of tree map viewer should be \"Selected\"", () => readingReads(page, "row source", el("tree map viewer"), "Selected"));
      await session.step(55, "And the \"rows shown\" reading of tree map viewer should be 42", () => readingIs(page, "rows shown", el("tree map viewer"), 42));
      await session.step(56, "And the \"leaves\" reading of tree map viewer should be 2", () => readingIs(page, "leaves", el("tree map viewer"), 2));
      await session.step(57, "And the \"selected rows of Asian\" reading of tree map viewer should be 15", () => readingIs(page, "selected rows of Asian", el("tree map viewer"), 15));
      await session.step(58, "And the \"selected rows of Black\" reading of tree map viewer should be 27", () => readingIs(page, "selected rows of Black", el("tree map viewer"), 27));
      await session.step(59, "And tree map viewer should have a \"leaf Asian\" area", () => hasArea(page, el("tree map viewer"), "leaf Asian"));
      await session.step(60, "And tree map viewer should not have a \"selection of Asian\" area", () => hasNoArea(page, el("tree map viewer"), "selection of Asian"));
      await session.step(61, "And tree map viewer should not have a \"selection of Black\" area", () => hasNoArea(page, el("tree map viewer"), "selection of Black"));
      await session.step(62, "When user sets \"rowSource\" property of tree map viewer to \"Filtered\"", () => setProperty(page, "rowSource", el("tree map viewer"), "Filtered"));
      await session.step(63, "Then the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
      await session.step(64, "And tree map viewer should have a \"selection of Asian\" area", () => hasArea(page, el("tree map viewer"), "selection of Asian"));
      await session.step(65, "When user clears the row selection", () => clearSelection(page));
      await session.step(66, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A bar chart with On Click Filter narrows the tree map to the clicked bar (public#935)", async () => {
      await session.step(69, "When user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","SEX"],["Value Aggr Type","count"],["On Click","Filter"]]), [["Split","SEX"],["Value Aggr Type","count"],["On Click","Filter"]]);
      await session.step(73, "Then the \"bars\" reading of bar chart viewer should be 2", () => readingIs(page, "bars", el("bar chart viewer"), 2));
      await session.step(74, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(75, "When user clicks on the \"bar M\" area of bar chart viewer", () => clickArea(page, "bar M", el("bar chart viewer")));
      await session.step(76, "Then the filter should pass exactly the rows where \"SEX\" is \"M\"", () => filterIsExactlyCategory(page, "SEX", "M"));
      await session.step(77, "And 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(78, "And the \"rows shown\" reading of tree map viewer should be 447", () => readingIs(page, "rows shown", el("tree map viewer"), 447));
      await session.step(79, "And the \"rows of Caucasian\" reading of tree map viewer should be 416", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 416));
      await session.step(80, "And the \"rows of Asian\" reading of tree map viewer should be 8", () => readingIs(page, "rows of Asian", el("tree map viewer"), 8));
      await session.step(81, "And the \"leaves\" reading of tree map viewer should be 4", () => readingIs(page, "leaves", el("tree map viewer"), 4));
      await session.step(82, "When user clicks on empty plot space of bar chart viewer", () => clickEmptySpace(page, el("bar chart viewer")));
      await session.step(83, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(84, "And the \"rows shown\" reading of tree map viewer should be 1000", () => readingIs(page, "rows shown", el("tree map viewer"), 1000));
      await session.step(85, "And the \"rows of Caucasian\" reading of tree map viewer should be 896", () => readingIs(page, "rows of Caucasian", el("tree map viewer"), 896));
      await session.step(86, "When user clicks on close icon of bar chart viewer", () => clickOn(page, el("close icon of bar chart viewer")));
      await session.step(87, "Then the open tableview should have 0 bar chart viewers", () => viewerCount(page, 0, "bar chart"));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Auto Layout hides the selectors on a small map and keeps them once it is off (public#2131)", async () => {
      await session.step(91, "Then the \"selection panel shown\" reading of tree map viewer should be \"true\"", () => readingReads(page, "selection panel shown", el("tree map viewer"), "true"));
      await session.step(92, "And tree map viewer should have a \"split selector 1\" area", () => hasArea(page, el("tree map viewer"), "split selector 1"));
      await session.step(93, "And tree map viewer should have a \"color selector\" area", () => hasArea(page, el("tree map viewer"), "color selector"));
      await session.step(94, "When user resizes tree map viewer to 180 by 180", () => resizeTo(page, el("tree map viewer"), 180, 180));
      await session.step(95, "Then the \"selection panel shown\" reading of tree map viewer should be \"false\"", () => readingReads(page, "selection panel shown", el("tree map viewer"), "false"));
      await session.step(96, "And tree map viewer should not have a \"split selector 1\" area", () => hasNoArea(page, el("tree map viewer"), "split selector 1"));
      await session.step(97, "And tree map viewer should not have a \"color selector\" area", () => hasNoArea(page, el("tree map viewer"), "color selector"));
      await session.step(98, "And the \"leaves\" reading of tree map viewer should be 4", () => readingIs(page, "leaves", el("tree map viewer"), 4));
      await session.step(99, "When user sets \"autoLayout\" property of tree map viewer to \"false\"", () => setProperty(page, "autoLayout", el("tree map viewer"), "false"));
      await session.step(100, "Then the \"selection panel shown\" reading of tree map viewer should be \"true\"", () => readingReads(page, "selection panel shown", el("tree map viewer"), "true"));
      await session.step(101, "And tree map viewer should have a \"split selector 1\" area", () => hasArea(page, el("tree map viewer"), "split selector 1"));
      await session.step(102, "And tree map viewer should have a \"color selector\" area", () => hasArea(page, el("tree map viewer"), "color selector"));
      await session.step(103, "When user sets \"autoLayout\" property of tree map viewer to \"true\"", () => setProperty(page, "autoLayout", el("tree map viewer"), "true"));
      await session.step(104, "Then the \"selection panel shown\" reading of tree map viewer should be \"false\"", () => readingReads(page, "selection panel shown", el("tree map viewer"), "false"));
      await session.step(105, "And tree map viewer should not have a \"split selector 1\" area", () => hasNoArea(page, el("tree map viewer"), "split selector 1"));
      await session.step(106, "When user restores the size of tree map viewer", () => restoreSize(page, el("tree map viewer")));
      await session.step(107, "Then the \"selection panel shown\" reading of tree map viewer should be \"true\"", () => readingReads(page, "selection panel shown", el("tree map viewer"), "true"));
      await session.step(108, "And tree map viewer should have a \"split selector 1\" area", () => hasArea(page, el("tree map viewer"), "split selector 1"));
      await session.step(109, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A table with no categorical column gets the categorical gate", async () => {
      await session.step(112, "Given user opens a table \"tree-map-numbers\" with:", () => openTableOf(page, "tree-map-numbers", [["x","y"],["1","2.5"],["2","3.5"],["3","4.5"]]), [["x","y"],["1","2.5"],["2","3.5"],["3","4.5"]]);
      await session.step(117, "And user adds a tree map viewer", () => addViewer(page, "tree map"));
      await session.step(118, "Then the \"error\" reading of tree map viewer should be \"Requires categorical columns\"", () => readingReads(page, "error", el("tree map viewer"), "Requires categorical columns"));
      await session.step(119, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
