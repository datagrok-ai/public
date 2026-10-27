/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-row-source-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart, public#2350, public#1650, public#1628, public#3709, public#1582, public#2485, public#3210]
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
import {clearSelection, selectWhereIs, selectWhereOneOf} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hoverArea, loadLayout, noErrors, notRepainted, painted, pointerAway, propertyShouldBe, readingBetween, readingIs, repainted, reportsNoError, saveLayout, setProperty, showsSelectedRows, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart regressions — row source, row groups, categories and the log axis", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-row-source-regressions.feature", import.meta.url);
  test("A layout with a categorical X and Row Source Selected loads with nothing selected (public#2350)", {tag: ["@viewers", "@realizes:viewers.line-chart", "@realizes:public#2350"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(24, "When user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","RACE"],["yColumnNames","AGE"],["rowSource","Selected"]]), [["xColumnName","RACE"],["yColumnNames","AGE"],["rowSource","Selected"]]);
    await session.step(28, "And user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
    await session.step(29, "Then the \"rows shown\" reading of line chart viewer should be 15", () => readingIs(page, "rows shown", el("line chart viewer"), 15));
    await session.step(30, "When user saves the layout of the current table view", () => saveLayout(page));
    await session.step(31, "And user clears the row selection", () => clearSelection(page));
    await session.step(32, "Then the \"rows shown\" reading of line chart viewer should be 0", () => readingIs(page, "rows shown", el("line chart viewer"), 0));
    await session.step(33, "When user loads the saved layout", () => loadLayout(page));
    await session.step(34, "Then line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await session.step(35, "And \"rowSource\" property of line chart viewer should be \"Selected\"", () => propertyShouldBe(page, "rowSource", el("line chart viewer"), "Selected"));
    await session.step(36, "And the \"rows shown\" reading of line chart viewer should be 0", () => readingIs(page, "rows shown", el("line chart viewer"), 0));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
    await session.step(38, "When user selects rows where \"RACE\" is \"Black\"", () => selectWhereIs(page, "RACE", "Black"));
    await session.step(39, "Then the \"rows shown\" reading of line chart viewer should be 27", () => readingIs(page, "rows shown", el("line chart viewer"), 27));
    await session.step(40, "And the \"x categories\" reading of line chart viewer should be 1", () => readingIs(page, "x categories", el("line chart viewer"), 1));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
  });
  test("A split chart under Row Source Selected logs nothing and draws exactly the selected rows (public#1650)", {tag: ["@viewers", "@realizes:viewers.line-chart", "@realizes:public#1650", "@realizes:public#1628"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(47, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(48, "When user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","Average Mass"],["yColumnNames","Chemical Space X"],["splitColumnNames","Series"]]), [["xColumnName","Average Mass"],["yColumnNames","Chemical Space X"],["splitColumnNames","Series"]]);
    await session.step(52, "Then the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
    await session.step(53, "When user sets \"rowSource\" property of line chart viewer to \"Selected\"", () => setProperty(page, "rowSource", el("line chart viewer"), "Selected"));
    await session.step(54, "Then the \"rows shown\" reading of line chart viewer should be 0", () => readingIs(page, "rows shown", el("line chart viewer"), 0));
    await session.step(55, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await session.step(56, "And no errors should have been logged", () => noErrors(page));
    await session.step(57, "When user selects rows where \"Series\" is \"Pyrrolidines\"", () => selectWhereIs(page, "Series", "Pyrrolidines"));
    await session.step(58, "Then line chart viewer should show the selected rows of its table", () => showsSelectedRows(page, el("line chart viewer")));
    await session.step(59, "And the \"rows shown\" reading of line chart viewer should be 21", () => readingIs(page, "rows shown", el("line chart viewer"), 21));
    await session.step(60, "And the \"lines\" reading of line chart viewer should be 1", () => readingIs(page, "lines", el("line chart viewer"), 1));
    await session.step(61, "And the \"markers drawn\" reading of line chart viewer should be 20", () => readingIs(page, "markers drawn", el("line chart viewer"), 20));
    await session.step(62, "When user selects rows where \"Series\" is one of \"Pyrrolidines, Aminopiperidines\"", () => selectWhereOneOf(page, "Series", "Pyrrolidines, Aminopiperidines"));
    await session.step(63, "Then the \"rows shown\" reading of line chart viewer should be 25", () => readingIs(page, "rows shown", el("line chart viewer"), 25));
    await session.step(64, "And the \"lines\" reading of line chart viewer should be 2", () => readingIs(page, "lines", el("line chart viewer"), 2));
    await session.step(65, "And the \"markers drawn\" reading of line chart viewer should be 24", () => readingIs(page, "markers drawn", el("line chart viewer"), 24));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
  });
  test("A row group hovered in a bar chart repaints the line chart only while Show Mouse Over Row Group is on (public#3709)", {tag: ["@viewers", "@realizes:viewers.line-chart", "@realizes:public#3709"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(71, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(72, "When user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","HEIGHT"]]), [["xColumnName","AGE"],["yColumnNames","HEIGHT"]]);
    await session.step(75, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
    await session.step(77, "And user takes a snapshot of line chart viewer", () => takeSnapshot(page, el("line chart viewer")));
    await session.step(78, "And user hovers over the \"bar Asian\" area of bar chart viewer", () => hoverArea(page, "bar Asian", el("bar chart viewer")));
    await session.step(79, "Then line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
    await session.step(80, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(81, "And user sets \"showMouseOverRowGroup\" property of line chart viewer to \"false\"", () => setProperty(page, "showMouseOverRowGroup", el("line chart viewer"), "false"));
    await session.step(82, "And user takes a snapshot of line chart viewer", () => takeSnapshot(page, el("line chart viewer")));
    await session.step(83, "And user hovers over the \"bar Asian\" area of bar chart viewer", () => hoverArea(page, "bar Asian", el("bar chart viewer")));
    await session.step(84, "Then bar chart viewer should have repainted", () => repainted(page, el("bar chart viewer")));
    await session.step(85, "And line chart viewer should not have repainted", () => notRepainted(page, el("line chart viewer")));
    await session.step(86, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
    await session.step(87, "And user sets \"showMouseOverRowGroup\" property of line chart viewer to \"true\"", () => setProperty(page, "showMouseOverRowGroup", el("line chart viewer"), "true"));
    await session.step(88, "Then no errors should have been logged", () => noErrors(page));
  });
  test("A molecule column on X draws a category per structure (public#1582)", {tag: ["@viewers", "@realizes:viewers.line-chart", "@realizes:public#1582"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(93, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(94, "When user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","Core"],["yColumnNames","Chemical Space X"]]), [["xColumnName","Core"],["yColumnNames","Chemical Space X"]]);
    await session.step(97, "Then the \"x categories\" reading of line chart viewer should be 6", () => readingIs(page, "x categories", el("line chart viewer"), 6));
    await session.step(98, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await session.step(99, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
    await session.step(100, "And no errors should have been logged", () => noErrors(page));
  });
  test("A logarithmic Y axis fits the data (public#2485)", {tag: ["@viewers", "@realizes:viewers.line-chart", "@realizes:public#2485"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(105, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(106, "When user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","Num Heavy Atoms"],["yColumnNames","Whole blood assay 1"]]), [["xColumnName","Num Heavy Atoms"],["yColumnNames","Whole blood assay 1"]]);
    await session.step(109, "And user sets \"yAxisType\" property of line chart viewer to \"logarithmic\"", () => setProperty(page, "yAxisType", el("line chart viewer"), "logarithmic"));
    await session.step(110, "Then the \"y axis max of \\\"Whole blood assay 1\\\"\" reading of line chart viewer should be between 31.5 and 37.8", () => readingBetween(page, "y axis max of \"Whole blood assay 1\"", el("line chart viewer"), 31.5, 37.8));
    await session.step(111, "And the \"y axis min of \\\"Whole blood assay 1\\\"\" reading of line chart viewer should be between 0.04 and 0.14", () => readingBetween(page, "y axis min of \"Whole blood assay 1\"", el("line chart viewer"), 0.04, 0.14));
    await session.step(112, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await session.step(113, "And no errors should have been logged", () => noErrors(page));
  });
  test("A split with a category per row keeps the chart working (public#3210)", {tag: ["@viewers", "@realizes:viewers.line-chart", "@realizes:public#3210"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(118, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(119, "When user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","STARTED"],["yColumnNames","AGE"],["splitColumnNames","USUBJID"]]), [["xColumnName","STARTED"],["yColumnNames","AGE"],["splitColumnNames","USUBJID"]]);
    await session.step(123, "Then the \"split columns\" reading of line chart viewer should be 1", () => readingIs(page, "split columns", el("line chart viewer"), 1));
    await session.step(124, "And the \"categories\" reading of line chart viewer should be 1000", () => readingIs(page, "categories", el("line chart viewer"), 1000));
    await session.step(125, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await session.step(126, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
    await session.step(127, "Then the \"split columns\" reading of line chart viewer should be 0", () => readingIs(page, "split columns", el("line chart viewer"), 0));
    await session.step(128, "And the \"categories\" reading of line chart viewer should be 1", () => readingIs(page, "categories", el("line chart viewer"), 1));
    await session.step(129, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
    await session.step(130, "And no errors should have been logged", () => noErrors(page));
  });
});
