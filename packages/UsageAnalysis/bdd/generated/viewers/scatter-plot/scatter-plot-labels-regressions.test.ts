/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-labels-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot]
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
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, filterPasses, noneSelected, removeColumn, selectFirstRows, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, dragAreaBy, hasArea, loadLayout, noBalloons, noErrors, propertyShouldBe, readingAtLeast, readingIs, readingNotAsRemembered, rememberReading, reportsNoError, saveLayout, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot regressions — labels and missing columns", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-labels-regressions.feature", import.meta.url);
  test("A label is dragged under Mouse Drag = Select (public#3120)", {tag: ["@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"],["Label Columns","AGE"],["Show Labels For","Selected"],["mouseDrag","Select"]]), [["X","WEIGHT"],["Y","HEIGHT"],["Label Columns","AGE"],["Show Labels For","Selected"],["mouseDrag","Select"]]);
    await session.step(24, "When user selects the first 1 rows", () => selectFirstRows(page, 1));
    await session.step(25, "Then scatter plot viewer should have a \"label of row 1\" area", () => hasArea(page, el("scatter plot viewer"), "label of row 1"));
    await session.step(26, "And the \"labels shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "labels shown", el("scatter plot viewer"), 1));
    await session.step(27, "When user drags the \"label of row 1\" area of scatter plot viewer by 60 pixels to the right", () => dragAreaBy(page, "label of row 1", el("scatter plot viewer"), 60, "right"));
    await session.step(28, "Then 1 row should be selected", () => selectedRowCount(page, 1));
    await session.step(29, "And the \"labels shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "labels shown", el("scatter plot viewer"), 0));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
  });
  test("A label is dragged when Show Labels For is All and a panel filter is set (public#3323)", {tag: ["@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(34, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(35, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"],["Label Columns","AGE"],["Show Labels For","All"],["mouseDrag","Select"]]), [["X","WEIGHT"],["Y","HEIGHT"],["Label Columns","AGE"],["Show Labels For","All"],["mouseDrag","Select"]]);
    await session.step(41, "When user adds a categorical filter on \"SEX\" keeping \"F\"", () => addCategoricalFilter(page, "SEX", "F"));
    await session.step(42, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
    await session.step(43, "And scatter plot viewer should have a \"label of row 1\" area", () => hasArea(page, el("scatter plot viewer"), "label of row 1"));
    await session.step(44, "When user remembers the \"labels shown\" reading of scatter plot viewer", () => rememberReading(page, "labels shown", el("scatter plot viewer")));
    await session.step(45, "And user drags the \"label of row 1\" area of scatter plot viewer by 60 pixels to the right", () => dragAreaBy(page, "label of row 1", el("scatter plot viewer"), 60, "right"));
    await session.step(46, "Then the \"labels shown\" reading of scatter plot viewer should not be as remembered", () => readingNotAsRemembered(page, "labels shown", el("scatter plot viewer")));
    await session.step(47, "And no rows should be selected", () => noneSelected(page));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
  });
  test("A molecule column as the label column draws labels (public#3041)", {tag: ["@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(52, "Given user opens spgi-100 dataset", () => openDataset(page, ds("spgi-100")));
    await session.step(53, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","Stereo Category"],["Y","Chemist 521"]]), [["X","Stereo Category"],["Y","Chemist 521"]]);
    await session.step(56, "Then the \"labels shown\" reading of scatter plot viewer should be 0", () => readingIs(page, "labels shown", el("scatter plot viewer"), 0));
    await session.step(57, "When user sets \"Label Columns\" property of scatter plot viewer to \"Core\"", () => setProperty(page, "Label Columns", el("scatter plot viewer"), "Core"));
    await session.step(58, "Then the \"labels shown\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "labels shown", el("scatter plot viewer"), 1));
    await session.step(59, "And scatter plot viewer should report no error", () => reportsNoError(page, el("scatter plot viewer")));
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
  });
  test("A plot whose logarithmic X column is removed lets the layout be saved (public#2867)", {tag: ["@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(64, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(65, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"],["X Axis Type","logarithmic"],["Y Axis Type","logarithmic"]]), [["X","WEIGHT"],["Y","HEIGHT"],["X Axis Type","logarithmic"],["Y Axis Type","logarithmic"]]);
    await session.step(70, "When user saves the layout of the current table view", () => saveLayout(page));
    await session.step(71, "And user removes \"WEIGHT\" column", () => removeColumn(page, "WEIGHT"));
    await session.step(72, "And user loads the saved layout", () => loadLayout(page));
    await session.step(73, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(74, "And \"X\" property of scatter plot viewer should be \"WEIGHT\"", () => propertyShouldBe(page, "X", el("scatter plot viewer"), "WEIGHT"));
    await session.step(75, "When user saves the layout of the current table view", () => saveLayout(page));
    await session.step(76, "And user loads the saved layout", () => loadLayout(page));
    await session.step(77, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(78, "And \"X\" property of scatter plot viewer should be \"WEIGHT\"", () => propertyShouldBe(page, "X", el("scatter plot viewer"), "WEIGHT"));
    await session.step(79, "And \"X Axis Type\" property of scatter plot viewer should be \"logarithmic\"", () => propertyShouldBe(page, "X Axis Type", el("scatter plot viewer"), "logarithmic"));
    await session.step(80, "And no errors should have been logged", () => noErrors(page));
    await session.step(81, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
