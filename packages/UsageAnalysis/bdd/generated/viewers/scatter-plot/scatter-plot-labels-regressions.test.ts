/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-labels-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, public#3120, public#3323]
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
import {addCategoricalFilter, filterPasses, noneSelected, selectFirstRows, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, dragAreaBy, hasArea, noErrors, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot regressions — dragging a label", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-labels-regressions.feature", import.meta.url);
  test("A label is dragged under Mouse Drag = Select (public#3120)", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:public#3120"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
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
  test("A label is dragged when Show Labels For is All and a panel filter is set (public#3323)", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:public#3323"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(36, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"],["Label Columns","AGE"],["Show Labels For","All"],["mouseDrag","Select"]]), [["X","WEIGHT"],["Y","HEIGHT"],["Label Columns","AGE"],["Show Labels For","All"],["mouseDrag","Select"]]);
    await session.step(42, "When user adds a categorical filter on \"SEX\" keeping \"F\"", () => addCategoricalFilter(page, "SEX", "F"));
    await session.step(43, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
    await session.step(44, "And scatter plot viewer should have a \"label of row 1\" area", () => hasArea(page, el("scatter plot viewer"), "label of row 1"));
    await session.step(45, "And the \"labels shown\" reading of scatter plot viewer should be 274", () => readingIs(page, "labels shown", el("scatter plot viewer"), 274));
    await session.step(46, "When user drags the \"label of row 1\" area of scatter plot viewer by 60 pixels to the right", () => dragAreaBy(page, "label of row 1", el("scatter plot viewer"), 60, "right"));
    await session.step(47, "Then the \"labels shown\" reading of scatter plot viewer should be 301", () => readingIs(page, "labels shown", el("scatter plot viewer"), 301));
    await session.step(48, "And no rows should be selected", () => noneSelected(page));
    await session.step(49, "And no errors should have been logged", () => noErrors(page));
  });
});
