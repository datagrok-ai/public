/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/nx/nx-regression-line.feature
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
import {selectWhereIs, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, closeContextMenu, hasArea, hasNoArea, noErrors, pickFromContextMenu, readingIs, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("NX — Tools > Show regression line on a scatter plot that shows only the selected rows", () => {
  const session = feature(test, "features/viewers/nx/nx-regression-line.feature", import.meta.url);
  test("The regression line is drawn on a scatter plot that shows only the selected rows", {tag: ["@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"]]), [["xColumnName","HEIGHT"],["yColumnName","WEIGHT"]]);
    await session.step(23, "And user selects rows where \"SEX\" is \"F\"", () => selectWhereIs(page, "SEX", "F"));
    await session.step(24, "Then 553 rows should be selected", () => selectedRowCount(page, 553));
    await session.step(25, "When user sets \"Row Source\" property of scatter plot viewer to \"FilteredSelected\"", () => setProperty(page, "Row Source", el("scatter plot viewer"), "FilteredSelected"));
    await session.step(26, "Then the \"rows shown\" reading of scatter plot viewer should be 499", () => readingIs(page, "rows shown", el("scatter plot viewer"), 499));
    await session.step(27, "And the \"regression lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "regression lines", el("scatter plot viewer"), 0));
    await session.step(28, "And scatter plot viewer should not have a \"regression line\" area", () => hasNoArea(page, el("scatter plot viewer"), "regression line"));
    await session.step(29, "When user picks \"Tools > Show Regression Line\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tools > Show Regression Line", el("scatter plot viewer")));
    await session.step(30, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(31, "Then the \"regression lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "regression lines", el("scatter plot viewer"), 1));
    await session.step(32, "And scatter plot viewer should have a \"regression line\" area", () => hasArea(page, el("scatter plot viewer"), "regression line"));
    await session.step(33, "And the \"rows shown\" reading of scatter plot viewer should be 499", () => readingIs(page, "rows shown", el("scatter plot viewer"), 499));
    await session.step(34, "And no errors should have been logged", () => noErrors(page));
  });
});
