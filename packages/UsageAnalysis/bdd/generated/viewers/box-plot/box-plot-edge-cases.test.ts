/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/box-plot/box-plot-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot, core#960]
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
import {hoverOver} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, deleteSelected, filterPasses, filterPassesAll, rowCount, selectAllRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, hasNoArea, noErrors, painted, pointerAway, reportsNoError, resizeTo, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Box plot edge cases", () => {
  const session = feature(test, "features/viewers/box-plot/box-plot-edge-cases.feature", import.meta.url);
  test("A box plot over a table with every row deleted logs no error", {tag: ["@viewers", "@realizes:viewers.box-plot", "@edge-cases", "@realizes:core#960"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","SEX"]]), [["Value","AGE"],["Category 1","SEX"]]);
    await session.step(18, "Then box plot viewer should be painted", () => painted(page, el("box plot viewer")));
    await session.step(22, "When user selects all rows", () => selectAllRows(page));
    await session.step(23, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(24, "Then the table should have 0 rows", () => rowCount(page, 0));
    await session.step(25, "And box plot viewer should report no error", () => reportsNoError(page, el("box plot viewer")));
    await session.step(26, "And no errors should have been logged", () => noErrors(page));
  });
  test("The pointer over a box plot shrunk below its axes logs no error", {tag: ["@viewers", "@realizes:viewers.box-plot", "@edge-cases", "@realizes:core#960"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","SEX"]]), [["Value","AGE"],["Category 1","SEX"]]);
    await session.step(18, "Then box plot viewer should be painted", () => painted(page, el("box plot viewer")));
    await session.step(30, "When user resizes box plot viewer to 170 by 150", () => resizeTo(page, el("box plot viewer"), 170, 150));
    await session.step(31, "Then box plot viewer should not have a \"x axis\" area", () => hasNoArea(page, el("box plot viewer"), "x axis"));
    await session.step(32, "When user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
    await session.step(33, "And user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
    await session.step(34, "And user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
    await session.step(35, "Then box plot viewer should be painted", () => painted(page, el("box plot viewer")));
    await session.step(36, "And no errors should have been logged", () => noErrors(page));
  });
  test("The group comparison goes with a single group left and comes back with two", {tag: ["@viewers", "@realizes:viewers.box-plot", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","SEX"]]), [["Value","AGE"],["Category 1","SEX"]]);
    await session.step(18, "Then box plot viewer should be painted", () => painted(page, el("box plot viewer")));
    await session.step(40, "When user sets \"Show Group Comparison\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Group Comparison", el("box plot viewer"), "true"));
    await session.step(41, "Then box plot viewer should have a \"group comparison\" area", () => hasArea(page, el("box plot viewer"), "group comparison"));
    await session.step(42, "When user adds a categorical filter on \"SEX\" keeping \"F\"", () => addCategoricalFilter(page, "SEX", "F"));
    await session.step(43, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
    await session.step(44, "And box plot viewer should not have a \"group comparison\" area", () => hasNoArea(page, el("box plot viewer"), "group comparison"));
    await session.step(45, "When user adds a categorical filter on \"SEX\" keeping \"F,M\"", () => addCategoricalFilter(page, "SEX", "F,M"));
    await session.step(46, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(47, "And box plot viewer should have a \"group comparison\" area", () => hasArea(page, el("box plot viewer"), "group comparison"));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
  });
});
