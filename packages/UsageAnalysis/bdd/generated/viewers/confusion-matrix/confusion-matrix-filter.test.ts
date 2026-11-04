/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/confusion-matrix/confusion-matrix-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.confusion-matrix, core#233]
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
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {addCategoricalFilter, filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, noErrors, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The confusion matrix filters the way the other viewers do", () => {
  const session = feature(test, "features/viewers/confusion-matrix/confusion-matrix-filter.feature", import.meta.url);
  test("The confusion matrix filters the way the other viewers do", {tag: ["@journey", "@viewers", "@realizes:viewers.confusion-matrix", "@realizes:core#233"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a confusion matrix viewer with:", () => addViewerWith(page, "confusion matrix", [["xColumnName","SEX"],["yColumnName","CONTROL"]]), [["xColumnName","SEX"],["yColumnName","CONTROL"]]);
    await run.scenario("Closing the matrix lets go of its filter", async () => {
      await session.step(18, "When user clicks on the \"cell M as false\" area of confusion matrix viewer", () => clickArea(page, "cell M as false", el("confusion matrix viewer")));
      await session.step(19, "Then 446 rows should pass the filter", () => filterPasses(page, 446));
      await session.step(20, "When user clicks on close icon of confusion matrix viewer", () => clickOn(page, el("close icon of confusion matrix viewer")));
      await session.step(21, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(22, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second click takes the filter off, another cell replaces it", async () => {
      await session.step(26, "When user adds a confusion matrix viewer with:", () => addViewerWith(page, "confusion matrix", [["xColumnName","SEX"],["yColumnName","CONTROL"]]), [["xColumnName","SEX"],["yColumnName","CONTROL"]]);
      await session.step(29, "And user clicks on the \"cell M as false\" area of confusion matrix viewer", () => clickArea(page, "cell M as false", el("confusion matrix viewer")));
      await session.step(30, "Then 446 rows should pass the filter", () => filterPasses(page, 446));
      await session.step(31, "When user clicks on the \"cell M as false\" area of confusion matrix viewer", () => clickArea(page, "cell M as false", el("confusion matrix viewer")));
      await session.step(32, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(33, "When user clicks on the \"cell M as false\" area of confusion matrix viewer", () => clickArea(page, "cell M as false", el("confusion matrix viewer")));
      await session.step(34, "Then 446 rows should pass the filter", () => filterPasses(page, 446));
      await session.step(35, "When user clicks on the \"cell F as false\" area of confusion matrix viewer", () => clickArea(page, "cell F as false", el("confusion matrix viewer")));
      await session.step(36, "Then 548 rows should pass the filter", () => filterPasses(page, 548));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Changing a column and resetting the filter let go of it", async () => {
      await session.step(41, "When user sets \"yColumnName\" property of confusion matrix viewer to \"DIS_POP\"", () => setProperty(page, "yColumnName", el("confusion matrix viewer"), "DIS_POP"));
      await session.step(42, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(43, "When user sets \"yColumnName\" property of confusion matrix viewer to \"CONTROL\"", () => setProperty(page, "yColumnName", el("confusion matrix viewer"), "CONTROL"));
      await session.step(44, "And user clicks on the \"cell M as false\" area of confusion matrix viewer", () => clickArea(page, "cell M as false", el("confusion matrix viewer")));
      await session.step(45, "Then 446 rows should pass the filter", () => filterPasses(page, 446));
      await session.step(46, "When user picks \"View > Reset Filter\" from the top menu", () => pickFromTopMenu(page, "View > Reset Filter"));
      await session.step(47, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(48, "When user clicks on the \"cell M as false\" area of confusion matrix viewer", () => clickArea(page, "cell M as false", el("confusion matrix viewer")));
      await session.step(49, "Then 446 rows should pass the filter", () => filterPasses(page, 446));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The matrix's filter and the filter panel apply together, whichever comes first", async () => {
      await session.step(54, "When user adds a categorical filter on \"RACE\" keeping \"Asian\"", () => addCategoricalFilter(page, "RACE", "Asian"));
      await session.step(55, "Then 8 rows should pass the filter", () => filterPasses(page, 8));
      await session.step(56, "When user picks \"View > Reset Filter\" from the top menu", () => pickFromTopMenu(page, "View > Reset Filter"));
      await session.step(57, "And user clicks on the \"cell M as false\" area of confusion matrix viewer", () => clickArea(page, "cell M as false", el("confusion matrix viewer")));
      await session.step(58, "And user adds a categorical filter on \"RACE\" keeping \"Asian\"", () => addCategoricalFilter(page, "RACE", "Asian"));
      await session.step(59, "Then 8 rows should pass the filter", () => filterPasses(page, 8));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
