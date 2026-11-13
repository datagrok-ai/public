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
import {clickOn, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {addCategoricalFilter, filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, hoverArea, noErrors, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The confusion matrix filters the way the other viewers do", () => {
  const session = feature(test, "features/viewers/confusion-matrix/confusion-matrix-filter.feature", import.meta.url);
  test("The confusion matrix filters the way the other viewers do", {tag: ["@journey", "@viewers", "@realizes:viewers.confusion-matrix", "@realizes:core#233"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a confusion matrix viewer with:", () => addViewerWith(page, "confusion matrix", [["xColumnName","SEX"],["yColumnName","CONTROL"]]), [["xColumnName","SEX"],["yColumnName","CONTROL"]]);
    await run.scenario("Closing the matrix lets go of its filter", async () => {
      await session.step(19, "When user hovers over the \"cell M as false\" area of confusion matrix viewer", () => hoverArea(page, "cell M as false", el("confusion matrix viewer")));
      await session.step(20, "Then tooltip should contain text \"Number of M predicted as false\"", () => shouldContainText(page, el("tooltip"), "Number of M predicted as false"));
      await session.step(21, "When user clicks on the \"cell M as false\" area of confusion matrix viewer", () => clickArea(page, "cell M as false", el("confusion matrix viewer")));
      await session.step(22, "Then 446 rows should pass the filter", () => filterPasses(page, 446));
      await session.step(23, "When user clicks on close icon of confusion matrix viewer", () => clickOn(page, el("close icon of confusion matrix viewer")));
      await session.step(24, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(25, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second click takes the filter off, another cell replaces it", async () => {
      await session.step(29, "When user adds a confusion matrix viewer with:", () => addViewerWith(page, "confusion matrix", [["xColumnName","SEX"],["yColumnName","CONTROL"]]), [["xColumnName","SEX"],["yColumnName","CONTROL"]]);
      await session.step(32, "And user clicks on the \"cell M as false\" area of confusion matrix viewer", () => clickArea(page, "cell M as false", el("confusion matrix viewer")));
      await session.step(33, "Then 446 rows should pass the filter", () => filterPasses(page, 446));
      await session.step(34, "When user clicks on the \"cell M as false\" area of confusion matrix viewer", () => clickArea(page, "cell M as false", el("confusion matrix viewer")));
      await session.step(35, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(36, "When user clicks on the \"cell M as false\" area of confusion matrix viewer", () => clickArea(page, "cell M as false", el("confusion matrix viewer")));
      await session.step(37, "Then 446 rows should pass the filter", () => filterPasses(page, 446));
      await session.step(38, "When user clicks on the \"cell F as false\" area of confusion matrix viewer", () => clickArea(page, "cell F as false", el("confusion matrix viewer")));
      await session.step(39, "Then 548 rows should pass the filter", () => filterPasses(page, 548));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Changing a column and resetting the filter let go of it", async () => {
      await session.step(44, "When user sets \"yColumnName\" property of confusion matrix viewer to \"DIS_POP\"", () => setProperty(page, "yColumnName", el("confusion matrix viewer"), "DIS_POP"));
      await session.step(45, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(46, "When user sets \"yColumnName\" property of confusion matrix viewer to \"CONTROL\"", () => setProperty(page, "yColumnName", el("confusion matrix viewer"), "CONTROL"));
      await session.step(47, "And user clicks on the \"cell M as false\" area of confusion matrix viewer", () => clickArea(page, "cell M as false", el("confusion matrix viewer")));
      await session.step(48, "Then 446 rows should pass the filter", () => filterPasses(page, 446));
      await session.step(49, "When user picks \"View > Reset Filter\" from the top menu", () => pickFromTopMenu(page, "View > Reset Filter"));
      await session.step(50, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(51, "When user clicks on the \"cell M as false\" area of confusion matrix viewer", () => clickArea(page, "cell M as false", el("confusion matrix viewer")));
      await session.step(52, "Then 446 rows should pass the filter", () => filterPasses(page, 446));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The matrix's filter and the filter panel apply together, whichever comes first", async () => {
      await session.step(57, "When user adds a categorical filter on \"RACE\" keeping \"Asian\"", () => addCategoricalFilter(page, "RACE", "Asian"));
      await session.step(58, "Then 8 rows should pass the filter", () => filterPasses(page, 8));
      await session.step(59, "When user picks \"View > Reset Filter\" from the top menu", () => pickFromTopMenu(page, "View > Reset Filter"));
      await session.step(60, "And user clicks on the \"cell M as false\" area of confusion matrix viewer", () => clickArea(page, "cell M as false", el("confusion matrix viewer")));
      await session.step(61, "And user adds a categorical filter on \"RACE\" keeping \"Asian\"", () => addCategoricalFilter(page, "RACE", "Asian"));
      await session.step(62, "Then 8 rows should pass the filter", () => filterPasses(page, 8));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
