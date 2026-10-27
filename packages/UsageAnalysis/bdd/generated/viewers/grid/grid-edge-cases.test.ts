/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid, core#1137]
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
import {clickOn, pressKey} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hasArea, noErrors, pickFromAreaContextMenu, readingReads, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid edge cases", () => {
  const session = feature(test, "features/viewers/grid/grid-edge-cases.feature", import.meta.url);
  test("An HTML summary column keeps rendering after a column it shows is removed", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(17, "When user picks \"Add > Summary Columns > Default HTML Form\" from the context menu of the \"cell 2 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Add > Summary Columns > Default HTML Form", "cell 2 of USUBJID", el("grid")));
    await session.step(18, "And user clicks on All link in \"Select columns...\" dialog", () => clickOn(page, el("All link in \"Select columns...\" dialog")));
    await session.step(19, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(20, "Then the \"cell type of html\" reading of grid should be \"html\"", () => readingReads(page, "cell type of html", el("grid"), "html"));
    await session.step(21, "When user picks \"Remove\" from the context menu of the \"header CONTROL\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header CONTROL", el("grid")));
    await session.step(22, "Then the table should not have a column \"CONTROL\"", () => hasNoColumn(page, "CONTROL"));
    await session.step(23, "When user clicks on the \"cell 2 of AGE\" area of grid", () => clickArea(page, "cell 2 of AGE", el("grid")));
    await session.step(24, "And user presses Control+End", () => pressKey(page, "Control+End"));
    await session.step(25, "Then grid should have a \"cell 1000 of html\" area", () => hasArea(page, el("grid"), "cell 1000 of html"));
    await session.step(26, "And the \"cell type of html\" reading of grid should be \"html\"", () => readingReads(page, "cell type of html", el("grid"), "html"));
    await session.step(27, "And no errors should have been logged", () => noErrors(page));
  });
  test("Unpin Column with Frozen Columns above the column count logs no error", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:36", "@realizes:core#1137"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(34, "When user sets \"Frozen Columns\" property of grid to \"20\"", () => setProperty(page, "Frozen Columns", el("grid"), "20"));
      await session.step(35, "And user picks \"Pin > Unpin Column\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Unpin Column", "header AGE", el("grid")));
      await session.step(36, "Then no errors should have been logged", () => noErrors(page));
    }, 36);
  });
});
