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
import {clickOn, pressKey, pressKeyIn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasNoColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, eventNotFired, hasArea, listenFor, noErrors, pickFromAreaContextMenu, readingReads, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid edge cases", () => {
  const session = feature(test, "features/viewers/grid/grid-edge-cases.feature", import.meta.url);
  test("An HTML summary column keeps rendering after a column it shows is removed", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(18, "When user picks \"Add > Summary Columns > Default HTML Form\" from the context menu of the \"cell 2 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Add > Summary Columns > Default HTML Form", "cell 2 of USUBJID", el("grid")));
    await session.step(19, "And user clicks on All link in \"Select columns...\" dialog", () => clickOn(page, el("All link in \"Select columns...\" dialog")));
    await session.step(20, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(21, "Then the \"cell type of html\" reading of grid should be \"html\"", () => readingReads(page, "cell type of html", el("grid"), "html"));
    await session.step(22, "When user picks \"Remove\" from the context menu of the \"header CONTROL\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header CONTROL", el("grid")));
    await session.step(23, "Then the table should not have a column \"CONTROL\"", () => hasNoColumn(page, "CONTROL"));
    await session.step(24, "When user clicks on the \"cell 2 of AGE\" area of grid", () => clickArea(page, "cell 2 of AGE", el("grid")));
    await session.step(25, "And user presses Control+End", () => pressKey(page, "Control+End"));
    await session.step(26, "Then grid should have a \"cell 1000 of html\" area", () => hasArea(page, el("grid"), "cell 1000 of html"));
    await session.step(27, "And the \"cell type of html\" reading of grid should be \"html\"", () => readingReads(page, "cell type of html", el("grid"), "html"));
    await session.step(28, "And no errors should have been logged", () => noErrors(page));
  });
  test("Unpin Column with Frozen Columns above the column count logs no error", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:37", "@realizes:core#1137"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(35, "When user sets \"Frozen Columns\" property of grid to \"20\"", () => setProperty(page, "Frozen Columns", el("grid"), "20"));
      await session.step(36, "And user picks \"Pin > Unpin Column\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Unpin Column", "header AGE", el("grid")));
      await session.step(37, "Then no errors should have been logged", () => noErrors(page));
    }, 37);
  });
  test("Text typed into a numeric cell is not reported as an edit", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@ticket-pending", "@known-failure:51"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(44, "Given user listens for \"d4-grid-cell-value-edited\" event on grid", () => listenFor(page, "d4-grid-cell-value-edited", el("grid")));
      await session.step(45, "When user double-clicks on the \"cell 4 of AGE\" area of grid", () => doubleClickArea(page, "cell 4 of AGE", el("grid")));
      await session.step(46, "And user presses Control+A in cell editor", () => pressKeyIn(page, "Control+A", el("cell editor")));
      await session.step(47, "And user types \"abc\" into cell editor", () => typeInto(page, "abc", el("cell editor")));
      await session.step(48, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(49, "Then cell editor should be hidden", () => shouldBe(page, el("cell editor"), "hidden"));
      await session.step(50, "And the value of \"AGE\" column in row 4 should be \"45\"", () => valueInRow(page, "AGE", 4, "45"));
      await session.step(51, "And \"d4-grid-cell-value-edited\" event should not have fired on grid", () => eventNotFired(page, "d4-grid-cell-value-edited", el("grid")));
    }, 51);
  });
});
