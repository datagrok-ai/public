/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scripts/script-run-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.scripts, core#1155]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, shouldBe, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {dialogCloses, scriptOnServer, scriptsView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noTableViewOpen} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Script run edge cases: Run... with no table to give", () => {
  const session = feature(test, "features/scripts/script-run-edge-cases.feature", import.meta.url);
  test("Run... with no table open keeps OK disabled", {tag: ["@serial", "@edge-cases", "@realizes:views.scripts", "@realizes:core#1155"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And a script \"BddScriptNoTable{time}\" is on the server:", () => scriptOnServer(page, session.text("BddScriptNoTable{time}"), "//language: javascript\n//input: dataframe table\n//output: int count\ncount = table.rowCount;"));
    await session.step(26, "And no table view is open", () => noTableViewOpen(page));
    await session.step(27, "And user opens the Scripts view", () => scriptsView(page));
    await session.step(31, "When user clears gallery search", () => clearField(page, el("gallery search")));
    await session.step(32, "And user types \"BddScriptNoTable{time}\" into gallery search", () => typeInto(page, session.text("BddScriptNoTable{time}"), el("gallery search")));
    await session.step(33, "And user picks \"Run...\" from the context menu of \"BddScriptNoTable{time}\" link in gallery", () => pickFromContextMenu(page, "Run...", el(session.text("\"BddScriptNoTable{time}\" link in gallery"))));
    await session.step(34, "Then \"BddScriptNoTable{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BddScriptNoTable{time}\" dialog")), "visible"));
    await session.step(35, "And Table input in \"BddScriptNoTable{time}\" dialog should have value \"\"", () => shouldHaveValue(page, el(session.text("Table input in \"BddScriptNoTable{time}\" dialog")), ""));
    await session.step(36, "And OK button in \"BddScriptNoTable{time}\" dialog should be disabled", () => shouldBe(page, el(session.text("OK button in \"BddScriptNoTable{time}\" dialog")), "disabled"));
    await session.step(37, "When user clicks on CANCEL button in \"BddScriptNoTable{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BddScriptNoTable{time}\" dialog"))));
    await session.step(38, "Then the \"BddScriptNoTable{time}\" dialog should close", () => dialogCloses(page, session.text("BddScriptNoTable{time}")));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
    await session.step(40, "When user clears gallery search", () => clearField(page, el("gallery search")));
  });
});
