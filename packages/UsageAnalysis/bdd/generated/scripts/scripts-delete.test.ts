/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scripts/scripts-delete.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.scripts]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {dialogCloses, scriptOnServer, scriptsOnServer, scriptsView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuLists, noBalloons, noErrors, openContextMenu, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Deleting a script", () => {
  const session = feature(test, "features/scripts/scripts-delete.feature", import.meta.url);
  test("Deleting a script", {tag: ["@journey", "@serial", "@realizes:views.scripts"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And a script \"BddScriptDelete{time}\" is on the server:", () => scriptOnServer(page, session.text("BddScriptDelete{time}"), "#language: r\n#input: dataframe table\n#output: int count\ncount <- nrow(table) * ncol(table)"));
    await session.step(23, "And user opens the Scripts view", () => scriptsView(page));
    await run.scenario("The context menu offers Delete", async () => {
      await session.step(26, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(27, "And user types \"BddScriptDelete{time}\" into gallery search", () => typeInto(page, session.text("BddScriptDelete{time}"), el("gallery search")));
      await session.step(28, "Then \"BddScriptDelete{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BddScriptDelete{time}\" link in gallery")), "visible"));
      await session.step(29, "When user opens the context menu of \"BddScriptDelete{time}\" link in gallery", () => openContextMenu(page, el(session.text("\"BddScriptDelete{time}\" link in gallery"))));
      await session.step(30, "Then the open menu should list \"Delete\"", () => menuLists(page, "Delete"));
      await session.step(31, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(32, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("CANCEL keeps the script", async () => {
      await session.step(35, "When user picks \"Delete\" from the context menu of \"BddScriptDelete{time}\" link in gallery", () => pickFromContextMenu(page, "Delete", el(session.text("\"BddScriptDelete{time}\" link in gallery"))));
      await session.step(36, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(37, "And \"Are you sure?\" dialog should contain text \"Delete script \\\"BddScriptDelete{time}\\\"?\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete script \"BddScriptDelete{time}\"?")));
      await session.step(38, "When user clicks on CANCEL button in \"Are you sure?\" dialog", () => clickOn(page, el("CANCEL button in \"Are you sure?\" dialog")));
      await session.step(39, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(40, "And \"BddScriptDelete{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BddScriptDelete{time}\" link in gallery")), "visible"));
      await session.step(41, "And 1 script named \"BddScriptDelete{time}\" should be on the server", () => scriptsOnServer(page, 1, session.text("BddScriptDelete{time}")));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("YES deletes the script", async () => {
      await session.step(45, "When user picks \"Delete\" from the context menu of \"BddScriptDelete{time}\" link in gallery", () => pickFromContextMenu(page, "Delete", el(session.text("\"BddScriptDelete{time}\" link in gallery"))));
      await session.step(46, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(47, "When user clicks on YES button in \"Are you sure?\" dialog", () => clickOn(page, el("YES button in \"Are you sure?\" dialog")));
      await session.step(48, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(49, "And \"BddScriptDelete{time}\" link in gallery should be absent", () => shouldBe(page, el(session.text("\"BddScriptDelete{time}\" link in gallery")), "absent"));
      await session.step(50, "And 0 scripts named \"BddScriptDelete{time}\" should be on the server", () => scriptsOnServer(page, 0, session.text("BddScriptDelete{time}")));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
      await session.step(52, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(53, "When user clears gallery search", () => clearField(page, el("gallery search")));
    });
    run.finish();
  });
});
