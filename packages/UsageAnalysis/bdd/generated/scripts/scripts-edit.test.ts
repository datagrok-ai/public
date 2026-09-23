/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scripts/scripts-edit.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.scripts]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {saveScript} from '../../bindings/scripts.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {appendToEditor, clearField, doubleClickOn, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeCurrentView, scriptContains, scriptOnServer, scriptsView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Editing a script", () => {
  const session = feature(test, "features/scripts/scripts-edit.feature", import.meta.url);
  test("Editing a script", {tag: ["@journey", "@serial", "@realizes:views.scripts"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And a script \"BddScriptEdit{time}\" is on the server:", () => scriptOnServer(page, session.text("BddScriptEdit{time}"), "#language: r\n#input: dataframe table\n#output: int count\ncount <- nrow(table) * ncol(table)"));
    await session.step(25, "And user opens the Scripts view", () => scriptsView(page));
    await run.scenario("A double-click opens the script in the editor, with nothing to save", async () => {
      await session.step(28, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(29, "And user types \"BddScriptEdit{time}\" into gallery search", () => typeInto(page, session.text("BddScriptEdit{time}"), el("gallery search")));
      await session.step(30, "And user double-clicks on \"BddScriptEdit{time}\" link in gallery", () => doubleClickOn(page, el(session.text("\"BddScriptEdit{time}\" link in gallery"))));
      await session.step(31, "Then the \"BddScriptEdit{time}\" view should be current", () => viewIsCurrent(page, session.text("BddScriptEdit{time}")));
      await session.step(32, "And code editor should contain the text \"count <- nrow(table) * ncol(table)\"", () => shouldContainText(page, el("code editor"), "count <- nrow(table) * ncol(table)"));
      await session.step(33, "And Save button should be disabled", () => shouldBe(page, el("Save button"), "disabled"));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An appended line is saved to the server", async () => {
      await session.step(37, "When user appends \"newParam = \\\"test\\\"\" to code editor", () => appendToEditor(page, "newParam = \"test\"", el("code editor")));
      await session.step(38, "Then Save button should be enabled", () => shouldBe(page, el("Save button"), "enabled"));
      await session.step(39, "When user saves the script", () => saveScript(page));
      await session.step(40, "Then an info balloon containing \"Script saved.\" should have been shown", () => infoBalloonText(page, "Script saved."));
      await session.step(41, "And the script \"BddScriptEdit{time}\" on the server should contain \"newParam = \\\"test\\\"\"", () => scriptContains(page, session.text("BddScriptEdit{time}"), "newParam = \"test\""));
      await session.step(42, "And Save button should be disabled", () => shouldBe(page, el("Save button"), "disabled"));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
      await session.step(44, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Opened again, the editor shows the saved line", async () => {
      await session.step(47, "When user closes the current view", () => closeCurrentView(page));
      await session.step(48, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
      await session.step(49, "When user double-clicks on \"BddScriptEdit{time}\" link in gallery", () => doubleClickOn(page, el(session.text("\"BddScriptEdit{time}\" link in gallery"))));
      await session.step(50, "Then the \"BddScriptEdit{time}\" view should be current", () => viewIsCurrent(page, session.text("BddScriptEdit{time}")));
      await session.step(51, "And code editor should contain the text \"newParam = \\\"test\\\"\"", () => shouldContainText(page, el("code editor"), "newParam = \"test\""));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
      await session.step(53, "When user closes the current view", () => closeCurrentView(page));
      await session.step(54, "And user clears gallery search", () => clearField(page, el("gallery search")));
    });
    run.finish();
  });
});
