/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/data/curve-json-columns.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/curves.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded, pasteInto, pressKey, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnsExactly, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {currentColumnIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, closeAllViews, contextPanelOpen, dialogCloses, openDataset, openTableOf, packageInstalled, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Curve JSON typed into a formula, a series without a name, and a broken curve cell", () => {
  const session = feature(test, "features/data/curve-json-columns.feature", import.meta.url);
  test("Curve JSON typed into a formula, a series without a name, and a broken curve cell", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await run.scenario("A JSON formula adds its column and Add New Column opens again", async () => {
      await session.step(20, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(21, "When user picks \"Edit > Add New Column...\" from the top menu", () => pickFromTopMenu(page, "Edit > Add New Column..."));
      await session.step(22, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(23, "When user pastes \"{}\" into code editor in \"Add New Column\" dialog", () => pasteInto(page, "{}", el("code editor in \"Add New Column\" dialog")));
      await session.step(24, "And user types \"jsoncol\" into column name input", () => typeInto(page, "jsoncol", el("column name input")));
      await session.step(25, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(26, "Then the \"Add New Column\" dialog should close", () => dialogCloses(page, "Add New Column"));
      await session.step(27, "And the table should have a column \"jsoncol\"", () => hasColumn(page, "jsoncol"));
      await session.step(28, "And no errors should have been logged", () => noErrors(page));
      await session.step(29, "When user picks \"Edit > Add New Column...\" from the top menu", () => pickFromTopMenu(page, "Edit > Add New Column..."));
      await session.step(30, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
      await session.step(31, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(32, "Then the \"Add New Column\" dialog should close", () => dialogCloses(page, "Add New Column"));
    });
    await run.scenario("The statistic of a series without a name becomes a column", async () => {
      await session.step(35, "Given user closes all views", () => closeAllViews(page));
      await session.step(36, "And user opens a table \"unnamed-series\" with:", () => openTableOf(page, "unnamed-series", [["id","curve"],["1","{\"series\":[{\"points\":[{\"x\":0.1,\"y\":1},{\"x\":1,\"y\":5},{\"x\":10,\"y\":9}]}]}"],["2","{\"series\":[{\"points\":[{\"x\":0.1,\"y\":2},{\"x\":1,\"y\":6},{\"x\":10,\"y\":8}]}]}"]]), [["id","curve"],["1","{\"series\":[{\"points\":[{\"x\":0.1,\"y\":1},{\"x\":1,\"y\":5},{\"x\":10,\"y\":9}]}]}"],["2","{\"series\":[{\"points\":[{\"x\":0.1,\"y\":2},{\"x\":1,\"y\":6},{\"x\":10,\"y\":8}]}]}"]]);
      await session.step(40, "And the context panel is open", () => contextPanelOpen(page));
      await session.step(41, "When user clicks on the \"cell 1 of curve\" area of grid", () => clickArea(page, "cell 1 of curve", el("grid")));
      await session.step(42, "Then the current column should be \"curve\"", () => currentColumnIs(page, "curve"));
      await session.step(43, "Given \"Fit\" pane in context panel is expanded", () => isExpanded(page, el("\"Fit\" pane in context panel")));
      await session.step(44, "When user clicks on first \"Calculate rSquared for the whole column\" icon in context panel", () => clickOn(page, el("first \"Calculate rSquared for the whole column\" icon in context panel")));
      await session.step(45, "Then the table should have the columns \"id, curve, curve 1 rSquared\"", () => columnsExactly(page, "id, curve, curve 1 rSquared"));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
      await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A click on a truncated curve cell logs no error", async () => {
      await session.step(50, "Given user closes all views", () => closeAllViews(page));
      await session.step(51, "And user opens a table \"broken-curves\" with:", () => openTableOf(page, "broken-curves", [["id","curve"],["1","{\"series\":[{\"points\":[{\"x\":0.1,\"y\":1},{\"x\":1,\"y\":5},{\"x\":10,\"y\":9}]}]}"],["2","{\"series\":[{\"points\":[{\"x\":0.1"],["3","{\"series\":[{\"points\":[{\"x\":0.1,\"y\":2},{\"x\":1,\"y\":6},{\"x\":10,\"y\":8}]}]}"]]), [["id","curve"],["1","{\"series\":[{\"points\":[{\"x\":0.1,\"y\":1},{\"x\":1,\"y\":5},{\"x\":10,\"y\":9}]}]}"],["2","{\"series\":[{\"points\":[{\"x\":0.1"],["3","{\"series\":[{\"points\":[{\"x\":0.1,\"y\":2},{\"x\":1,\"y\":6},{\"x\":10,\"y\":8}]}]}"]]);
      await session.step(56, "Then \"curve\" column should have semantic type \"fit\"", () => columnSemType(page, "curve", "fit"));
      await session.step(57, "When user clicks on the \"cell 2 of curve\" area of grid", () => clickArea(page, "cell 2 of curve", el("grid")));
      await session.step(58, "Then the current column should be \"curve\"", () => currentColumnIs(page, "curve"));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
      await session.step(60, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
