/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/transform/convert-notation-once.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [GROK-17964]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, selectIn, shouldBe, shouldNotBe, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {columnIsCurrentObject} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, contextPanelOpen, contextPanelShows, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Convert Notation is offered once in a molecule column's Actions", () => {
  const session = feature(test, "features/transform/convert-notation-once.feature", import.meta.url);
  test("Convert Notation is offered once in a molecule column's Actions", {tag: ["@journey", "@realizes:GROK-17964"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(10, "And user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
    await session.step(11, "And the context panel is open", () => contextPanelOpen(page));
    await run.scenario("The action is listed once for the molecule column", async () => {
      await session.step(14, "Given the \"canonical_smiles\" column is the current object", () => columnIsCurrentObject(page, "canonical_smiles"));
      await session.step(15, "Then the context panel should show \"canonical_smiles\"", () => contextPanelShows(page, "canonical_smiles"));
      await session.step(16, "When user expands Actions accordion header in context panel", () => expand(page, el("Actions accordion header in context panel")));
      await session.step(17, "Then there should be 1 visible \"Convert Notation...\" action in context panel", () => visibleCount(page, 1, el("\"Convert Notation...\" action in context panel")));
      await session.step(18, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Cancelling the action's dialog leaves one action", async () => {
      await session.step(21, "When user clicks on \"Convert Notation...\" action in context panel", () => clickOn(page, el("\"Convert Notation...\" action in context panel")));
      await session.step(22, "Then \"Convert Notation\" dialog should be visible", () => shouldBe(page, el("\"Convert Notation\" dialog"), "visible"));
      await session.step(23, "And Overwrite input in \"Convert Notation\" dialog should not be checked", () => shouldNotBe(page, el("Overwrite input in \"Convert Notation\" dialog"), "checked"));
      await session.step(24, "And Join input in \"Convert Notation\" dialog should be checked", () => shouldBe(page, el("Join input in \"Convert Notation\" dialog"), "checked"));
      await session.step(25, "And Kekulize input in \"Convert Notation\" dialog should not be checked", () => shouldNotBe(page, el("Kekulize input in \"Convert Notation\" dialog"), "checked"));
      await session.step(26, "When user clicks on CANCEL button in \"Convert Notation\" dialog", () => clickOn(page, el("CANCEL button in \"Convert Notation\" dialog")));
      await session.step(27, "And user clicks on the \"header canonical_smiles\" area of grid", () => clickArea(page, "header canonical_smiles", el("grid")));
      await session.step(28, "Then the context panel should show \"canonical_smiles\"", () => contextPanelShows(page, "canonical_smiles"));
      await session.step(29, "When user expands Actions accordion header in context panel", () => expand(page, el("Actions accordion header in context panel")));
      await session.step(30, "Then there should be 1 visible \"Convert Notation...\" action in context panel", () => visibleCount(page, 1, el("\"Convert Notation...\" action in context panel")));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Running a conversion from the action leaves one action", async () => {
      await session.step(34, "When user clicks on \"Convert Notation...\" action in context panel", () => clickOn(page, el("\"Convert Notation...\" action in context panel")));
      await session.step(35, "And user selects \"molblock\" in \"Target Notation\" input in \"Convert Notation\" dialog", () => selectIn(page, "molblock", el("\"Target Notation\" input in \"Convert Notation\" dialog")));
      await session.step(36, "And user clicks on OK button in \"Convert Notation\" dialog", () => clickOn(page, el("OK button in \"Convert Notation\" dialog")));
      await session.step(37, "Then \"Convert Notation\" dialog should be hidden", () => shouldBe(page, el("\"Convert Notation\" dialog"), "hidden"));
      await session.step(38, "And the table should have a column \"canonical_smiles_molblock\"", () => hasColumn(page, "canonical_smiles_molblock"));
      await session.step(39, "Given the \"canonical_smiles\" column is the current object", () => columnIsCurrentObject(page, "canonical_smiles"));
      await session.step(40, "Then the context panel should show \"canonical_smiles\"", () => contextPanelShows(page, "canonical_smiles"));
      await session.step(41, "When user expands Actions accordion header in context panel", () => expand(page, el("Actions accordion header in context panel")));
      await session.step(42, "Then there should be 1 visible \"Convert Notation...\" action in context panel", () => visibleCount(page, 1, el("\"Convert Notation...\" action in context panel")));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Opening and cancelling the dialog twice more leaves one action", async () => {
      await session.step(46, "When user clicks on \"Convert Notation...\" action in context panel", () => clickOn(page, el("\"Convert Notation...\" action in context panel")));
      await session.step(47, "And user clicks on CANCEL button in \"Convert Notation\" dialog", () => clickOn(page, el("CANCEL button in \"Convert Notation\" dialog")));
      await session.step(48, "And user clicks on \"Convert Notation...\" action in context panel", () => clickOn(page, el("\"Convert Notation...\" action in context panel")));
      await session.step(49, "And user clicks on CANCEL button in \"Convert Notation\" dialog", () => clickOn(page, el("CANCEL button in \"Convert Notation\" dialog")));
      await session.step(50, "And user clicks on the \"header canonical_smiles\" area of grid", () => clickArea(page, "header canonical_smiles", el("grid")));
      await session.step(51, "Then the context panel should show \"canonical_smiles\"", () => contextPanelShows(page, "canonical_smiles"));
      await session.step(52, "When user expands Actions accordion header in context panel", () => expand(page, el("Actions accordion header in context panel")));
      await session.step(53, "Then there should be 1 visible \"Convert Notation...\" action in context panel", () => visibleCount(page, 1, el("\"Convert Notation...\" action in context panel")));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
