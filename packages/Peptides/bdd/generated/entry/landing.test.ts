/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/entry/landing.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openLanding, peptidesInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeCurrentView, switchView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Peptides landing view", () => {
  const session = feature(test, "features/entry/landing.feature", import.meta.url);
  test("Peptides landing view", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(6, "Given user is logged in", () => loggedIn(page));
    await session.step(7, "And the Peptides package is initialized", () => peptidesInitialized(page));
    await session.step(8, "And user opens the Peptides landing view", () => openLanding(page));
    await session.step(9, "Then no errors should have been logged", () => noErrors(page));
    await session.step(10, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await run.scenario("The landing view offers three demos and folds the side panels", async () => {
      await session.step(13, "Then the \"Peptides\" view should be current", () => viewIsCurrent(page, "Peptides"));
      await session.step(14, "And \"Simple demo\" button should be visible", () => shouldBe(page, el("\"Simple demo\" button"), "visible"));
      await session.step(15, "And \"Complex demo\" button should be visible", () => shouldBe(page, el("\"Complex demo\" button"), "visible"));
      await session.step(16, "And \"HELM demo\" button should be visible", () => shouldBe(page, el("\"HELM demo\" button"), "visible"));
      await session.step(17, "And toolbox should be hidden", () => shouldBe(page, el("toolbox"), "hidden"));
      await session.step(18, "And context panel should be hidden", () => shouldBe(page, el("context panel"), "hidden"));
      await session.step(19, "And help panel should be hidden", () => shouldBe(page, el("help panel"), "hidden"));
      await session.step(20, "Then no errors should have been logged", () => noErrors(page));
      await session.step(21, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A demo button opens the dataset in its original notation [demo=Simple, rows=647, column=AlignedSequence, notation=fasta]", async () => {
      await session.step(24, "When user clicks on \"Simple demo\" button", () => clickOn(page, el("\"Simple demo\" button")));
      await session.step(25, "Then the \"PeptidesView\" view should be current", () => viewIsCurrent(page, "PeptidesView"));
      await session.step(26, "And the table should have 647 rows", () => rowCount(page, 647));
      await session.step(27, "And \"AlignedSequence\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "AlignedSequence", "Macromolecule"));
      await session.step(28, "And \"AlignedSequence\" column should have units \"fasta\"", () => columnUnits(page, "AlignedSequence", "fasta"));
      await session.step(29, "And context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(30, "When user closes the current view", () => closeCurrentView(page));
      await session.step(31, "And user switches to the \"Peptides\" view", () => switchView(page, "Peptides"));
      await session.step(32, "Then no errors should have been logged", () => noErrors(page));
      await session.step(33, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A demo button opens the dataset in its original notation [demo=Complex, rows=540, column=MSA, notation=separator]", async () => {
      await session.step(24, "When user clicks on \"Complex demo\" button", () => clickOn(page, el("\"Complex demo\" button")));
      await session.step(25, "Then the \"PeptidesView\" view should be current", () => viewIsCurrent(page, "PeptidesView"));
      await session.step(26, "And the table should have 540 rows", () => rowCount(page, 540));
      await session.step(27, "And \"MSA\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "MSA", "Macromolecule"));
      await session.step(28, "And \"MSA\" column should have units \"separator\"", () => columnUnits(page, "MSA", "separator"));
      await session.step(29, "And context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(30, "When user closes the current view", () => closeCurrentView(page));
      await session.step(31, "And user switches to the \"Peptides\" view", () => switchView(page, "Peptides"));
      await session.step(32, "Then no errors should have been logged", () => noErrors(page));
      await session.step(33, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A demo button opens the dataset in its original notation [demo=HELM, rows=334, column=HELM, notation=helm]", async () => {
      await session.step(24, "When user clicks on \"HELM demo\" button", () => clickOn(page, el("\"HELM demo\" button")));
      await session.step(25, "Then the \"PeptidesView\" view should be current", () => viewIsCurrent(page, "PeptidesView"));
      await session.step(26, "And the table should have 334 rows", () => rowCount(page, 334));
      await session.step(27, "And \"HELM\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "HELM", "Macromolecule"));
      await session.step(28, "And \"HELM\" column should have units \"helm\"", () => columnUnits(page, "HELM", "helm"));
      await session.step(29, "And context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(30, "When user closes the current view", () => closeCurrentView(page));
      await session.step(31, "And user switches to the \"Peptides\" view", () => switchView(page, "Peptides"));
      await session.step(32, "Then no errors should have been logged", () => noErrors(page));
      await session.step(33, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
