/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/composition.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.analyze.composition]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, isExpanded, shouldBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {clearSelection, noneSelected, onlyStartingWithSelected, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDatasetRows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, clickArea, hasArea, noBalloons, noErrors, painted, propertyShouldBe, readingAtLeast, readingHigher, readingIs, repaintedBy, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Composition analysis", () => {
  const session = feature(test, "features/analyze/composition.feature", import.meta.url);
  test("Composition analysis", {tag: ["@journey", "@realizes:bio.analyze.composition"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens filter_FASTA dataset keeping the first 9 rows", () => openDatasetRows(page, ds("filter_FASTA"), 9));
    await session.step(14, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(15, "Then \"fasta\" column should have units \"fasta\"", () => columnUnits(page, "fasta", "fasta"));
    await run.scenario("The command docks a WebLogo bound to the sequence column, with no dialog", async () => {
      await session.step(18, "When user picks \"Bio > Analyze > Composition\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Composition"));
      await session.step(19, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(20, "And WebLogo viewer should be visible", () => shouldBe(page, el("WebLogo viewer"), "visible"));
      await session.step(21, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(22, "And WebLogo viewer should be bound to table \"filter_FASTA\"", () => boundTable(page, el("WebLogo viewer"), "filter_FASTA"));
      await session.step(23, "And \"Sequence Column Name\" property of WebLogo viewer should be \"fasta\"", () => propertyShouldBe(page, "Sequence Column Name", el("WebLogo viewer"), "fasta"));
      await session.step(24, "And WebLogo viewer should be painted", () => painted(page, el("WebLogo viewer")));
      await session.step(25, "And WebLogo viewer should have a \"position 1\" area", () => hasArea(page, el("WebLogo viewer"), "position 1"));
      await session.step(26, "And WebLogo viewer should have a \"monomer M at position 1\" area", () => hasArea(page, el("WebLogo viewer"), "monomer M at position 1"));
      await session.step(27, "And the \"positions shown\" reading of WebLogo viewer should be at least 30", () => readingAtLeast(page, "positions shown", el("WebLogo viewer"), 30));
      await session.step(28, "And the \"rows shown\" reading of WebLogo viewer should be 9", () => readingIs(page, "rows shown", el("WebLogo viewer"), 9));
      await session.step(29, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A click on a glyph selects the rows with that monomer at that position", async () => {
      await session.step(32, "Given no rows should be selected", () => noneSelected(page));
      await session.step(33, "When user clicks on the \"monomer M at position 1\" area of WebLogo viewer", () => clickArea(page, "monomer M at position 1", el("WebLogo viewer")));
      await session.step(34, "Then some rows should be selected", () => someSelected(page));
      await session.step(35, "And the \"rows selected\" reading of WebLogo viewer should be higher than before", () => readingHigher(page, "rows selected", el("WebLogo viewer")));
      await session.step(36, "And only rows where \"fasta\" starts with \"M\" should be selected", () => onlyStartingWithSelected(page, "fasta", "M"));
      await session.step(37, "When user clears the row selection", () => clearSelection(page));
      await session.step(38, "Then no rows should be selected", () => noneSelected(page));
    });
    await run.scenario("The gear opens the viewer's properties in the context panel", async () => {
      await session.step(41, "When user clicks on settings icon of WebLogo viewer", () => clickOn(page, el("settings icon of WebLogo viewer")));
      await session.step(42, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(43, "And \"Show Position Labels\" property in context panel should be present", () => shouldBe(page, el("\"Show Position Labels\" property in context panel"), "present"));
      await session.step(44, "Given \"Layout\" category in context panel is expanded", () => isExpanded(page, el("\"Layout\" category in context panel")));
      await session.step(45, "Then \"Show Position Labels\" property in context panel should be visible", () => shouldBe(page, el("\"Show Position Labels\" property in context panel"), "visible"));
      await session.step(46, "And \"Show Position Labels\" property of WebLogo viewer should be \"true\"", () => propertyShouldBe(page, "Show Position Labels", el("WebLogo viewer"), "true"));
    });
    await run.scenario("Switching a property off in the context panel repaints the logo", async () => {
      await session.step(49, "When user takes a snapshot of WebLogo viewer", () => takeSnapshot(page, el("WebLogo viewer")));
      await session.step(50, "And user unchecks \"Show Position Labels\" property in context panel", () => uncheck(page, el("\"Show Position Labels\" property in context panel")));
      await session.step(51, "Then \"Show Position Labels\" property of WebLogo viewer should be \"false\"", () => propertyShouldBe(page, "Show Position Labels", el("WebLogo viewer"), "false"));
      await session.step(52, "And WebLogo viewer should have repainted by at least 100 pixels", () => repaintedBy(page, el("WebLogo viewer"), 100));
      await session.step(53, "When user takes a snapshot of WebLogo viewer", () => takeSnapshot(page, el("WebLogo viewer")));
      await session.step(54, "And user checks \"Show Position Labels\" property in context panel", () => check(page, el("\"Show Position Labels\" property in context panel")));
      await session.step(55, "Then \"Show Position Labels\" property of WebLogo viewer should be \"true\"", () => propertyShouldBe(page, "Show Position Labels", el("WebLogo viewer"), "true"));
      await session.step(56, "And WebLogo viewer should have repainted by at least 100 pixels", () => repaintedBy(page, el("WebLogo viewer"), 100));
      await session.step(57, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
