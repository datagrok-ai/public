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
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {clearSelection, noneSelected, onlyStartingWithSelected, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDatasetRows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, clickArea, hasArea, noBalloons, noErrors, painted, propertyShouldBe, readingAtLeast, readingHigher, readingIs, repainted, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Composition analysis", () => {
  const session = feature(test, "features/analyze/composition.feature", import.meta.url);
  test("Composition analysis", {tag: ["@journey", "@realizes:bio.analyze.composition"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens filter_FASTA dataset keeping the first 9 rows", () => openDatasetRows(page, ds("filter_FASTA"), 9));
    await session.step(11, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(12, "Then \"fasta\" column should have units \"fasta\"", () => columnUnits(page, "fasta", "fasta"));
    await run.scenario("The command docks a WebLogo bound to the sequence column, with no dialog", async () => {
      await session.step(15, "When user picks \"Bio > Analyze > Composition\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Composition"));
      await session.step(16, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(17, "And WebLogo viewer should be visible", () => shouldBe(page, el("WebLogo viewer"), "visible"));
      await session.step(18, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(19, "And WebLogo viewer should be bound to table \"filter_FASTA\"", () => boundTable(page, el("WebLogo viewer"), "filter_FASTA"));
      await session.step(20, "And \"Sequence Column Name\" property of WebLogo viewer should be \"fasta\"", () => propertyShouldBe(page, "Sequence Column Name", el("WebLogo viewer"), "fasta"));
      await session.step(21, "And WebLogo viewer should be painted", () => painted(page, el("WebLogo viewer")));
      await session.step(22, "And WebLogo viewer should have a \"position 1\" area", () => hasArea(page, el("WebLogo viewer"), "position 1"));
      await session.step(23, "And WebLogo viewer should have a \"monomer M at position 1\" area", () => hasArea(page, el("WebLogo viewer"), "monomer M at position 1"));
      await session.step(24, "And the \"positions shown\" reading of WebLogo viewer should be at least 30", () => readingAtLeast(page, "positions shown", el("WebLogo viewer"), 30));
      await session.step(25, "And the \"rows shown\" reading of WebLogo viewer should be 9", () => readingIs(page, "rows shown", el("WebLogo viewer"), 9));
      await session.step(26, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A click on a glyph selects the rows with that monomer at that position", async () => {
      await session.step(29, "Given no rows should be selected", () => noneSelected(page));
      await session.step(30, "When user clicks on the \"monomer M at position 1\" area of WebLogo viewer", () => clickArea(page, "monomer M at position 1", el("WebLogo viewer")));
      await session.step(31, "Then some rows should be selected", () => someSelected(page));
      await session.step(32, "And the \"rows selected\" reading of WebLogo viewer should be higher than before", () => readingHigher(page, "rows selected", el("WebLogo viewer")));
      await session.step(33, "And only rows where \"fasta\" starts with \"M\" should be selected", () => onlyStartingWithSelected(page, "fasta", "M"));
      await session.step(34, "When user clears the row selection", () => clearSelection(page));
      await session.step(35, "Then no rows should be selected", () => noneSelected(page));
    });
    await run.scenario("The gear opens the viewer's properties in the context panel", async () => {
      await session.step(38, "When user clicks on settings icon of WebLogo viewer", () => clickOn(page, el("settings icon of WebLogo viewer")));
      await session.step(39, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(40, "And \"Show Position Labels\" property in context panel should be present", () => shouldBe(page, el("\"Show Position Labels\" property in context panel"), "present"));
      await session.step(41, "When user clicks on \"Layout\" category in context panel", () => clickOn(page, el("\"Layout\" category in context panel")));
      await session.step(42, "Then \"Show Position Labels\" property in context panel should be visible", () => shouldBe(page, el("\"Show Position Labels\" property in context panel"), "visible"));
      await session.step(43, "And \"Show Position Labels\" property of WebLogo viewer should be \"true\"", () => propertyShouldBe(page, "Show Position Labels", el("WebLogo viewer"), "true"));
    });
    await run.scenario("A property change repaints the logo", async () => {
      await session.step(46, "When user sets \"Show Position Labels\" property of WebLogo viewer to \"false\"", () => setProperty(page, "Show Position Labels", el("WebLogo viewer"), "false"));
      await session.step(47, "Then \"Show Position Labels\" property of WebLogo viewer should be \"false\"", () => propertyShouldBe(page, "Show Position Labels", el("WebLogo viewer"), "false"));
      await session.step(48, "And WebLogo viewer should have repainted", () => repainted(page, el("WebLogo viewer")));
      await session.step(49, "When user sets \"Show Position Labels\" property of WebLogo viewer to \"true\"", () => setProperty(page, "Show Position Labels", el("WebLogo viewer"), "true"));
      await session.step(50, "Then WebLogo viewer should have repainted", () => repainted(page, el("WebLogo viewer")));
      await session.step(51, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
