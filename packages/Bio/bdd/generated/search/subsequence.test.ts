/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/search/subsequence.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.search.subsequence]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '../../bindings/monomer-form.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, hoverOver, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterIsExactlyContains, filterPanelCount, filterPanelHas, filterPasses, filterPassesAll, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {panelHasNoCardOfType} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {noBalloons, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickFromViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Subsequence search on the filter panel", () => {
  const session = feature(test, "features/search/subsequence.feature", import.meta.url);
  test("Subsequence search on the filter panel", {tag: ["@journey", "@realizes:bio.search.subsequence"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens filter_FASTA dataset", () => openDataset(page, ds("filter_FASTA")));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(17, "Then the table should have 14 rows", () => rowCount(page, 14));
    await session.step(18, "And \"fasta\" column should have units \"fasta\"", () => columnUnits(page, "fasta", "fasta"));
    await run.scenario("The command adds a filter bound to the sequence column", async () => {
      await session.step(21, "When user picks \"Bio > Search > Subsequence Search ...\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Subsequence Search ..."));
      await session.step(22, "Then filters viewer should be visible", () => shouldBe(page, el("filters viewer"), "visible"));
      await session.step(23, "And \"Substructure\" input in filters viewer should be visible", () => shouldBe(page, el("\"Substructure\" input in filters viewer"), "visible"));
      await session.step(24, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
      await session.step(25, "And the filter panel should have a filter on \"fasta\" column", () => filterPanelHas(page, "fasta"));
    });
    await run.scenario("A subsequence one row contains keeps that row alone", async () => {
      await session.step(28, "When user enters \"RTDEVSNHTHDKPTLTWFEEIFEEYHSP\" into \"Substructure\" input in filters viewer", () => enterInto(page, "RTDEVSNHTHDKPTLTWFEEIFEEYHSP", el("\"Substructure\" input in filters viewer")));
      await session.step(29, "Then 1 row should pass the filter", () => filterPasses(page, 1));
      await session.step(30, "And the filter should pass exactly the rows where \"fasta\" contains \"RTDEVSNHTHDKPTLTWFEEIFEEYHSP\"", () => filterIsExactlyContains(page, "fasta", "RTDEVSNHTHDKPTLTWFEEIFEEYHSP"));
      await session.step(31, "And the table should have 14 rows", () => rowCount(page, 14));
    });
    await run.scenario("Reset restores every row and empties the query", async () => {
      await session.step(34, "When user clicks on \"arrow rotate left\" icon in filters viewer", () => clickOn(page, el("\"arrow rotate left\" icon in filters viewer")));
      await session.step(35, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(36, "And \"Substructure\" input in filters viewer should have value \"\"", () => shouldHaveValue(page, el("\"Substructure\" input in filters viewer"), ""));
      await session.step(37, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Add Filter on the panel puts the card back on the column picked", async () => {
      await session.step(41, "When user hovers over \"fasta\" filter card", () => hoverOver(page, el("\"fasta\" filter card")));
      await session.step(42, "And user clicks on close of \"fasta\" filter card", () => clickOn(page, el("close of \"fasta\" filter card")));
      await session.step(43, "Then the filter panel should have no \"Bio:bioSubstructureFilter\" filter card", () => panelHasNoCardOfType(page, "Bio:bioSubstructureFilter"));
      await session.step(44, "When user picks \"Add Filter | Bio Substructure Filter...\" from the viewer menu of filter panel", () => pickFromViewerMenu(page, "Add Filter | Bio Substructure Filter...", el("filter panel")));
      await session.step(45, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(46, "When user clicks on All label in \"Select columns...\" dialog", () => clickOn(page, el("All label in \"Select columns...\" dialog")));
      await session.step(47, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(48, "Then \"fasta\" filter card should be visible", () => shouldBe(page, el("\"fasta\" filter card"), "visible"));
      await session.step(49, "And the \"type of fasta\" reading of filter panel should be \"Bio:bioSubstructureFilter\"", () => readingReads(page, "type of fasta", el("filter panel"), "Bio:bioSubstructureFilter"));
      await session.step(50, "When user enters \"RTDEVSNHTHDKPTLTWFEEIFEEYHSP\" into \"Substructure\" input in filters viewer", () => enterInto(page, "RTDEVSNHTHDKPTLTWFEEIFEEYHSP", el("\"Substructure\" input in filters viewer")));
      await session.step(51, "Then 1 row should pass the filter", () => filterPasses(page, 1));
      await session.step(52, "When user clicks on \"arrow rotate left\" icon in filters viewer", () => clickOn(page, el("\"arrow rotate left\" icon in filters viewer")));
      await session.step(53, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
