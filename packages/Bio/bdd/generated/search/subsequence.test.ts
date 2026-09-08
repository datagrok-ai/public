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
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterIsExactlyContains, filterPanelCount, filterPanelHas, filterPasses, filterPassesAll, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Subsequence search on the filter panel", () => {
  const session = feature(test, "features/search/subsequence.feature", import.meta.url);
  test("Subsequence search on the filter panel", {tag: ["@journey", "@realizes:bio.search.subsequence"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And user opens filter_FASTA dataset", () => openDataset(page, ds("filter_FASTA")));
    await session.step(10, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(11, "Then the table should have 14 rows", () => rowCount(page, 14));
    await session.step(12, "And \"fasta\" column should have units \"fasta\"", () => columnUnits(page, "fasta", "fasta"));
    await run.scenario("The command adds a filter bound to the sequence column", async () => {
      await session.step(15, "When user picks \"Bio > Search > Subsequence Search ...\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Subsequence Search ..."));
      await session.step(16, "Then filters viewer should be visible", () => shouldBe(page, el("filters viewer"), "visible"));
      await session.step(17, "And \"Substructure\" input in filters viewer should be visible", () => shouldBe(page, el("\"Substructure\" input in filters viewer"), "visible"));
      await session.step(18, "And the filter panel should have 1 filter", () => filterPanelCount(page, 1));
      await session.step(19, "And the filter panel should have a filter on \"fasta\" column", () => filterPanelHas(page, "fasta"));
    });
    await run.scenario("A subsequence one row contains keeps that row alone", async () => {
      await session.step(22, "When user enters \"RTDEVSNHTHDKPTLTWFEEIFEEYHSP\" into \"Substructure\" input in filters viewer", () => enterInto(page, "RTDEVSNHTHDKPTLTWFEEIFEEYHSP", el("\"Substructure\" input in filters viewer")));
      await session.step(23, "Then 1 row should pass the filter", () => filterPasses(page, 1));
      await session.step(24, "And the filter should pass exactly the rows where \"fasta\" contains \"RTDEVSNHTHDKPTLTWFEEIFEEYHSP\"", () => filterIsExactlyContains(page, "fasta", "RTDEVSNHTHDKPTLTWFEEIFEEYHSP"));
      await session.step(25, "And the table should have 14 rows", () => rowCount(page, 14));
    });
    await run.scenario("Reset restores every row and empties the query", async () => {
      await session.step(28, "When user clicks on \"arrow rotate left\" icon in filters viewer", () => clickOn(page, el("\"arrow rotate left\" icon in filters viewer")));
      await session.step(29, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(30, "And \"Substructure\" input in filters viewer should have value \"\"", () => shouldHaveValue(page, el("\"Substructure\" input in filters viewer"), ""));
      await session.step(31, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
