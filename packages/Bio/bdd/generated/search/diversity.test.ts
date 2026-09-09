/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/search/diversity.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.search.diversity]
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
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, noBalloons, noErrors, painted, readingAtLeast, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Diversity search", () => {
  const session = feature(test, "features/search/diversity.feature", import.meta.url);
  test("Diversity search", {tag: ["@journey", "@realizes:bio.search.diversity"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And user opens filter_FASTA dataset", () => openDataset(page, ds("filter_FASTA")));
    await session.step(10, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(11, "Then \"fasta\" column should have units \"fasta\"", () => columnUnits(page, "fasta", "fasta"));
    await run.scenario("The command docks the viewer with a varied subset", async () => {
      await session.step(14, "When user picks \"Bio > Search > Diversity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Diversity Search"));
      await session.step(15, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(16, "And \"Sequence Diversity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Diversity Search\" viewer"), "visible"));
      await session.step(17, "And \"Sequence Diversity Search\" viewer should be bound to table \"filter_FASTA\"", () => boundTable(page, el("\"Sequence Diversity Search\" viewer"), "filter_FASTA"));
      await session.step(18, "And the \"source column\" reading of \"Sequence Diversity Search\" viewer should be \"fasta\"", () => readingReads(page, "source column", el("\"Sequence Diversity Search\" viewer"), "fasta"));
      await session.step(19, "And the \"subset size\" reading of \"Sequence Diversity Search\" viewer should be 10", () => readingIs(page, "subset size", el("\"Sequence Diversity Search\" viewer"), 10));
      await session.step(20, "And the \"distinct sequences\" reading of \"Sequence Diversity Search\" viewer should be at least 2", () => readingAtLeast(page, "distinct sequences", el("\"Sequence Diversity Search\" viewer"), 2));
      await session.step(21, "And \"Sequence Diversity Search\" viewer should be painted", () => painted(page, el("\"Sequence Diversity Search\" viewer")));
      await session.step(22, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Closing the viewer removes it", async () => {
      await session.step(25, "When user clicks on close icon of \"Sequence Diversity Search\" viewer", () => clickOn(page, el("close icon of \"Sequence Diversity Search\" viewer")));
      await session.step(26, "Then \"Sequence Diversity Search\" viewer should be hidden", () => shouldBe(page, el("\"Sequence Diversity Search\" viewer"), "hidden"));
    });
    await run.scenario("On a HELM table the subset is HELM", async () => {
      await session.step(29, "Given user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
      await session.step(30, "Then \"HELM string\" column should have units \"helm\"", () => columnUnits(page, "HELM string", "helm"));
      await session.step(31, "When user picks \"Bio > Search > Diversity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Diversity Search"));
      await session.step(32, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(33, "And \"Sequence Diversity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Diversity Search\" viewer"), "visible"));
      await session.step(34, "And \"Sequence Diversity Search\" viewer should be bound to table \"filter_HELM\"", () => boundTable(page, el("\"Sequence Diversity Search\" viewer"), "filter_HELM"));
      await session.step(35, "And the \"source column\" reading of \"Sequence Diversity Search\" viewer should be \"HELM string\"", () => readingReads(page, "source column", el("\"Sequence Diversity Search\" viewer"), "HELM string"));
      await session.step(36, "And the \"subset size\" reading of \"Sequence Diversity Search\" viewer should be 4", () => readingIs(page, "subset size", el("\"Sequence Diversity Search\" viewer"), 4));
      await session.step(37, "And the \"distinct sequences\" reading of \"Sequence Diversity Search\" viewer should be 4", () => readingIs(page, "distinct sequences", el("\"Sequence Diversity Search\" viewer"), 4));
      await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
