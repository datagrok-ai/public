/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/search/diversity.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.diversity-search]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {cardsFromFilteredRows} from '../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {enterInto, expand, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterBetween} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu, readingIs, readingNotAsRemembered, readingReads, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Diversity Search viewer and its properties", () => {
  const session = feature(test, "features/search/diversity.feature", import.meta.url);
  test("The Diversity Search viewer and its properties", {tag: ["@journey", "@realizes:chem.cp.diversity-search"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(12, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await run.scenario("The viewer shows a diverse set of 12 molecules", async () => {
      await session.step(15, "When user picks \"Chem > Search > Diversity Search...\" from the top menu", () => pickFromTopMenu(page, "Chem > Search > Diversity Search..."));
      await session.step(16, "Then Chem Diversity Search viewer should be visible", () => shouldBe(page, el("Chem Diversity Search viewer"), "visible"));
      await session.step(17, "And the \"cards\" reading of Chem Diversity Search viewer should be 12", () => readingIs(page, "cards", el("Chem Diversity Search viewer"), 12));
      await session.step(18, "And the \"header\" reading of Chem Diversity Search viewer should be \"Tanimoto, Morgan\"", () => readingReads(page, "header", el("Chem Diversity Search viewer"), "Tanimoto, Morgan"));
      await session.step(19, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Another metric picks another set, and the limit sets the number of cards", async () => {
      await session.step(22, "When user remembers the \"card row set\" reading of Chem Diversity Search viewer", () => rememberReading(page, "card row set", el("Chem Diversity Search viewer")));
      await session.step(23, "And user picks \"Properties...\" from the context menu of Chem Diversity Search viewer", () => pickFromContextMenu(page, "Properties...", el("Chem Diversity Search viewer")));
      await session.step(24, "And user expands Misc category", () => expand(page, el("Misc category")));
      await session.step(25, "And user selects \"Cosine\" in \"Distance Metric\" property", () => selectIn(page, "Cosine", el("\"Distance Metric\" property")));
      await session.step(26, "Then the \"header\" reading of Chem Diversity Search viewer should be \"Cosine, Morgan\"", () => readingReads(page, "header", el("Chem Diversity Search viewer"), "Cosine, Morgan"));
      await session.step(27, "And the \"card row set\" reading of Chem Diversity Search viewer should not be as remembered", () => readingNotAsRemembered(page, "card row set", el("Chem Diversity Search viewer")));
      await session.step(28, "When user enters \"6\" into Limit property", () => enterInto(page, "6", el("Limit property")));
      await session.step(29, "Then the \"cards\" reading of Chem Diversity Search viewer should be 6", () => readingIs(page, "cards", el("Chem Diversity Search viewer"), 6));
      await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Another fingerprint picks another set, and Size sets every card", async () => {
      await session.step(34, "When user remembers the \"card row set\" reading of Chem Diversity Search viewer", () => rememberReading(page, "card row set", el("Chem Diversity Search viewer")));
      await session.step(35, "And user selects \"MACCS\" in Fingerprint property", () => selectIn(page, "MACCS", el("Fingerprint property")));
      await session.step(36, "Then the \"header\" reading of Chem Diversity Search viewer should be \"Cosine, MACCS\"", () => readingReads(page, "header", el("Chem Diversity Search viewer"), "Cosine, MACCS"));
      await session.step(37, "And the \"card row set\" reading of Chem Diversity Search viewer should not be as remembered", () => readingNotAsRemembered(page, "card row set", el("Chem Diversity Search viewer")));
      await session.step(38, "When user selects \"large\" in Size property", () => selectIn(page, "large", el("Size property")));
      await session.step(39, "Then the \"card sizes\" reading of Chem Diversity Search viewer should be \"300x150\"", () => readingReads(page, "card sizes", el("Chem Diversity Search viewer"), "300x150"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Row Source Filtered keeps to the rows that pass the filter", async () => {
      await session.step(43, "When user filters rows where \"NumAromaticRings\" is between 0 and 1", () => filterBetween(page, "NumAromaticRings", 0, 1));
      await session.step(44, "And user selects \"Filtered\" in \"Row Source\" property", () => selectIn(page, "Filtered", el("\"Row Source\" property")));
      await session.step(45, "Then every card of Chem Diversity Search viewer should show a row that passes the filter", () => cardsFromFilteredRows(page, el("Chem Diversity Search viewer")));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
