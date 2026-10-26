/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/search/similarity-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.similarity-search, public#4288]
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {cardsFromFilteredRows} from '../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {expand, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterBetween} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu, readingAsRemembered, readingIs, readingNotAsRemembered, readingReads, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Similarity Search regressions from closed issues", () => {
  const session = feature(test, "features/search/similarity-regressions.feature", import.meta.url);
  test("Similarity Search regressions from closed issues", {tag: ["@journey", "@realizes:chem.cp.similarity-search", "@realizes:public#4288"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(13, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await run.scenario("The Similarity Search viewer follows the filter, and Row Source All gives back every row (public#4288)", async () => {
      await session.step(17, "When user picks \"Chem > Search > Similarity Search...\" from the top menu", () => pickFromTopMenu(page, "Chem > Search > Similarity Search..."));
      await session.step(18, "Then Chem Similarity Search viewer should be visible", () => shouldBe(page, el("Chem Similarity Search viewer"), "visible"));
      await session.step(19, "And the \"cards\" reading of Chem Similarity Search viewer should be 12", () => readingIs(page, "cards", el("Chem Similarity Search viewer"), 12));
      await session.step(20, "And the \"row source\" reading of Chem Similarity Search viewer should be \"Filtered\"", () => readingReads(page, "row source", el("Chem Similarity Search viewer"), "Filtered"));
      await session.step(21, "When user remembers the \"card row set\" reading of Chem Similarity Search viewer", () => rememberReading(page, "card row set", el("Chem Similarity Search viewer")));
      await session.step(22, "And user filters rows where \"NumAromaticRings\" is between 3 and 7", () => filterBetween(page, "NumAromaticRings", 3, 7));
      await session.step(23, "Then the \"card row set\" reading of Chem Similarity Search viewer should not be as remembered", () => readingNotAsRemembered(page, "card row set", el("Chem Similarity Search viewer")));
      await session.step(24, "And the \"cards\" reading of Chem Similarity Search viewer should be 12", () => readingIs(page, "cards", el("Chem Similarity Search viewer"), 12));
      await session.step(25, "And every card of Chem Similarity Search viewer should show a row that passes the filter", () => cardsFromFilteredRows(page, el("Chem Similarity Search viewer")));
      await session.step(26, "When user picks \"Properties...\" from the context menu of Chem Similarity Search viewer", () => pickFromContextMenu(page, "Properties...", el("Chem Similarity Search viewer")));
      await session.step(27, "And user expands Misc category", () => expand(page, el("Misc category")));
      await session.step(28, "And user selects \"All\" in \"Row Source\" property", () => selectIn(page, "All", el("\"Row Source\" property")));
      await session.step(29, "Then the \"card row set\" reading of Chem Similarity Search viewer should be as remembered", () => readingAsRemembered(page, "card row set", el("Chem Similarity Search viewer")));
      await session.step(30, "When user selects \"Filtered\" in \"Row Source\" property", () => selectIn(page, "Filtered", el("\"Row Source\" property")));
      await session.step(31, "Then the \"card row set\" reading of Chem Similarity Search viewer should not be as remembered", () => readingNotAsRemembered(page, "card row set", el("Chem Similarity Search viewer")));
      await session.step(32, "And every card of Chem Similarity Search viewer should show a row that passes the filter", () => cardsFromFilteredRows(page, el("Chem Similarity Search viewer")));
      await session.step(33, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
