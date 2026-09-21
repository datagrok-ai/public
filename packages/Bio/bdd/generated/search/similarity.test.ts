/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/search/similarity.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.search.similarity]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, makeLastRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, noBalloons, noErrors, painted, readingDiffers, readingIs, readingReads, readingSame, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Similarity search", () => {
  const session = feature(test, "features/search/similarity.feature", import.meta.url);
  test("Similarity search", {tag: ["@journey", "@realizes:bio.search.similarity"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And user opens filter_FASTA dataset", () => openDataset(page, ds("filter_FASTA")));
    await session.step(10, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(11, "Then \"fasta\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "fasta", "Macromolecule"));
    await run.scenario("The command docks the viewer with a full neighbour list", async () => {
      await session.step(14, "When user picks \"Bio > Search > Similarity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Similarity Search"));
      await session.step(15, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(16, "And \"Sequence Similarity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Similarity Search\" viewer"), "visible"));
      await session.step(17, "And \"Sequence Similarity Search\" viewer should be bound to table \"filter_FASTA\"", () => boundTable(page, el("\"Sequence Similarity Search\" viewer"), "filter_FASTA"));
      await session.step(18, "And the \"source column\" reading of \"Sequence Similarity Search\" viewer should be \"fasta\"", () => readingReads(page, "source column", el("\"Sequence Similarity Search\" viewer"), "fasta"));
      await session.step(19, "And the \"limit\" reading of \"Sequence Similarity Search\" viewer should be 10", () => readingIs(page, "limit", el("\"Sequence Similarity Search\" viewer"), 10));
      await session.step(20, "And the \"target row\" reading of \"Sequence Similarity Search\" viewer should be 0", () => readingIs(page, "target row", el("\"Sequence Similarity Search\" viewer"), 0));
      await session.step(21, "And the \"neighbours\" reading of \"Sequence Similarity Search\" viewer should be 11", () => readingIs(page, "neighbours", el("\"Sequence Similarity Search\" viewer"), 11));
      await session.step(22, "And \"Sequence Similarity Search\" viewer should be painted", () => painted(page, el("\"Sequence Similarity Search\" viewer")));
      await session.step(23, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Another current row is another query", async () => {
      await session.step(26, "When user takes a snapshot of \"Sequence Similarity Search\" viewer", () => takeSnapshot(page, el("\"Sequence Similarity Search\" viewer")));
      await session.step(27, "And user makes the last row current", () => makeLastRowCurrent(page));
      await session.step(28, "Then the \"target row\" reading of \"Sequence Similarity Search\" viewer should be 13", () => readingIs(page, "target row", el("\"Sequence Similarity Search\" viewer"), 13));
      await session.step(29, "And the \"neighbour set\" reading of \"Sequence Similarity Search\" viewer should differ from before", () => readingDiffers(page, "neighbour set", el("\"Sequence Similarity Search\" viewer")));
      await session.step(30, "And the \"neighbours\" reading of \"Sequence Similarity Search\" viewer should be the same as before", () => readingSame(page, "neighbours", el("\"Sequence Similarity Search\" viewer")));
      await session.step(31, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
