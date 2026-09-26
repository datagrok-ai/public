/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/empty-current-row.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.int.empty-input-on-row-viewers]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {makeRowCurrent, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount, setCell} from '@datagrok-libraries/bdd/bindings/platform/data';
import {dialogCloses, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The row analyses on a table with an empty sequence", () => {
  const session = feature(test, "features/analyze/empty-current-row.feature", import.meta.url);
  test("The row analyses on a table with an empty sequence", {tag: ["@journey", "@realizes:bio.int.empty-input-on-row-viewers"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await run.scenario("Similarity Search takes an empty current row as its target", async () => {
      await session.step(19, "Given user opens filter_FASTA dataset", () => openDataset(page, ds("filter_FASTA")));
      await session.step(20, "When user sets \"fasta\" column in row 1 to \"\"", () => setCell(page, "fasta", 1, ""));
      await session.step(21, "And user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(22, "Then the value of \"fasta\" column in row 1 should be \"\"", () => valueInRow(page, "fasta", 1, ""));
      await session.step(23, "When user picks \"Bio > Search > Similarity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Similarity Search"));
      await session.step(24, "Then \"Sequence Similarity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Similarity Search\" viewer"), "visible"));
      await session.step(25, "And the \"target row\" reading of \"Sequence Similarity Search\" viewer should be 0", () => readingIs(page, "target row", el("\"Sequence Similarity Search\" viewer"), 0));
      await session.step(26, "And the \"neighbours\" reading of \"Sequence Similarity Search\" viewer should be 11", () => readingIs(page, "neighbours", el("\"Sequence Similarity Search\" viewer"), 11));
      await session.step(27, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(28, "And the table should have 14 rows", () => rowCount(page, 14));
      await session.step(29, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Diversity Search counts the empty sequence among the column's values", async () => {
      await session.step(32, "When user picks \"Bio > Search > Diversity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Diversity Search"));
      await session.step(33, "Then \"Sequence Diversity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Diversity Search\" viewer"), "visible"));
      await session.step(34, "And the \"subset size\" reading of \"Sequence Diversity Search\" viewer should be 10", () => readingIs(page, "subset size", el("\"Sequence Diversity Search\" viewer"), 10));
      await session.step(35, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(36, "And the table should have 14 rows", () => rowCount(page, 14));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Activity Cliffs runs over a column with an empty sequence", async () => {
      await session.step(40, "Given user opens FASTA_sample dataset", () => openDataset(page, ds("FASTA_sample")));
      await session.step(41, "When user sets \"Sequence\" column in row 1 to \"\"", () => setCell(page, "Sequence", 1, ""));
      await session.step(42, "And user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(43, "And user picks \"Bio > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Activity Cliffs..."));
      await session.step(44, "And user selects \"Activity\" in Activities input in \"Sequence Activity Cliffs\" dialog", () => selectIn(page, "Activity", el("Activities input in \"Sequence Activity Cliffs\" dialog")));
      await session.step(45, "And user clicks on OK button in \"Sequence Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Sequence Activity Cliffs\" dialog")));
      await session.step(46, "Then the \"Sequence Activity Cliffs\" dialog should close", () => dialogCloses(page, "Sequence Activity Cliffs"));
      await session.step(47, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(48, "And the \"cliffs\" reading of scatter plot viewer should be 2", () => readingIs(page, "cliffs", el("scatter plot viewer"), 2));
      await session.step(49, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(50, "And the table should have 64 rows", () => rowCount(page, 64));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
