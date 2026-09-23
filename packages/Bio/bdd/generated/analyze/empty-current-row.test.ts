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
import {commandCompleted, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount, setCell} from '@datagrok-libraries/bdd/bindings/platform/data';
import {dialogCloses, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorOrWarningBalloonMatching, noErrors, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The current-row analyses on an empty sequence", () => {
  const session = feature(test, "features/analyze/empty-current-row.feature", import.meta.url);
  test("The current-row analyses on an empty sequence", {tag: ["@journey", "@realizes:bio.int.empty-input-on-row-viewers", "@known-failure", "@GROK-16111"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the Bio package is initialized", () => bioInitialized(page));
    await run.scenario("Similarity Search on an empty current row keeps the table", async () => {
      await session.step(20, "Given user opens filter_FASTA dataset", () => openDataset(page, ds("filter_FASTA")));
      await session.step(21, "When user sets \"fasta\" column in row 1 to \"\"", () => setCell(page, "fasta", 1, ""));
      await session.step(22, "And user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(23, "Then the value of \"fasta\" column in row 1 should be \"\"", () => valueInRow(page, "fasta", 1, ""));
      await session.step(24, "When user picks \"Bio > Search > Similarity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Similarity Search"));
      await session.step(25, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(26, "And \"Sequence Similarity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Similarity Search\" viewer"), "visible"));
      await session.step(27, "And the \"target row\" reading of \"Sequence Similarity Search\" viewer should be 0", () => readingIs(page, "target row", el("\"Sequence Similarity Search\" viewer"), 0));
      await session.step(28, "And the table should have 14 rows", () => rowCount(page, 14));
      await session.step(29, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Similarity Search rejects an empty current row with a balloon", async () => {
      await session.step(34, "When user picks \"Bio > Search > Similarity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Similarity Search"));
      await session.step(35, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(36, "And an error or warning balloon matching \"empty|missing|no sequence\" should have been shown", () => errorOrWarningBalloonMatching(page, "empty|missing|no sequence"));
    }, {knownFailure: true});
    await run.scenario("Diversity Search on an empty current row keeps the table", async () => {
      await session.step(39, "When user picks \"Bio > Search > Diversity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Diversity Search"));
      await session.step(40, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(41, "And \"Sequence Diversity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Diversity Search\" viewer"), "visible"));
      await session.step(42, "And the table should have 14 rows", () => rowCount(page, 14));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Diversity Search rejects an empty current row with a balloon", async () => {
      await session.step(47, "When user picks \"Bio > Search > Diversity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Diversity Search"));
      await session.step(48, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(49, "And an error or warning balloon matching \"empty|missing|no sequence\" should have been shown", () => errorOrWarningBalloonMatching(page, "empty|missing|no sequence"));
    }, {knownFailure: true});
    await run.scenario("Activity Cliffs with an empty current row keeps the table", async () => {
      await session.step(52, "Given user opens FASTA_sample dataset", () => openDataset(page, ds("FASTA_sample")));
      await session.step(53, "When user sets \"Sequence\" column in row 1 to \"\"", () => setCell(page, "Sequence", 1, ""));
      await session.step(54, "And user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(55, "And user picks \"Bio > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Activity Cliffs..."));
      await session.step(56, "And user selects \"Activity\" in Activities input in \"Sequence Activity Cliffs\" dialog", () => selectIn(page, "Activity", el("Activities input in \"Sequence Activity Cliffs\" dialog")));
      await session.step(57, "And user clicks on OK button in \"Sequence Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Sequence Activity Cliffs\" dialog")));
      await session.step(58, "Then the \"Sequence Activity Cliffs\" dialog should close", () => dialogCloses(page, "Sequence Activity Cliffs"));
      await session.step(59, "And the top menu command should have completed", () => commandCompleted(page));
      await session.step(60, "And the table should have 64 rows", () => rowCount(page, 64));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Activity Cliffs rejects an empty current row with a balloon", async () => {
      await session.step(67, "When user picks \"Bio > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Activity Cliffs..."));
      await session.step(68, "And user selects \"Activity\" in Activities input in \"Sequence Activity Cliffs\" dialog", () => selectIn(page, "Activity", el("Activities input in \"Sequence Activity Cliffs\" dialog")));
      await session.step(69, "And user clicks on OK button in \"Sequence Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Sequence Activity Cliffs\" dialog")));
      await session.step(70, "Then an error or warning balloon matching \"empty|missing|no sequence\" should have been shown", () => errorOrWarningBalloonMatching(page, "empty|missing|no sequence"));
    }, {knownFailure: true});
    run.finish();
  });
});
