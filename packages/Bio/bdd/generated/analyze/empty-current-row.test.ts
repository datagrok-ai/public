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
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("The current-row analyses on an empty sequence", () => {
  const session = feature(test, "features/analyze/empty-current-row.feature", import.meta.url);
  test("Similarity Search on an empty current row keeps the table", {tag: ["@realizes:bio.int.empty-input-on-row-viewers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(19, "Given user opens filter_FASTA dataset", () => openDataset(page, ds("filter_FASTA")));
    await session.step(20, "When user sets \"fasta\" column in row 1 to \"\"", () => setCell(page, "fasta", 1, ""));
    await session.step(21, "And user makes row 1 current", () => makeRowCurrent(page, 1));
    await session.step(22, "Then the value of \"fasta\" column in row 1 should be \"\"", () => valueInRow(page, "fasta", 1, ""));
    await session.step(23, "When user picks \"Bio > Search > Similarity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Similarity Search"));
    await session.step(24, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(25, "And \"Sequence Similarity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Similarity Search\" viewer"), "visible"));
    await session.step(26, "And the \"target row\" reading of \"Sequence Similarity Search\" viewer should be 0", () => readingIs(page, "target row", el("\"Sequence Similarity Search\" viewer"), 0));
    await session.step(27, "And the table should have 14 rows", () => rowCount(page, 14));
    await session.step(28, "And no errors should have been logged", () => noErrors(page));
  });
  test("Similarity Search rejects an empty current row with a balloon", {tag: ["@realizes:bio.int.empty-input-on-row-viewers", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await knownFailure(async () => {
      await session.step(33, "Given user opens filter_FASTA dataset", () => openDataset(page, ds("filter_FASTA")));
      await session.step(34, "When user sets \"fasta\" column in row 1 to \"\"", () => setCell(page, "fasta", 1, ""));
      await session.step(35, "And user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(36, "And user picks \"Bio > Search > Similarity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Similarity Search"));
      await session.step(37, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(38, "And an error or warning balloon matching \"empty|missing|null|no sequence\" should have been shown", () => errorOrWarningBalloonMatching(page, "empty|missing|null|no sequence"));
    });
  });
  test("Diversity Search on an empty current row keeps the table", {tag: ["@realizes:bio.int.empty-input-on-row-viewers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(41, "Given user opens filter_FASTA dataset", () => openDataset(page, ds("filter_FASTA")));
    await session.step(42, "When user sets \"fasta\" column in row 1 to \"\"", () => setCell(page, "fasta", 1, ""));
    await session.step(43, "And user makes row 1 current", () => makeRowCurrent(page, 1));
    await session.step(44, "And user picks \"Bio > Search > Diversity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Diversity Search"));
    await session.step(45, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(46, "And \"Sequence Diversity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Diversity Search\" viewer"), "visible"));
    await session.step(47, "And the table should have 14 rows", () => rowCount(page, 14));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
  });
  test("Diversity Search rejects an empty current row with a balloon", {tag: ["@realizes:bio.int.empty-input-on-row-viewers", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await knownFailure(async () => {
      await session.step(53, "Given user opens filter_FASTA dataset", () => openDataset(page, ds("filter_FASTA")));
      await session.step(54, "When user sets \"fasta\" column in row 1 to \"\"", () => setCell(page, "fasta", 1, ""));
      await session.step(55, "And user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(56, "And user picks \"Bio > Search > Diversity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Diversity Search"));
      await session.step(57, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(58, "And an error or warning balloon matching \"empty|missing|null|no sequence\" should have been shown", () => errorOrWarningBalloonMatching(page, "empty|missing|null|no sequence"));
    });
  });
  test("Activity Cliffs with an empty current row keeps the table", {tag: ["@realizes:bio.int.empty-input-on-row-viewers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(61, "Given user opens FASTA_sample dataset", () => openDataset(page, ds("FASTA_sample")));
    await session.step(62, "When user sets \"Sequence\" column in row 1 to \"\"", () => setCell(page, "Sequence", 1, ""));
    await session.step(63, "And user makes row 1 current", () => makeRowCurrent(page, 1));
    await session.step(64, "And user picks \"Bio > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Activity Cliffs..."));
    await session.step(65, "And user selects \"Activity\" in Activities input in \"Sequence Activity Cliffs\" dialog", () => selectIn(page, "Activity", el("Activities input in \"Sequence Activity Cliffs\" dialog")));
    await session.step(66, "And user clicks on OK button in \"Sequence Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Sequence Activity Cliffs\" dialog")));
    await session.step(67, "Then the \"Sequence Activity Cliffs\" dialog should close", () => dialogCloses(page, "Sequence Activity Cliffs"));
    await session.step(68, "And the top menu command should have completed", () => commandCompleted(page));
    await session.step(69, "And the table should have 64 rows", () => rowCount(page, 64));
    await session.step(70, "And no errors should have been logged", () => noErrors(page));
  });
  test("Activity Cliffs rejects an empty current row with a balloon", {tag: ["@realizes:bio.int.empty-input-on-row-viewers", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await knownFailure(async () => {
      await session.step(77, "Given user opens FASTA_sample dataset", () => openDataset(page, ds("FASTA_sample")));
      await session.step(78, "When user sets \"Sequence\" column in row 1 to \"\"", () => setCell(page, "Sequence", 1, ""));
      await session.step(79, "And user makes row 1 current", () => makeRowCurrent(page, 1));
      await session.step(80, "And user picks \"Bio > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Activity Cliffs..."));
      await session.step(81, "And user selects \"Activity\" in Activities input in \"Sequence Activity Cliffs\" dialog", () => selectIn(page, "Activity", el("Activities input in \"Sequence Activity Cliffs\" dialog")));
      await session.step(82, "And user clicks on OK button in \"Sequence Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Sequence Activity Cliffs\" dialog")));
      await session.step(83, "Then an error or warning balloon matching \"empty|missing|null|no sequence\" should have been shown", () => errorOrWarningBalloonMatching(page, "empty|missing|null|no sequence"));
    });
  });
});
