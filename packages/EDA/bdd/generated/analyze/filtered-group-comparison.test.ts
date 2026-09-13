/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/filtered-group-comparison.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [ml.menu.analyze.group-comparison.control-comparisons]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {commandCompleted, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPasses, filterTo, tableColumnComplete, tableOpen, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, noBalloons, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Group comparison of a filtered table", () => {
  const session = feature(test, "features/analyze/filtered-group-comparison.feature", import.meta.url);
  test("Group comparison of a filtered table", {tag: ["@journey", "@eda", "@realizes:ml.menu.analyze.group-comparison.control-comparisons", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await run.scenario("Running control comparisons preserves the filter and produces the comparison table", async () => {
      await session.step(16, "When user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
      await session.step(17, "Then 3243 rows should pass the filter", () => filterPasses(page, 3243));
      await session.step(18, "When user picks \"ML > Analyze > Group Comparison > Control Comparisons...\" from the top menu", () => pickFromTopMenu(page, "ML > Analyze > Group Comparison > Control Comparisons..."));
      await session.step(19, "Then \"Control comparisons\" dialog should be visible", () => shouldBe(page, el("\"Control comparisons\" dialog"), "visible"));
      await session.step(20, "And 3243 rows should pass the filter", () => filterPasses(page, 3243));
      await session.step(21, "When user clicks on Run button in \"Control comparisons\" dialog", () => clickOn(page, el("Run button in \"Control comparisons\" dialog")));
      await session.step(22, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(23, "And \"Control comparisons\" dialog should be hidden", () => shouldBe(page, el("\"Control comparisons\" dialog"), "hidden"));
      await session.step(24, "And 3243 rows should pass the filter", () => filterPasses(page, 3243));
      await session.step(25, "And table \"Control comparisons result\" should be open", () => tableOpen(page, "Control comparisons result"));
      await session.step(26, "And table \"Control comparisons result\" should have 3 rows", () => tableRows(page, "Control comparisons result", 3));
      await session.step(27, "And table \"Control comparisons result\" should have no missing values in \"n\" column", () => tableColumnComplete(page, "Control comparisons result", "n"));
      await session.step(28, "And second grid viewer should be bound to table \"Control comparisons result\"", () => boundTable(page, el("second grid viewer"), "Control comparisons result"));
      await session.step(29, "And the \"text of cell 1 of Group\" reading of second grid viewer should be \"Black\"", () => readingReads(page, "text of cell 1 of Group", el("second grid viewer"), "Black"));
      await session.step(30, "And the \"text of cell 2 of Group\" reading of second grid viewer should be \"Caucasian\"", () => readingReads(page, "text of cell 2 of Group", el("second grid viewer"), "Caucasian"));
      await session.step(31, "And the \"text of cell 3 of Group\" reading of second grid viewer should be \"Other\"", () => readingReads(page, "text of cell 3 of Group", el("second grid viewer"), "Other"));
      await session.step(32, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(33, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The comparison sizes count only the women", async () => {
      await session.step(37, "Then the \"text of cell 1 of n\" reading of second grid viewer should be \"104\"", () => readingReads(page, "text of cell 1 of n", el("second grid viewer"), "104"));
      await session.step(38, "And the \"text of cell 2 of n\" reading of second grid viewer should be \"2823\"", () => readingReads(page, "text of cell 2 of n", el("second grid viewer"), "2823"));
      await session.step(39, "And the \"text of cell 3 of n\" reading of second grid viewer should be \"279\"", () => readingReads(page, "text of cell 3 of n", el("second grid viewer"), "279"));
    }, {knownFailure: true});
    run.finish();
  });
});
