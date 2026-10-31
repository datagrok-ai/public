/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/activity-cliffs.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.activity-cliffs, public#1512]
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
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe, shouldContainText, shouldHaveValue, shouldNotBe, switchOff, switchOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {commandCompleted, newColumnMatching, newestMatchingFilled, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPasses, filterPassesAll, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, switchTableView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Activity Cliffs over molecules and their activity", () => {
  const session = feature(test, "features/analyze/activity-cliffs.feature", import.meta.url);
  test("Activity Cliffs over molecules and their activity", {tag: ["@journey", "@realizes:chem.cp.activity-cliffs", "@realizes:public#1512"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens activity-cliffs dataset", () => openDataset(page, ds("activity-cliffs")));
    await run.scenario("The dialog opens on the molecule column with its defaults", async () => {
      await session.step(22, "When user picks \"Chem > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Activity Cliffs..."));
      await session.step(23, "Then \"Activity Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Activity Cliffs\" dialog"), "visible"));
      await session.step(24, "And Column input in \"Activity Cliffs\" dialog should contain text \"smiles\"", () => shouldContainText(page, el("Column input in \"Activity Cliffs\" dialog"), "smiles"));
      await session.step(25, "And Activities input in \"Activity Cliffs\" dialog should contain text \"Activity\"", () => shouldContainText(page, el("Activities input in \"Activity Cliffs\" dialog"), "Activity"));
      await session.step(26, "And Method input in \"Activity Cliffs\" dialog should have value \"UMAP\"", () => shouldHaveValue(page, el("Method input in \"Activity Cliffs\" dialog"), "UMAP"));
      await session.step(27, "And Similarity input in \"Activity Cliffs\" dialog should have value \"Tanimoto\"", () => shouldHaveValue(page, el("Similarity input in \"Activity Cliffs\" dialog"), "Tanimoto"));
      await session.step(28, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The run embeds the molecules, plots them and counts the cliffs", async () => {
      await session.step(31, "When user clicks on OK button in \"Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Activity Cliffs\" dialog")));
      await session.step(32, "Then a new column matching \"^Embed_X_\" should have been added", () => newColumnMatching(page, "^Embed_X_"));
      await session.step(33, "And a new column matching \"^Embed_Y_\" should have been added", () => newColumnMatching(page, "^Embed_Y_"));
      await session.step(34, "And the newest column matching \"^Embed_X_\" should have no missing values", () => newestMatchingFilled(page, "^Embed_X_"));
      await session.step(35, "And a new column matching \"^sali_\" should have been added", () => newColumnMatching(page, "^sali_"));
      await session.step(36, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(37, "And the \"cliffs\" reading of scatter plot viewer should be 2", () => readingIs(page, "cliffs", el("scatter plot viewer"), 2));
      await session.step(38, "And the \"only cliffs\" reading of scatter plot viewer should be \"false\"", () => readingReads(page, "only cliffs", el("scatter plot viewer"), "false"));
      await session.step(39, "And the table should have 29 rows", () => rowCount(page, 29));
      await session.step(40, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The cliff link opens the panel that lists the cliffs", async () => {
      await session.step(44, "When user clicks on \"cliffs\" button", () => clickOn(page, el("\"cliffs\" button")));
      await session.step(45, "Then \"Activity cliffs\" dock panel should be visible", () => shouldBe(page, el("\"Activity cliffs\" dock panel"), "visible"));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Activity cliffs panel stays with the view of its analysis (public#1512)", async () => {
      await session.step(50, "Given user opens cars dataset", () => openDataset(page, ds("cars")));
      await session.step(51, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
      await session.step(52, "And \"Activity cliffs\" dock panel should not be visible", () => shouldNotBe(page, el("\"Activity cliffs\" dock panel"), "visible"));
      await session.step(53, "Given user switches to the \"activity_cliffs_test\" table view", () => switchTableView(page, "activity_cliffs_test"));
      await session.step(54, "Then the \"activity_cliffs_test\" view should be current", () => viewIsCurrent(page, "activity_cliffs_test"));
      await session.step(55, "And \"Activity cliffs\" dock panel should be visible", () => shouldBe(page, el("\"Activity cliffs\" dock panel"), "visible"));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show only cliffs keeps the rows that take part in one", async () => {
      await session.step(59, "When user switches on \"Show only cliffs\" input", () => switchOn(page, el("\"Show only cliffs\" input")));
      await session.step(60, "Then the \"only cliffs\" reading of scatter plot viewer should be \"true\"", () => readingReads(page, "only cliffs", el("scatter plot viewer"), "true"));
      await session.step(61, "And 4 rows should pass the filter", () => filterPasses(page, 4));
      await session.step(62, "When user switches off \"Show only cliffs\" input", () => switchOff(page, el("\"Show only cliffs\" input")));
      await session.step(63, "Then the \"only cliffs\" reading of scatter plot viewer should be \"false\"", () => readingReads(page, "only cliffs", el("scatter plot viewer"), "false"));
      await session.step(64, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A stricter similarity cutoff finds fewer", async () => {
      await session.step(68, "Given user opens activity-cliffs dataset", () => openDataset(page, ds("activity-cliffs")));
      await session.step(69, "When user picks \"Chem > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Activity Cliffs..."));
      await session.step(70, "And user enters \"95\" into \"Similarity cutoff\" input in \"Activity Cliffs\" dialog", () => enterInto(page, "95", el("\"Similarity cutoff\" input in \"Activity Cliffs\" dialog")));
      await session.step(71, "And user clicks on OK button in \"Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Activity Cliffs\" dialog")));
      await session.step(72, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(73, "And a new column matching \"^sali_\" should have been added", () => newColumnMatching(page, "^sali_"));
      await session.step(74, "And the \"cliffs\" reading of scatter plot viewer should be 1", () => readingIs(page, "cliffs", el("scatter plot viewer"), 1));
      await session.step(75, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
