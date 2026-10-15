/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/r-groups.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.r-group-analysis]
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
import {check, clickOn, finishedUpdating, shouldBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn, hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnMatching, newColumnNamed, noNewColumn, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {autostartsCompleted, openDataset, openTableOf, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorBalloonText, noErrors, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("R-Groups Analysis with MCS, Replace latest and no core", () => {
  const session = feature(test, "features/analyze/r-groups.feature", import.meta.url);
  test("R-Groups Analysis with MCS, Replace latest and no core", {tag: ["@journey", "@realizes:chem.cp.r-group-analysis"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(22, "And the package autostarts have completed", () => autostartsCompleted(page));
    await run.scenario("MCS over molecules that are all the same finds no R-groups", async () => {
      await session.step(25, "Given user opens a table \"same_mols\" with:", () => openTableOf(page, "same_mols", [["id","smiles"],["1","c1ccccc1"],["2","c1ccccc1"],["3","c1ccccc1"],["4","c1ccccc1"],["5","c1ccccc1"]]), [["id","smiles"],["1","c1ccccc1"],["2","c1ccccc1"],["3","c1ccccc1"],["4","c1ccccc1"],["5","c1ccccc1"]]);
      await session.step(32, "When user picks \"Chem > Analyze > R-Groups Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > R-Groups Analysis..."));
      await session.step(33, "Then \"R-Groups Analysis\" dialog should be visible", () => shouldBe(page, el("\"R-Groups Analysis\" dialog"), "visible"));
      await session.step(34, "And \"Visual analysis\" input in \"R-Groups Analysis\" dialog should be checked", () => shouldBe(page, el("\"Visual analysis\" input in \"R-Groups Analysis\" dialog"), "checked"));
      await session.step(36, "When user clicks on R-Groups settings icon", () => clickOn(page, el("R-Groups settings icon")));
      await session.step(37, "And user unchecks \"Only match at R groups\" input in \"R-Groups Analysis\" dialog", () => uncheck(page, el("\"Only match at R groups\" input in \"R-Groups Analysis\" dialog")));
      await session.step(38, "When user clicks on MCS button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("MCS button in \"R-Groups Analysis\" dialog")));
      await session.step(39, "And \"R-Groups Analysis\" dialog should have finished updating", () => finishedUpdating(page, el("\"R-Groups Analysis\" dialog")));
      await session.step(40, "And user clicks on OK button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("OK button in \"R-Groups Analysis\" dialog")));
      await session.step(41, "Then an error balloon containing \"No R-Groups were found\" should have been shown", () => errorBalloonText(page, "No R-Groups were found"));
      await session.step(42, "And the table should not have a column \"R1\"", () => hasNoColumn(page, "R1"));
      await session.step(43, "And the table should have 5 rows", () => rowCount(page, 5));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("MCS over one series adds R-group columns and a trellis plot", async () => {
      await session.step(47, "Given user opens sar-small dataset", () => openDataset(page, ds("sar-small")));
      await session.step(48, "When user picks \"Chem > Analyze > R-Groups Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > R-Groups Analysis..."));
      await session.step(49, "And user clicks on MCS button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("MCS button in \"R-Groups Analysis\" dialog")));
      await session.step(50, "And \"R-Groups Analysis\" dialog should have finished updating", () => finishedUpdating(page, el("\"R-Groups Analysis\" dialog")));
      await session.step(51, "And user clicks on OK button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("OK button in \"R-Groups Analysis\" dialog")));
      await session.step(52, "Then a new column matching \"^R1\" should have been added", () => newColumnMatching(page, "^R1"));
      await session.step(53, "And trellis plot viewer should be visible", () => shouldBe(page, el("trellis plot viewer"), "visible"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Replace latest off the second run adds its own columns beside the first", async () => {
      await session.step(57, "When user picks \"Chem > Analyze > R-Groups Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > R-Groups Analysis..."));
      await session.step(58, "Then \"Replace latest\" input in \"R-Groups Analysis\" dialog should be visible", () => shouldBe(page, el("\"Replace latest\" input in \"R-Groups Analysis\" dialog"), "visible"));
      await session.step(59, "When user clicks on MCS button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("MCS button in \"R-Groups Analysis\" dialog")));
      await session.step(60, "And \"R-Groups Analysis\" dialog should have finished updating", () => finishedUpdating(page, el("\"R-Groups Analysis\" dialog")));
      await session.step(61, "And user unchecks \"Replace latest\" input in \"R-Groups Analysis\" dialog", () => uncheck(page, el("\"Replace latest\" input in \"R-Groups Analysis\" dialog")));
      await session.step(62, "And user clicks on OK button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("OK button in \"R-Groups Analysis\" dialog")));
      await session.step(63, "Then a new column \"R1_1\" should have been added", () => newColumnNamed(page, "R1_1"));
      await session.step(64, "And the open tableview should have 2 trellis plot viewers", () => viewerCount(page, 2, "trellis plot"));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Replace latest on the third run takes the latest set away", async () => {
      await session.step(68, "When user picks \"Chem > Analyze > R-Groups Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > R-Groups Analysis..."));
      await session.step(69, "And user clicks on MCS button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("MCS button in \"R-Groups Analysis\" dialog")));
      await session.step(70, "And \"R-Groups Analysis\" dialog should have finished updating", () => finishedUpdating(page, el("\"R-Groups Analysis\" dialog")));
      await session.step(71, "And user checks \"Replace latest\" input in \"R-Groups Analysis\" dialog", () => check(page, el("\"Replace latest\" input in \"R-Groups Analysis\" dialog")));
      await session.step(72, "And user watches the task bar", () => watchTaskBar(page));
      await session.step(73, "And user clicks on OK button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("OK button in \"R-Groups Analysis\" dialog")));
      await session.step(74, "Then the task bar should have finished \"R-Group analysis running\"", () => taskBarFinished(page, "R-Group analysis running"));
      await session.step(75, "And the open tableview should have 2 trellis plot viewers", () => viewerCount(page, 2, "trellis plot"));
      await session.step(76, "And the table should have a column \"R1_1\"", () => hasColumn(page, "R1_1"));
      await session.step(77, "And no new column should have been added", () => noNewColumn(page));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A run without a core says so and keeps the latest results", async () => {
      await session.step(81, "When user picks \"Chem > Analyze > R-Groups Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > R-Groups Analysis..."));
      await session.step(82, "And user clicks on OK button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("OK button in \"R-Groups Analysis\" dialog")));
      await session.step(83, "Then an error balloon containing \"No core was provided\" should have been shown", () => errorBalloonText(page, "No core was provided"));
      await session.step(84, "And the table should have a column \"R1_1\"", () => hasColumn(page, "R1_1"));
      await session.step(85, "And the open tableview should have 2 trellis plot viewers", () => viewerCount(page, 2, "trellis plot"));
    });
    run.finish();
  });
});
