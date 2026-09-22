/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/r-groups.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.r-group-analysis]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {canvasColors} from '@datagrok-libraries/bdd/bindings/common/pixels';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, shouldBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnMatching, noNewColumn, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorBalloonText, noErrors, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("R-Groups Analysis with MCS, Replace latest and no core", () => {
  const session = feature(test, "features/analyze/r-groups.feature", import.meta.url);
  test("R-Groups Analysis with MCS, Replace latest and no core", {tag: ["@journey", "@realizes:chem.cp.r-group-analysis"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await run.scenario("MCS over molecules that are all the same finds no R-groups", async () => {
      await session.step(15, "Given user opens a table \"same_mols\" with:", () => openTableOf(page, "same_mols", [["id","smiles"],["1","c1ccccc1"],["2","c1ccccc1"],["3","c1ccccc1"],["4","c1ccccc1"],["5","c1ccccc1"]]));
      await session.step(22, "When user picks \"Chem > Analyze > R-Groups Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > R-Groups Analysis..."));
      await session.step(23, "Then \"R-Groups Analysis\" dialog should be visible", () => shouldBe(page, el("\"R-Groups Analysis\" dialog"), "visible"));
      await session.step(24, "And \"Visual analysis\" input in \"R-Groups Analysis\" dialog should be checked", () => shouldBe(page, el("\"Visual analysis\" input in \"R-Groups Analysis\" dialog"), "checked"));
      await session.step(25, "When user clicks on MCS button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("MCS button in \"R-Groups Analysis\" dialog")));
      await session.step(26, "And user clicks on OK button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("OK button in \"R-Groups Analysis\" dialog")));
      await session.step(27, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(28, "And no new column should have been added", () => noNewColumn(page));
      await session.step(29, "And the table should have 5 rows", () => rowCount(page, 5));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("MCS over one series adds R-group columns and a trellis plot", async () => {
      await session.step(33, "Given user opens sar-small dataset", () => openDataset(page, ds("sar-small")));
      await session.step(34, "When user picks \"Chem > Analyze > R-Groups Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > R-Groups Analysis..."));
      await session.step(35, "And user clicks on MCS button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("MCS button in \"R-Groups Analysis\" dialog")));
      await session.step(36, "And user clicks on OK button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("OK button in \"R-Groups Analysis\" dialog")));
      await session.step(37, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(38, "And a new column matching \"^R1\" should have been added", () => newColumnMatching(page, "^R1"));
      await session.step(39, "And trellis plot viewer should be visible", () => shouldBe(page, el("trellis plot viewer"), "visible"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Replace latest off the second run adds its own columns beside the first", async () => {
      await session.step(43, "When user picks \"Chem > Analyze > R-Groups Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > R-Groups Analysis..."));
      await session.step(44, "Then \"Replace latest\" input in \"R-Groups Analysis\" dialog should be visible", () => shouldBe(page, el("\"Replace latest\" input in \"R-Groups Analysis\" dialog"), "visible"));
      await session.step(45, "When user clicks on MCS button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("MCS button in \"R-Groups Analysis\" dialog")));
      await session.step(46, "And user unchecks \"Replace latest\" input in \"R-Groups Analysis\" dialog", () => uncheck(page, el("\"Replace latest\" input in \"R-Groups Analysis\" dialog")));
      await session.step(47, "And user clicks on OK button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("OK button in \"R-Groups Analysis\" dialog")));
      await session.step(48, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(49, "And a new column matching \"^R1\" should have been added", () => newColumnMatching(page, "^R1"));
      await session.step(50, "And the open tableview should have 2 trellis plot viewers", () => viewerCount(page, 2, "trellis plot"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Replace latest on the third run takes the latest set away", async () => {
      await session.step(54, "When user picks \"Chem > Analyze > R-Groups Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > R-Groups Analysis..."));
      await session.step(55, "And user clicks on MCS button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("MCS button in \"R-Groups Analysis\" dialog")));
      await session.step(56, "And user checks \"Replace latest\" input in \"R-Groups Analysis\" dialog", () => check(page, el("\"Replace latest\" input in \"R-Groups Analysis\" dialog")));
      await session.step(57, "And user clicks on OK button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("OK button in \"R-Groups Analysis\" dialog")));
      await session.step(58, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(59, "And the open tableview should have 2 trellis plot viewers", () => viewerCount(page, 2, "trellis plot"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A run without a core says so and keeps the results of the run before it", async () => {
      await session.step(64, "When user picks \"Chem > Analyze > R-Groups Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > R-Groups Analysis..."));
      await session.step(65, "And user clicks on OK button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("OK button in \"R-Groups Analysis\" dialog")));
      await session.step(66, "Then an error balloon containing \"No core was provided\" should have been shown", () => errorBalloonText(page, "No core was provided"));
      await session.step(67, "And the table should have a column \"R1\"", () => hasColumn(page, "R1"));
    });
    await run.scenario("MCS over a table of unrelated molecules finds no R-groups", async () => {
      await session.step(70, "Given user opens smiles dataset", () => openDataset(page, ds("smiles")));
      await session.step(71, "When user picks \"Chem > Analyze > R-Groups Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > R-Groups Analysis..."));
      await session.step(72, "And user clicks on MCS button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("MCS button in \"R-Groups Analysis\" dialog")));
      await session.step(73, "Then the canvases of \"R-Groups Analysis\" dialog should be painted in at least 2 colors", () => canvasColors(page, el("\"R-Groups Analysis\" dialog"), 2));
      await session.step(74, "When user clicks on OK button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("OK button in \"R-Groups Analysis\" dialog")));
      await session.step(75, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(76, "And an error balloon containing \"No R-Groups were found\" should have been shown", () => errorBalloonText(page, "No R-Groups were found"));
      await session.step(77, "And no new column should have been added", () => noNewColumn(page));
      await session.step(78, "And the table should have 1000 rows", () => rowCount(page, 1000));
    });
    run.finish();
  });
});
