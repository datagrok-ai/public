/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/empty-input.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.int.empty-input-analyses]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, finishedUpdating, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs, setColumnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnMatching, newestMatchingFilled, noNewColumn, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {currentColumnIs, filterPassesAll, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, dialogCloses, openTableOf, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, errorBalloonText, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("R-Groups Analysis and Chemical Space on an all-empty molecule column", () => {
  const session = feature(test, "features/analyze/empty-input.feature", import.meta.url);
  test("R-Groups Analysis and Chemical Space on an all-empty molecule column", {tag: ["@journey", "@realizes:chem.int.empty-input-analyses"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(11, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(12, "And user opens a table \"empty_mols\" with:", () => openTableOf(page, "empty_mols", [["id","structure"],["1",""],["2",""],["3",""],["4",""],["5",""],["6",""],["7",""],["8",""],["9",""],["10",""]]), [["id","structure"],["1",""],["2",""],["3",""],["4",""],["5",""],["6",""],["7",""],["8",""],["9",""],["10",""]]);
    await session.step(24, "And user sets the semantic type of \"structure\" column to \"Molecule\"", () => setColumnSemType(page, "structure", "Molecule"));
    await run.scenario("R-Groups Analysis with MCS says it has no core and adds nothing", async () => {
      await session.step(27, "When user picks \"Chem > Analyze > R-Groups Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > R-Groups Analysis..."));
      await session.step(28, "Then \"R-Groups Analysis\" dialog should be visible", () => shouldBe(page, el("\"R-Groups Analysis\" dialog"), "visible"));
      await session.step(29, "When user clicks on MCS button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("MCS button in \"R-Groups Analysis\" dialog")));
      await session.step(30, "And \"R-Groups Analysis\" dialog should have finished updating", () => finishedUpdating(page, el("\"R-Groups Analysis\" dialog")));
      await session.step(31, "And user clicks on OK button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("OK button in \"R-Groups Analysis\" dialog")));
      await session.step(32, "Then an error balloon containing \"No core was provided\" should have been shown", () => errorBalloonText(page, "No core was provided"));
      await session.step(33, "And no new column should have been added", () => noNewColumn(page));
      await session.step(34, "And the table should have 10 rows", () => rowCount(page, 10));
      await session.step(35, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(36, "When user clicks on the \"cell 4 of structure\" area of grid", () => clickArea(page, "cell 4 of structure", el("grid")));
      await session.step(37, "Then row 4 should be current", () => currentRowIs(page, 4));
      await session.step(38, "And the current column should be \"structure\"", () => currentColumnIs(page, "structure"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Chemical Space embeds the empty column and plots it", async () => {
      await session.step(42, "When user picks \"Chem > Analyze > Chemical Space...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Chemical Space..."));
      await session.step(43, "Then \"Chem Space\" dialog should be visible", () => shouldBe(page, el("\"Chem Space\" dialog"), "visible"));
      await session.step(44, "And Column input in \"Chem Space\" dialog should contain text \"structure\"", () => shouldContainText(page, el("Column input in \"Chem Space\" dialog"), "structure"));
      await session.step(45, "When user clicks on OK button in \"Chem Space\" dialog", () => clickOn(page, el("OK button in \"Chem Space\" dialog")));
      await session.step(46, "Then the \"Chem Space\" dialog should close", () => dialogCloses(page, "Chem Space"));
      await session.step(47, "And the top menu command should have completed", () => commandCompleted(page));
      await session.step(48, "And a new column matching \"^Embed_X_\" should have been added", () => newColumnMatching(page, "^Embed_X_"));
      await session.step(49, "And a new column matching \"^Embed_Y_\" should have been added", () => newColumnMatching(page, "^Embed_Y_"));
      await session.step(50, "And the newest column matching \"^Embed_X_\" should have no missing values", () => newestMatchingFilled(page, "^Embed_X_"));
      await session.step(51, "And a new column matching \"^Cluster \" should have been added", () => newColumnMatching(page, "^Cluster "));
      await session.step(52, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(53, "And the table should have 10 rows", () => rowCount(page, 10));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
