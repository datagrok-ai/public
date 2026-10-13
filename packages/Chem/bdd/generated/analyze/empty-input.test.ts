/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/empty-input.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.int.empty-input-analyses]
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
import {clickOn, finishedUpdating, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs, setColumnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {noNewColumn, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {currentColumnIs, filterPassesAll, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openTableOf, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, errorBalloonText, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("R-Groups Analysis on an all-empty molecule column", () => {
  const session = feature(test, "features/analyze/empty-input.feature", import.meta.url);
  test("R-Groups Analysis on an all-empty molecule column", {tag: ["@journey", "@realizes:chem.int.empty-input-analyses"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(13, "And user opens a table \"empty_mols\" with:", () => openTableOf(page, "empty_mols", [["id","structure"],["1",""],["2",""],["3",""],["4",""],["5",""],["6",""],["7",""],["8",""],["9",""],["10",""]]), [["id","structure"],["1",""],["2",""],["3",""],["4",""],["5",""],["6",""],["7",""],["8",""],["9",""],["10",""]]);
    await session.step(25, "And user sets the semantic type of \"structure\" column to \"Molecule\"", () => setColumnSemType(page, "structure", "Molecule"));
    await run.scenario("R-Groups Analysis with MCS says it has no core and adds nothing", async () => {
      await session.step(28, "When user picks \"Chem > Analyze > R-Groups Analysis...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > R-Groups Analysis..."));
      await session.step(29, "Then \"R-Groups Analysis\" dialog should be visible", () => shouldBe(page, el("\"R-Groups Analysis\" dialog"), "visible"));
      await session.step(30, "When user clicks on MCS button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("MCS button in \"R-Groups Analysis\" dialog")));
      await session.step(31, "And \"R-Groups Analysis\" dialog should have finished updating", () => finishedUpdating(page, el("\"R-Groups Analysis\" dialog")));
      await session.step(32, "And user clicks on OK button in \"R-Groups Analysis\" dialog", () => clickOn(page, el("OK button in \"R-Groups Analysis\" dialog")));
      await session.step(33, "Then an error balloon containing \"No core was provided\" should have been shown", () => errorBalloonText(page, "No core was provided"));
      await session.step(34, "And no new column should have been added", () => noNewColumn(page));
      await session.step(35, "And the table should have 10 rows", () => rowCount(page, 10));
      await session.step(36, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(37, "When user clicks on the \"cell 4 of structure\" area of grid", () => clickArea(page, "cell 4 of structure", el("grid")));
      await session.step(38, "Then row 4 should be current", () => currentRowIs(page, 4));
      await session.step(39, "And the current column should be \"structure\"", () => currentColumnIs(page, "structure"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
