/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/mpo/coverage-score.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
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
import {shippedMpoProfile} from '../../bindings/mpo.js';
import {below} from '@datagrok-libraries/bdd/bindings/common/geometry';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, selectIn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, columnType, displayedInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, dialogCloses, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("MPO Score counts the measured profile parameters of a row with Coverage Score", () => {
  const session = feature(test, "features/mpo/coverage-score.feature", import.meta.url);
  test("MPO Score counts the measured profile parameters of a row with Coverage Score", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And the shipped MPO profile \"ADME Test\" is on the server", () => shippedMpoProfile(page, "ADME Test"));
    await session.step(21, "And user opens a table \"adme coverage\" with:", () => openTableOf(page, "adme coverage", [["smiles","Caco2","Lipophilicity","Solubility"],["c1ccccc1","-6","2","0"],["CCO","","2","0"],["c1ccncc1","","3",""],["CC(=O)O","","",""],["CCN","<-5","2","0"]]), [["smiles","Caco2","Lipophilicity","Solubility"],["c1ccccc1","-6","2","0"],["CCO","","2","0"],["c1ccncc1","","3",""],["CC(=O)O","","",""],["CCN","<-5","2","0"]]);
    await run.scenario("Coverage Score is off by default and a run without it adds no coverage column", async () => {
      await session.step(30, "When user picks \"Chem > Calculate > MPO Score...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > MPO Score..."));
      await session.step(31, "Then \"MPO Score\" dialog should be visible", () => shouldBe(page, el("\"MPO Score\" dialog"), "visible"));
      await session.step(32, "And \"Coverage Score\" input in \"MPO Score\" dialog should be below \"Pareto front\" input in \"MPO Score\" dialog", () => below(page, el("\"Coverage Score\" input in \"MPO Score\" dialog"), el("\"Pareto front\" input in \"MPO Score\" dialog")));
      await session.step(33, "And \"Coverage Score\" input in \"MPO Score\" dialog should have value \"false\"", () => shouldHaveValue(page, el("\"Coverage Score\" input in \"MPO Score\" dialog"), "false"));
      await session.step(34, "When user selects \"✓ ADME Test\" in Profile input in \"MPO Score\" dialog", () => selectIn(page, "✓ ADME Test", el("Profile input in \"MPO Score\" dialog")));
      await session.step(35, "And user clicks on OK button in \"MPO Score\" dialog", () => clickOn(page, el("OK button in \"MPO Score\" dialog")));
      await session.step(36, "Then the \"MPO Score\" dialog should close", () => dialogCloses(page, "MPO Score"));
      await session.step(37, "And a new column \"MPO ADME Test\" should have been added", () => newColumnNamed(page, "MPO ADME Test"));
      await session.step(38, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Checked, the run adds the coverage column with the count of measured parameters per row", async () => {
      await session.step(42, "When user picks \"Chem > Calculate > MPO Score...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > MPO Score..."));
      await session.step(43, "And user selects \"✓ ADME Test\" in Profile input in \"MPO Score\" dialog", () => selectIn(page, "✓ ADME Test", el("Profile input in \"MPO Score\" dialog")));
      await session.step(44, "And user checks \"Coverage Score\" input in \"MPO Score\" dialog", () => check(page, el("\"Coverage Score\" input in \"MPO Score\" dialog")));
      await session.step(45, "And user clicks on OK button in \"MPO Score\" dialog", () => clickOn(page, el("OK button in \"MPO Score\" dialog")));
      await session.step(46, "Then the \"MPO Score\" dialog should close", () => dialogCloses(page, "MPO Score"));
      await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(48, "And a new column \"MPO ADME Test coverage\" should have been added", () => newColumnNamed(page, "MPO ADME Test coverage"));
      await session.step(49, "And \"MPO ADME Test coverage\" column should have type \"int\"", () => columnType(page, "MPO ADME Test coverage", "int"));
      await session.step(50, "And the \"MPO ADME Test coverage\" cell of row 1 should be displayed as \"3\"", () => displayedInRow(page, "MPO ADME Test coverage", 1, "3"));
      await session.step(51, "And the \"MPO ADME Test coverage\" cell of row 2 should be displayed as \"2\"", () => displayedInRow(page, "MPO ADME Test coverage", 2, "2"));
      await session.step(52, "And the \"MPO ADME Test coverage\" cell of row 3 should be displayed as \"1\"", () => displayedInRow(page, "MPO ADME Test coverage", 3, "1"));
      await session.step(53, "And the \"MPO ADME Test coverage\" cell of row 4 should be displayed as \"0\"", () => displayedInRow(page, "MPO ADME Test coverage", 4, "0"));
      await session.step(54, "And the \"MPO ADME Test coverage\" cell of row 5 should be displayed as \"3\"", () => displayedInRow(page, "MPO ADME Test coverage", 5, "3"));
      await session.step(55, "And the \"MPO ADME Test\" cell of row 2 should be displayed as \"\"", () => displayedInRow(page, "MPO ADME Test", 2, ""));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second run with Coverage Score replaces both columns instead of adding copies", async () => {
      await session.step(59, "When user picks \"Chem > Calculate > MPO Score...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > MPO Score..."));
      await session.step(60, "And user selects \"✓ ADME Test\" in Profile input in \"MPO Score\" dialog", () => selectIn(page, "✓ ADME Test", el("Profile input in \"MPO Score\" dialog")));
      await session.step(61, "And user checks \"Coverage Score\" input in \"MPO Score\" dialog", () => check(page, el("\"Coverage Score\" input in \"MPO Score\" dialog")));
      await session.step(62, "And user clicks on OK button in \"MPO Score\" dialog", () => clickOn(page, el("OK button in \"MPO Score\" dialog")));
      await session.step(63, "Then the \"MPO Score\" dialog should close", () => dialogCloses(page, "MPO Score"));
      await session.step(64, "And the table should have 6 columns", () => columnCount(page, 6));
      await session.step(65, "And the \"MPO ADME Test coverage\" cell of row 4 should be displayed as \"0\"", () => displayedInRow(page, "MPO ADME Test coverage", 4, "0"));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
