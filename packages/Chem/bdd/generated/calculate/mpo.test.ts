/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/calculate/mpo.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.mpo-profile-crud, public#4403]
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
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, selectIn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnType, everyValueBetween, everyValueMatches} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, dialogCloses, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("MPO Score over a profile's own properties", () => {
  const session = feature(test, "features/calculate/mpo.feature", import.meta.url);
  test("MPO Score over a profile's own properties", {tag: ["@journey", "@realizes:chem.cp.mpo-profile-crud", "@realizes:public#4403"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(24, "And the shipped MPO profile \"ADME Test\" is on the server", () => shippedMpoProfile(page, "ADME Test"));
    await run.scenario("The dialog offers the ADME Test profile, which maps the table's three columns", async () => {
      await session.step(27, "Given user opens a table \"adme\" with:", () => openTableOf(page, "adme", [["smiles","Caco2","Lipophilicity","Solubility"],["c1ccccc1","-6","2","-4"],["CCO","-7","3","-5"],["c1ccncc1","-5","1","-3"],["CC(=O)O","-4","4","-6"]]), [["smiles","Caco2","Lipophilicity","Solubility"],["c1ccccc1","-6","2","-4"],["CCO","-7","3","-5"],["c1ccncc1","-5","1","-3"],["CC(=O)O","-4","4","-6"]]);
      await session.step(33, "When user picks \"Chem > Calculate > MPO Score...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > MPO Score..."));
      await session.step(34, "Then \"MPO Score\" dialog should be visible", () => shouldBe(page, el("\"MPO Score\" dialog"), "visible"));
      await session.step(35, "When user selects \"✓ ADME Test\" in Profile input in \"MPO Score\" dialog", () => selectIn(page, "✓ ADME Test", el("Profile input in \"MPO Score\" dialog")));
      await session.step(36, "Then Aggregation input in \"MPO Score\" dialog should have value \"Average\"", () => shouldHaveValue(page, el("Aggregation input in \"MPO Score\" dialog"), "Average"));
      await session.step(37, "And \"Coverage Score\" input in \"MPO Score\" dialog should have value \"false\"", () => shouldHaveValue(page, el("\"Coverage Score\" input in \"MPO Score\" dialog"), "false"));
      await session.step(38, "And OK button in \"MPO Score\" dialog should be enabled", () => shouldBe(page, el("OK button in \"MPO Score\" dialog"), "enabled"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The run appends the profile's score between 0 and 1 for every row", async () => {
      await session.step(42, "When user clicks on OK button in \"MPO Score\" dialog", () => clickOn(page, el("OK button in \"MPO Score\" dialog")));
      await session.step(43, "Then a new column \"MPO ADME Test\" should have been added", () => newColumnNamed(page, "MPO ADME Test"));
      await session.step(44, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(45, "And \"MPO ADME Test\" column should have no missing values", () => columnComplete(page, "MPO ADME Test"));
      await session.step(46, "And every value of \"MPO ADME Test\" column should lie between 0 and 1", () => everyValueBetween(page, "MPO ADME Test", 0, 1));
      await session.step(47, "And the table should have 4 rows", () => rowCount(page, 4));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Coverage Score adds the count of the profile's measured properties, and a second run replaces the score column", async () => {
      await session.step(52, "When user picks \"Chem > Calculate > MPO Score...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > MPO Score..."));
      await session.step(53, "And user selects \"✓ ADME Test\" in Profile input in \"MPO Score\" dialog", () => selectIn(page, "✓ ADME Test", el("Profile input in \"MPO Score\" dialog")));
      await session.step(54, "And user checks \"Coverage Score\" input in \"MPO Score\" dialog", () => check(page, el("\"Coverage Score\" input in \"MPO Score\" dialog")));
      await session.step(55, "And user clicks on OK button in \"MPO Score\" dialog", () => clickOn(page, el("OK button in \"MPO Score\" dialog")));
      await session.step(56, "Then the \"MPO Score\" dialog should close", () => dialogCloses(page, "MPO Score"));
      await session.step(57, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(58, "And a new column \"MPO ADME Test coverage\" should have been added", () => newColumnNamed(page, "MPO ADME Test coverage"));
      await session.step(59, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(60, "And \"MPO ADME Test coverage\" column should have type \"int\"", () => columnType(page, "MPO ADME Test coverage", "int"));
      await session.step(61, "And every value of \"MPO ADME Test coverage\" column should match \"^3$\"", () => everyValueMatches(page, "MPO ADME Test coverage", "^3$"));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
