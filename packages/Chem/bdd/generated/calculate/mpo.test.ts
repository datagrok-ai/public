/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/calculate/mpo.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.mpo-profile-crud]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {shippedMpoProfile} from '../../bindings/mpo.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, everyValueBetween} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("MPO Score over a profile's own properties", () => {
  const session = feature(test, "features/calculate/mpo.feature", import.meta.url);
  test("MPO Score over a profile's own properties", {tag: ["@journey", "@realizes:chem.cp.mpo-profile-crud"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And the shipped MPO profile \"ADME Test\" is on the server", () => shippedMpoProfile(page, "ADME Test"));
    await run.scenario("The dialog offers the ADME Test profile, which maps the table's three columns", async () => {
      await session.step(23, "Given user opens a table \"adme\" with:", () => openTableOf(page, "adme", [["smiles","Caco2","Lipophilicity","Solubility"],["c1ccccc1","-6","2","-4"],["CCO","-7","3","-5"],["c1ccncc1","-5","1","-3"],["CC(=O)O","-4","4","-6"]]), [["smiles","Caco2","Lipophilicity","Solubility"],["c1ccccc1","-6","2","-4"],["CCO","-7","3","-5"],["c1ccncc1","-5","1","-3"],["CC(=O)O","-4","4","-6"]]);
      await session.step(29, "When user picks \"Chem > Calculate > MPO Score...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > MPO Score..."));
      await session.step(30, "Then \"MPO Score\" dialog should be visible", () => shouldBe(page, el("\"MPO Score\" dialog"), "visible"));
      await session.step(31, "When user selects \"✓ ADME Test\" in Profile input in \"MPO Score\" dialog", () => selectIn(page, "✓ ADME Test", el("Profile input in \"MPO Score\" dialog")));
      await session.step(32, "Then Aggregation input in \"MPO Score\" dialog should have value \"Average\"", () => shouldHaveValue(page, el("Aggregation input in \"MPO Score\" dialog"), "Average"));
      await session.step(33, "And OK button in \"MPO Score\" dialog should be enabled", () => shouldBe(page, el("OK button in \"MPO Score\" dialog"), "enabled"));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The run appends the profile's score between 0 and 1 for every row", async () => {
      await session.step(37, "When user clicks on OK button in \"MPO Score\" dialog", () => clickOn(page, el("OK button in \"MPO Score\" dialog")));
      await session.step(38, "Then a new column \"MPO ADME Test\" should have been added", () => newColumnNamed(page, "MPO ADME Test"));
      await session.step(39, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(40, "And \"MPO ADME Test\" column should have no missing values", () => columnComplete(page, "MPO ADME Test"));
      await session.step(41, "And every value of \"MPO ADME Test\" column should lie between 0 and 1", () => everyValueBetween(page, "MPO ADME Test", 0, 1));
      await session.step(42, "And the table should have 4 rows", () => rowCount(page, 4));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
