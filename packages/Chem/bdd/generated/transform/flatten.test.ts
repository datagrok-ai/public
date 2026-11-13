/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/transform/flatten.feature
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
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {everyValueMatches} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {detectTypes} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Flatten Molecules keeps the E/Z of double bonds unless asked", () => {
  const session = feature(test, "features/transform/flatten.feature", import.meta.url);
  test("Flatten Molecules keeps the E/Z of double bonds unless asked", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(12, "And user opens a table \"alkenes\" with:", () => openTableOf(page, "alkenes", [["smiles"],["C/C=C/[C@H](N)C"],["F/C=C\\F"],["C1CC/C=C/CCC1"]]), [["smiles"],["C/C=C/[C@H](N)C"],["F/C=C\\F"],["C1CC/C=C/CCC1"]]);
    await session.step(17, "And the semantic types of the current table are detected", () => detectTypes(page));
    await run.scenario("Flattening keeps the E/Z of double bonds by default", async () => {
      await session.step(20, "When user picks \"Chem > Transform > Flatten Molecules...\" from the top menu", () => pickFromTopMenu(page, "Chem > Transform > Flatten Molecules..."));
      await session.step(21, "Then \"Flatten E/Z\" input in \"Flatten Molecules\" dialog should not be checked", () => shouldNotBe(page, el("\"Flatten E/Z\" input in \"Flatten Molecules\" dialog"), "checked"));
      await session.step(22, "When user clicks on OK button in \"Flatten Molecules\" dialog", () => clickOn(page, el("OK button in \"Flatten Molecules\" dialog")));
      await session.step(23, "Then a new column \"smiles_flat\" should have been added", () => newColumnNamed(page, "smiles_flat"));
      await session.step(24, "And every value of \"smiles_flat\" column should match \"^[^@]*[/\\\\][^@]*$\"", () => everyValueMatches(page, "smiles_flat", "^[^@]*[/\\\\][^@]*$"));
      await session.step(25, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Flatten E/Z flattens the double bonds too", async () => {
      await session.step(28, "When user picks \"Chem > Transform > Flatten Molecules...\" from the top menu", () => pickFromTopMenu(page, "Chem > Transform > Flatten Molecules..."));
      await session.step(29, "And user checks \"Flatten E/Z\" input in \"Flatten Molecules\" dialog", () => check(page, el("\"Flatten E/Z\" input in \"Flatten Molecules\" dialog")));
      await session.step(30, "And user clicks on OK button in \"Flatten Molecules\" dialog", () => clickOn(page, el("OK button in \"Flatten Molecules\" dialog")));
      await session.step(31, "Then a new column \"smiles_flat (2)\" should have been added", () => newColumnNamed(page, "smiles_flat (2)"));
      await session.step(32, "And every value of \"smiles_flat (2)\" column should match \"^[^@/\\\\]*$\"", () => everyValueMatches(page, "smiles_flat (2)", "^[^@/\\\\]*$"));
      await session.step(33, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
