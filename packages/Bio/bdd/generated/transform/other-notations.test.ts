/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/transform/other-notations.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.transform.convert-notation]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '../../bindings/monomer-form.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Converting an aligned column to HELM", () => {
  const session = feature(test, "features/transform/other-notations.feature", import.meta.url);
  test("An aligned column of multi-letter monomers converts to HELM", {tag: ["@realizes:bio.transform.convert-notation"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(22, "Given user opens filter_MSA dataset", () => openDataset(page, ds("filter_MSA")));
    await session.step(23, "When user picks \"Bio > Transform > Convert Sequence Notation...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Convert Sequence Notation..."));
    await session.step(24, "Then \"Convert Sequence Notation\" dialog should be visible", () => shouldBe(page, el("\"Convert Sequence Notation\" dialog"), "visible"));
    await session.step(25, "And \"Convert Sequence Notation\" dialog should contain text \"Current notation: separator\"", () => shouldContainText(page, el("\"Convert Sequence Notation\" dialog"), "Current notation: separator"));
    await session.step(26, "When user selects \"helm\" in \"Convert to\" input in \"Convert Sequence Notation\" dialog", () => selectIn(page, "helm", el("\"Convert to\" input in \"Convert Sequence Notation\" dialog")));
    await session.step(27, "And user clicks on OK button in \"Convert Sequence Notation\" dialog", () => clickOn(page, el("OK button in \"Convert Sequence Notation\" dialog")));
    await session.step(28, "Then 1 new column should have been added", () => newColumnsCount(page, 1));
    await session.step(29, "And a new column \"helm(MSA)\" should have been added", () => newColumnNamed(page, "helm(MSA)"));
    await session.step(30, "And \"helm(MSA)\" column should have units \"helm\"", () => columnUnits(page, "helm(MSA)", "helm"));
    await session.step(31, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(32, "And no errors should have been logged", () => noErrors(page));
  });
});
