/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/other-notations.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.analyze.composition]
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
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {noneSelected, onlyStartingWithSelected, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, clickArea, hasArea, noBalloons, noErrors, painted, propertyShouldBe, readingAtLeast, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Composition on HELM and MSA columns", () => {
  const session = feature(test, "features/analyze/other-notations.feature", import.meta.url);
  test("Composition docks a WebLogo over a HELM column", {tag: ["@realizes:bio.analyze.composition"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(17, "Given user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
    await session.step(18, "When user picks \"Bio > Analyze > Composition\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Composition"));
    await session.step(19, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(20, "And WebLogo viewer should be visible", () => shouldBe(page, el("WebLogo viewer"), "visible"));
    await session.step(21, "And WebLogo viewer should be bound to table \"filter_HELM\"", () => boundTable(page, el("WebLogo viewer"), "filter_HELM"));
    await session.step(22, "And \"Sequence Column Name\" property of WebLogo viewer should be \"HELM string\"", () => propertyShouldBe(page, "Sequence Column Name", el("WebLogo viewer"), "HELM string"));
    await session.step(23, "And WebLogo viewer should be painted", () => painted(page, el("WebLogo viewer")));
    await session.step(24, "And WebLogo viewer should have a \"position 1\" area", () => hasArea(page, el("WebLogo viewer"), "position 1"));
    await session.step(25, "And the \"rows shown\" reading of WebLogo viewer should be 4", () => readingIs(page, "rows shown", el("WebLogo viewer"), 4));
    await session.step(26, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(27, "And no errors should have been logged", () => noErrors(page));
  });
  test("Composition on an aligned column selects the rows of a multi-letter monomer", {tag: ["@realizes:bio.analyze.composition"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(30, "Given user opens filter_MSA dataset", () => openDataset(page, ds("filter_MSA")));
    await session.step(31, "Then \"MSA\" column should have units \"separator\"", () => columnUnits(page, "MSA", "separator"));
    await session.step(32, "When user picks \"Bio > Analyze > Composition\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Composition"));
    await session.step(33, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(34, "And WebLogo viewer should be visible", () => shouldBe(page, el("WebLogo viewer"), "visible"));
    await session.step(35, "And \"Sequence Column Name\" property of WebLogo viewer should be \"MSA\"", () => propertyShouldBe(page, "Sequence Column Name", el("WebLogo viewer"), "MSA"));
    await session.step(36, "And WebLogo viewer should be painted", () => painted(page, el("WebLogo viewer")));
    await session.step(37, "And WebLogo viewer should have a \"monomer meI at position 1\" area", () => hasArea(page, el("WebLogo viewer"), "monomer meI at position 1"));
    await session.step(38, "And the \"positions shown\" reading of WebLogo viewer should be at least 15", () => readingAtLeast(page, "positions shown", el("WebLogo viewer"), 15));
    await session.step(39, "And no rows should be selected", () => noneSelected(page));
    await session.step(40, "When user clicks on the \"monomer meI at position 1\" area of WebLogo viewer", () => clickArea(page, "monomer meI at position 1", el("WebLogo viewer")));
    await session.step(41, "Then some rows should be selected", () => someSelected(page));
    await session.step(42, "And only rows where \"MSA\" starts with \"meI/\" should be selected", () => onlyStartingWithSelected(page, "MSA", "meI/"));
    await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
  });
});
