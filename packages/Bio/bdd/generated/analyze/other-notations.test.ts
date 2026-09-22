/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/other-notations.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.analyze.sequence-space, bio.analyze.activity-cliffs, bio.analyze.composition]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {cliffCount} from '../../bindings/bio-a.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnUnits, distinctValues} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnMatching, newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {noneSelected, onlyStartingWithSelected, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {dialogCloses, openDataset, openDatasetRows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, clickArea, hasArea, noBalloons, noErrors, painted, propertyShouldBe, readingAtLeast, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Analyze commands on HELM and MSA columns", () => {
  const session = feature(test, "features/analyze/other-notations.feature", import.meta.url);
  test("Sequence Space embeds a HELM column [notation=HELM, dataset=HELM_sample, column=HELM, units=helm]", {tag: ["@realizes:bio.analyze.sequence-space", "@realizes:bio.analyze.activity-cliffs", "@realizes:bio.analyze.composition"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(22, "Given user opens HELM_sample dataset keeping the first 100 rows", () => openDatasetRows(page, ds("HELM_sample"), 100));
    await session.step(23, "Then \"HELM\" column should have units \"helm\"", () => columnUnits(page, "HELM", "helm"));
    await session.step(24, "When user picks \"Bio > Analyze > Sequence Space...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Sequence Space..."));
    await session.step(25, "Then \"Sequence Space\" dialog should be visible", () => shouldBe(page, el("\"Sequence Space\" dialog"), "visible"));
    await session.step(26, "And editor of Column input in \"Sequence Space\" dialog should have text \"HELM\"", () => shouldHaveText(page, el("editor of Column input in \"Sequence Space\" dialog"), "HELM"));
    await session.step(27, "When user clicks on OK button in \"Sequence Space\" dialog", () => clickOn(page, el("OK button in \"Sequence Space\" dialog")));
    await session.step(28, "Then the \"Sequence Space\" dialog should close", () => dialogCloses(page, "Sequence Space"));
    await session.step(29, "And the top menu command should have completed", () => commandCompleted(page));
    await session.step(30, "And a new column \"Embed_X_1\" should have been added", () => newColumnNamed(page, "Embed_X_1"));
    await session.step(31, "And a new column \"Embed_Y_1\" should have been added", () => newColumnNamed(page, "Embed_Y_1"));
    await session.step(32, "And \"Embed_X_1\" column should have no missing values", () => columnComplete(page, "Embed_X_1"));
    await session.step(33, "And \"Embed_X_1\" column should have at least 10 distinct values", () => distinctValues(page, "Embed_X_1", 10));
    await session.step(34, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(35, "And \"X\" property of scatter plot viewer should be \"Embed_X_1\"", () => propertyShouldBe(page, "X", el("scatter plot viewer"), "Embed_X_1"));
    await session.step(36, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
    await session.step(37, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
  test("Sequence Space embeds a MSA column [notation=MSA, dataset=MSA_sample, column=MSA, units=separator]", {tag: ["@realizes:bio.analyze.sequence-space", "@realizes:bio.analyze.activity-cliffs", "@realizes:bio.analyze.composition"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(22, "Given user opens MSA_sample dataset keeping the first 100 rows", () => openDatasetRows(page, ds("MSA_sample"), 100));
    await session.step(23, "Then \"MSA\" column should have units \"separator\"", () => columnUnits(page, "MSA", "separator"));
    await session.step(24, "When user picks \"Bio > Analyze > Sequence Space...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Sequence Space..."));
    await session.step(25, "Then \"Sequence Space\" dialog should be visible", () => shouldBe(page, el("\"Sequence Space\" dialog"), "visible"));
    await session.step(26, "And editor of Column input in \"Sequence Space\" dialog should have text \"MSA\"", () => shouldHaveText(page, el("editor of Column input in \"Sequence Space\" dialog"), "MSA"));
    await session.step(27, "When user clicks on OK button in \"Sequence Space\" dialog", () => clickOn(page, el("OK button in \"Sequence Space\" dialog")));
    await session.step(28, "Then the \"Sequence Space\" dialog should close", () => dialogCloses(page, "Sequence Space"));
    await session.step(29, "And the top menu command should have completed", () => commandCompleted(page));
    await session.step(30, "And a new column \"Embed_X_1\" should have been added", () => newColumnNamed(page, "Embed_X_1"));
    await session.step(31, "And a new column \"Embed_Y_1\" should have been added", () => newColumnNamed(page, "Embed_Y_1"));
    await session.step(32, "And \"Embed_X_1\" column should have no missing values", () => columnComplete(page, "Embed_X_1"));
    await session.step(33, "And \"Embed_X_1\" column should have at least 10 distinct values", () => distinctValues(page, "Embed_X_1", 10));
    await session.step(34, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(35, "And \"X\" property of scatter plot viewer should be \"Embed_X_1\"", () => propertyShouldBe(page, "X", el("scatter plot viewer"), "Embed_X_1"));
    await session.step(36, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
    await session.step(37, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
  test("Activity Cliffs finds cliffs in a HELM column [notation=HELM, dataset=HELM_sample, column=HELM]", {tag: ["@realizes:bio.analyze.sequence-space", "@realizes:bio.analyze.activity-cliffs", "@realizes:bio.analyze.composition"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(46, "Given user opens HELM_sample dataset keeping the first 100 rows", () => openDatasetRows(page, ds("HELM_sample"), 100));
    await session.step(47, "When user picks \"Bio > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Activity Cliffs..."));
    await session.step(48, "Then \"Sequence Activity Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Sequence Activity Cliffs\" dialog"), "visible"));
    await session.step(49, "And editor of Column input in \"Sequence Activity Cliffs\" dialog should have text \"HELM\"", () => shouldHaveText(page, el("editor of Column input in \"Sequence Activity Cliffs\" dialog"), "HELM"));
    await session.step(50, "When user selects \"Activity\" in Activities input in \"Sequence Activity Cliffs\" dialog", () => selectIn(page, "Activity", el("Activities input in \"Sequence Activity Cliffs\" dialog")));
    await session.step(51, "And user clicks on OK button in \"Sequence Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Sequence Activity Cliffs\" dialog")));
    await session.step(52, "Then the \"Sequence Activity Cliffs\" dialog should close", () => dialogCloses(page, "Sequence Activity Cliffs"));
    await session.step(53, "And the top menu command should have completed", () => commandCompleted(page));
    await session.step(54, "And a new column \"Embed_X_1\" should have been added", () => newColumnNamed(page, "Embed_X_1"));
    await session.step(55, "And a new column matching \"sali|SALI\" should have been added", () => newColumnMatching(page, "sali|SALI"));
    await session.step(56, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(57, "And title of scatter plot viewer should have text \"Activity cliffs\"", () => shouldHaveText(page, el("title of scatter plot viewer"), "Activity cliffs"));
    await session.step(58, "And the activity cliffs plot should report at least 1 cliff", () => cliffCount(page, 1));
    await session.step(59, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
    await session.step(60, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(61, "And no errors should have been logged", () => noErrors(page));
  });
  test("Activity Cliffs finds cliffs in a MSA column [notation=MSA, dataset=MSA_sample, column=MSA]", {tag: ["@realizes:bio.analyze.sequence-space", "@realizes:bio.analyze.activity-cliffs", "@realizes:bio.analyze.composition"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(46, "Given user opens MSA_sample dataset keeping the first 100 rows", () => openDatasetRows(page, ds("MSA_sample"), 100));
    await session.step(47, "When user picks \"Bio > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Activity Cliffs..."));
    await session.step(48, "Then \"Sequence Activity Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Sequence Activity Cliffs\" dialog"), "visible"));
    await session.step(49, "And editor of Column input in \"Sequence Activity Cliffs\" dialog should have text \"MSA\"", () => shouldHaveText(page, el("editor of Column input in \"Sequence Activity Cliffs\" dialog"), "MSA"));
    await session.step(50, "When user selects \"Activity\" in Activities input in \"Sequence Activity Cliffs\" dialog", () => selectIn(page, "Activity", el("Activities input in \"Sequence Activity Cliffs\" dialog")));
    await session.step(51, "And user clicks on OK button in \"Sequence Activity Cliffs\" dialog", () => clickOn(page, el("OK button in \"Sequence Activity Cliffs\" dialog")));
    await session.step(52, "Then the \"Sequence Activity Cliffs\" dialog should close", () => dialogCloses(page, "Sequence Activity Cliffs"));
    await session.step(53, "And the top menu command should have completed", () => commandCompleted(page));
    await session.step(54, "And a new column \"Embed_X_1\" should have been added", () => newColumnNamed(page, "Embed_X_1"));
    await session.step(55, "And a new column matching \"sali|SALI\" should have been added", () => newColumnMatching(page, "sali|SALI"));
    await session.step(56, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(57, "And title of scatter plot viewer should have text \"Activity cliffs\"", () => shouldHaveText(page, el("title of scatter plot viewer"), "Activity cliffs"));
    await session.step(58, "And the activity cliffs plot should report at least 1 cliff", () => cliffCount(page, 1));
    await session.step(59, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
    await session.step(60, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(61, "And no errors should have been logged", () => noErrors(page));
  });
  test("Composition docks a WebLogo over a HELM column", {tag: ["@realizes:bio.analyze.sequence-space", "@realizes:bio.analyze.activity-cliffs", "@realizes:bio.analyze.composition"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(69, "Given user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
    await session.step(70, "When user picks \"Bio > Analyze > Composition\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Composition"));
    await session.step(71, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(72, "And WebLogo viewer should be visible", () => shouldBe(page, el("WebLogo viewer"), "visible"));
    await session.step(73, "And WebLogo viewer should be bound to table \"filter_HELM\"", () => boundTable(page, el("WebLogo viewer"), "filter_HELM"));
    await session.step(74, "And \"Sequence Column Name\" property of WebLogo viewer should be \"HELM string\"", () => propertyShouldBe(page, "Sequence Column Name", el("WebLogo viewer"), "HELM string"));
    await session.step(75, "And WebLogo viewer should be painted", () => painted(page, el("WebLogo viewer")));
    await session.step(76, "And WebLogo viewer should have a \"position 1\" area", () => hasArea(page, el("WebLogo viewer"), "position 1"));
    await session.step(77, "And the \"rows shown\" reading of WebLogo viewer should be 4", () => readingIs(page, "rows shown", el("WebLogo viewer"), 4));
    await session.step(78, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(79, "And no errors should have been logged", () => noErrors(page));
  });
  test("Composition on an aligned column selects the rows of a multi-letter monomer", {tag: ["@realizes:bio.analyze.sequence-space", "@realizes:bio.analyze.activity-cliffs", "@realizes:bio.analyze.composition"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(82, "Given user opens filter_MSA dataset", () => openDataset(page, ds("filter_MSA")));
    await session.step(83, "Then \"MSA\" column should have units \"separator\"", () => columnUnits(page, "MSA", "separator"));
    await session.step(84, "When user picks \"Bio > Analyze > Composition\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Composition"));
    await session.step(85, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(86, "And WebLogo viewer should be visible", () => shouldBe(page, el("WebLogo viewer"), "visible"));
    await session.step(87, "And \"Sequence Column Name\" property of WebLogo viewer should be \"MSA\"", () => propertyShouldBe(page, "Sequence Column Name", el("WebLogo viewer"), "MSA"));
    await session.step(88, "And WebLogo viewer should be painted", () => painted(page, el("WebLogo viewer")));
    await session.step(89, "And WebLogo viewer should have a \"monomer meI at position 1\" area", () => hasArea(page, el("WebLogo viewer"), "monomer meI at position 1"));
    await session.step(90, "And the \"positions shown\" reading of WebLogo viewer should be at least 15", () => readingAtLeast(page, "positions shown", el("WebLogo viewer"), 15));
    await session.step(91, "And no rows should be selected", () => noneSelected(page));
    await session.step(92, "When user clicks on the \"monomer meI at position 1\" area of WebLogo viewer", () => clickArea(page, "monomer meI at position 1", el("WebLogo viewer")));
    await session.step(93, "Then some rows should be selected", () => someSelected(page));
    await session.step(94, "And only rows where \"MSA\" starts with \"meI/\" should be selected", () => onlyStartingWithSelected(page, "MSA", "meI/"));
    await session.step(95, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(96, "And no errors should have been logged", () => noErrors(page));
  });
});
