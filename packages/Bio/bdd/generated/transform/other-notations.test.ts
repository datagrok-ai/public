/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/transform/other-notations.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.calculate.extract-region, bio.transform.convert-notation, bio.transform.split-to-monomers]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, columnTag, columnUnits, everyValueMatches, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, newColumnsCount, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Transforming HELM and aligned columns", () => {
  const session = feature(test, "features/transform/other-notations.feature", import.meta.url);
  test("Extract Region on an aligned column keeps its separator notation", {tag: ["@realizes:bio.calculate.extract-region", "@realizes:bio.transform.convert-notation", "@realizes:bio.transform.split-to-monomers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(19, "Given user opens filter_MSA dataset", () => openDataset(page, ds("filter_MSA")));
    await session.step(20, "When user picks \"Bio > Calculate > Extract Region...\" from the top menu", () => pickFromTopMenu(page, "Bio > Calculate > Extract Region..."));
    await session.step(21, "Then \"Get Sequence Region\" dialog should be visible", () => shouldBe(page, el("\"Get Sequence Region\" dialog"), "visible"));
    await session.step(22, "And Start input in \"Get Sequence Region\" dialog should have value \"1\"", () => shouldHaveValue(page, el("Start input in \"Get Sequence Region\" dialog"), "1"));
    await session.step(23, "And End input in \"Get Sequence Region\" dialog should have value \"17\"", () => shouldHaveValue(page, el("End input in \"Get Sequence Region\" dialog"), "17"));
    await session.step(24, "When user selects \"3\" in Start input in \"Get Sequence Region\" dialog", () => selectIn(page, "3", el("Start input in \"Get Sequence Region\" dialog")));
    await session.step(25, "And user selects \"6\" in End input in \"Get Sequence Region\" dialog", () => selectIn(page, "6", el("End input in \"Get Sequence Region\" dialog")));
    await session.step(26, "And user enters \"region 3-6\" into \"Column name\" input in \"Get Sequence Region\" dialog", () => enterInto(page, "region 3-6", el("\"Column name\" input in \"Get Sequence Region\" dialog")));
    await session.step(27, "And user clicks on OK button in \"Get Sequence Region\" dialog", () => clickOn(page, el("OK button in \"Get Sequence Region\" dialog")));
    await session.step(28, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(29, "And 1 new column should have been added", () => newColumnsCount(page, 1));
    await session.step(30, "And \"region 3-6\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "region 3-6", "Macromolecule"));
    await session.step(31, "And \"region 3-6\" column should have units \"separator\"", () => columnUnits(page, "region 3-6", "separator"));
    await session.step(32, "And \"region 3-6\" column should have tag \"separator\" equal to \"/\"", () => columnTag(page, "region 3-6", "separator", "/"));
    await session.step(33, "And every value of \"region 3-6\" column should match \"^[^/]+(/[^/]+){3}$\"", () => everyValueMatches(page, "region 3-6", "^[^/]+(/[^/]+){3}$"));
    await session.step(34, "And the value of \"region 3-6\" column in row 1 should be \"Aca/N/T/dE\"", () => valueInRow(page, "region 3-6", 1, "Aca/N/T/dE"));
    await session.step(35, "And the value of \"region 3-6\" column in row 2 should be \"Aca/Cys_SEt/T/dK\"", () => valueInRow(page, "region 3-6", 2, "Aca/Cys_SEt/T/dK"));
    await session.step(36, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
  test("Extract Region on a HELM column cuts a HELM region", {tag: ["@realizes:bio.calculate.extract-region", "@realizes:bio.transform.convert-notation", "@realizes:bio.transform.split-to-monomers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(40, "Given user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
    await session.step(41, "When user picks \"Bio > Calculate > Extract Region...\" from the top menu", () => pickFromTopMenu(page, "Bio > Calculate > Extract Region..."));
    await session.step(42, "Then \"Get Sequence Region\" dialog should be visible", () => shouldBe(page, el("\"Get Sequence Region\" dialog"), "visible"));
    await session.step(43, "When user selects \"3\" in Start input in \"Get Sequence Region\" dialog", () => selectIn(page, "3", el("Start input in \"Get Sequence Region\" dialog")));
    await session.step(44, "And user selects \"6\" in End input in \"Get Sequence Region\" dialog", () => selectIn(page, "6", el("End input in \"Get Sequence Region\" dialog")));
    await session.step(45, "And user enters \"region 3-6\" into \"Column name\" input in \"Get Sequence Region\" dialog", () => enterInto(page, "region 3-6", el("\"Column name\" input in \"Get Sequence Region\" dialog")));
    await session.step(46, "And user clicks on OK button in \"Get Sequence Region\" dialog", () => clickOn(page, el("OK button in \"Get Sequence Region\" dialog")));
    await session.step(47, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(48, "And 1 new column should have been added", () => newColumnsCount(page, 1));
    await session.step(49, "And \"region 3-6\" column should have units \"helm\"", () => columnUnits(page, "region 3-6", "helm"));
    await session.step(50, "And every value of \"region 3-6\" column should match \"^PEPTIDE1\\{[^.{}]+(\\.[^.{}]+){3}\\}\\$\\$\\$\\$$\"", () => everyValueMatches(page, "region 3-6", "^PEPTIDE1\\{[^.{}]+(\\.[^.{}]+){3}\\}\\$\\$\\$\\$$"));
    await session.step(51, "And the value of \"region 3-6\" column in row 2 should be \"PEPTIDE1{P.Q.R.S}$$$$\"", () => valueInRow(page, "region 3-6", 2, "PEPTIDE1{P.Q.R.S}$$$$"));
    await session.step(52, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
  });
  test("An aligned column converts to HELM with its multi-letter monomers bracketed", {tag: ["@realizes:bio.calculate.extract-region", "@realizes:bio.transform.convert-notation", "@realizes:bio.transform.split-to-monomers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(57, "Given user opens filter_MSA dataset", () => openDataset(page, ds("filter_MSA")));
    await session.step(58, "When user picks \"Bio > Transform > Convert Sequence Notation...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Convert Sequence Notation..."));
    await session.step(59, "Then \"Convert Sequence Notation\" dialog should be visible", () => shouldBe(page, el("\"Convert Sequence Notation\" dialog"), "visible"));
    await session.step(60, "And \"Convert Sequence Notation\" dialog should contain text \"Current notation: separator\"", () => shouldContainText(page, el("\"Convert Sequence Notation\" dialog"), "Current notation: separator"));
    await session.step(61, "When user selects \"helm\" in \"Convert to\" input in \"Convert Sequence Notation\" dialog", () => selectIn(page, "helm", el("\"Convert to\" input in \"Convert Sequence Notation\" dialog")));
    await session.step(62, "And user clicks on OK button in \"Convert Sequence Notation\" dialog", () => clickOn(page, el("OK button in \"Convert Sequence Notation\" dialog")));
    await session.step(63, "Then 1 new column should have been added", () => newColumnsCount(page, 1));
    await session.step(64, "And a new column \"helm(MSA)\" should have been added", () => newColumnNamed(page, "helm(MSA)"));
    await session.step(65, "And \"helm(MSA)\" column should have units \"helm\"", () => columnUnits(page, "helm(MSA)", "helm"));
    await session.step(66, "And \"helm(MSA)\" column should have no missing values", () => columnComplete(page, "helm(MSA)"));
    await session.step(67, "And every value of \"helm(MSA)\" column should match \"^PEPTIDE1\\{(\\[[^\\]]+\\]|[A-Z*])(\\.(\\[[^\\]]+\\]|[A-Z*]))*\\}\\$\\$\\$\\$$\"", () => everyValueMatches(page, "helm(MSA)", "^PEPTIDE1\\{(\\[[^\\]]+\\]|[A-Z*])(\\.(\\[[^\\]]+\\]|[A-Z*]))*\\}\\$\\$\\$\\$$"));
    await session.step(68, "And the value of \"helm(MSA)\" column in row 1 should be \"PEPTIDE1{[meI].[hHis].[Aca].N.T.[dE].[Thr_PO3H2].[Aca].[D-Tyr_Et].[Tyr_ab-dehydroMe].[dV].E.N.[D-Orn].[D-aThr].*.[Phe_4Me]}$$$$\"", () => valueInRow(page, "helm(MSA)", 1, "PEPTIDE1{[meI].[hHis].[Aca].N.T.[dE].[Thr_PO3H2].[Aca].[D-Tyr_Et].[Tyr_ab-dehydroMe].[dV].E.N.[D-Orn].[D-aThr].*.[Phe_4Me]}$$$$"));
    await session.step(69, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(70, "And no errors should have been logged", () => noErrors(page));
  });
  test("Split to Monomers on a HELM column gives a column per position of the longest peptide", {tag: ["@realizes:bio.calculate.extract-region", "@realizes:bio.transform.convert-notation", "@realizes:bio.transform.split-to-monomers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(73, "Given user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
    await session.step(74, "When user picks \"Bio > Transform > Split to Monomers...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Split to Monomers..."));
    await session.step(75, "Then \"Split to Monomers\" dialog should be visible", () => shouldBe(page, el("\"Split to Monomers\" dialog"), "visible"));
    await session.step(76, "And editor of Sequence input in \"Split to Monomers\" dialog should have text \"HELM string\"", () => shouldHaveText(page, el("editor of Sequence input in \"Split to Monomers\" dialog"), "HELM string"));
    await session.step(77, "When user clicks on OK button in \"Split to Monomers\" dialog", () => clickOn(page, el("OK button in \"Split to Monomers\" dialog")));
    await session.step(78, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(79, "And 10 new columns should have been added", () => newColumnsCount(page, 10));
    await session.step(80, "And \"1\" column should have semantic type \"Monomer\"", () => columnSemType(page, "1", "Monomer"));
    await session.step(81, "And \"10\" column should have semantic type \"Monomer\"", () => columnSemType(page, "10", "Monomer"));
    await session.step(82, "And the value of \"1\" column in row 2 should be \"L\"", () => valueInRow(page, "1", 2, "L"));
    await session.step(83, "And the value of \"7\" column in row 2 should be \"T\"", () => valueInRow(page, "7", 2, "T"));
    await session.step(84, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(85, "And no errors should have been logged", () => noErrors(page));
  });
});
