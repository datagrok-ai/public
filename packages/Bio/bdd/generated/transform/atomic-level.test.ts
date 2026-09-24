/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/transform/atomic-level.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.transform.to-atomic-level]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {onlyLibrarySelected} from '../../bindings/monomer-libs.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, isExpanded, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, columnTag, columnUnits, everyValueContains, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnMatching, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, contextPanelOpen, contextPanelShows, dialogCloses, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("To Atomic Level", () => {
  const session = feature(test, "features/transform/atomic-level.feature", import.meta.url);
  test("To Atomic Level", {tag: ["@journey", "@serial", "@realizes:bio.transform.to-atomic-level"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And only \"HELMCoreLibrary.json\" monomer library is selected", () => onlyLibrarySelected(page, "HELMCoreLibrary.json"));
    await run.scenario("A fasta column becomes a molblock column", async () => {
      await session.step(22, "Given user opens filter_FASTA dataset", () => openDataset(page, ds("filter_FASTA")));
      await session.step(23, "Then \"fasta\" column should have units \"fasta\"", () => columnUnits(page, "fasta", "fasta"));
      await session.step(24, "When user picks \"Bio > Transform > To Atomic Level...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > To Atomic Level..."));
      await session.step(25, "Then \"To Atomic Level\" dialog should be visible", () => shouldBe(page, el("\"To Atomic Level\" dialog"), "visible"));
      await session.step(26, "And editor of Sequence input in \"To Atomic Level\" dialog should have text \"fasta\"", () => shouldHaveText(page, el("editor of Sequence input in \"To Atomic Level\" dialog"), "fasta"));
      await session.step(27, "And \"Non-linear\" checkbox in \"To Atomic Level\" dialog should be checked", () => shouldBe(page, el("\"Non-linear\" checkbox in \"To Atomic Level\" dialog"), "checked"));
      await session.step(28, "And \"Highlight monomers\" checkbox in \"To Atomic Level\" dialog should be unchecked", () => shouldBe(page, el("\"Highlight monomers\" checkbox in \"To Atomic Level\" dialog"), "unchecked"));
      await session.step(29, "When user clicks on OK button in \"To Atomic Level\" dialog", () => clickOn(page, el("OK button in \"To Atomic Level\" dialog")));
      await session.step(30, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(31, "And a new column matching \"^molfile\\(fasta\\)\" should have been added", () => newColumnMatching(page, "^molfile\\(fasta\\)"));
      await session.step(32, "And \"molfile(fasta)\" column should have semantic type \"Molecule\"", () => columnSemType(page, "molfile(fasta)", "Molecule"));
      await session.step(33, "And \"molfile(fasta)\" column should have units \"molblock\"", () => columnUnits(page, "molfile(fasta)", "molblock"));
      await session.step(34, "And every value of \"molfile(fasta)\" column should contain \"M  V30 BEGIN CTAB\"", () => everyValueContains(page, "molfile(fasta)", "M  V30 BEGIN CTAB"));
      await session.step(35, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A HELM column with cycles goes through the HELM converter", async () => {
      await session.step(38, "Given user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
      await session.step(39, "Then \"HELM string\" column should have units \"helm\"", () => columnUnits(page, "HELM string", "helm"));
      await session.step(40, "When user picks \"Bio > Transform > To Atomic Level...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > To Atomic Level..."));
      await session.step(41, "And user checks \"Highlight monomers\" checkbox in \"To Atomic Level\" dialog", () => check(page, el("\"Highlight monomers\" checkbox in \"To Atomic Level\" dialog")));
      await session.step(42, "And user clicks on OK button in \"To Atomic Level\" dialog", () => clickOn(page, el("OK button in \"To Atomic Level\" dialog")));
      await session.step(43, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(44, "And a new column matching \"^molfile\\(HELM string\\)\" should have been added", () => newColumnMatching(page, "^molfile\\(HELM string\\)"));
      await session.step(45, "And \"molfile(HELM string)\" column should have semantic type \"Molecule\"", () => columnSemType(page, "molfile(HELM string)", "Molecule"));
      await session.step(46, "And \"molfile(HELM string)\" column should have units \"molblock\"", () => columnUnits(page, "molfile(HELM string)", "molblock"));
      await session.step(47, "And \"molfile(HELM string)\" column should have no missing values", () => columnComplete(page, "molfile(HELM string)"));
      await session.step(48, "And every value of \"molfile(HELM string)\" column should contain \"M  V30 BEGIN CTAB\"", () => everyValueContains(page, "molfile(HELM string)", "M  V30 BEGIN CTAB"));
      await session.step(49, "And \"molfile(HELM string)\" column should have tag \".sequence-src-highlight-monomers\" equal to \"true\"", () => columnTag(page, "molfile(HELM string)", ".sequence-src-highlight-monomers", "true"));
      await session.step(50, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("An aligned column of multi-letter monomers becomes molblocks", async () => {
      await session.step(53, "Given user opens filter_MSA dataset", () => openDataset(page, ds("filter_MSA")));
      await session.step(54, "Then \"MSA\" column should have units \"separator\"", () => columnUnits(page, "MSA", "separator"));
      await session.step(55, "When user picks \"Bio > Transform > To Atomic Level...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > To Atomic Level..."));
      await session.step(56, "Then \"To Atomic Level\" dialog should be visible", () => shouldBe(page, el("\"To Atomic Level\" dialog"), "visible"));
      await session.step(57, "And editor of Sequence input in \"To Atomic Level\" dialog should have text \"MSA\"", () => shouldHaveText(page, el("editor of Sequence input in \"To Atomic Level\" dialog"), "MSA"));
      await session.step(58, "When user clicks on OK button in \"To Atomic Level\" dialog", () => clickOn(page, el("OK button in \"To Atomic Level\" dialog")));
      await session.step(59, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(60, "And a new column matching \"^molfile\\(MSA\\)\" should have been added", () => newColumnMatching(page, "^molfile\\(MSA\\)"));
      await session.step(61, "And \"molfile(MSA)\" column should have semantic type \"Molecule\"", () => columnSemType(page, "molfile(MSA)", "Molecule"));
      await session.step(62, "And \"molfile(MSA)\" column should have no missing values", () => columnComplete(page, "molfile(MSA)"));
      await session.step(63, "And every value of \"molfile(MSA)\" column should contain \"M  V30 BEGIN CTAB\"", () => everyValueContains(page, "molfile(MSA)", "M  V30 BEGIN CTAB"));
      await session.step(64, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The column's own action opens the same conversion on that column", async () => {
      await session.step(67, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(68, "And user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
      await session.step(69, "When user clicks on the \"header HELM string\" area of grid", () => clickArea(page, "header HELM string", el("grid")));
      await session.step(70, "Then the context panel should show \"HELM string\"", () => contextPanelShows(page, "HELM string"));
      await session.step(71, "Given Actions pane in context panel is expanded", () => isExpanded(page, el("Actions pane in context panel")));
      await session.step(72, "When user clicks on \"To Atomic Level...\" label in Actions pane in context panel", () => clickOn(page, el("\"To Atomic Level...\" label in Actions pane in context panel")));
      await session.step(73, "Then \"To Atomic Level\" dialog should be visible", () => shouldBe(page, el("\"To Atomic Level\" dialog"), "visible"));
      await session.step(74, "And editor of Sequence input in \"To Atomic Level\" dialog should have text \"HELM string\"", () => shouldHaveText(page, el("editor of Sequence input in \"To Atomic Level\" dialog"), "HELM string"));
      await session.step(75, "When user clicks on OK button in \"To Atomic Level\" dialog", () => clickOn(page, el("OK button in \"To Atomic Level\" dialog")));
      await session.step(76, "Then the \"To Atomic Level\" dialog should close", () => dialogCloses(page, "To Atomic Level"));
      await session.step(77, "And the table should have a column \"molfile(HELM string)\"", () => hasColumn(page, "molfile(HELM string)"));
      await session.step(78, "And \"molfile(HELM string)\" column should have units \"molblock\"", () => columnUnits(page, "molfile(HELM string)", "molblock"));
      await session.step(79, "And \"molfile(HELM string)\" column should have no missing values", () => columnComplete(page, "molfile(HELM string)"));
      await session.step(80, "And every value of \"molfile(HELM string)\" column should contain \"M  V30 BEGIN CTAB\"", () => everyValueContains(page, "molfile(HELM string)", "M  V30 BEGIN CTAB"));
      await session.step(81, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
