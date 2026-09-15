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
import {bioInitialized, noIsotopeFlag} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, columnTag, columnUnits, everyValueContains} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnMatching, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {callWith, resultContains} from '@datagrok-libraries/bdd/bindings/platform/functions';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("To Atomic Level", () => {
  const session = feature(test, "features/transform/atomic-level.feature", import.meta.url);
  test("To Atomic Level", {tag: ["@journey", "@realizes:bio.transform.to-atomic-level"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(12, "And only \"HELMCoreLibrary.json\" monomer library is selected", () => onlyLibrarySelected(page, "HELMCoreLibrary.json"));
    await run.scenario("A fasta column becomes a molblock column", async () => {
      await session.step(15, "Given user opens filter_FASTA dataset", () => openDataset(page, ds("filter_FASTA")));
      await session.step(16, "When user picks \"Bio > Transform > To Atomic Level...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > To Atomic Level..."));
      await session.step(17, "Then \"To Atomic Level\" dialog should be visible", () => shouldBe(page, el("\"To Atomic Level\" dialog"), "visible"));
      await session.step(18, "And editor of Sequence input in \"To Atomic Level\" dialog should have text \"fasta\"", () => shouldHaveText(page, el("editor of Sequence input in \"To Atomic Level\" dialog"), "fasta"));
      await session.step(19, "And \"Non-linear\" checkbox in \"To Atomic Level\" dialog should be checked", () => shouldBe(page, el("\"Non-linear\" checkbox in \"To Atomic Level\" dialog"), "checked"));
      await session.step(20, "And \"Highlight monomers\" checkbox in \"To Atomic Level\" dialog should be unchecked", () => shouldBe(page, el("\"Highlight monomers\" checkbox in \"To Atomic Level\" dialog"), "unchecked"));
      await session.step(21, "When user clicks on OK button in \"To Atomic Level\" dialog", () => clickOn(page, el("OK button in \"To Atomic Level\" dialog")));
      await session.step(22, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(23, "And a new column matching \"^molfile\\(fasta\\)\" should have been added", () => newColumnMatching(page, "^molfile\\(fasta\\)"));
      await session.step(24, "And \"molfile(fasta)\" column should have semantic type \"Molecule\"", () => columnSemType(page, "molfile(fasta)", "Molecule"));
      await session.step(25, "And \"molfile(fasta)\" column should have units \"molblock\"", () => columnUnits(page, "molfile(fasta)", "molblock"));
      await session.step(26, "And every value of \"molfile(fasta)\" column should contain \"M  V30 BEGIN CTAB\"", () => everyValueContains(page, "molfile(fasta)", "M  V30 BEGIN CTAB"));
      await session.step(27, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A HELM column with cycles goes through the HELM converter", async () => {
      await session.step(30, "Given user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
      await session.step(31, "Then \"HELM string\" column should have units \"helm\"", () => columnUnits(page, "HELM string", "helm"));
      await session.step(32, "When user picks \"Bio > Transform > To Atomic Level...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > To Atomic Level..."));
      await session.step(33, "And user checks \"Highlight monomers\" checkbox in \"To Atomic Level\" dialog", () => check(page, el("\"Highlight monomers\" checkbox in \"To Atomic Level\" dialog")));
      await session.step(34, "And user clicks on OK button in \"To Atomic Level\" dialog", () => clickOn(page, el("OK button in \"To Atomic Level\" dialog")));
      await session.step(35, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(36, "And a new column matching \"^molfile\\(HELM string\\)\" should have been added", () => newColumnMatching(page, "^molfile\\(HELM string\\)"));
      await session.step(37, "And \"molfile(HELM string)\" column should have semantic type \"Molecule\"", () => columnSemType(page, "molfile(HELM string)", "Molecule"));
      await session.step(38, "And \"molfile(HELM string)\" column should have units \"molblock\"", () => columnUnits(page, "molfile(HELM string)", "molblock"));
      await session.step(39, "And \"molfile(HELM string)\" column should have no missing values", () => columnComplete(page, "molfile(HELM string)"));
      await session.step(40, "And every value of \"molfile(HELM string)\" column should contain \"M  V30 BEGIN CTAB\"", () => everyValueContains(page, "molfile(HELM string)", "M  V30 BEGIN CTAB"));
      await session.step(41, "And \"molfile(HELM string)\" column should have tag \".sequence-src-highlight-monomers\" equal to \"true\"", () => columnTag(page, "molfile(HELM string)", ".sequence-src-highlight-monomers", "true"));
      await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The single-sequence functions give clean V3000 molfiles", async () => {
      await session.step(45, "When user calls \"Bio:toAtomicLevelSingleSeq\" function with:", () => callWith(page, "Bio:toAtomicLevelSingleSeq", [["sequence","ACDEFGHIK"]]));
      await session.step(47, "Then the result should contain text \"V3000\"", () => resultContains(page, "V3000"));
      await session.step(48, "And the result should contain text \"M  V30 BEGIN CTAB\"", () => resultContains(page, "M  V30 BEGIN CTAB"));
      await session.step(49, "And the result should not carry an isotope flag on a heavy atom", () => noIsotopeFlag(page));
      await session.step(50, "When user calls \"Bio:seq2atomic\" function with:", () => callWith(page, "Bio:seq2atomic", [["seq","PEPTIDE1{A.C.D.E.F.G.H.I.K}$$$$V2.0"],["nonlinear","true"]]));
      await session.step(53, "Then the result should contain text \"V3000\"", () => resultContains(page, "V3000"));
      await session.step(54, "And the result should not carry an isotope flag on a heavy atom", () => noIsotopeFlag(page));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
