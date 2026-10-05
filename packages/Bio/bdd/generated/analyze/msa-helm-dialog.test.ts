/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/msa-helm-dialog.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.analyze.msa]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnType, columnUnits, distinctValues} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {noNewColumn, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {addCalculated} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProjectWithTable, ownProjectGone, projectsOnServer, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The MSA dialog on a HELM column, before the PepSeA engine runs", () => {
  const session = feature(test, "features/analyze/msa-helm-dialog.feature", import.meta.url);
  test("The cluster column and both PepSeA gap penalties are offered on a HELM column", {tag: ["@realizes:bio.analyze.msa"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(25, "Given user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
    await session.step(26, "When user adds a calculated column \"Clusters\" with formula \"Length(${HELM string}) % 2\"", () => addCalculated(page, "Clusters", "Length(${HELM string}) % 2"));
    await session.step(27, "Then \"Clusters\" column should have type \"int\"", () => columnType(page, "Clusters", "int"));
    await session.step(28, "And \"Clusters\" column should have at least 2 distinct values", () => distinctValues(page, "Clusters", 2));
    await session.step(29, "When user picks \"Bio > Analyze > MSA...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > MSA..."));
    await session.step(30, "Then MSA dialog should be visible", () => shouldBe(page, el("MSA dialog"), "visible"));
    await session.step(31, "And Engine input in MSA dialog should have value \"PepSeA\"", () => shouldHaveValue(page, el("Engine input in MSA dialog"), "PepSeA"));
    await session.step(32, "And Clusters input in MSA dialog should be visible", () => shouldBe(page, el("Clusters input in MSA dialog"), "visible"));
    await session.step(33, "When user selects \"Clusters\" in Clusters input in MSA dialog", () => selectIn(page, "Clusters", el("Clusters input in MSA dialog")));
    await session.step(34, "Then editor of Clusters input in MSA dialog should have text \"Clusters\"", () => shouldHaveText(page, el("editor of Clusters input in MSA dialog"), "Clusters"));
    await session.step(35, "And \"Gap Open\" input in MSA dialog should have value \"1.53\"", () => shouldHaveValue(page, el("\"Gap Open\" input in MSA dialog"), "1.53"));
    await session.step(36, "And \"Gap Extend\" input in MSA dialog should have value \"0\"", () => shouldHaveValue(page, el("\"Gap Extend\" input in MSA dialog"), "0"));
    await session.step(37, "When user clicks on \"Alignment parameters\" button in MSA dialog", () => clickOn(page, el("\"Alignment parameters\" button in MSA dialog")));
    await session.step(38, "Then \"Gap Open\" input in MSA dialog should be hidden", () => shouldBe(page, el("\"Gap Open\" input in MSA dialog"), "hidden"));
    await session.step(39, "And \"Gap Extend\" input in MSA dialog should be hidden", () => shouldBe(page, el("\"Gap Extend\" input in MSA dialog"), "hidden"));
    await session.step(40, "When user clicks on \"Alignment parameters\" button in MSA dialog", () => clickOn(page, el("\"Alignment parameters\" button in MSA dialog")));
    await session.step(41, "Then \"Gap Open\" input in MSA dialog should be visible", () => shouldBe(page, el("\"Gap Open\" input in MSA dialog"), "visible"));
    await session.step(42, "And \"Gap Extend\" input in MSA dialog should be visible", () => shouldBe(page, el("\"Gap Extend\" input in MSA dialog"), "visible"));
    await session.step(43, "When user clicks on CANCEL button in MSA dialog", () => clickOn(page, el("CANCEL button in MSA dialog")));
    await session.step(44, "Then MSA dialog should be hidden", () => shouldBe(page, el("MSA dialog"), "hidden"));
    await session.step(45, "And no new column should have been added", () => noNewColumn(page));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
  });
  test("A HELM table reopened from a project offers the same engine", {tag: ["@realizes:bio.analyze.msa"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(49, "Given the user's own project \"bdd-bio-msa-helm-{run}\" is removed now and at feature end", () => ownProjectGone(page, session.text("bdd-bio-msa-helm-{run}")));
    await session.step(50, "And user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
    await session.step(51, "When user saves the current view as project \"bdd-bio-msa-helm-{run}\"", () => saveAsProject(page, session.text("bdd-bio-msa-helm-{run}")));
    await session.step(52, "Then 1 project named \"bdd-bio-msa-helm-{run}\" should be on the server", () => projectsOnServer(page, 1, session.text("bdd-bio-msa-helm-{run}")));
    await session.step(53, "When user closes all views", () => closeAllViews(page));
    await session.step(54, "And user opens the \"bdd-bio-msa-helm-{run}\" project and waits for its table", () => openProjectWithTable(page, session.text("bdd-bio-msa-helm-{run}")));
    await session.step(55, "Then \"HELM string\" column should have units \"helm\"", () => columnUnits(page, "HELM string", "helm"));
    await session.step(56, "When user picks \"Bio > Analyze > MSA...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > MSA..."));
    await session.step(57, "Then MSA dialog should be visible", () => shouldBe(page, el("MSA dialog"), "visible"));
    await session.step(58, "And editor of Sequence input in MSA dialog should have text \"HELM string\"", () => shouldHaveText(page, el("editor of Sequence input in MSA dialog"), "HELM string"));
    await session.step(59, "And Engine input in MSA dialog should have value \"PepSeA\"", () => shouldHaveValue(page, el("Engine input in MSA dialog"), "PepSeA"));
    await session.step(60, "And Method input in MSA dialog should have value \"mafft --auto\"", () => shouldHaveValue(page, el("Method input in MSA dialog"), "mafft --auto"));
    await session.step(61, "When user clicks on CANCEL button in MSA dialog", () => clickOn(page, el("CANCEL button in MSA dialog")));
    await session.step(62, "Then MSA dialog should be hidden", () => shouldBe(page, el("MSA dialog"), "hidden"));
    await session.step(63, "And no new column should have been added", () => noNewColumn(page));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
  });
});
