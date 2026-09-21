/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/helm-to-structures-to-smiles.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {everyValueContains, everyValueMatches, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnMatching, newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("HELM sequences to structures, then to SMILES", () => {
  const session = feature(test, "features/guides/helm-to-structures-to-smiles.feature", import.meta.url);
  test("Build molecules from the HELM sequences, then write them out as SMILES", {tag: ["@guide", "@help:datagrok/solutions/domains/bio"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens filter_HELM dataset", () => openDataset(page, ds("filter_HELM")));
    await session.step(14, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(15, "When user picks \"Bio > Transform > To Atomic Level...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > To Atomic Level..."));
    await session.step(16, "And user clicks on OK button in \"To Atomic Level\" dialog", () => clickOn(page, el("OK button in \"To Atomic Level\" dialog")));
    await session.step(17, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(18, "And a new column matching \"^molfile\\(HELM string\\)\" should have been added", () => newColumnMatching(page, "^molfile\\(HELM string\\)"));
    await session.step(19, "And every value of \"molfile(HELM string)\" column should contain \"V30 BEGIN CTAB\"", () => everyValueContains(page, "molfile(HELM string)", "V30 BEGIN CTAB"));
    await session.step(20, "When user picks \"Chem > Transform > Convert Notation...\" from the top menu", () => pickFromTopMenu(page, "Chem > Transform > Convert Notation..."));
    await session.step(21, "Then \"Target Notation\" input in \"Convert Notation\" dialog should have value \"smiles\"", () => shouldHaveValue(page, el("\"Target Notation\" input in \"Convert Notation\" dialog"), "smiles"));
    await session.step(22, "When user clicks on OK button in \"Convert Notation\" dialog", () => clickOn(page, el("OK button in \"Convert Notation\" dialog")));
    await session.step(23, "Then the top menu command should have completed", () => commandCompleted(page));
    await session.step(24, "And a new column \"molfile(HELM string)_smiles\" should have been added", () => newColumnNamed(page, "molfile(HELM string)_smiles"));
    await session.step(25, "And every value of \"molfile(HELM string)_smiles\" column should match \"^[A-Za-z0-9@+()\\[\\]\\\\/%=#$.:-]+$\"", () => everyValueMatches(page, "molfile(HELM string)_smiles", "^[A-Za-z0-9@+()\\[\\]\\\\/%=#$.:-]+$"));
    await session.step(26, "When user picks \"Column Properties...\" from the context menu of the \"header molfile(HELM string)_smiles\" area of grid", () => pickFromAreaContextMenu(page, "Column Properties...", "header molfile(HELM string)_smiles", el("grid")));
    await session.step(27, "And user enters \"canonical_smiles\" into \"New name\" input in \"molfile(HELM string)_smiles\" dialog", () => enterInto(page, "canonical_smiles", el("\"New name\" input in \"molfile(HELM string)_smiles\" dialog")));
    await session.step(28, "And user clicks on OK button in \"molfile(HELM string)_smiles\" dialog", () => clickOn(page, el("OK button in \"molfile(HELM string)_smiles\" dialog")));
    await session.step(29, "Then the table should have a column \"canonical_smiles\"", () => hasColumn(page, "canonical_smiles"));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
