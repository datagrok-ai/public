/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/integration/bio-menu.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [helm.cell.helm]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {helmInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnTag, columnUnits, someValueContains} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {walkToColumn} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Bio menu on the HELM showcase leaves the Helm renderer intact", () => {
  const session = feature(test, "features/integration/bio-menu.feature", import.meta.url);
  test("The Bio menu on the HELM showcase leaves the Helm renderer intact", {tag: ["@journey", "@realizes:helm.cell.helm"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the Helm package is initialized", () => helmInitialized(page));
    await session.step(31, "And user opens helm-showcase dataset", () => openDataset(page, ds("helm-showcase")));
    await session.step(32, "Then the table should have 53 rows", () => rowCount(page, 53));
    await session.step(33, "And \"HELM\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "HELM", "Macromolecule"));
    await run.scenario("Scan Liabilities annotates the showcase column", async () => {
      await session.step(36, "When user picks \"Bio > Annotate > Scan Liabilities...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Scan Liabilities..."));
      await session.step(37, "Then \"Scan Sequence Liabilities\" dialog should be visible", () => shouldBe(page, el("\"Scan Sequence Liabilities\" dialog"), "visible"));
      await session.step(38, "When user clicks on OK button in \"Scan Sequence Liabilities\" dialog", () => clickOn(page, el("OK button in \"Scan Sequence Liabilities\" dialog")));
      await session.step(39, "Then a new column \"~HELM_annotations\" should have been added", () => newColumnNamed(page, "~HELM_annotations"));
      await session.step(40, "And some value of \"~HELM_annotations\" column should contain \"oxid-m\"", () => someValueContains(page, "~HELM_annotations", "oxid-m"));
      await session.step(41, "And some value of \"~HELM_annotations\" column should contain \"oxid-w\"", () => someValueContains(page, "~HELM_annotations", "oxid-w"));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
      await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("After the sweep the column is still a HELM column painted by Helm", async () => {
      await session.step(46, "When user moves the current cell of grid to the \"HELM\" column", () => walkToColumn(page, el("grid"), "HELM"));
      await session.step(47, "Then the \"current row\" reading of grid should be 1", () => readingIs(page, "current row", el("grid"), 1));
      await session.step(48, "And the table should have 53 rows", () => rowCount(page, 53));
      await session.step(49, "And \"HELM\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "HELM", "Macromolecule"));
      await session.step(50, "And \"HELM\" column should have units \"helm\"", () => columnUnits(page, "HELM", "helm"));
      await session.step(51, "And \"HELM\" column should have tag \"cell.renderer\" equal to \"helm\"", () => columnTag(page, "HELM", "cell.renderer", "helm"));
      await session.step(52, "And the \"cell type of HELM\" reading of grid should be \"helm\"", () => readingReads(page, "cell type of HELM", el("grid"), "helm"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
      await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
