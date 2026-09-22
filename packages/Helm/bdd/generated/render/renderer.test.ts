/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/render/renderer.feature
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
import {columnSemType, columnTag, columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaColors, hasNoArea, noBalloons, noErrors, readingReads, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("HELM cell renderer", () => {
  const session = feature(test, "features/render/renderer.feature", import.meta.url);
  test("The showcase HELM column is detected and painted in monomer colors", {tag: ["@realizes:helm.cell.helm"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the Helm package is initialized", () => helmInitialized(page));
    await session.step(19, "Given user opens helm-showcase dataset", () => openDataset(page, ds("helm-showcase")));
    await session.step(20, "Then the table should have 53 rows", () => rowCount(page, 53));
    await session.step(21, "And \"HELM\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "HELM", "Macromolecule"));
    await session.step(22, "And \"HELM\" column should have units \"helm\"", () => columnUnits(page, "HELM", "helm"));
    await session.step(23, "And \"HELM\" column should have tag \"quality\" equal to \"Macromolecule\"", () => columnTag(page, "HELM", "quality", "Macromolecule"));
    await session.step(24, "And \"HELM\" column should have tag \"cell.renderer\" equal to \"helm\"", () => columnTag(page, "HELM", "cell.renderer", "helm"));
    await session.step(25, "And the \"cell type of HELM\" reading of grid should be \"helm\"", () => readingReads(page, "cell type of HELM", el("grid"), "helm"));
    await session.step(26, "And the \"cell 2 of HELM\" area of grid should be painted in at least 3 colors", () => areaColors(page, "cell 2 of HELM", el("grid"), 3));
    await session.step(27, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(28, "And no errors should have been logged", () => noErrors(page));
  });
  test("The HELM sample repaints its cells after scrolling away and back", {tag: ["@realizes:helm.cell.helm"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the Helm package is initialized", () => helmInitialized(page));
    await session.step(31, "Given user opens HELM dataset", () => openDataset(page, ds("HELM")));
    await session.step(32, "Then \"HELM\" column should have units \"helm\"", () => columnUnits(page, "HELM", "helm"));
    await session.step(33, "And the \"cell type of HELM\" reading of grid should be \"helm\"", () => readingReads(page, "cell type of HELM", el("grid"), "helm"));
    await session.step(34, "And the \"cell 1 of HELM\" area of grid should be painted in at least 3 colors", () => areaColors(page, "cell 1 of HELM", el("grid"), 3));
    await session.step(35, "When user scrolls the mouse wheel down over the \"cell 1 of HELM\" area of grid", () => wheelOverArea(page, "down", "cell 1 of HELM", el("grid")));
    await session.step(36, "Then grid should not have a \"cell 1 of HELM\" area", () => hasNoArea(page, el("grid"), "cell 1 of HELM"));
    await session.step(37, "And the \"cell 10 of HELM\" area of grid should be painted in at least 3 colors", () => areaColors(page, "cell 10 of HELM", el("grid"), 3));
    await session.step(38, "When user scrolls the mouse wheel up over the \"cell 10 of HELM\" area of grid", () => wheelOverArea(page, "up", "cell 10 of HELM", el("grid")));
    await session.step(39, "Then the \"cell 1 of HELM\" area of grid should be painted in at least 3 colors", () => areaColors(page, "cell 1 of HELM", el("grid"), 3));
    await session.step(40, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
  });
});
