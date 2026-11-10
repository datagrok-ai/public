/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/03b-helm-showcase-and-editor.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {helmInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterToAnyOf} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {doubleClickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight, gridShowsColumns, hoverMonomer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · HELM rendering and editor", () => {
  const session = feature(test, "features/demos/03b-helm-showcase-and-editor.feature", import.meta.url);
  test("Any HELM drawn in the grid, edited in one editor", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And the Helm package is initialized", () => helmInitialized(page));
    await session.step(15, "And user opens helm-showcase dataset", () => openDataset(page, ds("helm-showcase")));
    await session.step(16, "Then \"HELM\" column should have units \"helm\"", () => columnUnits(page, "HELM", "helm"));
    await session.step(17, "Given user filters rows where \"Name\" is one of \"Disulfide-rich (insulin-like), Lasso peptide (macrolactam + tail), siRNA with ss/as strand annotations, Cloverleaf — four-arm tRNA, Peptide – SMCC – RNA conjugate\"", () => filterToAnyOf(page, "Name", "Disulfide-rich (insulin-like), Lasso peptide (macrolactam + tail), siRNA with ss/as strand annotations, Cloverleaf — four-arm tRNA, Peptide – SMCC – RNA conjugate"));
    await session.step(18, "And the grid shows only the columns \"Name, Category, HELM\"", () => gridShowsColumns(page, "Name, Category, HELM"));
    await session.step(19, "And the \"Name\" column of the grid is 230 px wide", () => gridColumnWidth(page, "Name", 230));
    await session.step(20, "And the \"HELM\" column of the grid is 900 px wide", () => gridColumnWidth(page, "HELM", 900));
    await session.step(21, "And the rows of the grid are 150 px tall", () => gridRowHeight(page, 150));
    await session.step(22, "When user hovers over monomer 20 of the \"cell 33 of HELM\" area of grid", () => hoverMonomer(page, 20, "cell 33 of HELM", el("grid")));
    await session.step(23, "And the video says \"One renderer for peptides, duplexes, folded RNA and conjugates\"", () => videoSays(page, "One renderer for peptides, duplexes, folded RNA and conjugates"));
    await session.step(24, "And user hovers over monomer 6 of the \"cell 14 of HELM\" area of grid", () => hoverMonomer(page, 6, "cell 14 of HELM", el("grid")));
    await session.step(25, "And user double-clicks on the \"cell 14 of HELM\" area of grid", () => doubleClickArea(page, "cell 14 of HELM", el("grid")));
    await session.step(26, "Then HELM editor should be visible", () => shouldBe(page, el("HELM editor"), "visible"));
    await session.step(27, "When user clicks on HELM tab", () => clickOn(page, el("HELM tab")));
    await session.step(28, "Then notation pane should contain text \"PEPTIDE1\"", () => shouldContainText(page, el("notation pane"), "PEPTIDE1"));
    await session.step(29, "When user clicks on Properties tab", () => clickOn(page, el("Properties tab")));
    await session.step(30, "Then formula field should be visible", () => shouldBe(page, el("formula field"), "visible"));
    await session.step(31, "When user clicks on \"Molecular Structure\" tab in HELM editor", () => clickOn(page, el("\"Molecular Structure\" tab in HELM editor")));
    await session.step(32, "And user clicks on CANCEL button in HELM editor", () => clickOn(page, el("CANCEL button in HELM editor")));
    await session.step(33, "Then HELM editor should be absent", () => shouldBe(page, el("HELM editor"), "absent"));
    await session.step(34, "And no errors should have been logged", () => noErrors(page));
  });
});
