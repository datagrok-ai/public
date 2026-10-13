/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/editor/palette.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [helm.cell-editor.molecule]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {clickEmptyCanvas, helmInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, rememberVisibleCount, shouldBe, shouldContainText, shouldHaveText, typeInto, visibleCount, visibleFewerThanRemembered} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnUnits, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {doubleClickArea, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The HELM editor's monomer palette", () => {
  const session = feature(test, "features/editor/palette.feature", import.meta.url);
  test("The HELM editor's monomer palette", {tag: ["@journey", "@realizes:helm.cell-editor.molecule"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the Helm package is initialized", () => helmInitialized(page));
    await session.step(14, "And user opens helm-showcase dataset", () => openDataset(page, ds("helm-showcase")));
    await session.step(15, "Then \"HELM\" column should have units \"helm\"", () => columnUnits(page, "HELM", "helm"));
    await session.step(16, "When user double-clicks on the \"cell 1 of HELM\" area of grid", () => doubleClickArea(page, "cell 1 of HELM", el("grid")));
    await session.step(17, "Then HELM editor should be visible", () => shouldBe(page, el("HELM editor"), "visible"));
    await run.scenario("The Peptides tab lists the monomers and the search narrows them to a match", async () => {
      await session.step(20, "When user clicks on Peptides palette tab", () => clickOn(page, el("Peptides palette tab")));
      await session.step(21, "Then G monomer tile should be visible", () => shouldBe(page, el("G monomer tile"), "visible"));
      await session.step(22, "And Aca monomer tile should be visible", () => shouldBe(page, el("Aca monomer tile"), "visible"));
      await session.step(25, "When user remembers the number of visible monomer tiles", () => rememberVisibleCount(page, el("monomer tiles")));
      await session.step(26, "And user types \"Aca\" into palette search", () => typeInto(page, "Aca", el("palette search")));
      await session.step(27, "Then there should be fewer visible monomer tiles than remembered", () => visibleFewerThanRemembered(page, el("monomer tiles")));
      await session.step(28, "And Aca monomer tile should be visible", () => shouldBe(page, el("Aca monomer tile"), "visible"));
      await session.step(29, "And G monomer tile should be hidden", () => shouldBe(page, el("G monomer tile"), "hidden"));
      await session.step(30, "When user clears palette search", () => clearField(page, el("palette search")));
      await session.step(31, "Then G monomer tile should be visible", () => shouldBe(page, el("G monomer tile"), "visible"));
      await session.step(32, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A tile arms its monomer and a canvas click places it", async () => {
      await session.step(35, "When user clicks on G monomer tile", () => clickOn(page, el("G monomer tile")));
      await session.step(36, "Then editor status should contain text \"Next add: G\"", () => shouldContainText(page, el("editor status"), "Next add: G"));
      await session.step(37, "And there should be 2 visible drawn monomers", () => visibleCount(page, 2, el("drawn monomers")));
      await session.step(38, "When user clicks on an empty spot of editor canvas", () => clickEmptyCanvas(page));
      await session.step(39, "Then there should be 3 visible drawn monomers", () => visibleCount(page, 3, el("drawn monomers")));
      await session.step(40, "When user clicks on HELM tab", () => clickOn(page, el("HELM tab")));
      await session.step(41, "Then notation pane should have text \"PEPTIDE1{A.C}|PEPTIDE2{G}$$$$V2.0\"", () => shouldHaveText(page, el("notation pane"), "PEPTIDE1{A.C}|PEPTIDE2{G}$$$$V2.0"));
      await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The RNA tab shows the triplet builder", async () => {
      await session.step(46, "When user clicks on RNA palette tab", () => clickOn(page, el("RNA palette tab")));
      await session.step(47, "Then RNA builder should be visible", () => shouldBe(page, el("RNA builder"), "visible"));
      await session.step(48, "And there should be 5 visible triplets", () => visibleCount(page, 5, el("triplets")));
      await session.step(49, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Favorites is empty until a tile is starred", async () => {
      await session.step(52, "When user clicks on Favorites palette tab", () => clickOn(page, el("Favorites palette tab")));
      await session.step(53, "Then favorites empty note should contain text \"No favorites yet\"", () => shouldContainText(page, el("favorites empty note"), "No favorites yet"));
      await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Cancel throws the placed monomer away", async () => {
      await session.step(57, "When user clicks on CANCEL button in HELM editor", () => clickOn(page, el("CANCEL button in HELM editor")));
      await session.step(58, "Then HELM editor should be absent", () => shouldBe(page, el("HELM editor"), "absent"));
      await session.step(59, "And the value of \"HELM\" column in row 1 should be \"PEPTIDE1{A.C}$$$$\"", () => valueInRow(page, "HELM", 1, "PEPTIDE1{A.C}$$$$"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
