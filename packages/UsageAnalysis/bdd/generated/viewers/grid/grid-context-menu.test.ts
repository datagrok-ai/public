/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-context-menu.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid, core#1200, core#1198, core#1195, core#1193]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, clipboardHas, enterInto, isExpanded, pressKey, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag, currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {currentRowValue, rowCount, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, clickAreaHolding, closeContextMenu, hasArea, hasNoArea, hoverArea, menuLists, noErrors, pickFromAreaContextMenu, propertyShouldBe, readingIs, rightClickArea, setProperty, showsRows, wheelOverAreaTimes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("The context menu acts on the right-clicked cell", () => {
  const session = feature(test, "features/viewers/grid/grid-context-menu.feature", import.meta.url);
  test("A right click below the current row makes the clicked row current and keeps the scroll", {tag: ["@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(22, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(23, "When user clicks on the \"cell 5 of AGE\" area of grid", () => clickArea(page, "cell 5 of AGE", el("grid")));
    await session.step(24, "Then row 5 should be current", () => currentRowIs(page, 5));
    await session.step(25, "When user scrolls the mouse wheel down 60 times over the \"cell 5 of AGE\" area of grid", () => wheelOverAreaTimes(page, "down", 60, "cell 5 of AGE", el("grid")));
    await session.step(26, "Then grid should have a \"cell 1000 of AGE\" area", () => hasArea(page, el("grid"), "cell 1000 of AGE"));
    await session.step(27, "And grid should not have a \"cell 5 of AGE\" area", () => hasNoArea(page, el("grid"), "cell 5 of AGE"));
    await session.step(28, "When user right-clicks on the \"cell 1000 of AGE\" area of grid", () => rightClickArea(page, "cell 1000 of AGE", el("grid")));
    await session.step(29, "Then the open menu should list \"Current Column\"", () => menuLists(page, "Current Column"));
    await session.step(30, "And row 1000 should be current", () => currentRowIs(page, 1000));
    await session.step(31, "And the \"current row\" reading of grid should be 1000", () => readingIs(page, "current row", el("grid"), 1000));
    await session.step(32, "And grid should have a \"cell 1000 of AGE\" area", () => hasArea(page, el("grid"), "cell 1000 of AGE"));
    await session.step(33, "And grid should not have a \"cell 5 of AGE\" area", () => hasNoArea(page, el("grid"), "cell 5 of AGE"));
    await session.step(34, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(35, "Then row 1000 should be current", () => currentRowIs(page, 1000));
    await session.step(36, "And grid should have a \"cell 1000 of AGE\" area", () => hasArea(page, el("grid"), "cell 1000 of AGE"));
    await session.step(37, "And grid should not have a \"cell 5 of AGE\" area", () => hasNoArea(page, el("grid"), "cell 5 of AGE"));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
  test("Copy as SMILES copies the right-clicked molecule, not the one that was current", {tag: ["@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "Given user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await session.step(42, "Then the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(43, "When user clicks on the \"cell 5 of canonical_smiles\" area of grid", () => clickArea(page, "cell 5 of canonical_smiles", el("grid")));
    await session.step(44, "Then row 5 should be current", () => currentRowIs(page, 5));
    await session.step(45, "And \"canonical_smiles\" of the current row should be \"FC(F)(F)c1ccc(OC2CCNCC2)cc1\"", () => currentRowValue(page, "canonical_smiles", "FC(F)(F)c1ccc(OC2CCNCC2)cc1"));
    await session.step(46, "When user scrolls the mouse wheel down 200 times over the \"cell 5 of canonical_smiles\" area of grid", () => wheelOverAreaTimes(page, "down", 200, "cell 5 of canonical_smiles", el("grid")));
    await session.step(47, "Then the \"row order\" reading of grid should contain \"1000\"", () => readingContains(page, "row order", el("grid"), "1000"));
    await session.step(48, "When user picks \"Current Value > Copy as SMILES\" from the context menu of the \"cell 1000 of canonical_smiles\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Copy as SMILES", "cell 1000 of canonical_smiles", el("grid")));
    await session.step(49, "Then the clipboard should have the text \"CC(C)CCOC(=O)c1ccccc1N\"", () => clipboardHas(page, "CC(C)CCOC(=O)c1ccccc1N"));
    await session.step(50, "And row 1000 should be current", () => currentRowIs(page, 1000));
    await session.step(51, "And \"canonical_smiles\" of the current row should be \"CC(C)CCOC(=O)c1ccccc1N\"", () => currentRowValue(page, "canonical_smiles", "CC(C)CCOC(=O)c1ccccc1N"));
    await session.step(52, "And the \"row order\" reading of grid should contain \"1000\"", () => readingContains(page, "row order", el("grid"), "1000"));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
  });
  test("Edit Helm... opens the right-clicked peptide, not the one that was current", {tag: ["@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(56, "Given user opens helm-peptides dataset", () => openDataset(page, ds("helm-peptides")));
    await session.step(57, "Then the table should have 540 rows", () => rowCount(page, 540));
    await session.step(58, "When user clicks on the \"cell 5 of HELM\" area of grid", () => clickArea(page, "cell 5 of HELM", el("grid")));
    await session.step(59, "Then row 5 should be current", () => currentRowIs(page, 5));
    await session.step(60, "When user scrolls the mouse wheel down 200 times over the \"cell 5 of HELM\" area of grid", () => wheelOverAreaTimes(page, "down", 200, "cell 5 of HELM", el("grid")));
    await session.step(61, "Then the \"row order\" reading of grid should contain \"540\"", () => readingContains(page, "row order", el("grid"), "540"));
    await session.step(62, "When user picks \"Current Value > Edit Helm...\" from the context menu of the \"cell 540 of HELM\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Edit Helm...", "cell 540 of HELM", el("grid")));
    await session.step(63, "Then HELM notation tab should be visible", () => shouldBe(page, el("HELM notation tab"), "visible"));
    await session.step(64, "When user clicks on HELM notation tab", () => clickOn(page, el("HELM notation tab")));
    await session.step(65, "Then HELM notation should contain the text \"PEPTIDE1{meI.hHis.Hcy.Q.T.W.Q.Phe_4NH2.D-Tyr_Et.Tyr_ab-dehydroMe.dV.E.N.N.meK}$$$$\"", () => shouldContainText(page, el("HELM notation"), "PEPTIDE1{meI.hHis.Hcy.Q.T.W.Q.Phe_4NH2.D-Tyr_Et.Tyr_ab-dehydroMe.dV.E.N.N.meK}$$$$"));
    await session.step(66, "And row 540 should be current", () => currentRowIs(page, 540));
    await session.step(67, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(68, "Then HELM notation tab should be hidden", () => shouldBe(page, el("HELM notation tab"), "hidden"));
    await session.step(69, "And no errors should have been logged", () => noErrors(page));
  });
  test("Row Height typed in the context panel as 0 is not taken [typed=0]", {tag: ["@viewers", "@realizes:viewers.grid", "@known-failure:78", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(72, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(73, "When user hovers over the \"header AGE\" area of grid", () => hoverArea(page, "header AGE", el("grid")));
      await session.step(74, "And user clicks on \"Edit properties (F4)\" icon in grid", () => clickOn(page, el("\"Edit properties (F4)\" icon in grid")));
      await session.step(75, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(76, "Given \"Rows\" category in context panel is expanded", () => isExpanded(page, el("\"Rows\" category in context panel")));
      await session.step(77, "When user enters \"0\" into \"Row Height\" property in context panel", () => enterInto(page, "0", el("\"Row Height\" property in context panel")));
      await session.step(78, "Then \"Row Height\" property of grid should be \"28\"", () => propertyShouldBe(page, "Row Height", el("grid"), "28"));
    }, 78);
  });
  test("Row Height typed in the context panel as 501 is not taken [typed=501]", {tag: ["@viewers", "@realizes:viewers.grid", "@known-failure:78", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(72, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(73, "When user hovers over the \"header AGE\" area of grid", () => hoverArea(page, "header AGE", el("grid")));
      await session.step(74, "And user clicks on \"Edit properties (F4)\" icon in grid", () => clickOn(page, el("\"Edit properties (F4)\" icon in grid")));
      await session.step(75, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(76, "Given \"Rows\" category in context panel is expanded", () => isExpanded(page, el("\"Rows\" category in context panel")));
      await session.step(77, "When user enters \"501\" into \"Row Height\" property in context panel", () => enterInto(page, "501", el("\"Row Height\" property in context panel")));
      await session.step(78, "Then \"Row Height\" property of grid should be \"28\"", () => propertyShouldBe(page, "Row Height", el("grid"), "28"));
    }, 78);
  });
  test("A malformed number format marks the format field invalid", {tag: ["@viewers", "@realizes:viewers.grid", "@known-failure:92", "@realizes:core#1198"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(88, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(89, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(90, "Then Format AGE dialog should be visible", () => shouldBe(page, el("Format AGE dialog"), "visible"));
      await session.step(91, "When user enters \"0.0.0\" into Custom input in Format AGE dialog", () => enterInto(page, "0.0.0", el("Custom input in Format AGE dialog")));
      await session.step(92, "Then Custom input in Format AGE dialog should be invalid", () => shouldBe(page, el("Custom input in Format AGE dialog"), "invalid"));
    }, 92);
  });
  test("A malformed number format is not taken and the old format stays", {tag: ["@viewers", "@realizes:viewers.grid", "@known-failure:101", "@realizes:core#1198"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(96, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(97, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(98, "Then Format AGE dialog should be visible", () => shouldBe(page, el("Format AGE dialog"), "visible"));
      await session.step(99, "When user enters \"0.0.0\" into Custom input in Format AGE dialog", () => enterInto(page, "0.0.0", el("Custom input in Format AGE dialog")));
      await session.step(100, "And user clicks on OK button in Format AGE dialog", () => clickOn(page, el("OK button in Format AGE dialog")));
      await session.step(101, "Then \"AGE\" column should have tag \"format\" equal to \"\"", () => columnTag(page, "AGE", "format", ""));
    }, 101);
  });
  test("A column group named with spaces only is not made", {tag: ["@viewers", "@realizes:viewers.grid", "@known-failure:114", "@realizes:core#1195"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(105, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(106, "When user clicks on the \"header AGE\" area of grid holding Control", () => clickAreaHolding(page, "header AGE", el("grid"), "Control"));
      await session.step(107, "And user clicks on the \"header SEX\" area of grid holding Control", () => clickAreaHolding(page, "header SEX", el("grid"), "Control"));
      await session.step(108, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(109, "And Actions accordion header in context panel is expanded", () => isExpanded(page, el("Actions accordion header in context panel")));
      await session.step(110, "When user clicks on \"Group columns...\" text in context panel", () => clickOn(page, el("\"Group columns...\" text in context panel")));
      await session.step(111, "Then Group input in dialog should be visible", () => shouldBe(page, el("Group input in dialog"), "visible"));
      await session.step(112, "When user types \"   \" into Group input in dialog", () => typeInto(page, "   ", el("Group input in dialog")));
      await session.step(113, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
      await session.step(114, "Then \"AGE\" column should have tag \"group\" equal to \"\"", () => columnTag(page, "AGE", "group", ""));
    }, 114);
  });
  test("A read-only grid keeps its rows on Shift+Delete", {tag: ["@viewers", "@realizes:viewers.grid", "@known-failure:123", "@realizes:core#1193"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(118, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(119, "When user sets \"Allow Edit\" property of grid to \"false\"", () => setProperty(page, "Allow Edit", el("grid"), "false"));
      await session.step(120, "And user clicks on the \"row header 5\" area of grid holding Control", () => clickAreaHolding(page, "row header 5", el("grid"), "Control"));
      await session.step(121, "Then 1 row should be selected", () => selectedRowCount(page, 1));
      await session.step(122, "When user presses Shift+Delete", () => pressKey(page, "Shift+Delete"));
      await session.step(123, "Then the table should have 1000 rows", () => rowCount(page, 1000));
    }, 123);
  });
});
