/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sequence-translator/oligo-nucleotide-grid.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sequencetranslator.oligo-renderer.convert-helm-to-oligo, sequencetranslator.cell.oligo-nucleotide, sequencetranslator.panel.oligo-nucleotide, sequencetranslator.panel.oligo-structures]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, clipboardHas, clipboardImage, doubleClickOn, expand, isExpanded, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, columnSemType, columnsEqual, currentRowIs, hasColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelOpen, contextPanelOpen, dialogCloses, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, closeContextMenu, doubleClickArea, hasArea, infoBalloonText, menuLists, noBalloons, noErrors, pickFromAreaContextMenu, rightClickArea, wheelOverAreaTimes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {scrollGridTo} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("OligoNucleotide duplex column: conversion, context panels and cell actions", () => {
  const session = feature(test, "features/sequence-translator/oligo-nucleotide-grid.feature", import.meta.url);
  test("OligoNucleotide duplex column: conversion, context panels and cell actions", {tag: ["@journey", "@realizes:sequencetranslator.oligo-renderer.convert-helm-to-oligo", "@realizes:sequencetranslator.cell.oligo-nucleotide", "@realizes:sequencetranslator.panel.oligo-nucleotide", "@realizes:sequencetranslator.panel.oligo-structures"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(31, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(32, "And Files---App-Data---SequenceTranslator tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---SequenceTranslator tree node inside browse tree")));
    await session.step(33, "And Files---App-Data---SequenceTranslator---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---SequenceTranslator---samples tree node inside browse tree")));
    await session.step(34, "When user double-clicks Files---App-Data---SequenceTranslator---samples---sirna-demo.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---SequenceTranslator---samples---sirna-demo.csv tree node inside browse tree")));
    await session.step(35, "Then the \"sirna-demo\" view should be current", () => viewIsCurrent(page, "sirna-demo"));
    await session.step(36, "And the table should have 44 rows", () => rowCount(page, 44));
    await session.step(37, "And \"oligo_helm\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "oligo_helm", "Macromolecule"));
    await session.step(38, "When user scrolls the grid to the \"oligo_helm\" column", () => scrollGridTo(page, "oligo_helm"));
    await session.step(39, "And user picks \"Oligo > Convert HELM to Oligo\" from the context menu of the \"cell 1 of oligo_helm\" area of grid", () => pickFromAreaContextMenu(page, "Oligo > Convert HELM to Oligo", "cell 1 of oligo_helm", el("grid")));
    await session.step(40, "Then the table should have a column \"oligo_helm (oligo)\"", () => hasColumn(page, "oligo_helm (oligo)"));
    await run.scenario("Convert HELM to Oligo appends an OligoNucleotide column and leaves the HELM column as it was", async () => {
      await session.step(43, "Then \"oligo_helm (oligo)\" column should have semantic type \"OligoNucleotide\"", () => columnSemType(page, "oligo_helm (oligo)", "OligoNucleotide"));
      await session.step(44, "And \"oligo_helm (oligo)\" column should hold the same values as \"oligo_helm\" column", () => columnsEqual(page, "oligo_helm (oligo)", "oligo_helm"));
      await session.step(45, "And the value of \"oligo_helm\" column in row 1 should be \"RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$\"", () => valueInRow(page, "oligo_helm", 1, "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$"));
      await session.step(46, "And \"oligo_helm\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "oligo_helm", "Macromolecule"));
      await session.step(47, "And the table should have 13 columns", () => columnCount(page, 13));
      await session.step(48, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Oligo-Nucleotide pane summarises the current duplex and follows the current cell", async () => {
      await session.step(52, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(53, "When user scrolls the grid to the \"oligo_helm (oligo)\" column", () => scrollGridTo(page, "oligo_helm (oligo)"));
      await session.step(54, "And user clicks on the \"cell 1 of oligo_helm (oligo)\" area of grid", () => clickArea(page, "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(55, "And user expands \"Oligo-Nucleotide\" accordion header in context panel", () => expand(page, el("\"Oligo-Nucleotide\" accordion header in context panel")));
      await session.step(56, "Then \"Sense length\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"19 nt\"", () => shouldContainText(page, el("\"Sense length\" table row in \"Oligo-Nucleotide\" pane in context panel"), "19 nt"));
      await session.step(57, "And \"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"19 nt\"", () => shouldContainText(page, el("\"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel"), "19 nt"));
      await session.step(58, "And \"Modifications used\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"2'-OMe ×38, PS ×8\"", () => shouldContainText(page, el("\"Modifications used\" table row in \"Oligo-Nucleotide\" pane in context panel"), "2'-OMe ×38, PS ×8"));
      await session.step(59, "And \"Conjugates\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"GalNAc-L3 linker ×1\"", () => shouldContainText(page, el("\"Conjugates\" table row in \"Oligo-Nucleotide\" pane in context panel"), "GalNAc-L3 linker ×1"));
      await session.step(60, "And \"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"19 bp, blunt (auto-aligned)\"", () => shouldContainText(page, el("\"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel"), "19 bp, blunt (auto-aligned)"));
      await session.step(61, "And \"Oligo-Nucleotide\" pane in context panel should contain the text \"2'-O-Methyl\"", () => shouldContainText(page, el("\"Oligo-Nucleotide\" pane in context panel"), "2'-O-Methyl"));
      await session.step(62, "And \"Oligo-Nucleotide\" pane in context panel should contain the text \"Phosphorothioate (linkage)\"", () => shouldContainText(page, el("\"Oligo-Nucleotide\" pane in context panel"), "Phosphorothioate (linkage)"));
      await session.step(63, "When user clicks on the \"cell 2 of oligo_helm (oligo)\" area of grid", () => clickArea(page, "cell 2 of oligo_helm (oligo)", el("grid")));
      await session.step(64, "Then \"Modifications used\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"2'-OMe ×20, 2'-F ×18, PS ×8\"", () => shouldContainText(page, el("\"Modifications used\" table row in \"Oligo-Nucleotide\" pane in context panel"), "2'-OMe ×20, 2'-F ×18, PS ×8"));
      await session.step(65, "And \"Conjugates\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"Cholesterol ×1, GalNAc-L3 linker ×1\"", () => shouldContainText(page, el("\"Conjugates\" table row in \"Oligo-Nucleotide\" pane in context panel"), "Cholesterol ×1, GalNAc-L3 linker ×1"));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A single-strand ASO has no Duplex line; the Duplex line reports overhangs and HELM base pairs", async () => {
      await session.step(69, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(70, "When user scrolls the grid to the \"oligo_helm (oligo)\" column", () => scrollGridTo(page, "oligo_helm (oligo)"));
      await session.step(71, "And user clicks on the \"cell 1 of oligo_helm (oligo)\" area of grid", () => clickArea(page, "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(72, "And user expands \"Oligo-Nucleotide\" accordion header in context panel", () => expand(page, el("\"Oligo-Nucleotide\" accordion header in context panel")));
      await session.step(73, "Then \"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel should be visible", () => shouldBe(page, el("\"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel"), "visible"));
      await session.step(74, "When user scrolls the mouse wheel down 40 times over the \"cell 1 of oligo_helm (oligo)\" area of grid", () => wheelOverAreaTimes(page, "down", 40, "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(75, "And user clicks on the \"cell 38 of oligo_helm (oligo)\" area of grid", () => clickArea(page, "cell 38 of oligo_helm (oligo)", el("grid")));
      await session.step(76, "Then row 38 should be current", () => currentRowIs(page, 38));
      await session.step(77, "And \"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"19 bp, overhangs: 3' antisense +2, 3' sense +2 (auto-aligned)\"", () => shouldContainText(page, el("\"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel"), "19 bp, overhangs: 3' antisense +2, 3' sense +2 (auto-aligned)"));
      await session.step(78, "When user clicks on the \"cell 40 of oligo_helm (oligo)\" area of grid", () => clickArea(page, "cell 40 of oligo_helm (oligo)", el("grid")));
      await session.step(79, "Then row 40 should be current", () => currentRowIs(page, 40));
      await session.step(80, "And \"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"19 bp, blunt (from HELM pairs)\"", () => shouldContainText(page, el("\"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel"), "19 bp, blunt (from HELM pairs)"));
      await session.step(81, "When user scrolls the mouse wheel up 1 times over the \"cell 40 of oligo_helm (oligo)\" area of grid", () => wheelOverAreaTimes(page, "up", 1, "cell 40 of oligo_helm (oligo)", el("grid")));
      await session.step(82, "And user clicks on the \"cell 34 of oligo_helm (oligo)\" area of grid", () => clickArea(page, "cell 34 of oligo_helm (oligo)", el("grid")));
      await session.step(83, "Then row 34 should be current", () => currentRowIs(page, 34));
      await session.step(84, "And \"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"single-strand\"", () => shouldContainText(page, el("\"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel"), "single-strand"));
      await session.step(85, "And \"Modifications used\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"LNA ×6, PS ×18\"", () => shouldContainText(page, el("\"Modifications used\" table row in \"Oligo-Nucleotide\" pane in context panel"), "LNA ×6, PS ×18"));
      await session.step(86, "And \"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel should be absent", () => shouldBe(page, el("\"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel"), "absent"));
      await session.step(87, "When user scrolls the mouse wheel up 40 times over the \"cell 34 of oligo_helm (oligo)\" area of grid", () => wheelOverAreaTimes(page, "up", 40, "cell 34 of oligo_helm (oligo)", el("grid")));
      await session.step(88, "Then grid should have a \"cell 1 of oligo_helm (oligo)\" area", () => hasArea(page, el("grid"), "cell 1 of oligo_helm (oligo)"));
      await session.step(89, "When user clicks on the \"cell 1 of oligo_helm (oligo)\" area of grid", () => clickArea(page, "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(90, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(91, "And \"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"19 nt\"", () => shouldContainText(page, el("\"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel"), "19 nt"));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Oligo Structures pane builds the sense and antisense structures", async () => {
      await session.step(95, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(96, "When user scrolls the grid to the \"oligo_helm (oligo)\" column", () => scrollGridTo(page, "oligo_helm (oligo)"));
      await session.step(97, "And user clicks on the \"cell 1 of oligo_helm (oligo)\" area of grid", () => clickArea(page, "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(98, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(99, "And \"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"19 bp, blunt (auto-aligned)\"", () => shouldContainText(page, el("\"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel"), "19 bp, blunt (auto-aligned)"));
      await session.step(100, "When user expands \"Oligo Structures\" accordion header in context panel", () => expand(page, el("\"Oligo Structures\" accordion header in context panel")));
      await session.step(101, "Then \"Sense\" accordion header in \"Oligo Structures\" pane in context panel should be visible", () => shouldBe(page, el("\"Sense\" accordion header in \"Oligo Structures\" pane in context panel"), "visible"));
      await session.step(102, "And \"Antisense\" accordion header in \"Oligo Structures\" pane in context panel should be visible", () => shouldBe(page, el("\"Antisense\" accordion header in \"Oligo Structures\" pane in context panel"), "visible"));
      await session.step(103, "When user expands \"Sense\" accordion header in \"Oligo Structures\" pane in context panel", () => expand(page, el("\"Sense\" accordion header in \"Oligo Structures\" pane in context panel")));
      await session.step(104, "Then \"Explore\" accordion header in \"Sense\" pane in context panel should be visible", () => shouldBe(page, el("\"Explore\" accordion header in \"Sense\" pane in context panel"), "visible"));
      await session.step(105, "When user expands \"Antisense\" accordion header in \"Oligo Structures\" pane in context panel", () => expand(page, el("\"Antisense\" accordion header in \"Oligo Structures\" pane in context panel")));
      await session.step(106, "Then \"Explore\" accordion header in \"Antisense\" pane in context panel should be visible", () => shouldBe(page, el("\"Explore\" accordion header in \"Antisense\" pane in context panel"), "visible"));
      await session.step(107, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Copy as HELM copies the duplex HELM of the right-clicked cell", async () => {
      await session.step(111, "When user scrolls the grid to the \"oligo_helm (oligo)\" column", () => scrollGridTo(page, "oligo_helm (oligo)"));
      await session.step(112, "And user right-clicks on the \"cell 1 of oligo_helm (oligo)\" area of grid", () => rightClickArea(page, "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(113, "Then the open menu should list \"Current Value > Copy as HELM\"", () => menuLists(page, "Current Value > Copy as HELM"));
      await session.step(114, "And the open menu should list \"Current Value > Copy as Image\"", () => menuLists(page, "Current Value > Copy as Image"));
      await session.step(115, "And the open menu should list \"Current Value > Edit HELM\"", () => menuLists(page, "Current Value > Edit HELM"));
      await session.step(116, "And the open menu should list \"Enumerate Oligos\"", () => menuLists(page, "Enumerate Oligos"));
      await session.step(117, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(118, "And user picks \"Current Value > Copy as HELM\" from the context menu of the \"cell 1 of oligo_helm (oligo)\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Copy as HELM", "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(119, "Then an info balloon containing \"HELM copied to clipboard\" should have been shown", () => infoBalloonText(page, "HELM copied to clipboard"));
      await session.step(120, "And the clipboard should have the text \"RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$\"", () => clipboardHas(page, "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$"));
      await session.step(121, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Copy as Image puts a PNG picture of the cell on the clipboard", async () => {
      await session.step(124, "When user scrolls the grid to the \"oligo_helm (oligo)\" column", () => scrollGridTo(page, "oligo_helm (oligo)"));
      await session.step(125, "And user picks \"Current Value > Copy as Image\" from the context menu of the \"cell 1 of oligo_helm (oligo)\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Copy as Image", "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(126, "Then an info balloon containing \"Image copied to clipboard\" should have been shown", () => infoBalloonText(page, "Image copied to clipboard"));
      await session.step(127, "And the clipboard should hold a PNG image of at least 10000 bytes", () => clipboardImage(page, 10000));
      await session.step(128, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(129, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Edit HELM opens the HELM editor on the duplex, and CANCEL leaves the cell as it was", async () => {
      await session.step(132, "When user scrolls the grid to the \"oligo_helm (oligo)\" column", () => scrollGridTo(page, "oligo_helm (oligo)"));
      await session.step(133, "And user picks \"Current Value > Edit HELM\" from the context menu of the \"cell 1 of oligo_helm (oligo)\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Edit HELM", "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(134, "Then HELM notation tab should be visible", () => shouldBe(page, el("HELM notation tab"), "visible"));
      await session.step(135, "When user clicks on HELM notation tab", () => clickOn(page, el("HELM notation tab")));
      await session.step(136, "Then HELM notation should contain the text \"RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p\"", () => shouldContainText(page, el("HELM notation"), "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p"));
      await session.step(137, "And HELM notation should contain the text \"|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p\"", () => shouldContainText(page, el("HELM notation"), "|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p"));
      await session.step(138, "When user clicks on CANCEL button", () => clickOn(page, el("CANCEL button")));
      await session.step(139, "Then HELM notation tab should be hidden", () => shouldBe(page, el("HELM notation tab"), "hidden"));
      await session.step(140, "And the value of \"oligo_helm (oligo)\" column in row 1 should be \"RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$\"", () => valueInRow(page, "oligo_helm (oligo)", 1, "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$"));
      await session.step(141, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A double click shows the duplex full screen in a dialog without buttons", async () => {
      await session.step(144, "When user scrolls the grid to the \"oligo_helm (oligo)\" column", () => scrollGridTo(page, "oligo_helm (oligo)"));
      await session.step(145, "And user double-clicks on the \"cell 1 of oligo_helm (oligo)\" area of grid", () => doubleClickArea(page, "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(146, "Then \"Oligonucleotide\" dialog should be visible", () => shouldBe(page, el("\"Oligonucleotide\" dialog"), "visible"));
      await session.step(147, "And CANCEL button in \"Oligonucleotide\" dialog should be hidden", () => shouldBe(page, el("CANCEL button in \"Oligonucleotide\" dialog"), "hidden"));
      await session.step(148, "And OK button in \"Oligonucleotide\" dialog should be hidden", () => shouldBe(page, el("OK button in \"Oligonucleotide\" dialog"), "hidden"));
      await session.step(149, "When user clicks on Close icon in \"Oligonucleotide\" dialog", () => clickOn(page, el("Close icon in \"Oligonucleotide\" dialog")));
      await session.step(150, "Then the \"Oligonucleotide\" dialog should close", () => dialogCloses(page, "Oligonucleotide"));
      await session.step(151, "And the \"sirna-demo\" view should be current", () => viewIsCurrent(page, "sirna-demo"));
      await session.step(152, "And the value of \"oligo_helm (oligo)\" column in row 1 should be \"RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$\"", () => valueInRow(page, "oligo_helm (oligo)", 1, "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$"));
      await session.step(153, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(154, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Enumerate Oligos opens the HELM enumeration on the cell, and CANCEL adds nothing", async () => {
      await session.step(157, "When user scrolls the grid to the \"oligo_helm (oligo)\" column", () => scrollGridTo(page, "oligo_helm (oligo)"));
      await session.step(158, "And user picks \"Enumerate Oligos\" from the context menu of the \"cell 1 of oligo_helm (oligo)\" area of grid", () => pickFromAreaContextMenu(page, "Enumerate Oligos", "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(159, "Then \"PolyTool Helm Enumeration\" dialog should be visible", () => shouldBe(page, el("\"PolyTool Helm Enumeration\" dialog"), "visible"));
      await session.step(160, "And \"PolyTool Helm Enumeration\" dialog should contain the text \"m1G2sp3m4A5sp6m7C8p9m10U11p12\"", () => shouldContainText(page, el("\"PolyTool Helm Enumeration\" dialog"), "m1G2sp3m4A5sp6m7C8p9m10U11p12"));
      await session.step(161, "And \"PolyTool Helm Enumeration\" dialog should contain the text \"L3\"", () => shouldContainText(page, el("\"PolyTool Helm Enumeration\" dialog"), "L3"));
      await session.step(162, "When user clicks on CANCEL button in \"PolyTool Helm Enumeration\" dialog", () => clickOn(page, el("CANCEL button in \"PolyTool Helm Enumeration\" dialog")));
      await session.step(163, "Then the \"PolyTool Helm Enumeration\" dialog should close", () => dialogCloses(page, "PolyTool Helm Enumeration"));
      await session.step(164, "And the table should have 13 columns", () => columnCount(page, 13));
      await session.step(165, "And the table should have 44 rows", () => rowCount(page, 44));
      await session.step(166, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(167, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Combine Sense+Antisense to Oligo pairs the sense and antisense columns into a duplex column", async () => {
      await session.step(170, "When user scrolls the grid to the \"sense_helm\" column", () => scrollGridTo(page, "sense_helm"));
      await session.step(171, "And user picks \"Oligo > Combine Sense+Antisense to Oligo...\" from the context menu of the \"cell 1 of sense_helm\" area of grid", () => pickFromAreaContextMenu(page, "Oligo > Combine Sense+Antisense to Oligo...", "cell 1 of sense_helm", el("grid")));
      await session.step(172, "Then \"Combine Sense + Antisense to Oligonucleotide\" dialog should be visible", () => shouldBe(page, el("\"Combine Sense + Antisense to Oligonucleotide\" dialog"), "visible"));
      await session.step(173, "And Sense input in \"Combine Sense + Antisense to Oligonucleotide\" dialog should contain the text \"sense_helm\"", () => shouldContainText(page, el("Sense input in \"Combine Sense + Antisense to Oligonucleotide\" dialog"), "sense_helm"));
      await session.step(174, "When user selects \"antisense_helm\" in Antisense input in \"Combine Sense + Antisense to Oligonucleotide\" dialog", () => selectIn(page, "antisense_helm", el("Antisense input in \"Combine Sense + Antisense to Oligonucleotide\" dialog")));
      await session.step(175, "And user clicks on OK button in \"Combine Sense + Antisense to Oligonucleotide\" dialog", () => clickOn(page, el("OK button in \"Combine Sense + Antisense to Oligonucleotide\" dialog")));
      await session.step(176, "Then the \"Combine Sense + Antisense to Oligonucleotide\" dialog should close", () => dialogCloses(page, "Combine Sense + Antisense to Oligonucleotide"));
      await session.step(177, "And the table should have a column \"sense_helm+antisense_helm (oligo)\"", () => hasColumn(page, "sense_helm+antisense_helm (oligo)"));
      await session.step(178, "And \"sense_helm+antisense_helm (oligo)\" column should have semantic type \"OligoNucleotide\"", () => columnSemType(page, "sense_helm+antisense_helm (oligo)", "OligoNucleotide"));
      await session.step(179, "And the table should have 14 columns", () => columnCount(page, 14));
      await session.step(180, "And the value of \"sense_helm+antisense_helm (oligo)\" column in row 1 should be \"RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$\"", () => valueInRow(page, "sense_helm+antisense_helm (oligo)", 1, "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$"));
      await session.step(181, "And the value of \"sense_helm+antisense_helm (oligo)\" column in row 34 should be \"RNA1{[lna](C)[sp].[lna](A)[sp].[lna](G)[sp].d(T)[sp].d(G)[sp].d(T)[sp].d(T)[sp].d(C)[sp].d(T)[sp].d(T)[sp].d(G)[sp].d(C)[sp].d(T)[sp].d(C)[sp].d(T)[sp].d(A)[sp].[lna](T)[sp].[lna](A)[sp].[lna](A)}$$$$\"", () => valueInRow(page, "sense_helm+antisense_helm (oligo)", 34, "RNA1{[lna](C)[sp].[lna](A)[sp].[lna](G)[sp].d(T)[sp].d(G)[sp].d(T)[sp].d(T)[sp].d(C)[sp].d(T)[sp].d(T)[sp].d(G)[sp].d(C)[sp].d(T)[sp].d(C)[sp].d(T)[sp].d(A)[sp].[lna](T)[sp].[lna](A)[sp].[lna](A)}$$$$"));
      await session.step(182, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(183, "When user scrolls the grid to the \"sense_helm+antisense_helm (oligo)\" column", () => scrollGridTo(page, "sense_helm+antisense_helm (oligo)"));
      await session.step(184, "And user clicks on the \"cell 1 of sense_helm+antisense_helm (oligo)\" area of grid", () => clickArea(page, "cell 1 of sense_helm+antisense_helm (oligo)", el("grid")));
      await session.step(185, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(186, "And \"Sense length\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"19 nt\"", () => shouldContainText(page, el("\"Sense length\" table row in \"Oligo-Nucleotide\" pane in context panel"), "19 nt"));
      await session.step(187, "And \"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"19 nt\"", () => shouldContainText(page, el("\"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel"), "19 nt"));
      await session.step(188, "When user scrolls the mouse wheel down 40 times over the \"cell 1 of sense_helm+antisense_helm (oligo)\" area of grid", () => wheelOverAreaTimes(page, "down", 40, "cell 1 of sense_helm+antisense_helm (oligo)", el("grid")));
      await session.step(189, "And user scrolls the mouse wheel up 1 times over the \"cell 40 of sense_helm+antisense_helm (oligo)\" area of grid", () => wheelOverAreaTimes(page, "up", 1, "cell 40 of sense_helm+antisense_helm (oligo)", el("grid")));
      await session.step(190, "And user clicks on the \"cell 34 of sense_helm+antisense_helm (oligo)\" area of grid", () => clickArea(page, "cell 34 of sense_helm+antisense_helm (oligo)", el("grid")));
      await session.step(191, "Then row 34 should be current", () => currentRowIs(page, 34));
      await session.step(192, "And \"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"single-strand\"", () => shouldContainText(page, el("\"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel"), "single-strand"));
      await session.step(193, "And \"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel should be absent", () => shouldBe(page, el("\"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel"), "absent"));
      await session.step(194, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(195, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
