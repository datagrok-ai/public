/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sequence-translator/oligo-nucleotide-grid.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sequencetranslator.oligo-renderer.convert-helm-to-oligo, sequencetranslator.cell.oligo-nucleotide, sequencetranslator.panel.oligo-nucleotide, sequencetranslator.panel.oligo-structures]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
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
import {cellIsCurrentObject, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelOpen, contextPanelOpen, dialogCloses, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, closeContextMenu, doubleClickArea, hasArea, infoBalloonText, menuLists, noBalloons, noErrors, pickFromAreaContextMenu, rightClickArea, wheelOverAreaTimes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {scrollGridTo} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("OligoNucleotide duplex column: conversion, context panels and cell actions", () => {
  const session = feature(test, "features/sequence-translator/oligo-nucleotide-grid.feature", import.meta.url);
  test("OligoNucleotide duplex column: conversion, context panels and cell actions", {tag: ["@journey", "@realizes:sequencetranslator.oligo-renderer.convert-helm-to-oligo", "@realizes:sequencetranslator.cell.oligo-nucleotide", "@realizes:sequencetranslator.panel.oligo-nucleotide", "@realizes:sequencetranslator.panel.oligo-structures"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(34, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(35, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(36, "And Files---App-Data---SequenceTranslator tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---SequenceTranslator tree node inside browse tree")));
    await session.step(37, "And Files---App-Data---SequenceTranslator---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---SequenceTranslator---samples tree node inside browse tree")));
    await session.step(38, "When user double-clicks Files---App-Data---SequenceTranslator---samples---sirna-demo.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---SequenceTranslator---samples---sirna-demo.csv tree node inside browse tree")));
    await session.step(39, "Then the \"sirna-demo\" view should be current", () => viewIsCurrent(page, "sirna-demo"));
    await session.step(40, "And the table should have 44 rows", () => rowCount(page, 44));
    await session.step(41, "And \"oligo_helm\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "oligo_helm", "Macromolecule"));
    await session.step(42, "When user scrolls the grid to the \"oligo_helm\" column", () => scrollGridTo(page, "oligo_helm"));
    await session.step(43, "And user picks \"Oligo > Convert HELM to Oligo\" from the context menu of the \"cell 1 of oligo_helm\" area of grid", () => pickFromAreaContextMenu(page, "Oligo > Convert HELM to Oligo", "cell 1 of oligo_helm", el("grid")));
    await session.step(44, "Then the table should have a column \"oligo_helm (oligo)\"", () => hasColumn(page, "oligo_helm (oligo)"));
    await run.scenario("Convert HELM to Oligo appends an OligoNucleotide column and leaves the HELM column as it was", async () => {
      await session.step(47, "Then \"oligo_helm (oligo)\" column should have semantic type \"OligoNucleotide\"", () => columnSemType(page, "oligo_helm (oligo)", "OligoNucleotide"));
      await session.step(48, "And \"oligo_helm (oligo)\" column should hold the same values as \"oligo_helm\" column", () => columnsEqual(page, "oligo_helm (oligo)", "oligo_helm"));
      await session.step(49, "And the value of \"oligo_helm\" column in row 1 should be \"RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$\"", () => valueInRow(page, "oligo_helm", 1, "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$"));
      await session.step(50, "And \"oligo_helm\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "oligo_helm", "Macromolecule"));
      await session.step(51, "And the table should have 13 columns", () => columnCount(page, 13));
      await session.step(52, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Oligo-Nucleotide pane summarises the current duplex and follows the current cell", async () => {
      await session.step(56, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(57, "When user scrolls the grid to the \"oligo_helm (oligo)\" column", () => scrollGridTo(page, "oligo_helm (oligo)"));
      await session.step(58, "And user clicks on the \"cell 1 of oligo_helm (oligo)\" area of grid", () => clickArea(page, "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(59, "And user expands \"Oligo-Nucleotide\" accordion header in context panel", () => expand(page, el("\"Oligo-Nucleotide\" accordion header in context panel")));
      await session.step(60, "Then \"Sense length\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"19 nt\"", () => shouldContainText(page, el("\"Sense length\" table row in \"Oligo-Nucleotide\" pane in context panel"), "19 nt"));
      await session.step(61, "And \"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"19 nt\"", () => shouldContainText(page, el("\"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel"), "19 nt"));
      await session.step(62, "And \"Modifications used\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"2'-OMe ×38, PS ×8\"", () => shouldContainText(page, el("\"Modifications used\" table row in \"Oligo-Nucleotide\" pane in context panel"), "2'-OMe ×38, PS ×8"));
      await session.step(63, "And \"Conjugates\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"GalNAc-L3 linker ×1\"", () => shouldContainText(page, el("\"Conjugates\" table row in \"Oligo-Nucleotide\" pane in context panel"), "GalNAc-L3 linker ×1"));
      await session.step(64, "And \"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"19 bp, blunt (auto-aligned)\"", () => shouldContainText(page, el("\"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel"), "19 bp, blunt (auto-aligned)"));
      await session.step(65, "And \"Oligo-Nucleotide\" pane in context panel should contain the text \"2'-O-Methyl\"", () => shouldContainText(page, el("\"Oligo-Nucleotide\" pane in context panel"), "2'-O-Methyl"));
      await session.step(66, "And \"Oligo-Nucleotide\" pane in context panel should contain the text \"Phosphorothioate (linkage)\"", () => shouldContainText(page, el("\"Oligo-Nucleotide\" pane in context panel"), "Phosphorothioate (linkage)"));
      await session.step(67, "When user clicks on the \"cell 2 of oligo_helm (oligo)\" area of grid", () => clickArea(page, "cell 2 of oligo_helm (oligo)", el("grid")));
      await session.step(68, "Then \"Modifications used\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"2'-OMe ×20, 2'-F ×18, PS ×8\"", () => shouldContainText(page, el("\"Modifications used\" table row in \"Oligo-Nucleotide\" pane in context panel"), "2'-OMe ×20, 2'-F ×18, PS ×8"));
      await session.step(69, "And \"Conjugates\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"Cholesterol ×1, GalNAc-L3 linker ×1\"", () => shouldContainText(page, el("\"Conjugates\" table row in \"Oligo-Nucleotide\" pane in context panel"), "Cholesterol ×1, GalNAc-L3 linker ×1"));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A single-strand ASO has no Duplex line", async () => {
      await session.step(73, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(74, "When user scrolls the grid to the \"oligo_helm (oligo)\" column", () => scrollGridTo(page, "oligo_helm (oligo)"));
      await session.step(75, "And user clicks on the \"cell 1 of oligo_helm (oligo)\" area of grid", () => clickArea(page, "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(76, "And user expands \"Oligo-Nucleotide\" accordion header in context panel", () => expand(page, el("\"Oligo-Nucleotide\" accordion header in context panel")));
      await session.step(77, "Then \"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel should be visible", () => shouldBe(page, el("\"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel"), "visible"));
      await session.step(78, "Given the \"oligo_helm (oligo)\" cell of row 34 is the current object", () => cellIsCurrentObject(page, "oligo_helm (oligo)", 34));
      await session.step(79, "Then row 34 should be current", () => currentRowIs(page, 34));
      await session.step(80, "And \"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"single-strand\"", () => shouldContainText(page, el("\"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel"), "single-strand"));
      await session.step(81, "And \"Modifications used\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"LNA ×6, PS ×18\"", () => shouldContainText(page, el("\"Modifications used\" table row in \"Oligo-Nucleotide\" pane in context panel"), "LNA ×6, PS ×18"));
      await session.step(82, "And \"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel should be absent", () => shouldBe(page, el("\"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel"), "absent"));
      await session.step(83, "Given the \"oligo_helm (oligo)\" cell of row 1 is the current object", () => cellIsCurrentObject(page, "oligo_helm (oligo)", 1));
      await session.step(84, "Then grid should have a \"cell 1 of oligo_helm (oligo)\" area", () => hasArea(page, el("grid"), "cell 1 of oligo_helm (oligo)"));
      await session.step(85, "And row 1 should be current", () => currentRowIs(page, 1));
      await session.step(86, "And \"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"19 nt\"", () => shouldContainText(page, el("\"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel"), "19 nt"));
      await session.step(87, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Oligo Structures pane builds the sense and antisense structures", async () => {
      await session.step(90, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(91, "When user scrolls the grid to the \"oligo_helm (oligo)\" column", () => scrollGridTo(page, "oligo_helm (oligo)"));
      await session.step(92, "And user clicks on the \"cell 1 of oligo_helm (oligo)\" area of grid", () => clickArea(page, "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(93, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(94, "And \"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"19 bp, blunt (auto-aligned)\"", () => shouldContainText(page, el("\"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel"), "19 bp, blunt (auto-aligned)"));
      await session.step(95, "When user expands \"Oligo Structures\" accordion header in context panel", () => expand(page, el("\"Oligo Structures\" accordion header in context panel")));
      await session.step(96, "Then \"Sense\" accordion header in \"Oligo Structures\" pane in context panel should be visible", () => shouldBe(page, el("\"Sense\" accordion header in \"Oligo Structures\" pane in context panel"), "visible"));
      await session.step(97, "And \"Antisense\" accordion header in \"Oligo Structures\" pane in context panel should be visible", () => shouldBe(page, el("\"Antisense\" accordion header in \"Oligo Structures\" pane in context panel"), "visible"));
      await session.step(98, "When user expands \"Sense\" accordion header in \"Oligo Structures\" pane in context panel", () => expand(page, el("\"Sense\" accordion header in \"Oligo Structures\" pane in context panel")));
      await session.step(99, "Then \"Explore\" accordion header in \"Sense\" pane in context panel should be visible", () => shouldBe(page, el("\"Explore\" accordion header in \"Sense\" pane in context panel"), "visible"));
      await session.step(100, "When user expands \"Antisense\" accordion header in \"Oligo Structures\" pane in context panel", () => expand(page, el("\"Antisense\" accordion header in \"Oligo Structures\" pane in context panel")));
      await session.step(101, "Then \"Explore\" accordion header in \"Antisense\" pane in context panel should be visible", () => shouldBe(page, el("\"Explore\" accordion header in \"Antisense\" pane in context panel"), "visible"));
      await session.step(102, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(103, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Copy as HELM copies the duplex HELM of the right-clicked cell", async () => {
      await session.step(106, "When user scrolls the grid to the \"oligo_helm (oligo)\" column", () => scrollGridTo(page, "oligo_helm (oligo)"));
      await session.step(107, "And user right-clicks on the \"cell 1 of oligo_helm (oligo)\" area of grid", () => rightClickArea(page, "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(108, "Then the open menu should list \"Current Value > Copy as HELM\"", () => menuLists(page, "Current Value > Copy as HELM"));
      await session.step(109, "And the open menu should list \"Current Value > Copy as Image\"", () => menuLists(page, "Current Value > Copy as Image"));
      await session.step(110, "And the open menu should list \"Current Value > Edit HELM\"", () => menuLists(page, "Current Value > Edit HELM"));
      await session.step(111, "And the open menu should list \"Enumerate Oligos\"", () => menuLists(page, "Enumerate Oligos"));
      await session.step(112, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(113, "And user picks \"Current Value > Copy as HELM\" from the context menu of the \"cell 1 of oligo_helm (oligo)\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Copy as HELM", "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(114, "Then an info balloon containing \"HELM copied to clipboard\" should have been shown", () => infoBalloonText(page, "HELM copied to clipboard"));
      await session.step(115, "And the clipboard should have the text \"RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$\"", () => clipboardHas(page, "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$"));
      await session.step(116, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Copy as Image puts a PNG picture of the cell on the clipboard", async () => {
      await session.step(119, "When user scrolls the grid to the \"oligo_helm (oligo)\" column", () => scrollGridTo(page, "oligo_helm (oligo)"));
      await session.step(120, "And user picks \"Current Value > Copy as Image\" from the context menu of the \"cell 1 of oligo_helm (oligo)\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Copy as Image", "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(121, "Then an info balloon containing \"Image copied to clipboard\" should have been shown", () => infoBalloonText(page, "Image copied to clipboard"));
      await session.step(122, "And the clipboard should hold a PNG image of at least 10000 bytes", () => clipboardImage(page, 10000));
      await session.step(123, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(124, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Edit HELM opens the HELM editor on the duplex, and CANCEL leaves the cell as it was", async () => {
      await session.step(127, "When user scrolls the grid to the \"oligo_helm (oligo)\" column", () => scrollGridTo(page, "oligo_helm (oligo)"));
      await session.step(128, "And user picks \"Current Value > Edit HELM\" from the context menu of the \"cell 1 of oligo_helm (oligo)\" area of grid", () => pickFromAreaContextMenu(page, "Current Value > Edit HELM", "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(129, "Then HELM notation tab should be visible", () => shouldBe(page, el("HELM notation tab"), "visible"));
      await session.step(130, "When user clicks on HELM notation tab", () => clickOn(page, el("HELM notation tab")));
      await session.step(131, "Then HELM notation should contain the text \"RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p\"", () => shouldContainText(page, el("HELM notation"), "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p"));
      await session.step(132, "And HELM notation should contain the text \"|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p\"", () => shouldContainText(page, el("HELM notation"), "|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p"));
      await session.step(133, "When user clicks on CANCEL button", () => clickOn(page, el("CANCEL button")));
      await session.step(134, "Then HELM notation tab should be hidden", () => shouldBe(page, el("HELM notation tab"), "hidden"));
      await session.step(135, "And the value of \"oligo_helm (oligo)\" column in row 1 should be \"RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$\"", () => valueInRow(page, "oligo_helm (oligo)", 1, "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$"));
      await session.step(136, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A double click shows the duplex full screen in a dialog without buttons", async () => {
      await session.step(139, "When user scrolls the grid to the \"oligo_helm (oligo)\" column", () => scrollGridTo(page, "oligo_helm (oligo)"));
      await session.step(140, "And user double-clicks on the \"cell 1 of oligo_helm (oligo)\" area of grid", () => doubleClickArea(page, "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(141, "Then \"Oligonucleotide\" dialog should be visible", () => shouldBe(page, el("\"Oligonucleotide\" dialog"), "visible"));
      await session.step(142, "And CANCEL button in \"Oligonucleotide\" dialog should be hidden", () => shouldBe(page, el("CANCEL button in \"Oligonucleotide\" dialog"), "hidden"));
      await session.step(143, "And OK button in \"Oligonucleotide\" dialog should be hidden", () => shouldBe(page, el("OK button in \"Oligonucleotide\" dialog"), "hidden"));
      await session.step(144, "When user clicks on Close icon in \"Oligonucleotide\" dialog", () => clickOn(page, el("Close icon in \"Oligonucleotide\" dialog")));
      await session.step(145, "Then the \"Oligonucleotide\" dialog should close", () => dialogCloses(page, "Oligonucleotide"));
      await session.step(146, "And the \"sirna-demo\" view should be current", () => viewIsCurrent(page, "sirna-demo"));
      await session.step(147, "And the value of \"oligo_helm (oligo)\" column in row 1 should be \"RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$\"", () => valueInRow(page, "oligo_helm (oligo)", 1, "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$"));
      await session.step(148, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(149, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Enumerate Oligos opens the HELM enumeration on the cell, and CANCEL adds nothing", async () => {
      await session.step(152, "When user scrolls the grid to the \"oligo_helm (oligo)\" column", () => scrollGridTo(page, "oligo_helm (oligo)"));
      await session.step(153, "And user picks \"Enumerate Oligos\" from the context menu of the \"cell 1 of oligo_helm (oligo)\" area of grid", () => pickFromAreaContextMenu(page, "Enumerate Oligos", "cell 1 of oligo_helm (oligo)", el("grid")));
      await session.step(154, "Then \"PolyTool Helm Enumeration\" dialog should be visible", () => shouldBe(page, el("\"PolyTool Helm Enumeration\" dialog"), "visible"));
      await session.step(155, "And \"PolyTool Helm Enumeration\" dialog should contain the text \"m1G2sp3m4A5sp6m7C8p9m10U11p12\"", () => shouldContainText(page, el("\"PolyTool Helm Enumeration\" dialog"), "m1G2sp3m4A5sp6m7C8p9m10U11p12"));
      await session.step(156, "And \"PolyTool Helm Enumeration\" dialog should contain the text \"L3\"", () => shouldContainText(page, el("\"PolyTool Helm Enumeration\" dialog"), "L3"));
      await session.step(157, "When user clicks on CANCEL button in \"PolyTool Helm Enumeration\" dialog", () => clickOn(page, el("CANCEL button in \"PolyTool Helm Enumeration\" dialog")));
      await session.step(158, "Then the \"PolyTool Helm Enumeration\" dialog should close", () => dialogCloses(page, "PolyTool Helm Enumeration"));
      await session.step(159, "And the table should have 13 columns", () => columnCount(page, 13));
      await session.step(160, "And the table should have 44 rows", () => rowCount(page, 44));
      await session.step(161, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(162, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Combine Sense+Antisense to Oligo pairs the sense and antisense columns into a duplex column", async () => {
      await session.step(165, "When user scrolls the grid to the \"sense_helm\" column", () => scrollGridTo(page, "sense_helm"));
      await session.step(166, "And user picks \"Oligo > Combine Sense+Antisense to Oligo...\" from the context menu of the \"cell 1 of sense_helm\" area of grid", () => pickFromAreaContextMenu(page, "Oligo > Combine Sense+Antisense to Oligo...", "cell 1 of sense_helm", el("grid")));
      await session.step(167, "Then \"Combine Sense + Antisense to Oligonucleotide\" dialog should be visible", () => shouldBe(page, el("\"Combine Sense + Antisense to Oligonucleotide\" dialog"), "visible"));
      await session.step(168, "And Sense input in \"Combine Sense + Antisense to Oligonucleotide\" dialog should contain the text \"sense_helm\"", () => shouldContainText(page, el("Sense input in \"Combine Sense + Antisense to Oligonucleotide\" dialog"), "sense_helm"));
      await session.step(169, "When user selects \"antisense_helm\" in Antisense input in \"Combine Sense + Antisense to Oligonucleotide\" dialog", () => selectIn(page, "antisense_helm", el("Antisense input in \"Combine Sense + Antisense to Oligonucleotide\" dialog")));
      await session.step(170, "And user clicks on OK button in \"Combine Sense + Antisense to Oligonucleotide\" dialog", () => clickOn(page, el("OK button in \"Combine Sense + Antisense to Oligonucleotide\" dialog")));
      await session.step(171, "Then the \"Combine Sense + Antisense to Oligonucleotide\" dialog should close", () => dialogCloses(page, "Combine Sense + Antisense to Oligonucleotide"));
      await session.step(172, "And the table should have a column \"sense_helm+antisense_helm (oligo)\"", () => hasColumn(page, "sense_helm+antisense_helm (oligo)"));
      await session.step(173, "And \"sense_helm+antisense_helm (oligo)\" column should have semantic type \"OligoNucleotide\"", () => columnSemType(page, "sense_helm+antisense_helm (oligo)", "OligoNucleotide"));
      await session.step(174, "And the table should have 14 columns", () => columnCount(page, 14));
      await session.step(175, "And the value of \"sense_helm+antisense_helm (oligo)\" column in row 1 should be \"RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$\"", () => valueInRow(page, "sense_helm+antisense_helm (oligo)", 1, "RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$"));
      await session.step(176, "And the value of \"sense_helm+antisense_helm (oligo)\" column in row 34 should be \"RNA1{[lna](C)[sp].[lna](A)[sp].[lna](G)[sp].d(T)[sp].d(G)[sp].d(T)[sp].d(T)[sp].d(C)[sp].d(T)[sp].d(T)[sp].d(G)[sp].d(C)[sp].d(T)[sp].d(C)[sp].d(T)[sp].d(A)[sp].[lna](T)[sp].[lna](A)[sp].[lna](A)}$$$$\"", () => valueInRow(page, "sense_helm+antisense_helm (oligo)", 34, "RNA1{[lna](C)[sp].[lna](A)[sp].[lna](G)[sp].d(T)[sp].d(G)[sp].d(T)[sp].d(T)[sp].d(C)[sp].d(T)[sp].d(T)[sp].d(G)[sp].d(C)[sp].d(T)[sp].d(C)[sp].d(T)[sp].d(A)[sp].[lna](T)[sp].[lna](A)[sp].[lna](A)}$$$$"));
      await session.step(177, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(178, "When user scrolls the grid to the \"sense_helm+antisense_helm (oligo)\" column", () => scrollGridTo(page, "sense_helm+antisense_helm (oligo)"));
      await session.step(179, "And user clicks on the \"cell 1 of sense_helm+antisense_helm (oligo)\" area of grid", () => clickArea(page, "cell 1 of sense_helm+antisense_helm (oligo)", el("grid")));
      await session.step(180, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(181, "And \"Sense length\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"19 nt\"", () => shouldContainText(page, el("\"Sense length\" table row in \"Oligo-Nucleotide\" pane in context panel"), "19 nt"));
      await session.step(182, "And \"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"19 nt\"", () => shouldContainText(page, el("\"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel"), "19 nt"));
      await session.step(183, "When user scrolls the mouse wheel down 40 times over the \"cell 1 of sense_helm+antisense_helm (oligo)\" area of grid", () => wheelOverAreaTimes(page, "down", 40, "cell 1 of sense_helm+antisense_helm (oligo)", el("grid")));
      await session.step(184, "And user scrolls the mouse wheel up 1 times over the \"cell 40 of sense_helm+antisense_helm (oligo)\" area of grid", () => wheelOverAreaTimes(page, "up", 1, "cell 40 of sense_helm+antisense_helm (oligo)", el("grid")));
      await session.step(185, "And user clicks on the \"cell 34 of sense_helm+antisense_helm (oligo)\" area of grid", () => clickArea(page, "cell 34 of sense_helm+antisense_helm (oligo)", el("grid")));
      await session.step(186, "Then row 34 should be current", () => currentRowIs(page, 34));
      await session.step(187, "And \"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel should contain the text \"single-strand\"", () => shouldContainText(page, el("\"Antisense length\" table row in \"Oligo-Nucleotide\" pane in context panel"), "single-strand"));
      await session.step(188, "And \"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel should be absent", () => shouldBe(page, el("\"Duplex\" table row in \"Oligo-Nucleotide\" pane in context panel"), "absent"));
      await session.step(189, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(190, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
