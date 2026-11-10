/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/17-monomer-libraries.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '../../bindings/monomer-form.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {monomerSketcherReady} from '../../bindings/library-files.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, pressKey, selectIn, shouldBe, shouldContainText, shouldHaveValue, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelClosed, openDataset, packageInstalled, simpleModeOff, sketcherIs, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bioinformatics · Monomer libraries", () => {
  const session = feature(test, "features/demos/17-monomer-libraries.feature", import.meta.url);
  test("Browse the monomer libraries", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(19, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(20, "And user opens HELM_sample dataset", () => openDataset(page, ds("HELM_sample")));
    await session.step(21, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(22, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(23, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(24, "When the video says \"Monomers: the building blocks of every sequence\"", () => videoSays(page, "Monomers: the building blocks of every sequence"));
    await session.step(25, "And user picks \"Bio > Manage > Monomers\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Monomers"));
    await session.step(26, "Then the \"Manage Monomers\" view should be current", () => viewIsCurrent(page, "Manage Monomers"));
    await session.step(27, "And the monomer sketcher of the Manage Monomers view should be ready", () => monomerSketcherReady(page));
    await session.step(28, "And \"R1\" column should have semantic type \"Molecule\"", () => columnSemType(page, "R1", "Molecule"));
    await session.step(29, "Given the grid shows only the columns \"Monomer, Symbol, Name, R1, R2\"", () => gridShowsColumns(page, "Monomer, Symbol, Name, R1, R2"));
    await session.step(30, "And the \"R1\" column of the grid is 190 px wide", () => gridColumnWidth(page, "R1", 190));
    await session.step(31, "And the \"R2\" column of the grid is 190 px wide", () => gridColumnWidth(page, "R2", 190));
    await session.step(32, "When user clicks on the \"cell 4 of Name\" area of grid", () => clickArea(page, "cell 4 of Name", el("grid")));
    await session.step(33, "Then \"Name\" input in monomer form should have value \"Glutamic acid\"", () => shouldHaveValue(page, el("\"Name\" input in monomer form"), "Glutamic acid"));
    await session.step(34, "And \"R3\" R-group chip in monomer form should contain the text \"R3 OH\"", () => shouldContainText(page, el("\"R3\" R-group chip in monomer form"), "R3 OH"));
    await session.step(35, "When user clicks on \"R3\" R-group chip in monomer form", () => clickOn(page, el("\"R3\" R-group chip in monomer form")));
    await session.step(36, "Then \"OH\" menu item should be checked", () => shouldBe(page, el("\"OH\" menu item"), "checked"));
    await session.step(37, "When user hovers over \"Azide\" menu item", () => hoverOver(page, el("\"Azide\" menu item")));
    await session.step(38, "And user presses Escape", () => pressKey(page, "Escape"));
    await session.step(39, "And user selects \"polytool-lib\" in \"Library\" input in monomer form", () => selectIn(page, "polytool-lib", el("\"Library\" input in monomer form")));
    await session.step(40, "Then the table should have 12 rows", () => rowCount(page, 12));
    await session.step(41, "And \"R1\" column should have semantic type \"Molecule\"", () => columnSemType(page, "R1", "Molecule"));
    await session.step(42, "Given the grid shows only the columns \"Monomer, Symbol, Name, R1, R2\"", () => gridShowsColumns(page, "Monomer, Symbol, Name, R1, R2"));
    await session.step(43, "And the \"R1\" column of the grid is 190 px wide", () => gridColumnWidth(page, "R1", 190));
    await session.step(44, "And the \"R2\" column of the grid is 190 px wide", () => gridColumnWidth(page, "R2", 190));
    await session.step(45, "When user clicks on the \"cell 4 of Name\" area of grid", () => clickArea(page, "cell 4 of Name", el("grid")));
    await session.step(46, "Then \"Name\" input in monomer form should have value \"Azido Glycine\"", () => shouldHaveValue(page, el("\"Name\" input in monomer form"), "Azido Glycine"));
    await session.step(47, "When user picks \"Bio > Manage > Monomer Libraries\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Monomer Libraries"));
    await session.step(48, "Then the \"Manage Monomer Libraries\" view should be current", () => viewIsCurrent(page, "Manage Monomer Libraries"));
    await session.step(49, "And \"HELMCoreLibrary.json\" checkbox should be checked", () => shouldBe(page, el("\"HELMCoreLibrary.json\" checkbox"), "checked"));
    await session.step(50, "When the video says \"Library files load together; clashing symbols are resolved here\"", () => videoSays(page, "Library files load together; clashing symbols are resolved here"));
    await session.step(51, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(52, "And no errors should have been logged", () => noErrors(page));
  });
});
