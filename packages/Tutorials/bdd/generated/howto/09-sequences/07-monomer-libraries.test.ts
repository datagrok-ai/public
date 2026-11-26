/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/09-sequences/07-monomer-libraries.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {bioReady, monomerEditorReady} from '../../../bindings/howto-sequences.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKey, shouldBe, shouldContainText, shouldHaveValue, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, setGridColumnFormat} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {contextPanelClosed, openDataset, packageInstalled, simpleModeOff, sketcherIs, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Sequences", () => {
  const session = feature(test, "features/howto/09-sequences/07-monomer-libraries.feature", import.meta.url);
  test("Monomers and monomer libraries", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(18, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(19, "And user opens HELM_sample dataset", () => openDataset(page, ds("HELM_sample")));
    await session.step(20, "And the Bio package is initialized", () => bioReady(page));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(23, "And user sets the format of \"Activity\" column in the grid to \"0.00\"", () => setGridColumnFormat(page, "Activity", "0.00"));
    await session.step(24, "When the video says \"Monomers: the building blocks of every sequence\"", () => videoSays(page, "Monomers: the building blocks of every sequence"));
    await session.step(26, "And user picks \"Bio > Manage > Monomers\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Monomers"), undefined, "Bio > Manage > Monomers");
    await session.step(27, "Then the \"Manage Monomers\" view should be current", () => viewIsCurrent(page, "Manage Monomers"));
    await session.step(28, "And the monomer sketcher of the Manage Monomers view should be ready", () => monomerEditorReady(page));
    await session.step(29, "And \"R1\" column should have semantic type \"Molecule\"", () => columnSemType(page, "R1", "Molecule"));
    await session.step(30, "When the video says \"A library as a table: structure, symbol, name, attachment points\"", () => videoSays(page, "A library as a table: structure, symbol, name, attachment points"));
    await session.step(32, "And user clicks on the \"cell 4 of Name\" area of grid", () => clickArea(page, "cell 4 of Name", el("grid")), undefined, "Click one: its structure and details on the right");
    await session.step(33, "Then \"Name\" input in monomer form should have value \"Glutamic acid\"", () => shouldHaveValue(page, el("\"Name\" input in monomer form"), "Glutamic acid"));
    await session.step(34, "And \"R3\" R-group chip in monomer form should contain the text \"R3 OH\"", () => shouldContainText(page, el("\"R3\" R-group chip in monomer form"), "R3 OH"));
    await session.step(36, "When user clicks on \"R3\" R-group chip in monomer form", () => clickOn(page, el("\"R3\" R-group chip in monomer form")), undefined, "R1, R2, R3: where it bonds; the cap when it does not");
    await session.step(37, "Then \"OH\" menu item should be checked", () => shouldBe(page, el("\"OH\" menu item"), "checked"));
    await session.step(38, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(40, "And user picks \"Bio > Manage > Monomer Libraries\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Monomer Libraries"), undefined, "Bio > Manage > Monomer Libraries: the library files in use");
    await session.step(41, "Then the \"Manage Monomer Libraries\" view should be current", () => viewIsCurrent(page, "Manage Monomer Libraries"));
    await session.step(42, "And \"HELMCoreLibrary.json\" checkbox should be checked", () => shouldBe(page, el("\"HELMCoreLibrary.json\" checkbox"), "checked"));
    await session.step(43, "When the video says \"Library files load together; the ticked ones are in use\"", () => videoSays(page, "Library files load together; the ticked ones are in use"));
    await session.step(44, "And the video says \"A symbol two files define differently: you pick one on the right\"", () => videoSays(page, "A symbol two files define differently: you pick one on the right"));
    await session.step(45, "Then no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
  });
});
