/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demos/19b-library-chemistry.feature
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
import {clickHelmMonomer} from '../../bindings/demo-polytool.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, pasteInto, selectIn, shouldBe, shouldContainText, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnsExactly, everyValueBetween, maxInRow, setGridColumnFormat} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {newColumnNamed, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {balloonsClosed, contextPanelClosed, openDataset, packageInstalled, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridColumnWidth, gridRowHeight} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("5 · Design a library · Chemistry of the library", () => {
  const session = feature(test, "features/demos/19b-library-chemistry.feature", import.meta.url);
  test("Rank the library by chemistry", {tag: ["@guide", "@demo"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And simple mode is off", () => simpleModeOff(page));
    await session.step(19, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(20, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(21, "And user opens HELM_sample dataset", () => openDataset(page, ds("HELM_sample")));
    await session.step(22, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(23, "And user clicks on the \"cell 1 of HELM\" area of grid", () => clickArea(page, "cell 1 of HELM", el("grid")));
    await session.step(24, "And user picks \"Bio > PolyTool > Enumerate HELM...\" from the top menu", () => pickFromTopMenu(page, "Bio > PolyTool > Enumerate HELM..."));
    await session.step(25, "And user clicks on monomer 4 of the HELM in the \"PolyTool Helm Enumeration\" dialog", () => clickHelmMonomer(page, 4, "PolyTool Helm Enumeration"));
    await session.step(26, "And user double-clicks on the \"cell 1 of Monomers\" area of grid in Placeholders input in \"PolyTool Helm Enumeration\" dialog", () => doubleClickArea(page, "cell 1 of Monomers", el("grid in Placeholders input in \"PolyTool Helm Enumeration\" dialog")));
    await session.step(27, "And the notification balloons are closed", () => balloonsClosed(page));
    await session.step(28, "And user pastes \"D, Q, S\" into Search input in \"Select Monomers\" dialog", () => pasteInto(page, "D, Q, S", el("Search input in \"Select Monomers\" dialog")));
    await session.step(29, "And user clicks on OK button in \"Select Monomers\" dialog", () => clickOn(page, el("OK button in \"Select Monomers\" dialog")));
    await session.step(30, "And user clicks on monomer 12 of the HELM in the \"PolyTool Helm Enumeration\" dialog", () => clickHelmMonomer(page, 12, "PolyTool Helm Enumeration"));
    await session.step(31, "And user double-clicks on the \"cell 2 of Monomers\" area of grid in Placeholders input in \"PolyTool Helm Enumeration\" dialog", () => doubleClickArea(page, "cell 2 of Monomers", el("grid in Placeholders input in \"PolyTool Helm Enumeration\" dialog")));
    await session.step(32, "And the notification balloons are closed", () => balloonsClosed(page));
    await session.step(33, "And user pastes \"K, R\" into Search input in \"Select Monomers\" dialog", () => pasteInto(page, "K, R", el("Search input in \"Select Monomers\" dialog")));
    await session.step(34, "And user clicks on OK button in \"Select Monomers\" dialog", () => clickOn(page, el("OK button in \"Select Monomers\" dialog")));
    await session.step(35, "And user selects \"matrix\" in \"Enumerator type\" input in \"PolyTool Helm Enumeration\" dialog", () => selectIn(page, "matrix", el("\"Enumerator type\" input in \"PolyTool Helm Enumeration\" dialog")));
    await session.step(36, "Then \"PolyTool Helm Enumeration\" dialog should contain the text \"6 sequences will be generated\"", () => shouldContainText(page, el("\"PolyTool Helm Enumeration\" dialog"), "6 sequences will be generated"));
    await session.step(37, "Given user checks \"To atomic level\" checkbox in \"PolyTool Helm Enumeration\" dialog", () => check(page, el("\"To atomic level\" checkbox in \"PolyTool Helm Enumeration\" dialog")));
    await session.step(38, "And user clicks on OK button in \"PolyTool Helm Enumeration\" dialog", () => clickOn(page, el("OK button in \"PolyTool Helm Enumeration\" dialog")));
    await session.step(39, "Then the table should have 6 rows", () => rowCount(page, 6));
    await session.step(40, "And the table should have the columns \"Enumerated, molfile(Enumerated)\"", () => columnsExactly(page, "Enumerated, molfile(Enumerated)"));
    await session.step(41, "And \"molfile(Enumerated)\" column should have semantic type \"Molecule\"", () => columnSemType(page, "molfile(Enumerated)", "Molecule"));
    await session.step(42, "Given the notification balloons are closed", () => balloonsClosed(page));
    await session.step(43, "And the \"Enumerated\" column of the grid is 600 px wide", () => gridColumnWidth(page, "Enumerated", 600));
    await session.step(44, "And the \"molfile(Enumerated)\" column of the grid is 420 px wide", () => gridColumnWidth(page, "molfile(Enumerated)", 420));
    await session.step(45, "And the rows of the grid are 110 px tall", () => gridRowHeight(page, 110));
    await session.step(46, "When user picks \"Chem > Calculate > Chemical Properties...\" from the top menu", () => pickFromTopMenu(page, "Chem > Calculate > Chemical Properties..."));
    await session.step(47, "Then \"Chemical Properties\" dialog should be visible", () => shouldBe(page, el("\"Chemical Properties\" dialog"), "visible"));
    await session.step(48, "And Molecules input in \"Chemical Properties\" dialog should contain text \"molfile(Enumerated)\"", () => shouldContainText(page, el("Molecules input in \"Chemical Properties\" dialog"), "molfile(Enumerated)"));
    await session.step(49, "When user checks \"Chemical Properties\" calculator in \"Chemical Properties\" dialog", () => check(page, el("\"Chemical Properties\" calculator in \"Chemical Properties\" dialog")));
    await session.step(50, "Then MW input in \"Chemical Properties\" dialog should be checked", () => shouldBe(page, el("MW input in \"Chemical Properties\" dialog"), "checked"));
    await session.step(51, "When user checks PSA input in \"Chemical Properties\" dialog", () => check(page, el("PSA input in \"Chemical Properties\" dialog")));
    await session.step(52, "And user clicks on OK button in \"Chemical Properties\" dialog", () => clickOn(page, el("OK button in \"Chemical Properties\" dialog")));
    await session.step(53, "Then \"Chemical Properties\" dialog should be hidden", () => shouldBe(page, el("\"Chemical Properties\" dialog"), "hidden"));
    await session.step(54, "And a new column \"MW\" should have been added", () => newColumnNamed(page, "MW"));
    await session.step(55, "And a new column \"PSA\" should have been added", () => newColumnNamed(page, "PSA"));
    await session.step(56, "And every value of \"PSA\" column should lie between 100 and 2000", () => everyValueBetween(page, "PSA", 100, 2000));
    await session.step(57, "And the table should have 6 rows", () => rowCount(page, 6));
    await session.step(58, "Given user sets the format of \"MW\" column in the grid to \"0.0\"", () => setGridColumnFormat(page, "MW", "0.0"));
    await session.step(59, "And user sets the format of \"PSA\" column in the grid to \"0.0\"", () => setGridColumnFormat(page, "PSA", "0.0"));
    await session.step(60, "And the \"MW\" column of the grid is 120 px wide", () => gridColumnWidth(page, "MW", 120));
    await session.step(61, "And the \"PSA\" column of the grid is 120 px wide", () => gridColumnWidth(page, "PSA", 120));
    await session.step(62, "Then \"PSA\" column should have its maximum in row 4", () => maxInRow(page, "PSA", 4));
    await session.step(63, "When user double-clicks on the \"header PSA\" area of grid", () => doubleClickArea(page, "header PSA", el("grid")));
    await session.step(64, "And user hovers over the \"cell 4 of PSA\" area of grid", () => hoverArea(page, "cell 4 of PSA", el("grid")));
    await session.step(65, "And the video says \"Measured like small molecules, ranked by polar surface area\"", () => videoSays(page, "Measured like small molecules, ranked by polar surface area"));
    await session.step(66, "Then no errors should have been logged", () => noErrors(page));
  });
});
