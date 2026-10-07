/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-pinned-molecules.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKeyIn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {everyValueContains, everyValueMatches, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, openTableOf, packageInstalled, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {detectTypes, readingIsMolecule, readingIsRowMoleculeOf, rowMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {clickArea, doubleClickArea, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridPins} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Molecule cells of a pinned column, edited in Crux", () => {
  const session = feature(test, "features/viewers/grid/grid-pinned-molecules.feature", import.meta.url);
  test("In a pinned SMILES column, a typed C1CCCCC1 and a drawn edit are written as SMILES, and the table keeps its rows", {tag: ["@crux", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(15, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "Given user opens a table \"pinned-smiles\" with:", () => openTableOf(page, "pinned-smiles", [["id","molecule"],["1","CCO"],["2","c1ccccc1"],["3","CC(=O)O"]]), [["id","molecule"],["1","CCO"],["2","c1ccccc1"],["3","CC(=O)O"]]);
    await session.step(25, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(26, "When user picks \"Pin > Pin Column\" from the context menu of the \"header molecule\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Column", "header molecule", el("grid")));
    await session.step(27, "Then the grid should pin the columns \"molecule\"", () => gridPins(page, "molecule"));
    await session.step(28, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(29, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"CCO\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "CCO"));
    await session.step(30, "When user types \"C1CCCCC1\" into molecule input of sketcher dialog", () => typeInto(page, "C1CCCCC1", el("molecule input of sketcher dialog")));
    await session.step(31, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(32, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"C1CCCCC1\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "C1CCCCC1"));
    await session.step(33, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(34, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
    await session.step(35, "And the value of \"molecule\" column in row 1 should be \"C1CCCCC1\"", () => valueInRow(page, "molecule", 1, "C1CCCCC1"));
    await session.step(36, "When user double-clicks on the \"cell 2 of molecule\" area of grid", () => doubleClickArea(page, "cell 2 of molecule", el("grid")));
    await session.step(37, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"c1ccccc1\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "c1ccccc1"));
    await session.step(38, "When user clicks on crux single bond tool", () => clickOn(page, el("crux single bond tool")));
    await session.step(39, "And user clicks on the \"atom 0\" area of crux sketcher widget", () => clickArea(page, "atom 0", el("crux sketcher widget")));
    await session.step(40, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"Cc1ccccc1\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "Cc1ccccc1"));
    await session.step(41, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(42, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
    await session.step(43, "And the molecule in row 2 of \"molecule\" column should be \"Cc1ccccc1\"", () => rowMolecule(page, 2, "molecule", "Cc1ccccc1"));
    await session.step(44, "And every value of \"molecule\" column should match \"^\\S+$\"", () => everyValueMatches(page, "molecule", "^\\S+$"));
    await session.step(45, "And the table should have 3 rows", () => rowCount(page, 3));
    await session.step(46, "And the grid should pin the columns \"molecule\"", () => gridPins(page, "molecule"));
  });
  test("In a pinned molblock column, a typed C1CCCCC1 and a drawn molecule are written as molblocks, and the table keeps its rows", {tag: ["@crux", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(15, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(16, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(50, "Given user opens spgi-100 dataset", () => openDataset(page, ds("spgi-100")));
    await session.step(51, "When user picks \"Pin > Pin Column\" from the context menu of the \"header Structure\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Column", "header Structure", el("grid")));
    await session.step(52, "Then the grid should pin the columns \"Structure\"", () => gridPins(page, "Structure"));
    await session.step(53, "When user double-clicks on the \"cell 1 of Structure\" area of grid", () => doubleClickArea(page, "cell 1 of Structure", el("grid")));
    await session.step(54, "Then the \"smiles\" reading of crux sketcher widget should be the molecule in row 1 of \"Structure\" column", () => readingIsRowMoleculeOf(page, "smiles", el("crux sketcher widget"), 1, "Structure"));
    await session.step(55, "When user types \"C1CCCCC1\" into molecule input of sketcher dialog", () => typeInto(page, "C1CCCCC1", el("molecule input of sketcher dialog")));
    await session.step(56, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(57, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(58, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
    await session.step(59, "And the molecule in row 1 of \"Structure\" column should be \"C1CCCCC1\"", () => rowMolecule(page, 1, "Structure", "C1CCCCC1"));
    await session.step(60, "When user double-clicks on the \"cell 2 of Structure\" area of grid", () => doubleClickArea(page, "cell 2 of Structure", el("grid")));
    await session.step(61, "Then the \"smiles\" reading of crux sketcher widget should be the molecule in row 2 of \"Structure\" column", () => readingIsRowMoleculeOf(page, "smiles", el("crux sketcher widget"), 2, "Structure"));
    await session.step(62, "When user clicks on crux clear button", () => clickOn(page, el("crux clear button")));
    await session.step(63, "And user clicks on crux benzene tool", () => clickOn(page, el("crux benzene tool")));
    await session.step(64, "And user clicks on crux canvas", () => clickOn(page, el("crux canvas")));
    await session.step(65, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"c1ccccc1\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "c1ccccc1"));
    await session.step(66, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(67, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
    await session.step(68, "And the molecule in row 2 of \"Structure\" column should be \"c1ccccc1\"", () => rowMolecule(page, 2, "Structure", "c1ccccc1"));
    await session.step(69, "And every value of \"Structure\" column should contain \"M  END\"", () => everyValueContains(page, "Structure", "M  END"));
    await session.step(70, "And the table should have 100 rows", () => rowCount(page, 100));
    await session.step(71, "And the grid should pin the columns \"Structure\"", () => gridPins(page, "Structure"));
  });
});
