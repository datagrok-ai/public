/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/07-molecules/01-molecules-in-a-table.feature
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
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, expand, isExpanded, pressKeyIn, shouldBe, shouldContainText, shouldNotContainText, typeInto, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, currentRowIs, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelOpen, contextPanelClosed, contextPanelShowsCurrentCell, packageInstalled, simpleModeOff, sketcherIs, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Molecules", () => {
  const session = feature(test, "features/howto/07-molecules/01-molecules-in-a-table.feature", import.meta.url);
  test("View and edit molecules in a table", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(21, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(22, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(23, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(24, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(27, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(29, "When user expands Files---Demo---chem tree node inside browse tree", () => expand(page, el("Files---Demo---chem tree node inside browse tree")), undefined, "Files > Demo > chem holds tables of molecules");
    await session.step(30, "And user double-clicks on Files---Demo---chem---smiles.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---chem---smiles.csv tree node inside browse tree")));
    await session.step(31, "Then the \"smiles\" view should be current", () => viewIsCurrent(page, "smiles"));
    await session.step(32, "And the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(33, "And \"canonical_smiles\" column should have semantic type \"Molecule\"", () => columnSemType(page, "canonical_smiles", "Molecule"));
    await session.step(34, "When the video says \"SMILES text is recognized and drawn as structures\"", () => videoSays(page, "SMILES text is recognized and drawn as structures"));
    await session.step(36, "And user clicks on the \"cell 2 of canonical_smiles\" area of grid", () => clickArea(page, "cell 2 of canonical_smiles", el("grid")), undefined, "Click a molecule to make it current");
    await session.step(37, "Then row 2 should be current", () => currentRowIs(page, 2));
    await session.step(39, "When user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")), undefined, "The context panel shows what is known about it");
    await session.step(40, "Then the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
    await session.step(42, "When user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")), undefined, "Chemistry > Properties: weight, LogP, PSA, formula...");
    await session.step(43, "And user expands Properties accordion header in context panel", () => expand(page, el("Properties accordion header in context panel")));
    await session.step(44, "Then context panel should contain text \"C22H22O4\"", () => shouldContainText(page, el("context panel"), "C22H22O4"));
    await session.step(46, "When user clicks on the \"cell 5 of canonical_smiles\" area of grid", () => clickArea(page, "cell 5 of canonical_smiles", el("grid")), undefined, "Click another molecule: the panel follows");
    await session.step(47, "Then row 5 should be current", () => currentRowIs(page, 5));
    await session.step(48, "And context panel should not contain text \"C22H22O4\"", () => shouldNotContainText(page, el("context panel"), "C22H22O4"));
    await session.step(50, "When user double-clicks on the \"cell 1 of canonical_smiles\" area of grid", () => doubleClickArea(page, "cell 1 of canonical_smiles", el("grid")), undefined, "Double-click a structure to open the sketcher");
    await session.step(51, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(52, "When the video says \"Draw with the tools on the left, or type a SMILES above\" at the top", () => videoSaysAt(page, "Draw with the tools on the left, or type a SMILES above", "top"));
    await session.step(53, "And the video says \"Type or paste a SMILES, then press Enter\" at the top", () => videoSaysAt(page, "Type or paste a SMILES, then press Enter", "top"));
    await session.step(54, "And user clicks on molecule input of sketcher dialog", () => clickOn(page, el("molecule input of sketcher dialog")));
    await session.step(55, "And user types \"CC(=O)Oc1ccccc1C(=O)O\" into molecule input of sketcher dialog", () => typeInto(page, "CC(=O)Oc1ccccc1C(=O)O", el("molecule input of sketcher dialog")));
    await session.step(56, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
    await session.step(57, "And the video says \"OK puts the molecule in the cell\" at the top", () => videoSaysAt(page, "OK puts the molecule in the cell", "top"));
    await session.step(58, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
    await session.step(59, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
    await session.step(60, "And the value of \"canonical_smiles\" column in row 1 should be \"CC(=O)Oc1ccccc1C(=O)O\"", () => valueInRow(page, "canonical_smiles", 1, "CC(=O)Oc1ccccc1C(=O)O"));
    await session.step(61, "And context panel should contain text \"C9H8O4\"", () => shouldContainText(page, el("context panel"), "C9H8O4"));
    await session.step(62, "When the video says \"Row 1 is now aspirin, and the panel shows its properties\"", () => videoSays(page, "Row 1 is now aspirin, and the panel shows its properties"));
    await session.step(63, "Then the table should have 1000 rows", () => rowCount(page, 1000));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
  });
});
