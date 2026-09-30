/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/biostructure-viewer/file-open-and-preview.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [biostructureviewer.preview.biostructure, biostructureviewer.import.pdb, biostructureviewer.import.xyz]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, close, doubleClickOn, isExpanded, pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsExactly, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, browsePanelOpen, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Opening and previewing structure files from the Files browser", () => {
  const session = feature(test, "features/biostructure-viewer/file-open-and-preview.feature", import.meta.url);
  test("PDB, mmCIF, CIF and PDBQT files preview with the Mol* engine under their own names (GROK-17654, GROK-18999)", {tag: ["@viewers", "@realizes:biostructureviewer.preview.biostructure", "@realizes:biostructureviewer.import.pdb", "@realizes:biostructureviewer.import.xyz"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(42, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(43, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(44, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(45, "And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node inside browse tree")));
    await session.step(46, "And Files---App-Data---BiostructureViewer---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer---samples tree node inside browse tree")));
    await session.step(49, "When user clicks on Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree", () => clickOn(page, el("Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree")));
    await session.step(50, "Then the \"1bdq.pdb\" view should be current", () => viewIsCurrent(page, "1bdq.pdb"));
    await session.step(51, "And \"Reset Camera\" button should be visible", () => shouldBe(page, el("\"Reset Camera\" button"), "visible"));
    await session.step(52, "When user clicks on Files---App-Data---BiostructureViewer---samples---1RQ9.mmcif tree node inside browse tree", () => clickOn(page, el("Files---App-Data---BiostructureViewer---samples---1RQ9.mmcif tree node inside browse tree")));
    await session.step(53, "Then the \"1RQ9.mmcif\" view should be current", () => viewIsCurrent(page, "1RQ9.mmcif"));
    await session.step(54, "And \"Reset Camera\" button should be visible", () => shouldBe(page, el("\"Reset Camera\" button"), "visible"));
    await session.step(55, "When user clicks on Files---App-Data---BiostructureViewer---samples---1rq9-assembly1.cif tree node inside browse tree", () => clickOn(page, el("Files---App-Data---BiostructureViewer---samples---1rq9-assembly1.cif tree node inside browse tree")));
    await session.step(56, "Then the \"1rq9-assembly1.cif\" view should be current", () => viewIsCurrent(page, "1rq9-assembly1.cif"));
    await session.step(57, "And \"Reset Camera\" button should be visible", () => shouldBe(page, el("\"Reset Camera\" button"), "visible"));
    await session.step(58, "When user clicks on Files---App-Data---BiostructureViewer---samples---1bdq.autodock-gpu.pdbqt tree node inside browse tree", () => clickOn(page, el("Files---App-Data---BiostructureViewer---samples---1bdq.autodock-gpu.pdbqt tree node inside browse tree")));
    await session.step(59, "Then the \"1bdq.autodock-gpu.pdbqt\" view should be current", () => viewIsCurrent(page, "1bdq.autodock-gpu.pdbqt"));
    await session.step(60, "And \"Reset Camera\" button should be visible", () => shouldBe(page, el("\"Reset Camera\" button"), "visible"));
    await session.step(61, "And no errors should have been logged", () => noErrors(page));
    await session.step(62, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A double-click on a PDB opens it in a view of its own with the Mol* engine (GROK-14442, GROK-16968)", {tag: ["@viewers", "@realizes:biostructureviewer.preview.biostructure", "@realizes:biostructureviewer.import.pdb", "@realizes:biostructureviewer.import.xyz"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(42, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(43, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(44, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(45, "And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node inside browse tree")));
    await session.step(46, "And Files---App-Data---BiostructureViewer---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer---samples tree node inside browse tree")));
    await session.step(65, "Given simple mode is off", () => simpleModeOff(page));
    await session.step(66, "When user clicks on Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree", () => clickOn(page, el("Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree")));
    await session.step(67, "Then the \"1bdq.pdb\" view should be current", () => viewIsCurrent(page, "1bdq.pdb"));
    await session.step(68, "And \"Reset Camera\" button should be visible", () => shouldBe(page, el("\"Reset Camera\" button"), "visible"));
    await session.step(69, "When user clicks on Files---App-Data---BiostructureViewer---samples---1rq9-assembly1.cif tree node inside browse tree", () => clickOn(page, el("Files---App-Data---BiostructureViewer---samples---1rq9-assembly1.cif tree node inside browse tree")));
    await session.step(70, "Then the \"1rq9-assembly1.cif\" view should be current", () => viewIsCurrent(page, "1rq9-assembly1.cif"));
    await session.step(71, "And \"1bdq.pdb\" view should be absent", () => shouldBe(page, el("\"1bdq.pdb\" view"), "absent"));
    await session.step(72, "When user double-clicks on Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree")));
    await session.step(73, "Then \"1bdq.pdb\" view should be visible", () => shouldBe(page, el("\"1bdq.pdb\" view"), "visible"));
    await session.step(74, "And the \"1bdq.pdb\" view should be current", () => viewIsCurrent(page, "1bdq.pdb"));
    await session.step(75, "And \"Reset Camera\" button should be visible", () => shouldBe(page, el("\"Reset Camera\" button"), "visible"));
    await session.step(76, "And grid should be absent", () => shouldBe(page, el("grid"), "absent"));
    await session.step(77, "And \"Open file\" dialog should be absent", () => shouldBe(page, el("\"Open file\" dialog"), "absent"));
    await session.step(78, "When user clicks on Files---App-Data---BiostructureViewer---samples---1rq9-assembly1.cif tree node inside browse tree", () => clickOn(page, el("Files---App-Data---BiostructureViewer---samples---1rq9-assembly1.cif tree node inside browse tree")));
    await session.step(79, "Then the \"1rq9-assembly1.cif\" view should be current", () => viewIsCurrent(page, "1rq9-assembly1.cif"));
    await session.step(80, "And \"1bdq.pdb\" view should be visible", () => shouldBe(page, el("\"1bdq.pdb\" view"), "visible"));
    await session.step(81, "When user clicks on \"1bdq.pdb\" view", () => clickOn(page, el("\"1bdq.pdb\" view")));
    await session.step(82, "Then the \"1bdq.pdb\" view should be current", () => viewIsCurrent(page, "1bdq.pdb"));
    await session.step(83, "And \"Reset Camera\" button should be visible", () => shouldBe(page, el("\"Reset Camera\" button"), "visible"));
    await session.step(84, "And no errors should have been logged", () => noErrors(page));
    await session.step(85, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A double-click on an mmCIF opens it in a view of its own", {tag: ["@viewers", "@realizes:biostructureviewer.preview.biostructure", "@realizes:biostructureviewer.import.pdb", "@realizes:biostructureviewer.import.xyz"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(42, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(43, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(44, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(45, "And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node inside browse tree")));
    await session.step(46, "And Files---App-Data---BiostructureViewer---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer---samples tree node inside browse tree")));
    await session.step(88, "Given simple mode is off", () => simpleModeOff(page));
    await session.step(89, "When user double-clicks on Files---App-Data---BiostructureViewer---samples---1RQ9.mmcif tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---BiostructureViewer---samples---1RQ9.mmcif tree node inside browse tree")));
    await session.step(90, "Then \"1RQ9.mmcif\" view should be visible", () => shouldBe(page, el("\"1RQ9.mmcif\" view"), "visible"));
    await session.step(91, "And the \"1RQ9.mmcif\" view should be current", () => viewIsCurrent(page, "1RQ9.mmcif"));
    await session.step(92, "And \"Reset Camera\" button should be visible", () => shouldBe(page, el("\"Reset Camera\" button"), "visible"));
    await session.step(93, "When user clicks on Files---App-Data---BiostructureViewer---samples---1rq9-assembly1.cif tree node inside browse tree", () => clickOn(page, el("Files---App-Data---BiostructureViewer---samples---1rq9-assembly1.cif tree node inside browse tree")));
    await session.step(94, "Then the \"1rq9-assembly1.cif\" view should be current", () => viewIsCurrent(page, "1rq9-assembly1.cif"));
    await session.step(95, "And \"1RQ9.mmcif\" view should be visible", () => shouldBe(page, el("\"1RQ9.mmcif\" view"), "visible"));
    await session.step(96, "And no errors should have been logged", () => noErrors(page));
    await session.step(97, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Closing an unrelated view while a structure preview is shown raises nothing (CLAUDE-33)", {tag: ["@viewers", "@realizes:biostructureviewer.preview.biostructure", "@realizes:biostructureviewer.import.pdb", "@realizes:biostructureviewer.import.xyz"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(42, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(43, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(44, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(45, "And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node inside browse tree")));
    await session.step(46, "And Files---App-Data---BiostructureViewer---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer---samples tree node inside browse tree")));
    await session.step(100, "Given simple mode is off", () => simpleModeOff(page));
    await session.step(101, "When user clicks on Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree", () => clickOn(page, el("Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree")));
    await session.step(102, "Then the \"1bdq.pdb\" view should be current", () => viewIsCurrent(page, "1bdq.pdb"));
    await session.step(103, "And \"Reset Camera\" button should be visible", () => shouldBe(page, el("\"Reset Camera\" button"), "visible"));
    await session.step(104, "Given Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(105, "When user double-clicks on Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(106, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(107, "When user closes demog view", () => close(page, el("demog view")));
    await session.step(108, "Then \"demog\" view should be absent", () => shouldBe(page, el("\"demog\" view"), "absent"));
    await session.step(109, "And no errors should have been logged", () => noErrors(page));
    await session.step(110, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(111, "Given the browse panel is open", () => browsePanelOpen(page));
    await session.step(112, "When user clicks on Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree", () => clickOn(page, el("Files---App-Data---BiostructureViewer---samples---1bdq.pdb tree node inside browse tree")));
    await session.step(113, "Then the \"1bdq.pdb\" view should be current", () => viewIsCurrent(page, "1bdq.pdb"));
    await session.step(114, "And \"Reset Camera\" button should be visible", () => shouldBe(page, el("\"Reset Camera\" button"), "visible"));
    await session.step(115, "And help panel should be absent", () => shouldBe(page, el("help panel"), "absent"));
    await session.step(116, "When user presses F1", () => pressKey(page, "F1"));
    await session.step(117, "Then help panel should be visible", () => shouldBe(page, el("help panel"), "visible"));
    await session.step(118, "When user presses F1", () => pressKey(page, "F1"));
    await session.step(119, "Then help panel should be absent", () => shouldBe(page, el("help panel"), "absent"));
    await session.step(120, "And the \"1bdq.pdb\" view should be current", () => viewIsCurrent(page, "1bdq.pdb"));
    await session.step(121, "And no errors should have been logged", () => noErrors(page));
    await session.step(122, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Open table residues from a PDB file's menu opens a residue table with an NGL viewer", {tag: ["@viewers", "@realizes:biostructureviewer.preview.biostructure", "@realizes:biostructureviewer.import.pdb", "@realizes:biostructureviewer.import.xyz"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(42, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(43, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(44, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(45, "And Files---App-Data---BiostructureViewer tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer tree node inside browse tree")));
    await session.step(46, "And Files---App-Data---BiostructureViewer---samples tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---BiostructureViewer---samples tree node inside browse tree")));
    await session.step(125, "When user clicks on Files---App-Data---BiostructureViewer---samples tree node inside browse tree", () => clickOn(page, el("Files---App-Data---BiostructureViewer---samples tree node inside browse tree")));
    await session.step(126, "Then 1bdq.pdb link in gallery should be visible", () => shouldBe(page, el("1bdq.pdb link in gallery"), "visible"));
    await session.step(127, "When user picks \"Open table residues\" from the context menu of 1bdq.pdb link in gallery", () => pickFromContextMenu(page, "Open table residues", el("1bdq.pdb link in gallery")));
    await session.step(128, "Then the \"Table\" table view should open with 198 rows", () => tableViewOpened(page, "Table", 198));
    await session.step(129, "And the table should have the columns \"code, compId, seqId, label, seq, frame\"", () => columnsExactly(page, "code, compId, seqId, label, seq, frame"));
    await session.step(130, "And the value of \"code\" column in row 1 should be \"P\"", () => valueInRow(page, "code", 1, "P"));
    await session.step(131, "And the value of \"seqId\" column in row 1 should be \"1\"", () => valueInRow(page, "seqId", 1, "1"));
    await session.step(132, "And NGL viewer should be visible", () => shouldBe(page, el("NGL viewer"), "visible"));
    await session.step(133, "And no errors should have been logged", () => noErrors(page));
    await session.step(134, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
