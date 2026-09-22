/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/panels/info-panels.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cell.molecule, chem.panel.chemistry.rendering, chem.panel.chemistry.highlight, chem.panel.biology.toxicity, chem.panel.biology.drug-likeness, chem.panel.structure.identifiers, chem.panel.structure.2d-structure, chem.panel.structure.3d-structure]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {canvasColors} from '@datagrok-libraries/bdd/bindings/common/pixels';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {expand, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {cellIsCurrentObject} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, contextPanelOpen, contextPanelShows, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Chemistry, Biology and Structure panes of the Chem context panel", () => {
  const session = feature(test, "features/panels/info-panels.feature", import.meta.url);
  test("The Chemistry, Biology and Structure panes of the Chem context panel", {tag: ["@journey", "@realizes:chem.cell.molecule", "@realizes:chem.panel.chemistry.rendering", "@realizes:chem.panel.chemistry.highlight", "@realizes:chem.panel.biology.toxicity", "@realizes:chem.panel.biology.drug-likeness", "@realizes:chem.panel.structure.identifiers", "@realizes:chem.panel.structure.2d-structure", "@realizes:chem.panel.structure.3d-structure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
    await session.step(19, "And the context panel is open", () => contextPanelOpen(page));
    await run.scenario("The molecule column's Chemistry group offers Rendering and Highlight", async () => {
      await session.step(23, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(24, "And user clicks on the \"header canonical_smiles\" area of grid", () => clickArea(page, "header canonical_smiles", el("grid")));
      await session.step(25, "Then the context panel should show \"canonical_smiles\"", () => contextPanelShows(page, "canonical_smiles"));
      await session.step(26, "When user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
      await session.step(27, "Then \"Rendering\" pane in context panel should be visible", () => shouldBe(page, el("\"Rendering\" pane in context panel"), "visible"));
      await session.step(28, "And \"Highlight\" pane in context panel should be visible", () => shouldBe(page, el("\"Highlight\" pane in context panel"), "visible"));
      await session.step(29, "And \"Descriptors\" pane in context panel should be absent", () => shouldBe(page, el("\"Descriptors\" pane in context panel"), "absent"));
      await session.step(30, "When user expands Rendering accordion header in context panel", () => expand(page, el("Rendering accordion header in context panel")));
      await session.step(31, "Then \"Scaffold column\" choice input in \"Rendering\" pane in context panel should be visible", () => shouldBe(page, el("\"Scaffold column\" choice input in \"Rendering\" pane in context panel"), "visible"));
      await session.step(32, "And \"Highlight scaffold\" checkbox in \"Rendering\" pane in context panel should be visible", () => shouldBe(page, el("\"Highlight scaffold\" checkbox in \"Rendering\" pane in context panel"), "visible"));
      await session.step(33, "And \"Filter type\" choice input in \"Rendering\" pane in context panel should be visible", () => shouldBe(page, el("\"Filter type\" choice input in \"Rendering\" pane in context panel"), "visible"));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Toxicity grades the four risks and Drug Likeness reports a score", async () => {
      await session.step(38, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(39, "And user expands Biology accordion header in context panel", () => expand(page, el("Biology accordion header in context panel")));
      await session.step(40, "And user expands Toxicity accordion header in context panel", () => expand(page, el("Toxicity accordion header in context panel")));
      await session.step(41, "Then \"Toxicity\" pane in \"Biology\" pane in context panel should contain the text \"Mutagenicity\"", () => shouldContainText(page, el("\"Toxicity\" pane in \"Biology\" pane in context panel"), "Mutagenicity"));
      await session.step(42, "And \"Toxicity\" pane in \"Biology\" pane in context panel should contain the text \"Tumorigenicity\"", () => shouldContainText(page, el("\"Toxicity\" pane in \"Biology\" pane in context panel"), "Tumorigenicity"));
      await session.step(43, "And \"Toxicity\" pane in \"Biology\" pane in context panel should contain the text \"Irritating effects\"", () => shouldContainText(page, el("\"Toxicity\" pane in \"Biology\" pane in context panel"), "Irritating effects"));
      await session.step(44, "And \"Toxicity\" pane in \"Biology\" pane in context panel should contain the text \"Reproductive effects\"", () => shouldContainText(page, el("\"Toxicity\" pane in \"Biology\" pane in context panel"), "Reproductive effects"));
      await session.step(45, "And \"Toxicity\" pane in \"Biology\" pane in context panel should not contain the text \"Could not analyze toxicity\"", () => shouldNotContainText(page, el("\"Toxicity\" pane in \"Biology\" pane in context panel"), "Could not analyze toxicity"));
      await session.step(46, "When user expands \"Drug Likeness\" accordion header in context panel", () => expand(page, el("\"Drug Likeness\" accordion header in context panel")));
      await session.step(47, "Then \"Drug Likeness\" pane in \"Biology\" pane in context panel should contain the text \"Score:\"", () => shouldContainText(page, el("\"Drug Likeness\" pane in \"Biology\" pane in context panel"), "Score:"));
      await session.step(48, "And \"Drug Likeness\" pane in \"Biology\" pane in context panel should not contain the text \"Could not asses drug likeness\"", () => shouldNotContainText(page, el("\"Drug Likeness\" pane in \"Biology\" pane in context panel"), "Could not asses drug likeness"));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Identifiers pane holds the molecule's own Smiles, Inchi and Inchi key", async () => {
      await session.step(53, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(54, "And user expands Structure accordion header in context panel", () => expand(page, el("Structure accordion header in context panel")));
      await session.step(55, "And user expands Identifiers accordion header in context panel", () => expand(page, el("Identifiers accordion header in context panel")));
      await session.step(56, "Then \"Identifiers\" pane in context panel should contain the text \"Smiles\"", () => shouldContainText(page, el("\"Identifiers\" pane in context panel"), "Smiles"));
      await session.step(57, "And \"Identifiers\" pane in context panel should contain the text \"Inchi\"", () => shouldContainText(page, el("\"Identifiers\" pane in context panel"), "Inchi"));
      await session.step(58, "And \"Identifiers\" pane in context panel should contain the text \"Inchi key\"", () => shouldContainText(page, el("\"Identifiers\" pane in context panel"), "Inchi key"));
      await session.step(59, "And \"Identifiers\" pane in context panel should contain the text \"InChI=1S/\"", () => shouldContainText(page, el("\"Identifiers\" pane in context panel"), "InChI=1S/"));
      await session.step(60, "And \"Identifiers\" pane in context panel should not contain the text \"Malformed molecule\"", () => shouldNotContainText(page, el("\"Identifiers\" pane in context panel"), "Malformed molecule"));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The 2D Structure pane draws the molecule of the current cell", async () => {
      await session.step(65, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(66, "And user expands Structure accordion header in context panel", () => expand(page, el("Structure accordion header in context panel")));
      await session.step(67, "And user expands \"2D Structure\" accordion header in context panel", () => expand(page, el("\"2D Structure\" accordion header in context panel")));
      await session.step(68, "Then the canvases of \"2D Structure\" pane in context panel should be painted in at least 1 colors", () => canvasColors(page, el("\"2D Structure\" pane in context panel"), 1));
      await session.step(69, "And \"2D Structure\" pane in context panel should not contain the text \"Molecule is possibly malformed\"", () => shouldNotContainText(page, el("\"2D Structure\" pane in context panel"), "Molecule is possibly malformed"));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The 3D Structure pane builds a ball-and-stick view of the molecule", async () => {
      await session.step(74, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(75, "And user expands Structure accordion header in context panel", () => expand(page, el("Structure accordion header in context panel")));
      await session.step(76, "And user expands \"3D Structure\" accordion header in context panel", () => expand(page, el("\"3D Structure\" accordion header in context panel")));
      await session.step(77, "Then \"3D Structure\" pane in context panel should not contain the text \"Molecule has no atoms or malformed\"", () => shouldNotContainText(page, el("\"3D Structure\" pane in context panel"), "Molecule has no atoms or malformed"));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A V2000 molblock cell gets the same three groups", async () => {
      await session.step(82, "Given user opens mol1K.sdf dataset", () => openDataset(page, ds("mol1K.sdf")));
      await session.step(83, "Then \"molecule\" column should have semantic type \"Molecule\"", () => columnSemType(page, "molecule", "Molecule"));
      await session.step(84, "And \"molecule\" column should have units \"molblock\"", () => columnUnits(page, "molecule", "molblock"));
      await session.step(85, "Given the \"molecule\" cell of row 1 is the current object", () => cellIsCurrentObject(page, "molecule", 1));
      await session.step(86, "Then Chemistry accordion header in context panel should be visible", () => shouldBe(page, el("Chemistry accordion header in context panel"), "visible"));
      await session.step(87, "When user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
      await session.step(88, "And user expands Properties accordion header in context panel", () => expand(page, el("Properties accordion header in context panel")));
      await session.step(89, "Then \"Properties\" pane in context panel should contain the text \"MW\"", () => shouldContainText(page, el("\"Properties\" pane in context panel"), "MW"));
      await session.step(90, "And \"Properties\" pane in context panel should not contain the text \"Molecule is possibly malformed\"", () => shouldNotContainText(page, el("\"Properties\" pane in context panel"), "Molecule is possibly malformed"));
      await session.step(91, "When user expands Structure accordion header in context panel", () => expand(page, el("Structure accordion header in context panel")));
      await session.step(92, "And user expands \"2D Structure\" accordion header in context panel", () => expand(page, el("\"2D Structure\" accordion header in context panel")));
      await session.step(93, "Then the canvases of \"2D Structure\" pane in context panel should be painted in at least 1 colors", () => canvasColors(page, el("\"2D Structure\" pane in context panel"), 1));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A V3000 molblock cell gets the same three groups", async () => {
      await session.step(97, "Given user opens ApprovedDrugs2015 dataset", () => openDataset(page, ds("ApprovedDrugs2015")));
      await session.step(98, "Then \"molecule\" column should have semantic type \"Molecule\"", () => columnSemType(page, "molecule", "Molecule"));
      await session.step(99, "Given the \"molecule\" cell of row 1 is the current object", () => cellIsCurrentObject(page, "molecule", 1));
      await session.step(100, "When user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
      await session.step(101, "And user expands Properties accordion header in context panel", () => expand(page, el("Properties accordion header in context panel")));
      await session.step(102, "Then \"Properties\" pane in context panel should contain the text \"MW\"", () => shouldContainText(page, el("\"Properties\" pane in context panel"), "MW"));
      await session.step(103, "And \"Properties\" pane in context panel should not contain the text \"Molecule is possibly malformed\"", () => shouldNotContainText(page, el("\"Properties\" pane in context panel"), "Molecule is possibly malformed"));
      await session.step(104, "When user expands Structure accordion header in context panel", () => expand(page, el("Structure accordion header in context panel")));
      await session.step(105, "And user expands \"2D Structure\" accordion header in context panel", () => expand(page, el("\"2D Structure\" accordion header in context panel")));
      await session.step(106, "Then the canvases of \"2D Structure\" pane in context panel should be painted in at least 1 colors", () => canvasColors(page, el("\"2D Structure\" pane in context panel"), 1));
      await session.step(107, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A SMARTS cell gets its properties and its picture", async () => {
      await session.step(110, "Given user opens ex-smarts dataset", () => openDataset(page, ds("ex-smarts")));
      await session.step(111, "Then \"SMARTS\" column should have semantic type \"Molecule\"", () => columnSemType(page, "SMARTS", "Molecule"));
      await session.step(112, "Given the \"SMARTS\" cell of row 1 is the current object", () => cellIsCurrentObject(page, "SMARTS", 1));
      await session.step(113, "When user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
      await session.step(114, "And user expands Properties accordion header in context panel", () => expand(page, el("Properties accordion header in context panel")));
      await session.step(115, "Then \"Properties\" pane in context panel should contain the text \"MW\"", () => shouldContainText(page, el("\"Properties\" pane in context panel"), "MW"));
      await session.step(116, "And \"Properties\" pane in context panel should not contain the text \"Molecule is possibly malformed\"", () => shouldNotContainText(page, el("\"Properties\" pane in context panel"), "Molecule is possibly malformed"));
      await session.step(117, "When user expands Structure accordion header in context panel", () => expand(page, el("Structure accordion header in context panel")));
      await session.step(118, "And user expands \"2D Structure\" accordion header in context panel", () => expand(page, el("\"2D Structure\" accordion header in context panel")));
      await session.step(119, "Then the canvases of \"2D Structure\" pane in context panel should be painted in at least 1 colors", () => canvasColors(page, el("\"2D Structure\" pane in context panel"), 1));
    });
    run.finish();
  });
});
