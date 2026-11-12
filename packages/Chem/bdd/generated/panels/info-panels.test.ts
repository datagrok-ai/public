/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/panels/info-panels.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cell.molecule, chem.panel.chemistry.rendering, chem.panel.chemistry.highlight, chem.panel.biology.toxicity, chem.panel.biology.drug-likeness, chem.panel.structure.2d-structure, chem.panel.biology.structural-alerts, chem.panel.biology.pharmacophores, chem.panel.chemistry.properties]
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {engineMethod} from '../../bindings/engine.js';
import {canvasColors} from '@datagrok-libraries/bdd/bindings/common/pixels';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, enterInto, expand, hoverOver, selectIn, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, currentRowIs, hasColumn, hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, contextPanelOpen, contextPanelShows, contextPanelShowsCurrentCell, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaColor, areaNotColor, clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Chemistry, Biology and Structure panes of the Chem context panel", () => {
  const session = feature(test, "features/panels/info-panels.feature", import.meta.url);
  test("The Chemistry, Biology and Structure panes of the Chem context panel", {tag: ["@journey", "@realizes:chem.cell.molecule", "@realizes:chem.panel.chemistry.rendering", "@realizes:chem.panel.chemistry.highlight", "@realizes:chem.panel.biology.toxicity", "@realizes:chem.panel.biology.drug-likeness", "@realizes:chem.panel.structure.2d-structure", "@realizes:chem.panel.biology.structural-alerts", "@realizes:chem.panel.biology.pharmacophores", "@realizes:chem.panel.chemistry.properties"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(30, "And user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
    await session.step(31, "And the context panel is open", () => contextPanelOpen(page));
    await run.scenario("The molecule column's Chemistry group offers Rendering and Highlight", async () => {
      await session.step(35, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(36, "And user clicks on the \"header canonical_smiles\" area of grid", () => clickArea(page, "header canonical_smiles", el("grid")));
      await session.step(37, "Then the context panel should show \"canonical_smiles\"", () => contextPanelShows(page, "canonical_smiles"));
      await session.step(38, "When user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
      await session.step(39, "Then \"Rendering\" pane in context panel should be visible", () => shouldBe(page, el("\"Rendering\" pane in context panel"), "visible"));
      await session.step(40, "And \"Highlight\" pane in context panel should be visible", () => shouldBe(page, el("\"Highlight\" pane in context panel"), "visible"));
      await session.step(41, "And \"Descriptors\" pane in context panel should be absent", () => shouldBe(page, el("\"Descriptors\" pane in context panel"), "absent"));
      await session.step(42, "When user expands Rendering accordion header in context panel", () => expand(page, el("Rendering accordion header in context panel")));
      await session.step(43, "Then \"Scaffold column\" choice input in \"Rendering\" pane in context panel should be visible", () => shouldBe(page, el("\"Scaffold column\" choice input in \"Rendering\" pane in context panel"), "visible"));
      await session.step(44, "And \"Highlight scaffold\" checkbox in \"Rendering\" pane in context panel should be visible", () => shouldBe(page, el("\"Highlight scaffold\" checkbox in \"Rendering\" pane in context panel"), "visible"));
      await session.step(45, "And \"Filter type\" choice input in \"Rendering\" pane in context panel should be visible", () => shouldBe(page, el("\"Filter type\" choice input in \"Rendering\" pane in context panel"), "visible"));
      await session.step(46, "And \"Mark unspecified stereo\" checkbox in \"Rendering\" pane in context panel should be checked", () => shouldBe(page, el("\"Mark unspecified stereo\" checkbox in \"Rendering\" pane in context panel"), "checked"));
      await session.step(47, "And \"Simplified stereo label\" checkbox in \"Rendering\" pane in context panel should be checked", () => shouldBe(page, el("\"Simplified stereo label\" checkbox in \"Rendering\" pane in context panel"), "checked"));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Toxicity grades the four risks and Drug Likeness reports a score", async () => {
      await session.step(52, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(53, "And user expands Biology accordion header in context panel", () => expand(page, el("Biology accordion header in context panel")));
      await session.step(54, "And user expands Toxicity accordion header in context panel", () => expand(page, el("Toxicity accordion header in context panel")));
      await session.step(55, "Then \"Toxicity\" pane in \"Biology\" pane in context panel should show the Chemical Engine's method", () => engineMethod(page, "Toxicity", "Biology"));
      await session.step(56, "And \"Toxicity\" pane in \"Biology\" pane in context panel should not contain the text \"Could not analyze toxicity\"", () => shouldNotContainText(page, el("\"Toxicity\" pane in \"Biology\" pane in context panel"), "Could not analyze toxicity"));
      await session.step(57, "When user expands \"Drug Likeness\" accordion header in context panel", () => expand(page, el("\"Drug Likeness\" accordion header in context panel")));
      await session.step(58, "Then \"Drug Likeness\" pane in \"Biology\" pane in context panel should show the Chemical Engine's method", () => engineMethod(page, "Drug Likeness", "Biology"));
      await session.step(59, "And \"Drug Likeness\" pane in \"Biology\" pane in context panel should not contain the text \"Could not asses drug likeness\"", () => shouldNotContainText(page, el("\"Drug Likeness\" pane in \"Biology\" pane in context panel"), "Could not asses drug likeness"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The 2D Structure pane draws the molecule of the current cell", async () => {
      await session.step(64, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(65, "And user expands Structure accordion header in context panel", () => expand(page, el("Structure accordion header in context panel")));
      await session.step(66, "And user expands \"2D Structure\" accordion header in context panel", () => expand(page, el("\"2D Structure\" accordion header in context panel")));
      await session.step(67, "Then the canvases of \"2D Structure\" pane in context panel should be painted in at least 1 colors", () => canvasColors(page, el("\"2D Structure\" pane in context panel"), 1));
      await session.step(68, "And \"2D Structure\" pane in context panel should not contain the text \"Molecule is possibly malformed\"", () => shouldNotContainText(page, el("\"2D Structure\" pane in context panel"), "Molecule is possibly malformed"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Structural Alerts names the alerts of the current molecule and follows the current row", async () => {
      await session.step(73, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(74, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(75, "And the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
      await session.step(76, "When user expands Biology accordion header in context panel", () => expand(page, el("Biology accordion header in context panel")));
      await session.step(77, "And user expands \"Structural Alerts\" accordion header in context panel", () => expand(page, el("\"Structural Alerts\" accordion header in context panel")));
      await session.step(78, "Then \"Structural Alerts\" pane in \"Biology\" pane in context panel should contain the text \"Ketone\"", () => shouldContainText(page, el("\"Structural Alerts\" pane in \"Biology\" pane in context panel"), "Ketone"));
      await session.step(79, "And \"Structural Alerts\" pane in \"Biology\" pane in context panel should contain the text \"aliphatic ester, not lactones\"", () => shouldContainText(page, el("\"Structural Alerts\" pane in \"Biology\" pane in context panel"), "aliphatic ester, not lactones"));
      await session.step(80, "And \"Calculate for whole dataset\" button in \"Structural Alerts\" pane in context panel should be visible", () => shouldBe(page, el("\"Calculate for whole dataset\" button in \"Structural Alerts\" pane in context panel"), "visible"));
      await session.step(81, "When user clicks on the \"cell 5 of canonical_smiles\" area of grid", () => clickArea(page, "cell 5 of canonical_smiles", el("grid")));
      await session.step(82, "Then row 5 should be current", () => currentRowIs(page, 5));
      await session.step(83, "And the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
      await session.step(84, "And \"Structural Alerts\" pane in \"Biology\" pane in context panel should contain the text \"No alerts\"", () => shouldContainText(page, el("\"Structural Alerts\" pane in \"Biology\" pane in context panel"), "No alerts"));
      await session.step(85, "And \"Structural Alerts\" pane in \"Biology\" pane in context panel should not contain the text \"Ketone\"", () => shouldNotContainText(page, el("\"Structural Alerts\" pane in \"Biology\" pane in context panel"), "Ketone"));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Pharmacophores counts the feature families of the current molecule", async () => {
      await session.step(90, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(91, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(92, "And the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
      await session.step(93, "When user expands Biology accordion header in context panel", () => expand(page, el("Biology accordion header in context panel")));
      await session.step(94, "And user expands Pharmacophores accordion header in context panel", () => expand(page, el("Pharmacophores accordion header in context panel")));
      await session.step(95, "Then \"Pharmacophores\" pane in \"Biology\" pane in context panel should contain the text \"Hydrophobic (5)\"", () => shouldContainText(page, el("\"Pharmacophores\" pane in \"Biology\" pane in context panel"), "Hydrophobic (5)"));
      await session.step(96, "And \"Pharmacophores\" pane in \"Biology\" pane in context panel should contain the text \"Aromatic (3)\"", () => shouldContainText(page, el("\"Pharmacophores\" pane in \"Biology\" pane in context panel"), "Aromatic (3)"));
      await session.step(97, "And \"Pharmacophores\" pane in \"Biology\" pane in context panel should not contain the text \"Donor\"", () => shouldNotContainText(page, el("\"Pharmacophores\" pane in \"Biology\" pane in context panel"), "Donor"));
      await session.step(98, "When user clicks on the \"cell 5 of canonical_smiles\" area of grid", () => clickArea(page, "cell 5 of canonical_smiles", el("grid")));
      await session.step(99, "Then row 5 should be current", () => currentRowIs(page, 5));
      await session.step(100, "And the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
      await session.step(101, "And \"Pharmacophores\" pane in \"Biology\" pane in context panel should contain the text \"Donor (1)\"", () => shouldContainText(page, el("\"Pharmacophores\" pane in \"Biology\" pane in context panel"), "Donor (1)"));
      await session.step(102, "And \"Pharmacophores\" pane in \"Biology\" pane in context panel should contain the text \"Positive (1)\"", () => shouldContainText(page, el("\"Pharmacophores\" pane in \"Biology\" pane in context panel"), "Positive (1)"));
      await session.step(103, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Properties lists the molecule's properties and adds one as a column of the whole table", async () => {
      await session.step(107, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(108, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(109, "And the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
      await session.step(110, "When user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
      await session.step(111, "And user expands Properties accordion header in context panel", () => expand(page, el("Properties accordion header in context panel")));
      await session.step(112, "Then \"Properties\" pane in \"Chemistry\" pane in context panel should contain the text \"Rotatable bonds\"", () => shouldContainText(page, el("\"Properties\" pane in \"Chemistry\" pane in context panel"), "Rotatable bonds"));
      await session.step(113, "And \"Properties\" pane in \"Chemistry\" pane in context panel should contain the text \"C24H23NO3\"", () => shouldContainText(page, el("\"Properties\" pane in \"Chemistry\" pane in context panel"), "C24H23NO3"));
      await session.step(114, "And \"Properties\" pane in \"Chemistry\" pane in context panel should not contain the text \"Could not analyze properties\"", () => shouldNotContainText(page, el("\"Properties\" pane in \"Chemistry\" pane in context panel"), "Could not analyze properties"));
      await session.step(115, "When user clicks on the \"cell 5 of canonical_smiles\" area of grid", () => clickArea(page, "cell 5 of canonical_smiles", el("grid")));
      await session.step(116, "Then row 5 should be current", () => currentRowIs(page, 5));
      await session.step(117, "And the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
      await session.step(118, "And \"Properties\" pane in \"Chemistry\" pane in context panel should contain the text \"C12H14\"", () => shouldContainText(page, el("\"Properties\" pane in \"Chemistry\" pane in context panel"), "C12H14"));
      await session.step(119, "And the table should not have a column \"MW\"", () => hasNoColumn(page, "MW"));
      await session.step(120, "When user hovers over \"Properties\" pane in \"Chemistry\" pane in context panel", () => hoverOver(page, el("\"Properties\" pane in \"Chemistry\" pane in context panel")));
      await session.step(121, "And user clicks on \"Calculate MW for the whole table\" icon in \"Properties\" pane in \"Chemistry\" pane in context panel", () => clickOn(page, el("\"Calculate MW for the whole table\" icon in \"Properties\" pane in \"Chemistry\" pane in context panel")));
      await session.step(122, "Then the table should have a column \"MW\"", () => hasColumn(page, "MW"));
      await session.step(123, "And \"MW\" column should have no missing values", () => columnComplete(page, "MW"));
      await session.step(124, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Rendering pane's Highlight color repaints the column's scaffold highlight", async () => {
      await session.step(128, "When user opens chembl-scaffolds dataset", () => openDataset(page, ds("chembl-scaffolds")));
      await session.step(129, "And user clicks on the \"header Smiles\" area of grid", () => clickArea(page, "header Smiles", el("grid")));
      await session.step(130, "And user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
      await session.step(131, "And user expands Rendering accordion header in context panel", () => expand(page, el("Rendering accordion header in context panel")));
      await session.step(132, "And user selects \"Scaffold\" in \"Scaffold column\" choice input in \"Rendering\" pane in context panel", () => selectIn(page, "Scaffold", el("\"Scaffold column\" choice input in \"Rendering\" pane in context panel")));
      await session.step(133, "And user checks \"Highlight scaffold\" checkbox in \"Rendering\" pane in context panel", () => check(page, el("\"Highlight scaffold\" checkbox in \"Rendering\" pane in context panel")));
      await session.step(134, "Then the \"cell 1 of Smiles\" area of grid should not contain the color \"#00e0e0\"", () => areaNotColor(page, "cell 1 of Smiles", el("grid"), "#00e0e0"));
      await session.step(135, "When user enters \"#00e0e0\" into \"Highlight color\" input in \"Rendering\" pane in context panel", () => enterInto(page, "#00e0e0", el("\"Highlight color\" input in \"Rendering\" pane in context panel")));
      await session.step(136, "Then the \"cell 1 of Smiles\" area of grid should contain the color \"#00e0e0\"", () => areaColor(page, "cell 1 of Smiles", el("grid"), "#00e0e0"));
      await session.step(137, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
