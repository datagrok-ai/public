/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/panels/info-panels.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cell.molecule, chem.panel.chemistry.rendering, chem.panel.chemistry.highlight, chem.panel.biology.toxicity, chem.panel.biology.drug-likeness, chem.panel.structure.2d-structure]
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
import {canvasColors} from '@datagrok-libraries/bdd/bindings/common/pixels';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {expand, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, contextPanelOpen, contextPanelShows, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Chemistry, Biology and Structure panes of the Chem context panel", () => {
  const session = feature(test, "features/panels/info-panels.feature", import.meta.url);
  test("The Chemistry, Biology and Structure panes of the Chem context panel", {tag: ["@journey", "@realizes:chem.cell.molecule", "@realizes:chem.panel.chemistry.rendering", "@realizes:chem.panel.chemistry.highlight", "@realizes:chem.panel.biology.toxicity", "@realizes:chem.panel.biology.drug-likeness", "@realizes:chem.panel.structure.2d-structure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(21, "And user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await run.scenario("The molecule column's Chemistry group offers Rendering and Highlight", async () => {
      await session.step(26, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(27, "And user clicks on the \"header canonical_smiles\" area of grid", () => clickArea(page, "header canonical_smiles", el("grid")));
      await session.step(28, "Then the context panel should show \"canonical_smiles\"", () => contextPanelShows(page, "canonical_smiles"));
      await session.step(29, "When user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
      await session.step(30, "Then \"Rendering\" pane in context panel should be visible", () => shouldBe(page, el("\"Rendering\" pane in context panel"), "visible"));
      await session.step(31, "And \"Highlight\" pane in context panel should be visible", () => shouldBe(page, el("\"Highlight\" pane in context panel"), "visible"));
      await session.step(32, "And \"Descriptors\" pane in context panel should be absent", () => shouldBe(page, el("\"Descriptors\" pane in context panel"), "absent"));
      await session.step(33, "When user expands Rendering accordion header in context panel", () => expand(page, el("Rendering accordion header in context panel")));
      await session.step(34, "Then \"Scaffold column\" choice input in \"Rendering\" pane in context panel should be visible", () => shouldBe(page, el("\"Scaffold column\" choice input in \"Rendering\" pane in context panel"), "visible"));
      await session.step(35, "And \"Highlight scaffold\" checkbox in \"Rendering\" pane in context panel should be visible", () => shouldBe(page, el("\"Highlight scaffold\" checkbox in \"Rendering\" pane in context panel"), "visible"));
      await session.step(36, "And \"Filter type\" choice input in \"Rendering\" pane in context panel should be visible", () => shouldBe(page, el("\"Filter type\" choice input in \"Rendering\" pane in context panel"), "visible"));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Toxicity grades the four risks and Drug Likeness reports a score", async () => {
      await session.step(41, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(42, "And user expands Biology accordion header in context panel", () => expand(page, el("Biology accordion header in context panel")));
      await session.step(43, "And user expands Toxicity accordion header in context panel", () => expand(page, el("Toxicity accordion header in context panel")));
      await session.step(44, "Then \"Toxicity\" pane in \"Biology\" pane in context panel should contain the text \"Mutagenicity\"", () => shouldContainText(page, el("\"Toxicity\" pane in \"Biology\" pane in context panel"), "Mutagenicity"));
      await session.step(45, "And \"Toxicity\" pane in \"Biology\" pane in context panel should contain the text \"Tumorigenicity\"", () => shouldContainText(page, el("\"Toxicity\" pane in \"Biology\" pane in context panel"), "Tumorigenicity"));
      await session.step(46, "And \"Toxicity\" pane in \"Biology\" pane in context panel should contain the text \"Irritating effects\"", () => shouldContainText(page, el("\"Toxicity\" pane in \"Biology\" pane in context panel"), "Irritating effects"));
      await session.step(47, "And \"Toxicity\" pane in \"Biology\" pane in context panel should contain the text \"Reproductive effects\"", () => shouldContainText(page, el("\"Toxicity\" pane in \"Biology\" pane in context panel"), "Reproductive effects"));
      await session.step(48, "And \"Toxicity\" pane in \"Biology\" pane in context panel should not contain the text \"Could not analyze toxicity\"", () => shouldNotContainText(page, el("\"Toxicity\" pane in \"Biology\" pane in context panel"), "Could not analyze toxicity"));
      await session.step(49, "When user expands \"Drug Likeness\" accordion header in context panel", () => expand(page, el("\"Drug Likeness\" accordion header in context panel")));
      await session.step(50, "Then \"Drug Likeness\" pane in \"Biology\" pane in context panel should contain the text \"Score:\"", () => shouldContainText(page, el("\"Drug Likeness\" pane in \"Biology\" pane in context panel"), "Score:"));
      await session.step(51, "And \"Drug Likeness\" pane in \"Biology\" pane in context panel should not contain the text \"Could not asses drug likeness\"", () => shouldNotContainText(page, el("\"Drug Likeness\" pane in \"Biology\" pane in context panel"), "Could not asses drug likeness"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The 2D Structure pane draws the molecule of the current cell", async () => {
      await session.step(56, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(57, "And user expands Structure accordion header in context panel", () => expand(page, el("Structure accordion header in context panel")));
      await session.step(58, "And user expands \"2D Structure\" accordion header in context panel", () => expand(page, el("\"2D Structure\" accordion header in context panel")));
      await session.step(59, "Then the canvases of \"2D Structure\" pane in context panel should be painted in at least 1 colors", () => canvasColors(page, el("\"2D Structure\" pane in context panel"), 1));
      await session.step(60, "And \"2D Structure\" pane in context panel should not contain the text \"Molecule is possibly malformed\"", () => shouldNotContainText(page, el("\"2D Structure\" pane in context panel"), "Molecule is possibly malformed"));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
