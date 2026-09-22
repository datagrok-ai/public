/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/panels/molecule-panels.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.panels-chemistry-mixture, chem.cp.panels-synthon-search]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, shouldBe, shouldContainText, shouldHaveValue, shouldNotContainText, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {cellIsCurrentObject, tableOpen} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, contextPanelOpen, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, readingAtLeast} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Mixture, synthon-search and highlight panes of the Chem context panel", () => {
  const session = feature(test, "features/panels/molecule-panels.feature", import.meta.url);
  test("Mixture, synthon-search and highlight panes of the Chem context panel", {tag: ["@journey", "@realizes:chem.cp.panels-chemistry-mixture", "@realizes:chem.cp.panels-synthon-search"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens test_mixtures dataset", () => openDataset(page, ds("test_mixtures")));
    await session.step(20, "And the context panel is open", () => contextPanelOpen(page));
    await run.scenario("A mixture cell's Chemistry group offers Mixture and MixtureTree", async () => {
      await session.step(23, "Then \"mixture\" column should have semantic type \"ChemicalMixture\"", () => columnSemType(page, "mixture", "ChemicalMixture"));
      await session.step(24, "Given the \"mixture\" cell of row 1 is the current object", () => cellIsCurrentObject(page, "mixture", 1));
      await session.step(25, "Then Chemistry accordion header in context panel should be visible", () => shouldBe(page, el("Chemistry accordion header in context panel"), "visible"));
      await session.step(26, "When user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
      await session.step(27, "Then \"Mixture\" pane in context panel should be visible", () => shouldBe(page, el("\"Mixture\" pane in context panel"), "visible"));
      await session.step(28, "And \"MixtureTree\" pane in context panel should be visible", () => shouldBe(page, el("\"MixtureTree\" pane in context panel"), "visible"));
      await session.step(29, "And \"Descriptors\" pane in context panel should be absent", () => shouldBe(page, el("\"Descriptors\" pane in context panel"), "absent"));
      await session.step(30, "And \"Gasteiger Partial Charges\" pane in context panel should be absent", () => shouldBe(page, el("\"Gasteiger Partial Charges\" pane in context panel"), "absent"));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Mixture pane draws the components of the mixture", async () => {
      await session.step(34, "Given the \"mixture\" cell of row 1 is the current object", () => cellIsCurrentObject(page, "mixture", 1));
      await session.step(35, "Then Chemistry accordion header in context panel should be visible", () => shouldBe(page, el("Chemistry accordion header in context panel"), "visible"));
      await session.step(36, "When user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
      await session.step(37, "And user expands Mixture accordion header in context panel", () => expand(page, el("Mixture accordion header in context panel")));
      await session.step(38, "Then grid in \"Mixture\" pane in context panel should be visible", () => shouldBe(page, el("grid in \"Mixture\" pane in context panel"), "visible"));
      await session.step(39, "And \"Mixture\" pane in context panel should be visible", () => shouldBe(page, el("\"Mixture\" pane in context panel"), "visible"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("MixtureTree names the mixfile version and gives each component its own pane", async () => {
      await session.step(43, "Given the \"mixture\" cell of row 1 is the current object", () => cellIsCurrentObject(page, "mixture", 1));
      await session.step(44, "Then Chemistry accordion header in context panel should be visible", () => shouldBe(page, el("Chemistry accordion header in context panel"), "visible"));
      await session.step(45, "When user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
      await session.step(46, "And user expands MixtureTree accordion header in context panel", () => expand(page, el("MixtureTree accordion header in context panel")));
      await session.step(47, "Then \"MixtureTree\" pane in context panel should contain the text \"mixfileVersion: 1\"", () => shouldContainText(page, el("\"MixtureTree\" pane in context panel"), "mixfileVersion: 1"));
      await session.step(48, "And \"t-butyllithium\" pane in \"MixtureTree\" pane in context panel should be visible", () => shouldBe(page, el("\"t-butyllithium\" pane in \"MixtureTree\" pane in context panel"), "visible"));
      await session.step(49, "And \"pentane\" pane in \"MixtureTree\" pane in context panel should be visible", () => shouldBe(page, el("\"pentane\" pane in \"MixtureTree\" pane in context panel"), "visible"));
      await session.step(50, "When user expands \"t-butyllithium\" accordion header in \"MixtureTree\" pane in context panel", () => expand(page, el("\"t-butyllithium\" accordion header in \"MixtureTree\" pane in context panel")));
      await session.step(51, "Then \"t-butyllithium\" pane in \"MixtureTree\" pane in context panel should contain the text \"quantity\"", () => shouldContainText(page, el("\"t-butyllithium\" pane in \"MixtureTree\" pane in context panel"), "quantity"));
      await session.step(52, "And \"t-butyllithium\" pane in \"MixtureTree\" pane in context panel should contain the text \"1.7\"", () => shouldContainText(page, el("\"t-butyllithium\" pane in \"MixtureTree\" pane in context panel"), "1.7"));
      await session.step(53, "And \"t-butyllithium\" pane in \"MixtureTree\" pane in context panel should contain the text \"mol/L\"", () => shouldContainText(page, el("\"t-butyllithium\" pane in \"MixtureTree\" pane in context panel"), "mol/L"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A three-component mixture gets three component panes", async () => {
      await session.step(57, "When user clicks on the \"cell 2 of mixture\" area of grid", () => clickArea(page, "cell 2 of mixture", el("grid")));
      await session.step(58, "And user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
      await session.step(59, "And user expands MixtureTree accordion header in context panel", () => expand(page, el("MixtureTree accordion header in context panel")));
      await session.step(60, "Then \"MixtureTree\" pane in context panel should contain the text \"mixfileVersion: 0.01\"", () => shouldContainText(page, el("\"MixtureTree\" pane in context panel"), "mixfileVersion: 0.01"));
      await session.step(61, "And \"phenol\" pane in \"MixtureTree\" pane in context panel should be visible", () => shouldBe(page, el("\"phenol\" pane in \"MixtureTree\" pane in context panel"), "visible"));
      await session.step(62, "And \"chloroform\" pane in \"MixtureTree\" pane in context panel should be visible", () => shouldBe(page, el("\"chloroform\" pane in \"MixtureTree\" pane in context panel"), "visible"));
      await session.step(63, "And \"isoamyl alcohol\" pane in \"MixtureTree\" pane in context panel should be visible", () => shouldBe(page, el("\"isoamyl alcohol\" pane in \"MixtureTree\" pane in context panel"), "visible"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Databases group offers the two Synthon Search panes with their controls", async () => {
      await session.step(68, "Given user opens smiles dataset", () => openDataset(page, ds("smiles")));
      await session.step(69, "When user clicks on the \"cell 2 of canonical_smiles\" area of grid", () => clickArea(page, "cell 2 of canonical_smiles", el("grid")));
      await session.step(70, "And user expands Databases accordion header in context panel", () => expand(page, el("Databases accordion header in context panel")));
      await session.step(71, "And user expands \"Synthon Search\" accordion header in \"Databases\" pane in context panel", () => expand(page, el("\"Synthon Search\" accordion header in \"Databases\" pane in context panel")));
      await session.step(72, "Then \"Substructure Search\" pane in \"Synthon Search\" pane in context panel should be visible", () => shouldBe(page, el("\"Substructure Search\" pane in \"Synthon Search\" pane in context panel"), "visible"));
      await session.step(73, "And \"Similarity Search\" pane in \"Synthon Search\" pane in context panel should be visible", () => shouldBe(page, el("\"Similarity Search\" pane in \"Synthon Search\" pane in context panel"), "visible"));
      await session.step(74, "When user expands \"Substructure Search\" accordion header in \"Synthon Search\" pane in context panel", () => expand(page, el("\"Substructure Search\" accordion header in \"Synthon Search\" pane in context panel")));
      await session.step(75, "Then \"Space\" choice input in \"Substructure Search\" pane in context panel should have the value \"Syntons_5567.csv\"", () => shouldHaveValue(page, el("\"Space\" choice input in \"Substructure Search\" pane in context panel"), "Syntons_5567.csv"));
      await session.step(76, "And \"Max hits\" number input in \"Substructure Search\" pane in context panel should have the value \"100\"", () => shouldHaveValue(page, el("\"Max hits\" number input in \"Substructure Search\" pane in context panel"), "100"));
      await session.step(77, "And \"Include synthons\" checkbox in \"Substructure Search\" pane in context panel should be unchecked", () => shouldBe(page, el("\"Include synthons\" checkbox in \"Substructure Search\" pane in context panel"), "unchecked"));
      await session.step(78, "And \"Substructure Search\" pane in \"Synthon Search\" pane in context panel should not contain the text \"No synthon spaces found in synthon-data/\"", () => shouldNotContainText(page, el("\"Substructure Search\" pane in \"Synthon Search\" pane in context panel"), "No synthon spaces found in synthon-data/"));
      await session.step(79, "When user expands \"Similarity Search\" accordion header in \"Synthon Search\" pane in context panel", () => expand(page, el("\"Similarity Search\" accordion header in \"Synthon Search\" pane in context panel")));
      await session.step(80, "Then \"Cutoff\" slider in \"Similarity Search\" pane in context panel should have the value \"0.5\"", () => shouldHaveValue(page, el("\"Cutoff\" slider in \"Similarity Search\" pane in context panel"), "0.5"));
      await session.step(81, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Substructure Search pane returns hits and offers them as a table", async () => {
      await session.step(85, "Given user opens smiles dataset", () => openDataset(page, ds("smiles")));
      await session.step(86, "When user clicks on the \"cell 2 of canonical_smiles\" area of grid", () => clickArea(page, "cell 2 of canonical_smiles", el("grid")));
      await session.step(87, "And user expands Databases accordion header in context panel", () => expand(page, el("Databases accordion header in context panel")));
      await session.step(88, "And user expands \"Synthon Search\" accordion header in \"Databases\" pane in context panel", () => expand(page, el("\"Synthon Search\" accordion header in \"Databases\" pane in context panel")));
      await session.step(89, "And user expands \"Substructure Search\" accordion header in \"Synthon Search\" pane in context panel", () => expand(page, el("\"Substructure Search\" accordion header in \"Synthon Search\" pane in context panel")));
      await session.step(90, "Then grid in \"Substructure Search\" pane in context panel should be visible", () => shouldBe(page, el("grid in \"Substructure Search\" pane in context panel"), "visible"));
      await session.step(91, "And the \"rows shown\" reading of grid in \"Substructure Search\" pane in context panel should be at least 1", () => readingAtLeast(page, "rows shown", el("grid in \"Substructure Search\" pane in context panel"), 1));
      await session.step(92, "And there should be 1 visible \"Open compounds as table\" icon in \"Substructure Search\" pane in context panel", () => visibleCount(page, 1, el("\"Open compounds as table\" icon in \"Substructure Search\" pane in context panel")));
      await session.step(93, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Similarity Search pane returns hits and opens them as a table", async () => {
      await session.step(97, "Given user opens smiles dataset", () => openDataset(page, ds("smiles")));
      await session.step(98, "When user clicks on the \"cell 2 of canonical_smiles\" area of grid", () => clickArea(page, "cell 2 of canonical_smiles", el("grid")));
      await session.step(99, "And user expands Databases accordion header in context panel", () => expand(page, el("Databases accordion header in context panel")));
      await session.step(100, "And user expands \"Synthon Search\" accordion header in \"Databases\" pane in context panel", () => expand(page, el("\"Synthon Search\" accordion header in \"Databases\" pane in context panel")));
      await session.step(101, "And user expands \"Similarity Search\" accordion header in \"Synthon Search\" pane in context panel", () => expand(page, el("\"Similarity Search\" accordion header in \"Synthon Search\" pane in context panel")));
      await session.step(102, "Then grid in \"Similarity Search\" pane in context panel should be visible", () => shouldBe(page, el("grid in \"Similarity Search\" pane in context panel"), "visible"));
      await session.step(103, "And the \"rows shown\" reading of grid in \"Similarity Search\" pane in context panel should be at least 1", () => readingAtLeast(page, "rows shown", el("grid in \"Similarity Search\" pane in context panel"), 1));
      await session.step(104, "And \"Cutoff\" slider in \"Similarity Search\" pane in context panel should have the value \"0.5\"", () => shouldHaveValue(page, el("\"Cutoff\" slider in \"Similarity Search\" pane in context panel"), "0.5"));
      await session.step(105, "When user clicks on \"Open compounds as table\" icon in \"Similarity Search\" pane in context panel", () => clickOn(page, el("\"Open compounds as table\" icon in \"Similarity Search\" pane in context panel")));
      await session.step(106, "Then table \"Synthon Similarity Search Results\" should be open", () => tableOpen(page, "Synthon Similarity Search Results"));
      await session.step(107, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
