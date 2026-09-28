/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/panels/molecule-panels.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.panels-chemistry-mixture]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {expand, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {cellIsCurrentObject} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, contextPanelOpen, openDataset, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Mixture and highlight panes of the Chem context panel", () => {
  const session = feature(test, "features/panels/molecule-panels.feature", import.meta.url);
  test("Mixture and highlight panes of the Chem context panel", {tag: ["@journey", "@realizes:chem.cp.panels-chemistry-mixture"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And user opens test_mixtures dataset", () => openDataset(page, ds("test_mixtures")));
    await session.step(19, "And the context panel is open", () => contextPanelOpen(page));
    await run.scenario("A mixture cell's Chemistry group offers Mixture and MixtureTree", async () => {
      await session.step(22, "Then \"mixture\" column should have semantic type \"ChemicalMixture\"", () => columnSemType(page, "mixture", "ChemicalMixture"));
      await session.step(23, "Given the \"mixture\" cell of row 1 is the current object", () => cellIsCurrentObject(page, "mixture", 1));
      await session.step(24, "Then Chemistry accordion header in context panel should be visible", () => shouldBe(page, el("Chemistry accordion header in context panel"), "visible"));
      await session.step(25, "When user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
      await session.step(26, "Then \"Mixture\" pane in context panel should be visible", () => shouldBe(page, el("\"Mixture\" pane in context panel"), "visible"));
      await session.step(27, "And \"MixtureTree\" pane in context panel should be visible", () => shouldBe(page, el("\"MixtureTree\" pane in context panel"), "visible"));
      await session.step(28, "And \"Descriptors\" pane in context panel should be absent", () => shouldBe(page, el("\"Descriptors\" pane in context panel"), "absent"));
      await session.step(29, "And \"Gasteiger Partial Charges\" pane in context panel should be absent", () => shouldBe(page, el("\"Gasteiger Partial Charges\" pane in context panel"), "absent"));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Mixture pane draws the components of the mixture", async () => {
      await session.step(33, "Given the \"mixture\" cell of row 1 is the current object", () => cellIsCurrentObject(page, "mixture", 1));
      await session.step(34, "Then Chemistry accordion header in context panel should be visible", () => shouldBe(page, el("Chemistry accordion header in context panel"), "visible"));
      await session.step(35, "When user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
      await session.step(36, "And user expands Mixture accordion header in context panel", () => expand(page, el("Mixture accordion header in context panel")));
      await session.step(37, "Then grid in \"Mixture\" pane in context panel should be visible", () => shouldBe(page, el("grid in \"Mixture\" pane in context panel"), "visible"));
      await session.step(38, "And \"Mixture\" pane in context panel should be visible", () => shouldBe(page, el("\"Mixture\" pane in context panel"), "visible"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("MixtureTree names the mixfile version and gives each component its own pane", async () => {
      await session.step(42, "Given the \"mixture\" cell of row 1 is the current object", () => cellIsCurrentObject(page, "mixture", 1));
      await session.step(43, "Then Chemistry accordion header in context panel should be visible", () => shouldBe(page, el("Chemistry accordion header in context panel"), "visible"));
      await session.step(44, "When user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
      await session.step(45, "And user expands MixtureTree accordion header in context panel", () => expand(page, el("MixtureTree accordion header in context panel")));
      await session.step(46, "Then \"MixtureTree\" pane in context panel should contain the text \"mixfileVersion: 1\"", () => shouldContainText(page, el("\"MixtureTree\" pane in context panel"), "mixfileVersion: 1"));
      await session.step(47, "And \"t-butyllithium\" pane in \"MixtureTree\" pane in context panel should be visible", () => shouldBe(page, el("\"t-butyllithium\" pane in \"MixtureTree\" pane in context panel"), "visible"));
      await session.step(48, "And \"pentane\" pane in \"MixtureTree\" pane in context panel should be visible", () => shouldBe(page, el("\"pentane\" pane in \"MixtureTree\" pane in context panel"), "visible"));
      await session.step(49, "When user expands \"t-butyllithium\" accordion header in \"MixtureTree\" pane in context panel", () => expand(page, el("\"t-butyllithium\" accordion header in \"MixtureTree\" pane in context panel")));
      await session.step(50, "Then \"t-butyllithium\" pane in \"MixtureTree\" pane in context panel should contain the text \"quantity\"", () => shouldContainText(page, el("\"t-butyllithium\" pane in \"MixtureTree\" pane in context panel"), "quantity"));
      await session.step(51, "And \"t-butyllithium\" pane in \"MixtureTree\" pane in context panel should contain the text \"1.7\"", () => shouldContainText(page, el("\"t-butyllithium\" pane in \"MixtureTree\" pane in context panel"), "1.7"));
      await session.step(52, "And \"t-butyllithium\" pane in \"MixtureTree\" pane in context panel should contain the text \"mol/L\"", () => shouldContainText(page, el("\"t-butyllithium\" pane in \"MixtureTree\" pane in context panel"), "mol/L"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A three-component mixture gets three component panes", async () => {
      await session.step(56, "When user clicks on the \"cell 2 of mixture\" area of grid", () => clickArea(page, "cell 2 of mixture", el("grid")));
      await session.step(57, "And user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
      await session.step(58, "And user expands MixtureTree accordion header in context panel", () => expand(page, el("MixtureTree accordion header in context panel")));
      await session.step(59, "Then \"MixtureTree\" pane in context panel should contain the text \"mixfileVersion: 0.01\"", () => shouldContainText(page, el("\"MixtureTree\" pane in context panel"), "mixfileVersion: 0.01"));
      await session.step(60, "And \"phenol\" pane in \"MixtureTree\" pane in context panel should be visible", () => shouldBe(page, el("\"phenol\" pane in \"MixtureTree\" pane in context panel"), "visible"));
      await session.step(61, "And \"chloroform\" pane in \"MixtureTree\" pane in context panel should be visible", () => shouldBe(page, el("\"chloroform\" pane in \"MixtureTree\" pane in context panel"), "visible"));
      await session.step(62, "And \"isoamyl alcohol\" pane in \"MixtureTree\" pane in context panel should be visible", () => shouldBe(page, el("\"isoamyl alcohol\" pane in \"MixtureTree\" pane in context panel"), "visible"));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
