/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/panels/synthon-search.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.panels-synthon-search]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, shouldBecomeVisibleWithin, shouldHaveValue, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {tableOpen} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, contextPanelOpen, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, readingAtLeast} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Synthon Search panes return hits", () => {
  const session = feature(test, "features/panels/synthon-search.feature", import.meta.url);
  test("The Synthon Search panes return hits", {tag: ["@journey", "@full-stand", "@realizes:chem.cp.panels-synthon-search"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await session.step(15, "And the context panel is open", () => contextPanelOpen(page));
    await run.scenario("The Substructure Search pane returns hits and offers them as a table", async () => {
      await session.step(18, "When user clicks on the \"cell 2 of canonical_smiles\" area of grid", () => clickArea(page, "cell 2 of canonical_smiles", el("grid")));
      await session.step(19, "And user expands Databases accordion header in context panel", () => expand(page, el("Databases accordion header in context panel")));
      await session.step(20, "And user expands \"Synthon Search\" accordion header in \"Databases\" pane in context panel", () => expand(page, el("\"Synthon Search\" accordion header in \"Databases\" pane in context panel")));
      await session.step(21, "And user expands \"Substructure Search\" accordion header in \"Synthon Search\" pane in context panel", () => expand(page, el("\"Substructure Search\" accordion header in \"Synthon Search\" pane in context panel")));
      await session.step(22, "Then grid in \"Substructure Search\" pane in context panel should become visible within 120 seconds", () => shouldBecomeVisibleWithin(page, el("grid in \"Substructure Search\" pane in context panel"), 120));
      await session.step(23, "And the \"rows shown\" reading of grid in \"Substructure Search\" pane in context panel should be at least 1", () => readingAtLeast(page, "rows shown", el("grid in \"Substructure Search\" pane in context panel"), 1));
      await session.step(24, "And there should be 1 visible \"Open compounds as table\" icon in \"Substructure Search\" pane in context panel", () => visibleCount(page, 1, el("\"Open compounds as table\" icon in \"Substructure Search\" pane in context panel")));
      await session.step(25, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Similarity Search pane returns hits and opens them as a table", async () => {
      await session.step(28, "When user clicks on the \"cell 2 of canonical_smiles\" area of grid", () => clickArea(page, "cell 2 of canonical_smiles", el("grid")));
      await session.step(29, "And user expands Databases accordion header in context panel", () => expand(page, el("Databases accordion header in context panel")));
      await session.step(30, "And user expands \"Synthon Search\" accordion header in \"Databases\" pane in context panel", () => expand(page, el("\"Synthon Search\" accordion header in \"Databases\" pane in context panel")));
      await session.step(31, "And user expands \"Similarity Search\" accordion header in \"Synthon Search\" pane in context panel", () => expand(page, el("\"Similarity Search\" accordion header in \"Synthon Search\" pane in context panel")));
      await session.step(32, "Then grid in \"Similarity Search\" pane in context panel should become visible within 120 seconds", () => shouldBecomeVisibleWithin(page, el("grid in \"Similarity Search\" pane in context panel"), 120));
      await session.step(33, "And the \"rows shown\" reading of grid in \"Similarity Search\" pane in context panel should be at least 1", () => readingAtLeast(page, "rows shown", el("grid in \"Similarity Search\" pane in context panel"), 1));
      await session.step(34, "And \"Cutoff\" slider in \"Similarity Search\" pane in context panel should have the value \"0.5\"", () => shouldHaveValue(page, el("\"Cutoff\" slider in \"Similarity Search\" pane in context panel"), "0.5"));
      await session.step(35, "When user clicks on \"Open compounds as table\" icon in \"Similarity Search\" pane in context panel", () => clickOn(page, el("\"Open compounds as table\" icon in \"Similarity Search\" pane in context panel")));
      await session.step(36, "Then table \"Synthon Similarity Search Results\" should be open", () => tableOpen(page, "Synthon Similarity Search Results"));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
