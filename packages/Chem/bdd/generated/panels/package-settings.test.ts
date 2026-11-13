/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/panels/package-settings.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
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
import {above, atMostWidePx, inside, notClipped} from '@datagrok-libraries/bdd/bindings/common/geometry';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, dragBy, enterInto, expand, isExpanded, scrollIntoView, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, contextPanelOpen, contextPanelShows, openDataset, packageSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Chem's settings in the Plugins view", () => {
  const session = feature(test, "features/panels/package-settings.feature", import.meta.url);
  test("Chem's settings in the Plugins view", {tag: ["@journey", "@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And the \"Chem\" package settings are put back at feature end", () => packageSettingsPutBack(page, "Chem"));
    await session.step(19, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(20, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await run.scenario("The settings are grouped under their categories, in package.json's order", async () => {
      await session.step(23, "When user clicks on Platform---Plugins tree node inside browse tree", () => clickOn(page, el("Platform---Plugins tree node inside browse tree")));
      await session.step(24, "And user clicks on \"Chem\" gallery card", () => clickOn(page, el("\"Chem\" gallery card")));
      await session.step(25, "Then the context panel should show \"Chem\"", () => contextPanelShows(page, "Chem"));
      await session.step(26, "When user expands Settings accordion header in context panel", () => expand(page, el("Settings accordion header in context panel")));
      await session.step(27, "Then first \"Chemical engine\" label in \"Settings\" pane in context panel should be above \"Sketcher\" input in \"Settings\" pane in context panel", () => above(page, el("first \"Chemical engine\" label in \"Settings\" pane in context panel"), el("\"Sketcher\" input in \"Settings\" pane in context panel")));
      await session.step(28, "And \"Chemical Engine\" input in \"Settings\" pane in context panel should be above \"Rendering\" label in \"Settings\" pane in context panel", () => above(page, el("\"Chemical Engine\" input in \"Settings\" pane in context panel"), el("\"Rendering\" label in \"Settings\" pane in context panel")));
      await session.step(29, "And \"Rendering\" label in \"Settings\" pane in context panel should be above \"Highlight color\" input in \"Settings\" pane in context panel", () => above(page, el("\"Rendering\" label in \"Settings\" pane in context panel"), el("\"Highlight color\" input in \"Settings\" pane in context panel")));
      await session.step(30, "And \"Simplified stereo label\" checkbox in \"Settings\" pane in context panel should be above \"Enumeration\" label in \"Settings\" pane in context panel", () => above(page, el("\"Simplified stereo label\" checkbox in \"Settings\" pane in context panel"), el("\"Enumeration\" label in \"Settings\" pane in context panel")));
      await session.step(31, "And \"Enumeration\" label in \"Settings\" pane in context panel should be above \"Default Reactions\" input in \"Settings\" pane in context panel", () => above(page, el("\"Enumeration\" label in \"Settings\" pane in context panel"), el("\"Default Reactions\" input in \"Settings\" pane in context panel")));
      await session.step(32, "And \"Mark unspecified stereo\" checkbox in \"Settings\" pane in context panel should be checked", () => shouldBe(page, el("\"Mark unspecified stereo\" checkbox in \"Settings\" pane in context panel"), "checked"));
      await session.step(33, "And \"Simplified stereo label\" checkbox in \"Settings\" pane in context panel should be checked", () => shouldBe(page, el("\"Simplified stereo label\" checkbox in \"Settings\" pane in context panel"), "checked"));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("SAVE stays in sight in a narrow context panel", async () => {
      await session.step(37, "When user drags context panel splitter by 250 and 0 pixels", () => dragBy(page, el("context panel splitter"), 250, 0));
      await session.step(38, "Then context panel should be at most 300 px wide", () => atMostWidePx(page, el("context panel"), 300));
      await session.step(39, "When user scrolls SAVE button in \"Settings\" pane in context panel into view", () => scrollIntoView(page, el("SAVE button in \"Settings\" pane in context panel")));
      await session.step(40, "Then SAVE button in \"Settings\" pane in context panel should be inside context panel", () => inside(page, el("SAVE button in \"Settings\" pane in context panel"), el("context panel")));
      await session.step(41, "And SAVE button in \"Settings\" pane in context panel should not be clipped", () => notClipped(page, el("SAVE button in \"Settings\" pane in context panel")));
      await session.step(42, "When user drags context panel splitter by -250 and 0 pixels", () => dragBy(page, el("context panel splitter"), -250, 0));
      await session.step(43, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A saved Highlight color is every molecule column's at once", async () => {
      await session.step(46, "When user enters \"#2080ff\" into \"Highlight color\" input in \"Settings\" pane in context panel", () => enterInto(page, "#2080ff", el("\"Highlight color\" input in \"Settings\" pane in context panel")));
      await session.step(47, "And user clicks on SAVE button in \"Settings\" pane in context panel", () => clickOn(page, el("SAVE button in \"Settings\" pane in context panel")));
      await session.step(48, "Then SAVE button in \"Settings\" pane in context panel should be disabled", () => shouldBe(page, el("SAVE button in \"Settings\" pane in context panel"), "disabled"));
      await session.step(49, "When user opens chembl-scaffolds dataset", () => openDataset(page, ds("chembl-scaffolds")));
      await session.step(50, "And user clicks on the \"header Smiles\" area of grid", () => clickArea(page, "header Smiles", el("grid")));
      await session.step(51, "And user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
      await session.step(52, "And user expands Rendering accordion header in context panel", () => expand(page, el("Rendering accordion header in context panel")));
      await session.step(53, "Then \"Highlight color\" input in \"Rendering\" pane in context panel should have value \"#2080ff\"", () => shouldHaveValue(page, el("\"Highlight color\" input in \"Rendering\" pane in context panel"), "#2080ff"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
