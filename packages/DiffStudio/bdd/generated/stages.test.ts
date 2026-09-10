/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/stages.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [diffstudio.model.acid-production]
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openLibraryModel} from '../bindings/diff-studio.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, hoverOver, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, repainted, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A staged model and its inputs", () => {
  const session = feature(test, "features/stages.feature", import.meta.url);
  test("A staged model and its inputs", {tag: ["@journey", "@diffstudio", "@realizes:diffstudio.model.acid-production"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens the \"Acid production\" model of the Diff Studio library", () => openLibraryModel(page, "Acid production"));
    await run.scenario("The model arrives with its inputs and its plots", async () => {
      await session.step(15, "Then the \"Acid production\" view should be current", () => viewIsCurrent(page, "Acid production"));
      await session.step(16, "And \"1-st stage\" input should be visible", () => shouldBe(page, el("\"1-st stage\" input"), "visible"));
      await session.step(17, "And Multiaxis tab should be visible", () => shouldBe(page, el("Multiaxis tab"), "visible"));
      await session.step(18, "And Facet tab should be visible", () => shouldBe(page, el("Facet tab"), "visible"));
    });
    await run.scenario("Changing a stage duration redraws the solution", async () => {
      await session.step(21, "When user clicks on Multiaxis tab", () => clickOn(page, el("Multiaxis tab")));
      await session.step(22, "And user takes a snapshot of line chart viewer", () => takeSnapshot(page, el("line chart viewer")));
      await session.step(23, "And user enters \"50\" into \"1-st stage\" input", () => enterInto(page, "50", el("\"1-st stage\" input")));
      await session.step(24, "Then \"1-st stage\" input should have value \"50\"", () => shouldHaveValue(page, el("\"1-st stage\" input"), "50"));
      await session.step(25, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
    });
    await run.scenario("The inputs explain themselves on hover", async () => {
      await session.step(28, "When user hovers over \"1-st stage\" input", () => hoverOver(page, el("\"1-st stage\" input")));
      await session.step(29, "Then tooltip should contain text \"Duration of the 1-st stage\"", () => shouldContainText(page, el("tooltip"), "Duration of the 1-st stage"));
      await session.step(30, "When user hovers over biomass input", () => hoverOver(page, el("biomass input")));
      await session.step(31, "Then tooltip should contain text \"Aspergillus niger biomass\"", () => shouldContainText(page, el("tooltip"), "Aspergillus niger biomass"));
      await session.step(32, "When user hovers over glucose input", () => hoverOver(page, el("glucose input")));
      await session.step(33, "Then tooltip should contain text \"Glucose\"", () => shouldContainText(page, el("tooltip"), "Glucose"));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
