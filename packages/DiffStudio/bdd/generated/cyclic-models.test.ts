/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/cyclic-models.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [diffstudio.model.pk-pd]
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openLibraryModel} from '../bindings/diff-studio.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldContainText, shouldNotHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, repainted, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A cyclic model and its dosing inputs", () => {
  const session = feature(test, "features/cyclic-models.feature", import.meta.url);
  test("A cyclic model and its dosing inputs", {tag: ["@journey", "@diffstudio", "@realizes:diffstudio.model.pk-pd"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens the \"PK-PD\" model of the Diff Studio library", () => openLibraryModel(page, "PK-PD"));
    await run.scenario("The model arrives with its inputs and its plots", async () => {
      await session.step(16, "Then the \"PK-PD\" view should be current", () => viewIsCurrent(page, "PK-PD"));
      await session.step(17, "And count input should be visible", () => shouldBe(page, el("count input"), "visible"));
      await session.step(18, "And Multiaxis tab should be visible", () => shouldBe(page, el("Multiaxis tab"), "visible"));
      await session.step(19, "And Facet tab should be visible", () => shouldBe(page, el("Facet tab"), "visible"));
    });
    await run.scenario("The clickers move Count and the solution follows", async () => {
      await session.step(22, "When user clicks on Multiaxis tab", () => clickOn(page, el("Multiaxis tab")));
      await session.step(23, "And user takes a snapshot of line chart viewer", () => takeSnapshot(page, el("line chart viewer")));
      await session.step(24, "And user hovers over count input", () => hoverOver(page, el("count input")));
      await session.step(25, "And user clicks on plus icon in count input", () => clickOn(page, el("plus icon in count input")));
      await session.step(26, "Then count input should not have value \"1\"", () => shouldNotHaveValue(page, el("count input"), "1"));
      await session.step(27, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
    });
    await run.scenario("The dosing inputs explain themselves on hover", async () => {
      await session.step(30, "When user hovers over begin input", () => hoverOver(page, el("begin input")));
      await session.step(31, "Then tooltip should contain text \"Begin of dosing interval\"", () => shouldContainText(page, el("tooltip"), "Begin of dosing interval"));
      await session.step(32, "When user hovers over end input", () => hoverOver(page, el("end input")));
      await session.step(33, "Then tooltip should contain text \"End of dosing interval\"", () => shouldContainText(page, el("tooltip"), "End of dosing interval"));
      await session.step(34, "When user hovers over step input", () => hoverOver(page, el("step input")));
      await session.step(35, "Then tooltip should contain text \"Time step of simulation\"", () => shouldContainText(page, el("tooltip"), "Time step of simulation"));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
