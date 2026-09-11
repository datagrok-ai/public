/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/files-and-sharing.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [diffstudio.app.diff-studio]
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {reopenModelAddress} from '../bindings/diff-studio.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, dragSliderTo, enterInto, hoverOver, isExpanded, shouldBe, shouldHaveValue, shouldHaveValueBetween} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, urlShouldContain} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, repainted, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A model file previewed from Browse", () => {
  const session = feature(test, "features/files-and-sharing.feature", import.meta.url);
  test("A model file previewed from Browse", {tag: ["@journey", "@diffstudio", "@realizes:diffstudio.app.diff-studio"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await run.scenario("The library folder is reachable from the Files tree", async () => {
      await session.step(18, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(19, "And \"Files > App Data\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"Files > App Data\" tree node inside browse tree")));
      await session.step(20, "And \"Files > App Data > DiffStudio\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"Files > App Data > DiffStudio\" tree node inside browse tree")));
      await session.step(21, "Then \"Files > App Data > DiffStudio > library\" tree node inside browse tree should be visible", () => shouldBe(page, el("\"Files > App Data > DiffStudio > library\" tree node inside browse tree"), "visible"));
    });
    await run.scenario("A model file opens as a preview with its inputs", async () => {
      await session.step(24, "When user clicks on \"Files > App Data > DiffStudio > library\" tree node inside browse tree", () => clickOn(page, el("\"Files > App Data > DiffStudio > library\" tree node inside browse tree")));
      await session.step(25, "And user clicks on pk.ivp link in gallery", () => clickOn(page, el("pk.ivp link in gallery")));
      await session.step(26, "Then step input should be visible", () => shouldBe(page, el("step input"), "visible"));
      await session.step(27, "And count input should be visible", () => shouldBe(page, el("count input"), "visible"));
      await session.step(28, "And Multiaxis tab should be absent", () => shouldBe(page, el("Multiaxis tab"), "absent"));
      await session.step(29, "And Facet tab should be absent", () => shouldBe(page, el("Facet tab"), "absent"));
    });
    await run.scenario("The slider sets Step, as a reader would set it", async () => {
      await session.step(32, "When user takes a snapshot of line chart viewer", () => takeSnapshot(page, el("line chart viewer")));
      await session.step(33, "And user drags the slider of step input to 0.1", () => dragSliderTo(page, el("step input"), 0.1));
      await session.step(34, "Then step input should have a value between 0.09 and 0.1", () => shouldHaveValueBetween(page, el("step input"), 0.09, 0.1));
      await session.step(35, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
    });
    await run.scenario("The clicker counts Count up to four", async () => {
      await session.step(38, "When user hovers over count input", () => hoverOver(page, el("count input")));
      await session.step(39, "And user clicks on plus icon in count input", () => clickOn(page, el("plus icon in count input")));
      await session.step(40, "And user clicks on plus icon in count input", () => clickOn(page, el("plus icon in count input")));
      await session.step(41, "And user clicks on plus icon in count input", () => clickOn(page, el("plus icon in count input")));
      await session.step(42, "Then count input should have value \"4\"", () => shouldHaveValue(page, el("count input"), "4"));
    });
    await run.scenario("The inputs take a typed value too", async () => {
      await session.step(45, "When user enters \"0.1\" into step input", () => enterInto(page, "0.1", el("step input")));
      await session.step(46, "Then step input should have value \"0.1\"", () => shouldHaveValue(page, el("step input"), "0.1"));
    });
    await run.scenario("The address carries the inputs, and loading it again brings them back", async () => {
      await session.step(49, "Then the page address should contain \"step\"", () => urlShouldContain(page, "step"));
      await session.step(50, "When user opens the model at the page address", () => reopenModelAddress(page));
      await session.step(51, "Then step input should have value \"0.10\"", () => shouldHaveValue(page, el("step input"), "0.10"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
