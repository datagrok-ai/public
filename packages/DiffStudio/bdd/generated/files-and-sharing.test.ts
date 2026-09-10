/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/files-and-sharing.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [diffstudio.app.diff-studio]
--- */
import {test} from '@playwright/test';
import '../bindings/diff-studio.js';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, hoverOver, isExpanded, shouldBe, shouldHaveValue, shouldNotHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, openAddress, urlShouldContain} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A model file previewed from Browse", () => {
  const session = feature(test, "features/files-and-sharing.feature", import.meta.url);
  test("A model file previewed from Browse", {tag: ["@journey", "@diffstudio", "@realizes:diffstudio.app.diff-studio"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the browse panel is open", () => browsePanelOpen(page));
    await run.scenario("The library folder is reachable from the Files tree", async () => {
      await session.step(15, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(16, "And \"Files > App Data\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"Files > App Data\" tree node inside browse tree")));
      await session.step(17, "And \"Files > App Data > DiffStudio\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"Files > App Data > DiffStudio\" tree node inside browse tree")));
      await session.step(18, "Then \"Files > App Data > DiffStudio > library\" tree node inside browse tree should be visible", () => shouldBe(page, el("\"Files > App Data > DiffStudio > library\" tree node inside browse tree"), "visible"));
    });
    await run.scenario("A model file opens as a preview with its inputs", async () => {
      await session.step(21, "When user clicks on \"Files > App Data > DiffStudio > library\" tree node inside browse tree", () => clickOn(page, el("\"Files > App Data > DiffStudio > library\" tree node inside browse tree")));
      await session.step(22, "And user clicks on pk.ivp link in gallery", () => clickOn(page, el("pk.ivp link in gallery")));
      await session.step(23, "Then step input should be visible", () => shouldBe(page, el("step input"), "visible"));
      await session.step(24, "And count input should be visible", () => shouldBe(page, el("count input"), "visible"));
      await session.step(25, "And Multiaxis tab should be absent", () => shouldBe(page, el("Multiaxis tab"), "absent"));
      await session.step(26, "And Facet tab should be absent", () => shouldBe(page, el("Facet tab"), "absent"));
    });
    await run.scenario("The inputs the preview brings can be set", async () => {
      await session.step(29, "When user enters \"0.1\" into step input", () => enterInto(page, "0.1", el("step input")));
      await session.step(30, "Then step input should have value \"0.1\"", () => shouldHaveValue(page, el("step input"), "0.1"));
      await session.step(31, "When user hovers over count input", () => hoverOver(page, el("count input")));
      await session.step(32, "And user clicks on plus icon in count input", () => clickOn(page, el("plus icon in count input")));
      await session.step(33, "Then count input should not have value \"1\"", () => shouldNotHaveValue(page, el("count input"), "1"));
    });
    await run.scenario("The address carries the inputs, and loading it again brings them back", async () => {
      await session.step(36, "Then the page address should contain \"step\"", () => urlShouldContain(page, "step"));
      await session.step(37, "When user opens the page address of the current view", () => openAddress(page));
      await session.step(38, "Then step input should have value \"0.10\"", () => shouldHaveValue(page, el("step input"), "0.10"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
