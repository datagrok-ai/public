/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/hub.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [diffstudio.app.diff-studio]
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp} from '../bindings/steps.js';
import {clickOn, doubleClickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The app's hub, and the Open model menu", () => {
  const session = feature(test, "features/hub.feature", import.meta.url);
  test("The app's hub, and the Open model menu", {tag: ["@journey", "@diffstudio", "@realizes:diffstudio.app.diff-studio"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await run.scenario("The app opens on its hub, and a library card opens a model", async () => {
      await session.step(9, "Given user opens the Diff Studio app", () => openApp(page));
      await session.step(10, "Then the \"Diff Studio\" view should be current", () => viewIsCurrent(page, "Diff Studio"));
      await session.step(11, "And Create button should be visible", () => shouldBe(page, el("Create button"), "visible"));
      await session.step(12, "When user double-clicks on Bioreactor hub card", () => doubleClickOn(page, el("Bioreactor hub card")));
      await session.step(13, "Then the \"Bioreactor\" view should be current", () => viewIsCurrent(page, "Bioreactor"));
      await session.step(14, "And \"Process mode\" input should be visible", () => shouldBe(page, el("\"Process mode\" input"), "visible"));
    });
    await run.scenario("The Open model icon switches to another model of the library", async () => {
      await session.step(17, "When user clicks on open model button", () => clickOn(page, el("open model button")));
      await session.step(18, "And user picks \"Library > PK-PD\" from the open menu", () => pickFromOpenMenu(page, "Library > PK-PD"));
      await session.step(19, "Then the \"PK-PD\" view should be current", () => viewIsCurrent(page, "PK-PD"));
      await session.step(20, "And dose input should be visible", () => shouldBe(page, el("dose input"), "visible"));
      await session.step(21, "And \"Process mode\" input should be absent", () => shouldBe(page, el("\"Process mode\" input"), "absent"));
      await session.step(22, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
