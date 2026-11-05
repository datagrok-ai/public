/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/model-hub/model-hub-tools.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, isExpanded, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {browsePanelOpen, closeCurrentView, currentViewShows, galleryCountLower, rememberGalleryCount, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Model Hub's gallery search and its Tools menu", () => {
  const session = feature(test, "features/model-hub/model-hub-tools.feature", import.meta.url);
  test("The Model Hub's gallery search and its Tools menu", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(24, "And Apps---Compute tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Compute tree node inside browse tree")));
    await session.step(25, "When user clicks on Apps---Compute---Model-Hub tree node inside browse tree", () => clickOn(page, el("Apps---Compute---Model-Hub tree node inside browse tree")));
    await session.step(26, "Then the \"Model Hub\" view should be current", () => viewIsCurrent(page, "Model Hub"));
    await run.scenario("The search narrows the gallery to the models whose name matches", async () => {
      await session.step(29, "Then Bioreactor link in gallery should be visible", () => shouldBe(page, el("Bioreactor link in gallery"), "visible"));
      await session.step(30, "And PK-PD link in gallery should be visible", () => shouldBe(page, el("PK-PD link in gallery"), "visible"));
      await session.step(31, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(32, "And user types \"Bioreactor\" into gallery search", () => typeInto(page, "Bioreactor", el("gallery search")));
      await session.step(33, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(34, "And Bioreactor link in gallery should be visible", () => shouldBe(page, el("Bioreactor link in gallery"), "visible"));
      await session.step(35, "And PK-PD link in gallery should be hidden", () => shouldBe(page, el("PK-PD link in gallery"), "hidden"));
      await session.step(36, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(37, "Then PK-PD link in gallery should be visible", () => shouldBe(page, el("PK-PD link in gallery"), "visible"));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Tools > Compare Runs... opens an empty Run Comparison view that asks for runs", async () => {
      await session.step(41, "When user picks \"Tools > Compare Runs...\" from the top menu", () => pickFromTopMenu(page, "Tools > Compare Runs..."));
      await session.step(42, "Then the \"Run Comparison\" view should be current", () => viewIsCurrent(page, "Run Comparison"));
      await session.step(43, "And the current view should show the text \"Select model...\"", () => currentViewShows(page, "Select model..."));
      await session.step(44, "And the current view should show the text \"Add at least two runs (or tables) to compare\"", () => currentViewShows(page, "Add at least two runs (or tables) to compare"));
      await session.step(45, "And \"Add table\" input should be visible", () => shouldBe(page, el("\"Add table\" input"), "visible"));
      await session.step(46, "When user closes the current view", () => closeCurrentView(page));
      await session.step(47, "Then the \"Model Hub\" view should be current", () => viewIsCurrent(page, "Model Hub"));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
