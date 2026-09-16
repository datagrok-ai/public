/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-view-modes.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Browsing mode and persistent views", () => {
  const session = feature(test, "features/browse/browse-view-modes.feature", import.meta.url);
  test("A single click replaces the view the previous single click opened", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(26, "When user clicks on Tutorials tree node inside browse tree", () => clickOn(page, el("Tutorials tree node inside browse tree")));
    await session.step(27, "Then the \"Tutorials\" view should be current", () => viewIsCurrent(page, "Tutorials"));
    await session.step(28, "And Tutorials view should be visible", () => shouldBe(page, el("Tutorials view"), "visible"));
    await session.step(29, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(30, "Then Projects view should be visible", () => shouldBe(page, el("Projects view"), "visible"));
    await session.step(31, "And Tutorials view should be absent", () => shouldBe(page, el("Tutorials view"), "absent"));
    await session.step(32, "And no errors should have been logged", () => noErrors(page));
    await session.step(33, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A double click keeps the view through the next single click", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(36, "When user double-clicks on Tutorials tree node inside browse tree", () => doubleClickOn(page, el("Tutorials tree node inside browse tree")));
    await session.step(37, "Then the \"Tutorials\" view should be current", () => viewIsCurrent(page, "Tutorials"));
    await session.step(38, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(39, "Then Projects view should be visible", () => shouldBe(page, el("Projects view"), "visible"));
    await session.step(40, "And Tutorials view should be visible", () => shouldBe(page, el("Tutorials view"), "visible"));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
    await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A new browsing session does not unpin what is already persistent", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(49, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(50, "When user double-clicks on Tutorials tree node inside browse tree", () => doubleClickOn(page, el("Tutorials tree node inside browse tree")));
    await session.step(51, "Then Tutorials view should be visible", () => shouldBe(page, el("Tutorials view"), "visible"));
    await session.step(52, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(53, "Then Projects view should be visible", () => shouldBe(page, el("Projects view"), "visible"));
    await session.step(54, "When user clicks on Files---Demo tree node inside browse tree", () => clickOn(page, el("Files---Demo tree node inside browse tree")));
    await session.step(55, "Then Demo view should be visible", () => shouldBe(page, el("Demo view"), "visible"));
    await session.step(56, "And Tutorials view should be visible", () => shouldBe(page, el("Tutorials view"), "visible"));
    await session.step(57, "And Projects view should be absent", () => shouldBe(page, el("Projects view"), "absent"));
    await session.step(58, "And no errors should have been logged", () => noErrors(page));
    await session.step(59, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
