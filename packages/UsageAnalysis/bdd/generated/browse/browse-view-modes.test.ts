/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-view-modes.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, simpleModeOff, toolboxPaneHidden, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Browsing mode and persistent views", () => {
  const session = feature(test, "features/browse/browse-view-modes.feature", import.meta.url);
  test("A single click replaces the view the previous single click opened", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(29, "When user clicks on Tutorials tree node inside browse tree", () => clickOn(page, el("Tutorials tree node inside browse tree")));
    await session.step(30, "Then the \"Tutorials\" view should be current", () => viewIsCurrent(page, "Tutorials"));
    await session.step(31, "And Tutorials view should be visible", () => shouldBe(page, el("Tutorials view"), "visible"));
    await session.step(32, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(33, "Then Projects view should be visible", () => shouldBe(page, el("Projects view"), "visible"));
    await session.step(34, "And Tutorials view should be absent", () => shouldBe(page, el("Tutorials view"), "absent"));
    await session.step(35, "And no errors should have been logged", () => noErrors(page));
    await session.step(36, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A double click keeps the view through the next single click", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(39, "When user double-clicks on Tutorials tree node inside browse tree", () => doubleClickOn(page, el("Tutorials tree node inside browse tree")));
    await session.step(40, "Then the \"Tutorials\" view should be current", () => viewIsCurrent(page, "Tutorials"));
    await session.step(41, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(42, "Then Projects view should be visible", () => shouldBe(page, el("Projects view"), "visible"));
    await session.step(43, "And Tutorials view should be visible", () => shouldBe(page, el("Tutorials view"), "visible"));
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
    await session.step(45, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A new browsing session does not unpin what is already persistent", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
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
  test("The pin on a preview tab keeps the view through the next single click", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(62, "Given simple mode is off", () => simpleModeOff(page));
    await session.step(63, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(64, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(65, "When user clicks on Files---Demo---demog.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(66, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(67, "And \"This is Browse preview. Click to keep it open\" icon should be visible", () => shouldBe(page, el("\"This is Browse preview. Click to keep it open\" icon"), "visible"));
    await session.step(68, "When user clicks on \"This is Browse preview. Click to keep it open\" icon", () => clickOn(page, el("\"This is Browse preview. Click to keep it open\" icon")));
    await session.step(69, "Then \"This is Browse preview. Click to keep it open\" icon should be absent", () => shouldBe(page, el("\"This is Browse preview. Click to keep it open\" icon"), "absent"));
    await session.step(71, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(72, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(73, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(74, "Then Projects view should be visible", () => shouldBe(page, el("Projects view"), "visible"));
    await session.step(75, "And demog view should be present", () => shouldBe(page, el("demog view"), "present"));
    await session.step(76, "And no errors should have been logged", () => noErrors(page));
    await session.step(77, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("An edit keeps the previewed view through the next single click", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(80, "Given simple mode is off", () => simpleModeOff(page));
    await session.step(81, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(82, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(83, "When user clicks on Files---Demo---cars.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---cars.csv tree node inside browse tree")));
    await session.step(84, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
    await session.step(85, "And \"This is Browse preview. Click to keep it open\" icon should be visible", () => shouldBe(page, el("\"This is Browse preview. Click to keep it open\" icon"), "visible"));
    await session.step(86, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(87, "When user clicks on scatter plot icon in toolbox", () => clickOn(page, el("scatter plot icon in toolbox")));
    await session.step(88, "Then the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
    await session.step(89, "And \"This is Browse preview. Click to keep it open\" icon should be absent", () => shouldBe(page, el("\"This is Browse preview. Click to keep it open\" icon"), "absent"));
    await session.step(90, "Given the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(91, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(92, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(93, "Then Projects view should be visible", () => shouldBe(page, el("Projects view"), "visible"));
    await session.step(94, "And cars view should be present", () => shouldBe(page, el("cars view"), "present"));
    await session.step(95, "And no errors should have been logged", () => noErrors(page));
    await session.step(96, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
