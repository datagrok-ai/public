/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-shell-modes.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, openDataset, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The shell's modes around the Browse panel", () => {
  const session = feature(test, "features/browse/browse-shell-modes.feature", import.meta.url);
  test("Presentation mode puts the panels away and its back link gives them back", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(20, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(22, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "Then the browse tree should be visible", () => shouldBe(page, el("the browse tree"), "visible"));
    await session.step(25, "When user clicks on \"Presentation mode\" status bar toggle", () => clickOn(page, el("\"Presentation mode\" status bar toggle")));
    await session.step(26, "Then the browse tree should be hidden", () => shouldBe(page, el("the browse tree"), "hidden"));
    await session.step(27, "And status bar should be hidden", () => shouldBe(page, el("status bar"), "hidden"));
    await session.step(28, "And grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(29, "When user clicks on \"back to design mode\" link", () => clickOn(page, el("\"back to design mode\" link")));
    await session.step(30, "Then status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
    await session.step(31, "And the browse tree should be visible", () => shouldBe(page, el("the browse tree"), "visible"));
    await session.step(32, "And the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(33, "And no errors should have been logged", () => noErrors(page));
    await session.step(34, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("F7 enters presentation mode and leaves it", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(37, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(38, "And status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
    await session.step(39, "When user presses F7", () => pressKey(page, "F7"));
    await session.step(40, "Then status bar should be hidden", () => shouldBe(page, el("status bar"), "hidden"));
    await session.step(41, "And \"back to design mode\" link should be visible", () => shouldBe(page, el("\"back to design mode\" link"), "visible"));
    await session.step(42, "When user presses F7", () => pressKey(page, "F7"));
    await session.step(43, "Then status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
    await session.step(44, "And \"back to design mode\" link should be absent", () => shouldBe(page, el("\"back to design mode\" link"), "absent"));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
  });
  test("The Tabs toggle hides the view tabs and shows them again", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(49, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(50, "Then \"demog-1000\" view tab should be visible", () => shouldBe(page, el("\"demog-1000\" view tab"), "visible"));
    await session.step(51, "When user clicks on \"Tabs\" status bar toggle", () => clickOn(page, el("\"Tabs\" status bar toggle")));
    await session.step(52, "Then \"demog-1000\" view tab should be hidden", () => shouldBe(page, el("\"demog-1000\" view tab"), "hidden"));
    await session.step(53, "When user clicks on \"Tabs\" status bar toggle", () => clickOn(page, el("\"Tabs\" status bar toggle")));
    await session.step(54, "Then \"demog-1000\" view tab should be visible", () => shouldBe(page, el("\"demog-1000\" view tab"), "visible"));
    await session.step(55, "And the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(56, "And no errors should have been logged", () => noErrors(page));
    await session.step(57, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
