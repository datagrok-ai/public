/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/01-getting-started/06-arrange-your-workspace.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clickViewTab, contextPanelClosed, dockViewBeside, openDataset, resizeViewShare, simpleModeOff, toolboxPaneHidden, viewDockedBeside, viewIsCurrent, viewTakesShare} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Basics", () => {
  const session = feature(test, "features/howto/01-getting-started/06-arrange-your-workspace.feature", import.meta.url);
  test("Arrange your workspace", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And simple mode is off", () => simpleModeOff(page));
    await session.step(21, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(22, "And user opens iris dataset", () => openDataset(page, ds("iris")));
    await session.step(23, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(24, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(26, "When user clicks on the tab of the \"cars\" view", () => clickViewTab(page, "cars"), undefined, "Each open table has its own tab: click one to switch");
    await session.step(27, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
    await session.step(29, "When user drags the tab of the \"iris\" view to the right side of the \"cars\" view", () => dockViewBeside(page, "iris", "right", "cars"), undefined, "Drag a tab to the side of a view to see both");
    await session.step(30, "Then the \"iris\" view should be docked right-of the \"cars\" view", () => viewDockedBeside(page, "iris", "right-of", "cars"));
    await session.step(32, "When user clicks on the tab of the \"cars\" view", () => clickViewTab(page, "cars"), undefined, "Click cars to show it on the left");
    await session.step(33, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
    await session.step(35, "When user resizes the \"cars\" view to 60% of the width it shares with the \"iris\" view", () => resizeViewShare(page, "cars", 60, "width", "iris"), undefined, "Drag the line between them to share the width");
    await session.step(36, "Then the \"cars\" view should take about 60% of the width it shares with the \"iris\" view", () => viewTakesShare(page, "cars", 60, "width", "iris"));
    await session.step(38, "When user clicks on \"Toolbox\" status bar toggle", () => clickOn(page, el("\"Toolbox\" status bar toggle")), undefined, "The toggles at the bottom right show and hide panels");
    await session.step(39, "Then toolbox should be visible", () => shouldBe(page, el("toolbox"), "visible"));
    await session.step(40, "And \"Toolbox\" status bar toggle should be selected", () => shouldBe(page, el("\"Toolbox\" status bar toggle"), "selected"));
    await session.step(42, "When user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")), undefined, "Context Panel, also F4");
    await session.step(43, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
    await session.step(45, "When user clicks on \"Console\" status bar toggle", () => clickOn(page, el("\"Console\" status bar toggle")), undefined, "Console: commands and the log");
    await session.step(46, "Then console should be visible", () => shouldBe(page, el("console"), "visible"));
    await session.step(48, "When user clicks on \"Console\" status bar toggle", () => clickOn(page, el("\"Console\" status bar toggle")), undefined, "Click again to hide it");
    await session.step(49, "Then console should be hidden", () => shouldBe(page, el("console"), "hidden"));
    await session.step(51, "When user clicks on \"Tabs\" status bar toggle", () => clickOn(page, el("\"Tabs\" status bar toggle")), undefined, "Tabs: hide the view tabs for more room");
    await session.step(52, "Then \"cars\" view tab should be hidden", () => shouldBe(page, el("\"cars\" view tab"), "hidden"));
    await session.step(53, "When user clicks on \"Tabs\" status bar toggle", () => clickOn(page, el("\"Tabs\" status bar toggle")));
    await session.step(54, "Then \"cars\" view tab should be visible", () => shouldBe(page, el("\"cars\" view tab"), "visible"));
    await session.step(56, "When user clicks on \"Presentation mode\" status bar toggle", () => clickOn(page, el("\"Presentation mode\" status bar toggle")), undefined, "Presentation mode: just the current view, also F7");
    await session.step(57, "Then status bar should be hidden", () => shouldBe(page, el("status bar"), "hidden"));
    await session.step(58, "And toolbox should be hidden", () => shouldBe(page, el("toolbox"), "hidden"));
    await session.step(60, "When user clicks on \"back to design mode\" link", () => clickOn(page, el("\"back to design mode\" link")), undefined, "This link brings the panels back");
    await session.step(61, "Then status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
    await session.step(62, "And toolbox should be visible", () => shouldBe(page, el("toolbox"), "visible"));
    await session.step(63, "And no errors should have been logged", () => noErrors(page));
  });
});
