/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/table-view/shortcuts.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
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
import {isExpanded, pressKey, shouldBe, shouldContainText, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {onlyOfSelected, selectWhereIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeConsole, closeCurrentView, contextPanelOpen, contextPanelShows, currentViewType, layoutsDeleted, openTableOf, switchView, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {addViewer, clickArea, noErrors, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The table view's keyboard shortcuts open, clear and save what the help says", () => {
  const session = feature(test, "features/table-view/shortcuts.feature", import.meta.url);
  test("The table view's keyboard shortcuts open, clear and save what the help says", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the layouts named \"bdd-keys-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("bdd-keys-{time}")));
    await session.step(25, "And user opens a table \"bdd-keys-{time}\" with:", () => openTableOf(page, session.text("bdd-keys-{time}"), [["name","num"],["a","1"],["b","2"],["c","3"]]), [["name","num"],["a","1"],["b","2"],["c","3"]]);
    await session.step(30, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(31, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await run.scenario("Alt+X hides the Toolbox and shows it again", async () => {
      await session.step(34, "Then toolbox should be visible", () => shouldBe(page, el("toolbox"), "visible"));
      await session.step(35, "When user presses Alt+X", () => pressKey(page, "Alt+X"));
      await session.step(36, "Then toolbox should be hidden", () => shouldBe(page, el("toolbox"), "hidden"));
      await session.step(37, "And \"Toolbox\" status bar toggle should not be selected", () => shouldNotBe(page, el("\"Toolbox\" status bar toggle"), "selected"));
      await session.step(38, "When user presses Alt+X", () => pressKey(page, "Alt+X"));
      await session.step(39, "Then toolbox should be visible", () => shouldBe(page, el("toolbox"), "visible"));
      await session.step(40, "And \"Toolbox\" status bar toggle should be selected", () => shouldBe(page, el("\"Toolbox\" status bar toggle"), "selected"));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Alt+V docks the Variables panel and closes it again", async () => {
      await session.step(44, "Then \"Variables\" dock panel should be absent", () => shouldBe(page, el("\"Variables\" dock panel"), "absent"));
      await session.step(45, "When user presses Alt+V", () => pressKey(page, "Alt+V"));
      await session.step(46, "Then \"Variables\" dock panel should be visible", () => shouldBe(page, el("\"Variables\" dock panel"), "visible"));
      await session.step(47, "And \"Variables\" status bar toggle should be selected", () => shouldBe(page, el("\"Variables\" status bar toggle"), "selected"));
      await session.step(48, "When user presses Alt+V", () => pressKey(page, "Alt+V"));
      await session.step(49, "Then \"Variables\" dock panel should be absent", () => shouldBe(page, el("\"Variables\" dock panel"), "absent"));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The tilde opens the Console with its command line", async () => {
      await session.step(53, "Then console should be absent", () => shouldBe(page, el("console"), "absent"));
      await session.step(54, "When user presses Backquote", () => pressKey(page, "Backquote"));
      await session.step(55, "Then console should be visible", () => shouldBe(page, el("console"), "visible"));
      await session.step(56, "And \"Console\" status bar toggle should be selected", () => shouldBe(page, el("\"Console\" status bar toggle"), "selected"));
      await session.step(57, "When user closes the console", () => closeConsole(page));
      await session.step(58, "Then console should be hidden", () => shouldBe(page, el("console"), "hidden"));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Alt+S opens the Settings", async () => {
      await session.step(62, "When user presses Alt+S", () => pressKey(page, "Alt+S"));
      await session.step(63, "Then the current view should be a settings view", () => currentViewType(page, "settings"));
      await session.step(64, "And the \"Settings - General\" view should be current", () => viewIsCurrent(page, "Settings - General"));
      await session.step(65, "When user closes the current view", () => closeCurrentView(page));
      await session.step(66, "Then the \"bdd-keys-{time}\" view should be current", () => viewIsCurrent(page, session.text("bdd-keys-{time}")));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl+Shift+V opens a second view on the same table", async () => {
      await session.step(70, "When user presses Control+Shift+V", () => pressKey(page, "Control+Shift+V"));
      await session.step(71, "Then the \"bdd-keys-{time} copy\" view should be current", () => viewIsCurrent(page, session.text("bdd-keys-{time} copy")));
      await session.step(72, "And the open table views should be exactly \"bdd-keys-{time}, bdd-keys-{time}\"", () => openTableViewsExactly(page, session.text("bdd-keys-{time}, bdd-keys-{time}")));
      await session.step(73, "When user selects rows where \"name\" is \"b\"", () => selectWhereIs(page, "name", "b"));
      await session.step(74, "And user switches to the \"bdd-keys-{time}\" view", () => switchView(page, session.text("bdd-keys-{time}")));
      await session.step(75, "Then only rows where \"name\" is \"b\" should be selected", () => onlyOfSelected(page, "name", "b"));
      await session.step(76, "When user switches to the \"bdd-keys-{time} copy\" view", () => switchView(page, session.text("bdd-keys-{time} copy")));
      await session.step(77, "And user closes the current view", () => closeCurrentView(page));
      await session.step(78, "Then the open table views should be exactly \"bdd-keys-{time}\"", () => openTableViewsExactly(page, session.text("bdd-keys-{time}")));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl+Shift+R leaves only the grid, and Ctrl+Z gives the viewers back", async () => {
      await session.step(82, "Given user adds a histogram viewer", () => addViewer(page, "histogram"));
      await session.step(83, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
      await session.step(84, "Then the open tableview should have 1 histogram viewer", () => viewerCount(page, 1, "histogram"));
      await session.step(85, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(86, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
      await session.step(87, "And user presses Control+Shift+R", () => pressKey(page, "Control+Shift+R"));
      await session.step(88, "Then the open tableview should have 0 histogram viewers", () => viewerCount(page, 0, "histogram"));
      await session.step(89, "And the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(90, "And grid should be visible", () => shouldBe(page, el("grid"), "visible"));
      await session.step(91, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(92, "Then the open tableview should have 1 histogram viewer", () => viewerCount(page, 1, "histogram"));
      await session.step(93, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(94, "When user presses Control+Shift+R", () => pressKey(page, "Control+Shift+R"));
      await session.step(95, "Then the open tableview should have 0 histogram viewers", () => viewerCount(page, 0, "histogram"));
      await session.step(96, "And the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl+S saves the layout named after the table and shows it in the context panel", async () => {
      await session.step(100, "When user clicks on the \"header num\" area of grid", () => clickArea(page, "header num", el("grid")));
      await session.step(101, "Then the context panel should show \"num\"", () => contextPanelShows(page, "num"));
      await session.step(102, "When user presses Control+S", () => pressKey(page, "Control+S"));
      await session.step(103, "Then the context panel should show \"bdd-keys-{time}\"", () => contextPanelShows(page, session.text("bdd-keys-{time}")));
      await session.step(104, "And context panel should contain text \"Layout saved\"", () => shouldContainText(page, el("context panel"), "Layout saved"));
      await session.step(105, "Given Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
      await session.step(106, "And \"bdd-keys-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"bdd-keys-{time}\" layout card")), "visible"));
      await session.step(107, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ctrl+L docks the layout and viewer suggestion panels", async () => {
      await session.step(110, "Then \"Layout suggestions\" dock panel should be absent", () => shouldBe(page, el("\"Layout suggestions\" dock panel"), "absent"));
      await session.step(111, "When user presses Control+L", () => pressKey(page, "Control+L"));
      await session.step(112, "Then \"Layout suggestions\" dock panel should be visible", () => shouldBe(page, el("\"Layout suggestions\" dock panel"), "visible"));
      await session.step(113, "And \"Viewer suggestions\" dock panel should be visible", () => shouldBe(page, el("\"Viewer suggestions\" dock panel"), "visible"));
      await session.step(114, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
