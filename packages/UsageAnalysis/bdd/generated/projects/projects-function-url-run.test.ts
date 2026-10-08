/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-function-url-run.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.functions, views.projects]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {appendToEditor, clearField, clickOn, clipboardContains, enterInto, hoverOver, isExpanded, pressKey, replaceCode, shouldBe, shouldContainText, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeCurrentView, contextPanelShows, currentViewType, dialogCloses, noScriptOnServer, saveScript, scriptsOnServer, toolboxPaneShown, urlShouldContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu, pickFromOpenMenu, pointerAway, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A function run straight from its URL", () => {
  const session = feature(test, "features/projects/projects-function-url-run.feature", import.meta.url);
  test("A function run straight from its URL", {tag: ["@journey", "@serial", "@realizes:views.functions", "@realizes:views.projects"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(27, "And no script named \"BDDUrlRunScript{time}\" is on the server", () => noScriptOnServer(page, session.text("BDDUrlRunScript{time}")));
    await run.scenario("The script is written in the Scripts view", async () => {
      await session.step(30, "Given Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
      await session.step(31, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
      await session.step(32, "When user clicks on Platform---Functions---Scripts tree node inside browse tree", () => clickOn(page, el("Platform---Functions---Scripts tree node inside browse tree")));
      await session.step(33, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
      await session.step(34, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(35, "And user picks \"JavaScript Script...\" from the open menu", () => pickFromOpenMenu(page, "JavaScript Script..."));
      await session.step(36, "Then the \"Template\" view should be current", () => viewIsCurrent(page, "Template"));
      await session.step(37, "When user replaces the code of code editor with \"//name: BDDUrlRunScript{time}\"", () => replaceCode(page, el("code editor"), session.text("//name: BDDUrlRunScript{time}")));
      await session.step(38, "And user appends \"//language: javascript\" to code editor", () => appendToEditor(page, "//language: javascript", el("code editor")));
      await session.step(39, "And user appends \"//input: int rows = 10\" to code editor", () => appendToEditor(page, "//input: int rows = 10", el("code editor")));
      await session.step(40, "And user appends \"//output: dataframe df\" to code editor", () => appendToEditor(page, "//output: dataframe df", el("code editor")));
      await session.step(41, "And user appends \"df = grok.data.demo.demog(rows);\" to code editor", () => appendToEditor(page, "df = grok.data.demo.demog(rows);", el("code editor")));
      await session.step(42, "And user saves the script", () => saveScript(page));
      await session.step(43, "Then 1 script named \"BDDUrlRunScript{time}\" should be on the server", () => scriptsOnServer(page, 1, session.text("BDDUrlRunScript{time}")));
    });
    await run.scenario("The script's Links... dialog gives its Grok name", async () => {
      await session.step(46, "When user closes the current view", () => closeCurrentView(page));
      await session.step(47, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
      await session.step(48, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(49, "And user types \"BDDUrlRunScript{time}\" into gallery search", () => typeInto(page, session.text("BDDUrlRunScript{time}"), el("gallery search")));
      await session.step(50, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(51, "And user clicks on \"BDDUrlRunScript{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDDUrlRunScript{time}\" link in gallery"))));
      await session.step(52, "Then the context panel should show \"BDDUrlRunScript{time}\"", () => contextPanelShows(page, session.text("BDDUrlRunScript{time}")));
      await session.step(53, "When user clicks on \"Links...\" link in Details pane in context panel", () => clickOn(page, el("\"Links...\" link in Details pane in context panel")));
      await session.step(54, "Then \"Links to BDDUrlRunScript{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Links to BDDUrlRunScript{time}\" dialog")), "visible"));
      await session.step(55, "And \"Grok name\" input in \"Links to BDDUrlRunScript{time}\" dialog should have value \"Admin:BDDUrlRunScript{time}\"", () => shouldHaveValue(page, el(session.text("\"Grok name\" input in \"Links to BDDUrlRunScript{time}\" dialog")), session.text("Admin:BDDUrlRunScript{time}")));
      await session.step(56, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(57, "Then the \"Links to BDDUrlRunScript{time}\" dialog should close", () => dialogCloses(page, session.text("Links to BDDUrlRunScript{time}")));
      await session.step(59, "When user clears gallery search", () => clearField(page, el("gallery search")));
    });
    await run.scenario("Without run the form opens with the URL's value", async () => {
      await session.step(62, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(63, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(64, "When user opens the address \"/func/Admin.BDDUrlRunScript{time}?rows=25\"", () => openAddress(page, session.text("/func/Admin.BDDUrlRunScript{time}?rows=25")));
      await session.step(65, "Then \"Rows\" input should have value \"25\"", () => shouldHaveValue(page, el("\"Rows\" input"), "25"));
      await session.step(66, "And RUN button should be visible", () => shouldBe(page, el("RUN button"), "visible"));
      await session.step(67, "And status bar should contain text \"Rows: 0\"", () => shouldContainText(page, el("status bar"), "Rows: 0"));
    });
    await run.scenario("The copy icon gives a run link for the current value", async () => {
      await session.step(70, "When user enters \"40\" into \"Rows\" input", () => enterInto(page, "40", el("\"Rows\" input")));
      await session.step(71, "And user hovers over copy icon", () => hoverOver(page, el("copy icon")));
      await session.step(72, "Then tooltip should contain text \"Copy a link that runs this function with the current parameters\"", () => shouldContainText(page, el("tooltip"), "Copy a link that runs this function with the current parameters"));
      await session.step(73, "When user clicks on copy icon", () => clickOn(page, el("copy icon")));
      await session.step(74, "Then the clipboard should contain text \"/func/Admin.BDDUrlRunScript{time}?rows=40&run=true\"", () => clipboardContains(page, session.text("/func/Admin.BDDUrlRunScript{time}?rows=40&run=true")));
    });
    await run.scenario("The run link runs at once and shows only the result", async () => {
      await session.step(77, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(78, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(79, "When user opens the address \"/func/Admin.BDDUrlRunScript{time}?rows=40&run=true\"", () => openAddress(page, session.text("/func/Admin.BDDUrlRunScript{time}?rows=40&run=true")));
      await session.step(80, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(81, "And the table should have 40 rows", () => rowCount(page, 40));
      await session.step(82, "And status bar should contain text \"Rows: 40\"", () => shouldContainText(page, el("status bar"), "Rows: 40"));
      await session.step(83, "And the page address should contain \"run=true\"", () => urlShouldContain(page, "run=true"));
      await session.step(84, "And RUN button should be absent", () => shouldBe(page, el("RUN button"), "absent"));
      await session.step(85, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(86, "Then \"Rows\" input in Source pane in toolbox should have value \"40\"", () => shouldHaveValue(page, el("\"Rows\" input in Source pane in toolbox"), "40"));
      await session.step(87, "And REFRESH button in Source pane in toolbox should be visible", () => shouldBe(page, el("REFRESH button in Source pane in toolbox"), "visible"));
    });
    await run.scenario("A value changed in the result view is refreshed and carried by the link", async () => {
      await session.step(90, "When user enters \"15\" into \"Rows\" input in Source pane in toolbox", () => enterInto(page, "15", el("\"Rows\" input in Source pane in toolbox")));
      await session.step(91, "And user clicks on REFRESH button in Source pane in toolbox", () => clickOn(page, el("REFRESH button in Source pane in toolbox")));
      await session.step(92, "Then the \"rows\" reading of grid should be 15", () => readingIs(page, "rows", el("grid"), 15));
      await session.step(93, "And status bar should contain text \"Rows: 15\"", () => shouldContainText(page, el("status bar"), "Rows: 15"));
      await session.step(94, "When user hovers over copy icon in Source pane in toolbox", () => hoverOver(page, el("copy icon in Source pane in toolbox")));
      await session.step(95, "Then tooltip should contain text \"?rows=15&run=true\"", () => shouldContainText(page, el("tooltip"), "?rows=15&run=true"));
    });
    await run.scenario("run=false keeps the form", async () => {
      await session.step(98, "When user moves the pointer away from Source pane in toolbox", () => pointerAway(page, el("Source pane in toolbox")));
      await session.step(99, "And user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(100, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(101, "When user opens the address \"/func/Admin.BDDUrlRunScript{time}?rows=25&run=false\"", () => openAddress(page, session.text("/func/Admin.BDDUrlRunScript{time}?rows=25&run=false")));
      await session.step(102, "Then \"Rows\" input should have value \"25\"", () => shouldHaveValue(page, el("\"Rows\" input"), "25"));
      await session.step(103, "And RUN button should be visible", () => shouldBe(page, el("RUN button"), "visible"));
      await session.step(104, "And status bar should contain text \"Rows: 0\"", () => shouldContainText(page, el("status bar"), "Rows: 0"));
    });
    run.finish();
  });
});
