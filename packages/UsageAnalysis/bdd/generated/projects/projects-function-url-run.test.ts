/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-function-url-run.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.functions, views.projects, core#446]
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
import {clearField, clickOn, clipboardContains, enterInto, hoverOver, pressKey, shouldBe, shouldContainText, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, currentViewType, dialogCloses, scriptOnServer, toolboxPaneShown, urlShouldContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu, pointerAway, readingIs, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A function run straight from its URL", () => {
  const session = feature(test, "features/projects/projects-function-url-run.feature", import.meta.url);
  test("A function run straight from its URL", {tag: ["@journey", "@realizes:views.functions", "@realizes:views.projects", "@realizes:core#446"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(27, "And a script \"BDDUrlRunScript{time}\" is on the server:", () => scriptOnServer(page, session.text("BDDUrlRunScript{time}"), "//language: javascript\n//input: int rows = 10\n//output: dataframe df\ndf = grok.data.demo.demog(rows);"));
    await session.step(34, "And a script \"BDDUrlRunRequired{time}\" is on the server:", () => scriptOnServer(page, session.text("BDDUrlRunRequired{time}"), "//language: javascript\n//input: int rows\n//input: string label = \"x\"\n//output: dataframe df\ndf = grok.data.demo.demog(rows);"));
    await run.scenario("The script's Links... dialog gives its Grok name", async () => {
      await session.step(44, "When user opens the address \"/scripts\"", () => openAddress(page, "/scripts"));
      await session.step(45, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
      await session.step(46, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(47, "And user types \"BDDUrlRunScript{time}\" into gallery search", () => typeInto(page, session.text("BDDUrlRunScript{time}"), el("gallery search")));
      await session.step(48, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(49, "And user clicks on \"BDDUrlRunScript{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDDUrlRunScript{time}\" link in gallery"))));
      await session.step(50, "Then the context panel should show \"BDDUrlRunScript{time}\"", () => contextPanelShows(page, session.text("BDDUrlRunScript{time}")));
      await session.step(51, "When user clicks on \"Links...\" link in Details pane in context panel", () => clickOn(page, el("\"Links...\" link in Details pane in context panel")));
      await session.step(52, "Then \"Links to BDDUrlRunScript{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Links to BDDUrlRunScript{time}\" dialog")), "visible"));
      await session.step(53, "And \"Grok name\" input in \"Links to BDDUrlRunScript{time}\" dialog should have value \"Admin:BDDUrlRunScript{time}\"", () => shouldHaveValue(page, el(session.text("\"Grok name\" input in \"Links to BDDUrlRunScript{time}\" dialog")), session.text("Admin:BDDUrlRunScript{time}")));
      await session.step(54, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(55, "Then the \"Links to BDDUrlRunScript{time}\" dialog should close", () => dialogCloses(page, session.text("Links to BDDUrlRunScript{time}")));
      await session.step(57, "When user clears gallery search", () => clearField(page, el("gallery search")));
    });
    await run.scenario("Without run the form opens with the URL's value", async () => {
      await session.step(60, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(61, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(62, "When user opens the address \"/func/Admin.BDDUrlRunScript{time}?rows=25\"", () => openAddress(page, session.text("/func/Admin.BDDUrlRunScript{time}?rows=25")));
      await session.step(63, "Then \"Rows\" input should have value \"25\"", () => shouldHaveValue(page, el("\"Rows\" input"), "25"));
      await session.step(64, "And RUN button should be visible", () => shouldBe(page, el("RUN button"), "visible"));
      await session.step(65, "And status bar should contain text \"Rows: 0\"", () => shouldContainText(page, el("status bar"), "Rows: 0"));
    });
    await run.scenario("The copy icon gives a run link for the current value", async () => {
      await session.step(68, "When user enters \"40\" into \"Rows\" input", () => enterInto(page, "40", el("\"Rows\" input")));
      await session.step(69, "And user hovers over copy icon", () => hoverOver(page, el("copy icon")));
      await session.step(70, "Then tooltip should contain text \"Copy a link that runs this function with the current parameters\"", () => shouldContainText(page, el("tooltip"), "Copy a link that runs this function with the current parameters"));
      await session.step(71, "When user clicks on copy icon", () => clickOn(page, el("copy icon")));
      await session.step(72, "Then the clipboard should contain text \"/func/Admin.BDDUrlRunScript{time}?rows=40&run=true\"", () => clipboardContains(page, session.text("/func/Admin.BDDUrlRunScript{time}?rows=40&run=true")));
    });
    await run.scenario("The run link runs at once and shows only the result", async () => {
      await session.step(75, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(76, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(77, "When user opens the address \"/func/Admin.BDDUrlRunScript{time}?rows=40&run=true\"", () => openAddress(page, session.text("/func/Admin.BDDUrlRunScript{time}?rows=40&run=true")));
      await session.step(78, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(79, "And the table should have 40 rows", () => rowCount(page, 40));
      await session.step(80, "And status bar should contain text \"Rows: 40\"", () => shouldContainText(page, el("status bar"), "Rows: 40"));
      await session.step(81, "And the page address should contain \"run=true\"", () => urlShouldContain(page, "run=true"));
      await session.step(82, "And RUN button should be absent", () => shouldBe(page, el("RUN button"), "absent"));
      await session.step(83, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(84, "Then \"Rows\" input in Source pane in toolbox should have value \"40\"", () => shouldHaveValue(page, el("\"Rows\" input in Source pane in toolbox"), "40"));
      await session.step(85, "And REFRESH button in Source pane in toolbox should be visible", () => shouldBe(page, el("REFRESH button in Source pane in toolbox"), "visible"));
    });
    await run.scenario("A value changed in the result view is refreshed and carried by the link", async () => {
      await session.step(88, "When user enters \"15\" into \"Rows\" input in Source pane in toolbox", () => enterInto(page, "15", el("\"Rows\" input in Source pane in toolbox")));
      await session.step(89, "And user clicks on REFRESH button in Source pane in toolbox", () => clickOn(page, el("REFRESH button in Source pane in toolbox")));
      await session.step(90, "Then the \"rows\" reading of grid should be 15", () => readingIs(page, "rows", el("grid"), 15));
      await session.step(91, "And status bar should contain text \"Rows: 15\"", () => shouldContainText(page, el("status bar"), "Rows: 15"));
      await session.step(92, "When user hovers over copy icon in Source pane in toolbox", () => hoverOver(page, el("copy icon in Source pane in toolbox")));
      await session.step(93, "Then tooltip should contain text \"?rows=15&run=true\"", () => shouldContainText(page, el("tooltip"), "?rows=15&run=true"));
    });
    await run.scenario("run=false keeps the form", async () => {
      await session.step(96, "When user moves the pointer away from Source pane in toolbox", () => pointerAway(page, el("Source pane in toolbox")));
      await session.step(97, "And user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(98, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(99, "When user opens the address \"/func/Admin.BDDUrlRunScript{time}?rows=25&run=false\"", () => openAddress(page, session.text("/func/Admin.BDDUrlRunScript{time}?rows=25&run=false")));
      await session.step(100, "Then \"Rows\" input should have value \"25\"", () => shouldHaveValue(page, el("\"Rows\" input"), "25"));
      await session.step(101, "And RUN button should be visible", () => shouldBe(page, el("RUN button"), "visible"));
      await session.step(102, "And status bar should contain text \"Rows: 0\"", () => shouldContainText(page, el("status bar"), "Rows: 0"));
    });
    await run.scenario("The run link without a required value falls back to the form and names it (core#446)", async () => {
      await session.step(106, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(107, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(108, "When user opens the address \"/func/Admin.BDDUrlRunRequired{time}?label=%22y%22&run=true\"", () => openAddress(page, session.text("/func/Admin.BDDUrlRunRequired{time}?label=%22y%22&run=true")));
      await session.step(109, "Then the \"BDDUrlRunRequired{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDUrlRunRequired{time}")));
      await session.step(110, "And a warning balloon containing 'Unable to run \"BDDUrlRunRequired{time}\": rows' should have been shown", () => warningBalloonText(page, session.text("Unable to run \"BDDUrlRunRequired{time}\": rows")));
      await session.step(111, "And \"Label\" input should have value \"y\"", () => shouldHaveValue(page, el("\"Label\" input"), "y"));
      await session.step(112, "And RUN button should be visible", () => shouldBe(page, el("RUN button"), "visible"));
      await session.step(113, "And status bar should contain text \"Rows: 0\"", () => shouldContainText(page, el("status bar"), "Rows: 0"));
      await session.step(114, "When user enters \"5\" into \"Rows\" input", () => enterInto(page, "5", el("\"Rows\" input")));
      await session.step(115, "And user clicks on RUN button", () => clickOn(page, el("RUN button")));
      await session.step(116, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(117, "And the table should have 5 rows", () => rowCount(page, 5));
    });
    run.finish();
  });
});
