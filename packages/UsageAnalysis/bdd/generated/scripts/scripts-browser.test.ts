/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scripts/scripts-browser.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.scripts]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, followingShouldBe, pressKeyIn, selectIn, shouldBe, shouldContainText, shouldHaveText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeCurrentView, consoleShows, contextPanelOpen, contextPanelShows, dialogCloses, galleryMode, loadTable, noteConsole, openDataset, paneCountAtLeast, pickSharingUser, scriptOnServer, scriptsView, sharingPaneLists, sharingPaneListsNot, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors, pickFromContextMenu, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Scripts view and a script's context panel", () => {
  const session = feature(test, "features/scripts/scripts-browser.feature", import.meta.url);
  test("The Scripts view and a script's context panel", {tag: ["@journey", "@full-stand", "@serial", "@realizes:views.scripts"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And a script \"BddScriptBrowser{time}\" is on the server:", () => scriptOnServer(page, session.text("BddScriptBrowser{time}"), "#language: r\n#input: dataframe table\n#output: int count\n#output: string newParam\ncount <- nrow(table) * ncol(table)\nnewParam <- \"test\""));
    await session.step(39, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(40, "And user opens the Scripts view", () => scriptsView(page));
    await run.scenario("The view modes and the order of the gallery", async () => {
      await session.step(43, "When user clicks on \"Switch to brief view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to brief view\" icon inside gallery toolbar")));
      await session.step(44, "Then the gallery should be in brief mode", () => galleryMode(page, "brief"));
      await session.step(45, "When user clicks on \"Switch to grid view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to grid view\" icon inside gallery toolbar")));
      await session.step(46, "Then the gallery should be in grid mode", () => galleryMode(page, "grid"));
      await session.step(47, "When user clicks on \"Switch to card view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to card view\" icon inside gallery toolbar")));
      await session.step(48, "Then the gallery should be in card mode", () => galleryMode(page, "card"));
      await session.step(49, "Then no errors should have been logged", () => noErrors(page));
      await session.step(50, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A search narrows the gallery to the script", async () => {
      await session.step(53, "When user types \"BddScriptBrowser{time}\" into gallery search", () => typeInto(page, session.text("BddScriptBrowser{time}"), el("gallery search")));
      await session.step(54, "Then gallery counter should have text \"1\"", () => shouldHaveText(page, el("gallery counter"), "1"));
      await session.step(55, "And \"BddScriptBrowser{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BddScriptBrowser{time}\" link in gallery")), "visible"));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The script fills the context panel", async () => {
      await session.step(59, "When user clicks on \"BddScriptBrowser{time}\" link in gallery", () => clickOn(page, el(session.text("\"BddScriptBrowser{time}\" link in gallery"))));
      await session.step(60, "Then the context panel should show \"BddScriptBrowser{time}\"", () => contextPanelShows(page, session.text("BddScriptBrowser{time}")));
      await session.step(61, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Details\" section in context panel"],["\"Script\" section in context panel"],["\"Run\" section in context panel"],["\"Sharing\" section in context panel"],["\"Chats\" section in context panel"],["\"Activity\" section in context panel"]]), [["\"Details\" section in context panel"],["\"Script\" section in context panel"],["\"Run\" section in context panel"],["\"Sharing\" section in context panel"],["\"Chats\" section in context panel"],["\"Activity\" section in context panel"]]);
      await session.step(68, "When user clicks on \"Details\" pane header in context panel", () => clickOn(page, el("\"Details\" pane header in context panel")));
      await session.step(69, "Then \"Details\" section in context panel should contain text \"count, newParam\"", () => shouldContainText(page, el("\"Details\" section in context panel"), "count, newParam"));
      await session.step(70, "And \"Details\" section in context panel should contain text \"table\"", () => shouldContainText(page, el("\"Details\" section in context panel"), "table"));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Activity counts the creation and the run", async () => {
      await session.step(74, "Given user opens cars dataset", () => openDataset(page, ds("cars")));
      await session.step(75, "And user opens the Scripts view", () => scriptsView(page));
      await session.step(76, "When user types \"BddScriptBrowser{time}\" into gallery search", () => typeInto(page, session.text("BddScriptBrowser{time}"), el("gallery search")));
      await session.step(77, "And user picks \"Run...\" from the context menu of \"BddScriptBrowser{time}\" link in gallery", () => pickFromContextMenu(page, "Run...", el(session.text("\"BddScriptBrowser{time}\" link in gallery"))));
      await session.step(78, "And user selects \"cars\" in Table input in \"BddScriptBrowser{time}\" dialog", () => selectIn(page, "cars", el(session.text("Table input in \"BddScriptBrowser{time}\" dialog"))));
      await session.step(79, "And user clicks on OK button in \"BddScriptBrowser{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"BddScriptBrowser{time}\" dialog"))));
      await session.step(80, "Then the \"BddScriptBrowser{time}\" dialog should close", () => dialogCloses(page, session.text("BddScriptBrowser{time}")));
      await session.step(82, "When user closes the current view", () => closeCurrentView(page));
      await session.step(83, "And user opens the Scripts view", () => scriptsView(page));
      await session.step(84, "And user types \"BddScriptBrowser{time}\" into gallery search", () => typeInto(page, session.text("BddScriptBrowser{time}"), el("gallery search")));
      await session.step(85, "And user clicks on \"BddScriptBrowser{time}\" link in gallery", () => clickOn(page, el(session.text("\"BddScriptBrowser{time}\" link in gallery"))));
      await session.step(86, "Then the context panel should show \"BddScriptBrowser{time}\"", () => contextPanelShows(page, session.text("BddScriptBrowser{time}")));
      await session.step(87, "And the \"Activity\" pane of the context panel should count at least 1", () => paneCountAtLeast(page, "Activity", 1));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
      await session.step(89, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Share... gives the second account access", async () => {
      await session.step(92, "Then the sharing pane should not list the sharing user", () => sharingPaneListsNot(page));
      await session.step(93, "When user picks \"Share...\" from the context menu of \"BddScriptBrowser{time}\" link in gallery", () => pickFromContextMenu(page, "Share...", el(session.text("\"BddScriptBrowser{time}\" link in gallery"))));
      await session.step(94, "Then \"Share BddScriptBrowser{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BddScriptBrowser{time}\" dialog")), "visible"));
      await session.step(96, "And \"Share BddScriptBrowser{time}\" dialog should contain text \"Full access\"", () => shouldContainText(page, el(session.text("\"Share BddScriptBrowser{time}\" dialog")), "Full access"));
      await session.step(97, "When user picks the sharing user in \"User, group, or email\" input in \"Share BddScriptBrowser{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BddScriptBrowser{time}\" dialog"))));
      await session.step(98, "And user clicks on OK button in \"Share BddScriptBrowser{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BddScriptBrowser{time}\" dialog"))));
      await session.step(99, "Then the \"Share BddScriptBrowser{time}\" dialog should close", () => dialogCloses(page, session.text("Share BddScriptBrowser{time}")));
      await session.step(100, "When user clicks on \"BddScriptBrowser{time}\" link in gallery", () => clickOn(page, el(session.text("\"BddScriptBrowser{time}\" link in gallery"))));
      await session.step(101, "Then the sharing pane should list the sharing user", () => sharingPaneLists(page));
      await session.step(102, "And no errors should have been logged", () => noErrors(page));
      await session.step(103, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A chat message stays on the script", async () => {
      await session.step(106, "When user clicks on \"Chats\" pane header in context panel", () => clickOn(page, el("\"Chats\" pane header in context panel")));
      await session.step(107, "When user types \"bdd chat {time}\" into chat input in context panel", () => typeInto(page, session.text("bdd chat {time}"), el("chat input in context panel")));
      await session.step(108, "And user presses Enter in chat input in context panel", () => pressKeyIn(page, "Enter", el("chat input in context panel")));
      await session.step(109, "Then \"Chats\" section in context panel should contain text \"bdd chat {time}\"", () => shouldContainText(page, el("\"Chats\" section in context panel"), session.text("bdd chat {time}")));
      await session.step(110, "And no errors should have been logged", () => noErrors(page));
      await session.step(111, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The ACF sample runs on TSLA's Close", async () => {
      await session.step(114, "Given the \"System:DemoFiles/TSLA.csv\" file is loaded as a table", () => loadTable(page, "System:DemoFiles/TSLA.csv"));
      await session.step(115, "And user opens the Scripts view", () => scriptsView(page));
      await session.step(116, "When user types \"ACF\" into gallery search", () => typeInto(page, "ACF", el("gallery search")));
      await session.step(117, "And user notes the console output", () => noteConsole(page));
      await session.step(118, "And user picks \"Run...\" from the context menu of \"ACF\" link in gallery", () => pickFromContextMenu(page, "Run...", el("\"ACF\" link in gallery")));
      await session.step(119, "Then \"ACF\" dialog should be visible", () => shouldBe(page, el("\"ACF\" dialog"), "visible"));
      await session.step(120, "When user selects \"TSLA\" in Data input in \"ACF\" dialog", () => selectIn(page, "TSLA", el("Data input in \"ACF\" dialog")));
      await session.step(121, "And user clicks on editor of Columns input in \"ACF\" dialog", () => clickOn(page, el("editor of Columns input in \"ACF\" dialog")));
      await session.step(122, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(123, "And the \"text of cell 5 of __name\" reading of grid viewer in \"Select columns...\" dialog should be \"Close\"", () => readingReads(page, "text of cell 5 of __name", el("grid viewer in \"Select columns...\" dialog"), "Close"));
      await session.step(124, "When user clicks on the \"cell 5 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 5 of x", el("grid viewer in \"Select columns...\" dialog")));
      await session.step(125, "Then \"Select columns...\" dialog should contain text \"1 checked\"", () => shouldContainText(page, el("\"Select columns...\" dialog"), "1 checked"));
      await session.step(126, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
      await session.step(127, "Then editor of Columns input in \"ACF\" dialog should contain text \"(1)\"", () => shouldContainText(page, el("editor of Columns input in \"ACF\" dialog"), "(1)"));
      await session.step(128, "When user clicks on OK button in \"ACF\" dialog", () => clickOn(page, el("OK button in \"ACF\" dialog")));
      await session.step(129, "Then the \"ACF\" dialog should close", () => dialogCloses(page, "ACF"));
      await session.step(130, "And the console should show \"ACF(\"", () => consoleShows(page, "ACF("));
      await session.step(131, "And no errors should have been logged", () => noErrors(page));
      await session.step(132, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Edit... opens the ACF sample in the editor", async () => {
      await session.step(135, "Given user opens the Scripts view", () => scriptsView(page));
      await session.step(136, "When user types \"ACF\" into gallery search", () => typeInto(page, "ACF", el("gallery search")));
      await session.step(137, "And user picks \"Edit...\" from the context menu of \"ACF\" link in gallery", () => pickFromContextMenu(page, "Edit...", el("\"ACF\" link in gallery")));
      await session.step(138, "Then the \"ACF\" view should be current", () => viewIsCurrent(page, "ACF"));
      await session.step(139, "And code editor should contain the text \"#name: ACF\"", () => shouldContainText(page, el("code editor"), "#name: ACF"));
      await session.step(140, "And no errors should have been logged", () => noErrors(page));
      await session.step(141, "When user closes the current view", () => closeCurrentView(page));
      await session.step(142, "And user opens the Scripts view", () => scriptsView(page));
    });
    run.finish();
  });
});
