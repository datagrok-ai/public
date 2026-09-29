/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-lifecycle-script.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, sharing.share-dialog]
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
import {saveScript} from '../../bindings/scripts.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {appendToEditor, clearField, clickOn, enterInto, isExpanded, replaceCode, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeCurrentView, contextPanelShows, currentViewType, dialogCloses, noProjectOnServer, noScriptOnServer, pickSharingUser, projectsOnServer, savedWithDataSync, scriptContains, scriptsOnServer, scriptsView, sharingPaneLists, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project built from the user's own script, shared, then the script renamed and broken", () => {
  const session = feature(test, "features/projects/projects-lifecycle-script.feature", import.meta.url);
  test("A project built from the user's own script, shared, then the script renamed and broken", {tag: ["@journey", "@serial", "@realizes:views.projects", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And no project named \"BDDLifeScriptProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDLifeScriptProj{time}")));
    await session.step(27, "And no script named \"BDDLifeScript{time}\" is on the server", () => noScriptOnServer(page, session.text("BDDLifeScript{time}")));
    await session.step(28, "And no script named \"BDDLifeScriptRenamed{time}\" is on the server", () => noScriptOnServer(page, session.text("BDDLifeScriptRenamed{time}")));
    await run.scenario("A JavaScript script is written in a new editor and saved", async () => {
      await session.step(31, "Given Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
      await session.step(32, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
      await session.step(33, "When user clicks on Platform---Functions---Scripts tree node inside browse tree", () => clickOn(page, el("Platform---Functions---Scripts tree node inside browse tree")));
      await session.step(34, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
      await session.step(35, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(36, "And user picks \"JavaScript Script...\" from the open menu", () => pickFromOpenMenu(page, "JavaScript Script..."));
      await session.step(37, "Then the \"Template\" view should be current", () => viewIsCurrent(page, "Template"));
      await session.step(38, "And code editor should contain the text \"Hello World\"", () => shouldContainText(page, el("code editor"), "Hello World"));
      await session.step(39, "When user replaces the code of code editor with \"//name: BDDLifeScript{time}\"", () => replaceCode(page, el("code editor"), session.text("//name: BDDLifeScript{time}")));
      await session.step(40, "And user appends \"//language: javascript\" to code editor", () => appendToEditor(page, "//language: javascript", el("code editor")));
      await session.step(41, "And user appends \"//output: dataframe df\" to code editor", () => appendToEditor(page, "//output: dataframe df", el("code editor")));
      await session.step(42, "And user appends \"df = await grok.data.getDemoTable('demog.csv');\" to code editor", () => appendToEditor(page, "df = await grok.data.getDemoTable('demog.csv');", el("code editor")));
      await session.step(43, "And user saves the script", () => saveScript(page));
      await session.step(44, "Then an info balloon containing \"Script saved.\" should have been shown", () => infoBalloonText(page, "Script saved."));
      await session.step(45, "And 1 script named \"BDDLifeScript{time}\" should be on the server", () => scriptsOnServer(page, 1, session.text("BDDLifeScript{time}")));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Run... from the Scripts gallery opens the script's table", async () => {
      await session.step(49, "When user closes the current view", () => closeCurrentView(page));
      await session.step(50, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
      await session.step(51, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(52, "And user types \"BDDLifeScript{time}\" into gallery search", () => typeInto(page, session.text("BDDLifeScript{time}"), el("gallery search")));
      await session.step(53, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(54, "And user picks \"Run...\" from the context menu of \"BDDLifeScript{time}\" link in gallery", () => pickFromContextMenu(page, "Run...", el(session.text("\"BDDLifeScript{time}\" link in gallery"))));
      await session.step(55, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(56, "And the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(57, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
      await session.step(59, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Saved with Data sync, the table's creation script calls the script", async () => {
      await session.step(62, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(63, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(64, "And \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "visible"));
      await session.step(65, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(66, "And \"Save project\" dialog should contain text \"Some tables require this script for data sync.\"", () => shouldContainText(page, el("\"Save project\" dialog"), "Some tables require this script for data sync."));
      await session.step(67, "When user clicks on \"Creation script\" button in \"demog\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog")));
      await session.step(68, "Then \"demog\" project table in \"Save project\" dialog should contain text \":BDDLifeScript{time}()\"", () => shouldContainText(page, el("\"demog\" project table in \"Save project\" dialog"), session.text(":BDDLifeScript{time}()")));
      await session.step(69, "When user enters \"BDDLifeScriptProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDLifeScriptProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(70, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(71, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(72, "And an info balloon containing 'Project \"BDDLifeScriptProj{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDLifeScriptProj{time}\" uploaded")));
      await session.step(73, "And 1 project named \"BDDLifeScriptProj{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDLifeScriptProj{time}")));
      await session.step(74, "And the \"demog\" table of the \"BDDLifeScriptProj{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog", session.text("BDDLifeScriptProj{time}")));
      await session.step(75, "And \"Share BDDLifeScriptProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDLifeScriptProj{time}\" dialog")), "visible"));
      await session.step(76, "When user clicks on CANCEL button in \"Share BDDLifeScriptProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDLifeScriptProj{time}\" dialog"))));
      await session.step(77, "Then the \"Share BDDLifeScriptProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDLifeScriptProj{time}")));
      await session.step(78, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(79, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
      await session.step(81, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Only the project is shared with the second account", async () => {
      await session.step(84, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(85, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(86, "And user enters \"BDDLifeScript\" into gallery search", () => enterInto(page, "BDDLifeScript", el("gallery search")));
      await session.step(87, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(88, "And user picks \"Share...\" from the context menu of BDDLifeScriptProj{time} gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("BDDLifeScriptProj{time} gallery card"))));
      await session.step(89, "Then \"Share BDDLifeScriptProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDLifeScriptProj{time}\" dialog")), "visible"));
      await session.step(90, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(91, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDLifeScriptProj{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDLifeScriptProj{time}\" dialog"))));
      await session.step(92, "And user clicks on OK button in \"Share BDDLifeScriptProj{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDLifeScriptProj{time}\" dialog"))));
      await session.step(93, "Then the \"Share BDDLifeScriptProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDLifeScriptProj{time}")));
      await session.step(94, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(95, "When user clicks on BDDLifeScriptProj{time} gallery card", () => clickOn(page, el(session.text("BDDLifeScriptProj{time} gallery card"))));
      await session.step(96, "Then the context panel should show \"BDDLifeScriptProj{time}\"", () => contextPanelShows(page, session.text("BDDLifeScriptProj{time}")));
      await session.step(97, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
    });
    await run.scenario("The script is renamed in its editor with a new body", async () => {
      await session.step(100, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(101, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(102, "And user opens the Scripts view", () => scriptsView(page));
      await session.step(103, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(104, "And user types \"BDDLifeScript{time}\" into gallery search", () => typeInto(page, session.text("BDDLifeScript{time}"), el("gallery search")));
      await session.step(105, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(106, "And user picks \"Edit...\" from the context menu of \"BDDLifeScript{time}\" link in gallery", () => pickFromContextMenu(page, "Edit...", el(session.text("\"BDDLifeScript{time}\" link in gallery"))));
      await session.step(107, "Then the \"BDDLifeScript{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDLifeScript{time}")));
      await session.step(108, "When user replaces the code of code editor with \"//name: BDDLifeScriptRenamed{time}\"", () => replaceCode(page, el("code editor"), session.text("//name: BDDLifeScriptRenamed{time}")));
      await session.step(109, "And user appends \"//language: javascript\" to code editor", () => appendToEditor(page, "//language: javascript", el("code editor")));
      await session.step(110, "And user appends \"//output: dataframe df\" to code editor", () => appendToEditor(page, "//output: dataframe df", el("code editor")));
      await session.step(111, "And user appends \"df = grok.data.demo.demog(100);\" to code editor", () => appendToEditor(page, "df = grok.data.demo.demog(100);", el("code editor")));
      await session.step(112, "And user saves the script", () => saveScript(page));
      await session.step(113, "Then 1 script named \"BDDLifeScriptRenamed{time}\" should be on the server", () => scriptsOnServer(page, 1, session.text("BDDLifeScriptRenamed{time}")));
      await session.step(114, "And 0 scripts named \"BDDLifeScript{time}\" should be on the server", () => scriptsOnServer(page, 0, session.text("BDDLifeScript{time}")));
      await session.step(115, "And the script \"BDDLifeScriptRenamed{time}\" on the server should contain \"grok.data.demo.demog(100)\"", () => scriptContains(page, session.text("BDDLifeScriptRenamed{time}"), "grok.data.demo.demog(100)"));
      await session.step(116, "And no errors should have been logged", () => noErrors(page));
      await session.step(117, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The script is broken with a throw before its df line", async () => {
      await session.step(120, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(121, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(122, "And user opens the Scripts view", () => scriptsView(page));
      await session.step(123, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(124, "And user types \"BDDLifeScriptRenamed{time}\" into gallery search", () => typeInto(page, session.text("BDDLifeScriptRenamed{time}"), el("gallery search")));
      await session.step(125, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(126, "And user picks \"Edit...\" from the context menu of \"BDDLifeScriptRenamed{time}\" link in gallery", () => pickFromContextMenu(page, "Edit...", el(session.text("\"BDDLifeScriptRenamed{time}\" link in gallery"))));
      await session.step(127, "Then the \"BDDLifeScriptRenamed{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDLifeScriptRenamed{time}")));
      await session.step(128, "When user replaces the code of code editor with \"//name: BDDLifeScriptRenamed{time}\"", () => replaceCode(page, el("code editor"), session.text("//name: BDDLifeScriptRenamed{time}")));
      await session.step(129, "And user appends \"//language: javascript\" to code editor", () => appendToEditor(page, "//language: javascript", el("code editor")));
      await session.step(130, "And user appends \"//output: dataframe df\" to code editor", () => appendToEditor(page, "//output: dataframe df", el("code editor")));
      await session.step(131, "And user appends \"throw new Error('intentional break');\" to code editor", () => appendToEditor(page, "throw new Error('intentional break');", el("code editor")));
      await session.step(132, "And user appends \"df = grok.data.demo.demog(100);\" to code editor", () => appendToEditor(page, "df = grok.data.demo.demog(100);", el("code editor")));
      await session.step(133, "And user saves the script", () => saveScript(page));
      await session.step(134, "Then the script \"BDDLifeScriptRenamed{time}\" on the server should contain \"throw new Error('intentional break');\"", () => scriptContains(page, session.text("BDDLifeScriptRenamed{time}"), "throw new Error('intentional break');"));
      await session.step(135, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Delete Project and the script's Delete remove both", async () => {
      await session.step(138, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(139, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(140, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(141, "And user enters \"BDDLifeScriptProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeScriptProj{time}"), el("gallery search")));
      await session.step(142, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(143, "And user picks \"Delete Project\" from the context menu of BDDLifeScriptProj{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDLifeScriptProj{time} gallery card"))));
      await session.step(144, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(145, "And \"Are you sure?\" dialog should contain text 'Delete project \"BDDLifeScriptProj{time}\"?'", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete project \"BDDLifeScriptProj{time}\"?")));
      await session.step(146, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(148, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(149, "And 0 projects named \"BDDLifeScriptProj{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDLifeScriptProj{time}")));
      await session.step(150, "Given user opens the Scripts view", () => scriptsView(page));
      await session.step(151, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(152, "And user types \"BDDLifeScriptRenamed{time}\" into gallery search", () => typeInto(page, session.text("BDDLifeScriptRenamed{time}"), el("gallery search")));
      await session.step(153, "And user picks \"Delete\" from the context menu of \"BDDLifeScriptRenamed{time}\" link in gallery", () => pickFromContextMenu(page, "Delete", el(session.text("\"BDDLifeScriptRenamed{time}\" link in gallery"))));
      await session.step(154, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(155, "And \"Are you sure?\" dialog should contain text 'Delete script \"BDDLifeScriptRenamed{time}\"?'", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete script \"BDDLifeScriptRenamed{time}\"?")));
      await session.step(156, "When user clicks on YES button in \"Are you sure?\" dialog", () => clickOn(page, el("YES button in \"Are you sure?\" dialog")));
      await session.step(157, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(158, "And 0 scripts named \"BDDLifeScriptRenamed{time}\" should be on the server", () => scriptsOnServer(page, 0, session.text("BDDLifeScriptRenamed{time}")));
      await session.step(159, "And \"BDDLifeScriptRenamed{time}\" link in gallery should be absent", () => shouldBe(page, el(session.text("\"BDDLifeScriptRenamed{time}\" link in gallery")), "absent"));
      await session.step(160, "And no errors should have been logged", () => noErrors(page));
      await session.step(161, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(163, "When user clears gallery search", () => clearField(page, el("gallery search")));
    });
    run.finish();
  });
});
