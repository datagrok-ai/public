/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-lifecycle-script.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, sharing.share-dialog]
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
import {loggedIn, reloadPage} from '@datagrok-libraries/bdd/bindings/common/session';
import {appendToEditor, clearField, clickOn, doubleClickOn, enterInto, replaceCode, shouldBe, shouldContainText, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeCurrentView, contextPanelOpen, contextPanelShows, currentViewType, dialogCloses, doubleClickOpensProject, noProjectOnServer, noScriptOnServer, pickSharingUser, projectsOnServer, reloadedByDataSync, runningAccountSignedIn, saveScript, savedWithDataSync, scriptContains, scriptsOnServer, sharingPaneLists, sharingUserSignedIn, signInAsSelf, signInAsSharingUser, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project built from the user's own script, shared, then the script renamed and broken", () => {
  const session = feature(test, "features/projects/projects-lifecycle-script.feature", import.meta.url);
  test("A project built from the user's own script, shared, then the script renamed and broken", {tag: ["@journey", "@realizes:views.projects", "@realizes:sharing.share-dialog"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And no project named \"BDDLifeScriptProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDLifeScriptProj{time}")));
    await session.step(27, "And no script named \"BDDLifeScript{time}\" is on the server", () => noScriptOnServer(page, session.text("BDDLifeScript{time}")));
    await session.step(28, "And no script named \"BDDLifeScriptRenamed{time}\" is on the server", () => noScriptOnServer(page, session.text("BDDLifeScriptRenamed{time}")));
    await run.scenario("A JavaScript script is written in a new editor and saved", async () => {
      await session.step(31, "When user opens the address \"/scripts\"", () => openAddress(page, "/scripts"));
      await session.step(32, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
      await session.step(33, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(34, "And user picks \"JavaScript Script...\" from the open menu", () => pickFromOpenMenu(page, "JavaScript Script..."));
      await session.step(35, "Then the \"Template\" view should be current", () => viewIsCurrent(page, "Template"));
      await session.step(36, "And code editor should contain the text \"Hello World\"", () => shouldContainText(page, el("code editor"), "Hello World"));
      await session.step(37, "When user replaces the code of code editor with \"//name: BDDLifeScript{time}\"", () => replaceCode(page, el("code editor"), session.text("//name: BDDLifeScript{time}")));
      await session.step(38, "And user appends \"//language: javascript\" to code editor", () => appendToEditor(page, "//language: javascript", el("code editor")));
      await session.step(39, "And user appends \"//output: dataframe df\" to code editor", () => appendToEditor(page, "//output: dataframe df", el("code editor")));
      await session.step(40, "And user appends \"df = await grok.data.getDemoTable('demog.csv');\" to code editor", () => appendToEditor(page, "df = await grok.data.getDemoTable('demog.csv');", el("code editor")));
      await session.step(41, "And user saves the script", () => saveScript(page));
      await session.step(42, "Then an info balloon containing \"Script saved.\" should have been shown", () => infoBalloonText(page, "Script saved."));
      await session.step(43, "And 1 script named \"BDDLifeScript{time}\" should be on the server", () => scriptsOnServer(page, 1, session.text("BDDLifeScript{time}")));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Run... from the Scripts gallery opens the script's table", async () => {
      await session.step(47, "When user closes the current view", () => closeCurrentView(page));
      await session.step(48, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
      await session.step(49, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(50, "And user types \"BDDLifeScript{time}\" into gallery search", () => typeInto(page, session.text("BDDLifeScript{time}"), el("gallery search")));
      await session.step(51, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(52, "And user picks \"Run...\" from the context menu of \"BDDLifeScript{time}\" link in gallery", () => pickFromContextMenu(page, "Run...", el(session.text("\"BDDLifeScript{time}\" link in gallery"))));
      await session.step(53, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(54, "And the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(55, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
      await session.step(57, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Saved with Data sync, the table's creation script calls the script", async () => {
      await session.step(60, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(61, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(62, "And \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "visible"));
      await session.step(63, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(64, "And \"Save project\" dialog should contain text \"Some tables require this script for data sync.\"", () => shouldContainText(page, el("\"Save project\" dialog"), "Some tables require this script for data sync."));
      await session.step(65, "When user clicks on \"Creation script\" button in \"demog\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog")));
      await session.step(66, "Then \"demog\" project table in \"Save project\" dialog should contain text \":BDDLifeScript{time}()\"", () => shouldContainText(page, el("\"demog\" project table in \"Save project\" dialog"), session.text(":BDDLifeScript{time}()")));
      await session.step(67, "When user enters \"BDDLifeScriptProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDLifeScriptProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(68, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(69, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(70, "And an info balloon containing 'Project \"BDDLifeScriptProj{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDLifeScriptProj{time}\" uploaded")));
      await session.step(71, "And 1 project named \"BDDLifeScriptProj{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDLifeScriptProj{time}")));
      await session.step(72, "And the \"demog\" table of the \"BDDLifeScriptProj{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog", session.text("BDDLifeScriptProj{time}")));
      await session.step(73, "And \"Share BDDLifeScriptProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDLifeScriptProj{time}\" dialog")), "visible"));
      await session.step(74, "When user clicks on CANCEL button in \"Share BDDLifeScriptProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDLifeScriptProj{time}\" dialog"))));
      await session.step(75, "Then the \"Share BDDLifeScriptProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDLifeScriptProj{time}")));
      await session.step(76, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(77, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
      await session.step(79, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Only the project is shared with the second account", async () => {
      await session.step(82, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(83, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(84, "And user enters \"BDDLifeScript\" into gallery search", () => enterInto(page, "BDDLifeScript", el("gallery search")));
      await session.step(85, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(86, "And user picks \"Share...\" from the context menu of BDDLifeScriptProj{time} gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("BDDLifeScriptProj{time} gallery card"))));
      await session.step(87, "Then \"Share BDDLifeScriptProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDLifeScriptProj{time}\" dialog")), "visible"));
      await session.step(88, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(89, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDLifeScriptProj{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDLifeScriptProj{time}\" dialog"))));
      await session.step(90, "And user unchecks \"Send notifications\" input in \"Share BDDLifeScriptProj{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDDLifeScriptProj{time}\" dialog"))));
      await session.step(91, "And user clicks on OK button in \"Share BDDLifeScriptProj{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDLifeScriptProj{time}\" dialog"))));
      await session.step(92, "Then the \"Share BDDLifeScriptProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDLifeScriptProj{time}")));
      await session.step(93, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(94, "When user clicks on BDDLifeScriptProj{time} gallery card", () => clickOn(page, el(session.text("BDDLifeScriptProj{time} gallery card"))));
      await session.step(95, "Then the context panel should show \"BDDLifeScriptProj{time}\"", () => contextPanelShows(page, session.text("BDDLifeScriptProj{time}")));
      await session.step(96, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
    });
    await run.scenario("The second account gets the data of the shared project (GROK-19403)", async () => {
      await session.step(99, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(100, "And user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(101, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(102, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(103, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(104, "And user enters \"BDDLifeScriptProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeScriptProj{time}"), el("gallery search")));
      await session.step(105, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(106, "And user double-clicks on BDDLifeScriptProj{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDLifeScriptProj{time} gallery card"))));
      await session.step(107, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(108, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(109, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(110, "And \"Data loading error\" dialog should be absent", () => shouldBe(page, el("\"Data loading error\" dialog"), "absent"));
      await session.step(111, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(112, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(113, "And user signs in as themselves again", () => signInAsSelf(page));
      await session.step(114, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    });
    await run.scenario("The script is renamed in its editor with a new body", async () => {
      await session.step(117, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(118, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(119, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(120, "When user opens the address \"/scripts\"", () => openAddress(page, "/scripts"));
      await session.step(121, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
      await session.step(122, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(123, "And user types \"BDDLifeScript{time}\" into gallery search", () => typeInto(page, session.text("BDDLifeScript{time}"), el("gallery search")));
      await session.step(124, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(125, "And user picks \"Edit...\" from the context menu of \"BDDLifeScript{time}\" link in gallery", () => pickFromContextMenu(page, "Edit...", el(session.text("\"BDDLifeScript{time}\" link in gallery"))));
      await session.step(126, "Then the \"BDDLifeScript{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDLifeScript{time}")));
      await session.step(127, "When user replaces the code of code editor with \"//name: BDDLifeScriptRenamed{time}\"", () => replaceCode(page, el("code editor"), session.text("//name: BDDLifeScriptRenamed{time}")));
      await session.step(128, "And user appends \"//language: javascript\" to code editor", () => appendToEditor(page, "//language: javascript", el("code editor")));
      await session.step(129, "And user appends \"//output: dataframe df\" to code editor", () => appendToEditor(page, "//output: dataframe df", el("code editor")));
      await session.step(130, "And user appends \"df = grok.data.demo.demog(100);\" to code editor", () => appendToEditor(page, "df = grok.data.demo.demog(100);", el("code editor")));
      await session.step(131, "And user saves the script", () => saveScript(page));
      await session.step(132, "Then 1 script named \"BDDLifeScriptRenamed{time}\" should be on the server", () => scriptsOnServer(page, 1, session.text("BDDLifeScriptRenamed{time}")));
      await session.step(133, "And 0 scripts named \"BDDLifeScript{time}\" should be on the server", () => scriptsOnServer(page, 0, session.text("BDDLifeScript{time}")));
      await session.step(134, "And the script \"BDDLifeScriptRenamed{time}\" on the server should contain \"grok.data.demo.demog(100)\"", () => scriptContains(page, session.text("BDDLifeScriptRenamed{time}"), "grok.data.demo.demog(100)"));
      await session.step(135, "And no errors should have been logged", () => noErrors(page));
      await session.step(136, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("After a reload, the project runs the renamed script's new body", async () => {
      await session.step(139, "When user reloads the page", () => reloadPage(page));
      await session.step(140, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(141, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(142, "And user enters \"BDDLifeScriptProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeScriptProj{time}"), el("gallery search")));
      await session.step(143, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(144, "And user double-clicks on BDDLifeScriptProj{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDLifeScriptProj{time} gallery card"))));
      await session.step(145, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(146, "And the table should have 100 rows", () => rowCount(page, 100));
      await session.step(147, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(148, "And \"Data loading error\" dialog should be absent", () => shouldBe(page, el("\"Data loading error\" dialog"), "absent"));
      await session.step(149, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The script is broken with a throw before its df line", async () => {
      await session.step(152, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(153, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(154, "When user opens the address \"/scripts\"", () => openAddress(page, "/scripts"));
      await session.step(155, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
      await session.step(156, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(157, "And user types \"BDDLifeScriptRenamed{time}\" into gallery search", () => typeInto(page, session.text("BDDLifeScriptRenamed{time}"), el("gallery search")));
      await session.step(158, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(159, "And user picks \"Edit...\" from the context menu of \"BDDLifeScriptRenamed{time}\" link in gallery", () => pickFromContextMenu(page, "Edit...", el(session.text("\"BDDLifeScriptRenamed{time}\" link in gallery"))));
      await session.step(160, "Then the \"BDDLifeScriptRenamed{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDLifeScriptRenamed{time}")));
      await session.step(161, "When user replaces the code of code editor with \"//name: BDDLifeScriptRenamed{time}\"", () => replaceCode(page, el("code editor"), session.text("//name: BDDLifeScriptRenamed{time}")));
      await session.step(162, "And user appends \"//language: javascript\" to code editor", () => appendToEditor(page, "//language: javascript", el("code editor")));
      await session.step(163, "And user appends \"//output: dataframe df\" to code editor", () => appendToEditor(page, "//output: dataframe df", el("code editor")));
      await session.step(164, "And user appends \"throw new Error('intentional break');\" to code editor", () => appendToEditor(page, "throw new Error('intentional break');", el("code editor")));
      await session.step(165, "And user appends \"df = grok.data.demo.demog(100);\" to code editor", () => appendToEditor(page, "df = grok.data.demo.demog(100);", el("code editor")));
      await session.step(166, "And user saves the script", () => saveScript(page));
      await session.step(167, "Then the script \"BDDLifeScriptRenamed{time}\" on the server should contain \"throw new Error('intentional break');\"", () => scriptContains(page, session.text("BDDLifeScriptRenamed{time}"), "throw new Error('intentional break');"));
      await session.step(168, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The second account is told to ask the owner and gets no EDIT SCRIPT... (GROK-19728)", async () => {
      await session.step(171, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(172, "And user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(173, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(174, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(175, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(176, "And user enters \"BDDLifeScriptProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeScriptProj{time}"), el("gallery search")));
      await session.step(177, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(178, "And user double-clicks on BDDLifeScriptProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDLifeScriptProj{time} gallery card"))));
      await session.step(179, "Then \"Data loading error\" dialog should be visible", () => shouldBe(page, el("\"Data loading error\" dialog"), "visible"));
      await session.step(180, "And \"Data loading error\" dialog should contain text \"Ask the project owner to fix the script\"", () => shouldContainText(page, el("\"Data loading error\" dialog"), "Ask the project owner to fix the script"));
      await session.step(181, "And \"OPEN ANYWAY\" button in \"Data loading error\" dialog should be visible", () => shouldBe(page, el("\"OPEN ANYWAY\" button in \"Data loading error\" dialog"), "visible"));
      await session.step(182, "And \"CLOSE PROJECT\" button in \"Data loading error\" dialog should be visible", () => shouldBe(page, el("\"CLOSE PROJECT\" button in \"Data loading error\" dialog"), "visible"));
      await session.step(183, "And \"EDIT SCRIPT...\" button in \"Data loading error\" dialog should be absent", () => shouldBe(page, el("\"EDIT SCRIPT...\" button in \"Data loading error\" dialog"), "absent"));
      await session.step(184, "When user clicks on \"CLOSE PROJECT\" button in \"Data loading error\" dialog", () => clickOn(page, el("\"CLOSE PROJECT\" button in \"Data loading error\" dialog")));
      await session.step(185, "Then the \"Data loading error\" dialog should close", () => dialogCloses(page, "Data loading error"));
      await session.step(186, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(187, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    });
    await run.scenario("The owner is offered OPEN ANYWAY, EDIT SCRIPT... and CLOSE PROJECT", async () => {
      await session.step(190, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(191, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(192, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(193, "And user enters \"BDDLifeScriptProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeScriptProj{time}"), el("gallery search")));
      await session.step(194, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(195, "And user double-clicks on BDDLifeScriptProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDLifeScriptProj{time} gallery card"))));
      await session.step(196, "Then \"Data loading error\" dialog should be visible", () => shouldBe(page, el("\"Data loading error\" dialog"), "visible"));
      await session.step(197, "And \"Data loading error\" dialog should contain text \"could not load some of its data\"", () => shouldContainText(page, el("\"Data loading error\" dialog"), "could not load some of its data"));
      await session.step(198, "And \"Data loading error\" dialog should contain text \"Error: intentional break\"", () => shouldContainText(page, el("\"Data loading error\" dialog"), "Error: intentional break"));
      await session.step(199, "And \"OPEN ANYWAY\" button in \"Data loading error\" dialog should be visible", () => shouldBe(page, el("\"OPEN ANYWAY\" button in \"Data loading error\" dialog"), "visible"));
      await session.step(200, "And \"EDIT SCRIPT...\" button in \"Data loading error\" dialog should be visible", () => shouldBe(page, el("\"EDIT SCRIPT...\" button in \"Data loading error\" dialog"), "visible"));
      await session.step(201, "And \"CLOSE PROJECT\" button in \"Data loading error\" dialog should be visible", () => shouldBe(page, el("\"CLOSE PROJECT\" button in \"Data loading error\" dialog"), "visible"));
      await session.step(202, "When user clicks on \"CLOSE PROJECT\" button in \"Data loading error\" dialog", () => clickOn(page, el("\"CLOSE PROJECT\" button in \"Data loading error\" dialog")));
      await session.step(203, "Then the \"Data loading error\" dialog should close", () => dialogCloses(page, "Data loading error"));
    });
    run.finish();
  });
});
