/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-lifecycle-spaces.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, views.space, sharing.share-dialog, GROK-21025]
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
import {clickOn, doubleClickOn, dragTo, enterInto, isExpanded, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldOffer, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, contextPanelShows, dialogCloses, noProjectOnServer, noSpaceOnServer, pickSharingUser, projectsOnServer, reloadedByDataSync, runningAccountSignedIn, savedWithDataSync, sharingPaneLists, sharingUserSignedIn, signInAsSelf, signInAsSharingUser, spacesOnServer, urlShouldContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project built from a file in a space, shared, and reopened after the space is renamed", () => {
  const session = feature(test, "features/projects/projects-lifecycle-spaces.feature", import.meta.url);
  test("A project built from a file in a space, shared, and reopened after the space is renamed", {tag: ["@journey", "@serial", "@realizes:views.projects", "@realizes:views.space", "@realizes:sharing.share-dialog", "@known-failure", "@realizes:GROK-21025"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And no project named \"BDDLifeSpaceProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDLifeSpaceProj{time}")));
    await session.step(31, "And no space named \"BDDLifeSpace{time}, BDDLifeSpaceRen{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDLifeSpace{time}, BDDLifeSpaceRen{time}")));
    await run.scenario("A space is created", async () => {
      await session.step(34, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(35, "Then Create Space dialog should be visible", () => shouldBe(page, el("Create Space dialog"), "visible"));
      await session.step(36, "When user enters \"BDDLifeSpace{time}\" into Name input in Create Space dialog", () => enterInto(page, session.text("BDDLifeSpace{time}"), el("Name input in Create Space dialog")));
      await session.step(37, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(38, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(39, "And 1 space named \"BDDLifeSpace{time}\" should be on the server", () => spacesOnServer(page, 1, session.text("BDDLifeSpace{time}")));
      await session.step(40, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(41, "Then BDDLifeSpace{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("BDDLifeSpace{time} tree node inside browse tree")), "visible"));
    });
    await run.scenario("The demo file is copied into the space", async () => {
      await session.step(44, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(45, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(46, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
      await session.step(47, "When user drags demog.csv link in gallery to BDDLifeSpace{time} tree node inside browse tree", () => dragTo(page, el("demog.csv link in gallery"), el(session.text("BDDLifeSpace{time} tree node inside browse tree"))));
      await session.step(48, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(49, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
      await session.step(50, "And choice input in Move entity dialog should offer \"Move, Link, Copy\"", () => shouldOffer(page, el("choice input in Move entity dialog"), "Move, Link, Copy"));
      await session.step(51, "When user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(52, "Then choice input in Move entity dialog should have value \"Copy\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Copy"));
      await session.step(53, "When user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(54, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(55, "When user double-clicks on BDDLifeSpace{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDLifeSpace{time} tree node inside browse tree"))));
      await session.step(56, "Then the \"BDDLifeSpace{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDLifeSpace{time}")));
      await session.step(57, "And demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
      await session.step(58, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(59, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
      await session.step(60, "And demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
    });
    await run.scenario("The file opens from the space", async () => {
      await session.step(63, "When user double-clicks on BDDLifeSpace{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDLifeSpace{time} tree node inside browse tree"))));
      await session.step(64, "Then the \"BDDLifeSpace{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDLifeSpace{time}")));
      await session.step(65, "When user double-clicks on demog.csv link in gallery", () => doubleClickOn(page, el("demog.csv link in gallery")));
      await session.step(66, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(67, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(68, "And the page address should contain \"/file/BDDLifeSpace{time}.Files/demog.csv\"", () => urlShouldContain(page, session.text("/file/BDDLifeSpace{time}.Files/demog.csv")));
    });
    await run.scenario("The file from the space is saved as a project with Data sync", async () => {
      await session.step(71, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(72, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(73, "And \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "visible"));
      await session.step(74, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(75, "When user enters \"BDDLifeSpaceProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDLifeSpaceProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(76, "And user clicks on \"Creation script\" button in \"demog\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog")));
      await session.step(77, "Then \"demog\" project table in \"Save project\" dialog should contain text 'OpenFile(\"BDDLifeSpace{time}:Files/demog.csv\")'", () => shouldContainText(page, el("\"demog\" project table in \"Save project\" dialog"), session.text("OpenFile(\"BDDLifeSpace{time}:Files/demog.csv\")")));
      await session.step(78, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(79, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(80, "And an info balloon containing 'Project \"BDDLifeSpaceProj{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDLifeSpaceProj{time}\" uploaded")));
      await session.step(81, "And 1 project named \"BDDLifeSpaceProj{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDLifeSpaceProj{time}")));
      await session.step(82, "And the \"demog\" table of the \"BDDLifeSpaceProj{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog", session.text("BDDLifeSpaceProj{time}")));
      await session.step(83, "And \"Share BDDLifeSpaceProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDLifeSpaceProj{time}\" dialog")), "visible"));
      await session.step(84, "When user clicks on CANCEL button in \"Share BDDLifeSpaceProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDLifeSpaceProj{time}\" dialog"))));
      await session.step(85, "Then the \"Share BDDLifeSpaceProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDLifeSpaceProj{time}")));
      await session.step(86, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(87, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Only the project is shared with the second account", async () => {
      await session.step(91, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(92, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(93, "And user enters \"BDDLifeSpaceProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeSpaceProj{time}"), el("gallery search")));
      await session.step(94, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(95, "And user picks \"Share...\" from the context menu of BDDLifeSpaceProj{time} gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("BDDLifeSpaceProj{time} gallery card"))));
      await session.step(96, "Then \"Share BDDLifeSpaceProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDLifeSpaceProj{time}\" dialog")), "visible"));
      await session.step(97, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(98, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDLifeSpaceProj{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDLifeSpaceProj{time}\" dialog"))));
      await session.step(99, "And user unchecks \"Send notifications\" input in \"Share BDDLifeSpaceProj{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDDLifeSpaceProj{time}\" dialog"))));
      await session.step(100, "And user clicks on OK button in \"Share BDDLifeSpaceProj{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDLifeSpaceProj{time}\" dialog"))));
      await session.step(101, "Then the \"Share BDDLifeSpaceProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDLifeSpaceProj{time}")));
      await session.step(102, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(103, "When user clicks on BDDLifeSpaceProj{time} gallery card", () => clickOn(page, el(session.text("BDDLifeSpaceProj{time} gallery card"))));
      await session.step(104, "Then the context panel should show \"BDDLifeSpaceProj{time}\"", () => contextPanelShows(page, session.text("BDDLifeSpaceProj{time}")));
      await session.step(105, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
    });
    await run.scenario("The second account opens the shared project (GROK-18345)", async () => {
      await session.step(108, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(109, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(110, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(111, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(112, "And user enters \"BDDLifeSpaceProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeSpaceProj{time}"), el("gallery search")));
      await session.step(113, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(114, "And user double-clicks on BDDLifeSpaceProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDLifeSpaceProj{time} gallery card"))));
      await session.step(115, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(116, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(117, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(118, "And \"Data loading error\" dialog should be absent", () => shouldBe(page, el("\"Data loading error\" dialog"), "absent"));
      await session.step(119, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(120, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(121, "And user signs in as themselves again", () => signInAsSelf(page));
      await session.step(122, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    });
    await run.scenario("The owner renames the space", async () => {
      await session.step(125, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(126, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(127, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(128, "When user picks \"Rename...\" from the context menu of BDDLifeSpace{time} tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("BDDLifeSpace{time} tree node inside browse tree"))));
      await session.step(129, "Then Rename project dialog should be visible", () => shouldBe(page, el("Rename project dialog"), "visible"));
      await session.step(130, "And Name input in Rename project dialog should have value \"BDDLifeSpace{time}\"", () => shouldHaveValue(page, el("Name input in Rename project dialog"), session.text("BDDLifeSpace{time}")));
      await session.step(131, "When user enters \"BDDLifeSpaceRen{time}\" into Name input in Rename project dialog", () => enterInto(page, session.text("BDDLifeSpaceRen{time}"), el("Name input in Rename project dialog")));
      await session.step(132, "And user clicks on OK button in Rename project dialog", () => clickOn(page, el("OK button in Rename project dialog")));
      await session.step(133, "Then the \"Rename project\" dialog should close", () => dialogCloses(page, "Rename project"));
      await session.step(134, "And 1 space named \"BDDLifeSpaceRen{time}\" should be on the server", () => spacesOnServer(page, 1, session.text("BDDLifeSpaceRen{time}")));
      await session.step(135, "And 0 spaces named \"BDDLifeSpace{time}\" should be on the server", () => spacesOnServer(page, 0, session.text("BDDLifeSpace{time}")));
      await session.step(136, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(137, "And user enters \"BDDLifeSpaceProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeSpaceProj{time}"), el("gallery search")));
      await session.step(138, "And user clicks on Refresh icon in gallery toolbar", () => clickOn(page, el("Refresh icon in gallery toolbar")));
      await session.step(139, "Then BDDLifeSpaceProj{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDLifeSpaceProj{time} gallery card")), "visible"));
    });
    await run.scenario("After the space rename, the owner opens the project with its rows", async () => {
      await session.step(143, "When user double-clicks on BDDLifeSpaceProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDLifeSpaceProj{time} gallery card"))));
      await session.step(144, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(145, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(146, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
    }, {knownFailure: true});
    await run.scenario("The second account signs in after the space rename", async () => {
      await session.step(149, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(150, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(151, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(152, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(153, "And user enters \"BDDLifeSpaceProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeSpaceProj{time}"), el("gallery search")));
      await session.step(154, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(155, "Then BDDLifeSpaceProj{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDLifeSpaceProj{time} gallery card")), "visible"));
    });
    await run.scenario("After the space rename, the second account opens the project with its rows", async () => {
      await session.step(159, "When user double-clicks on BDDLifeSpaceProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDLifeSpaceProj{time} gallery card"))));
      await session.step(160, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(161, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(162, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
    }, {knownFailure: true});
    await run.scenario("The owner deletes the project and the space", async () => {
      await session.step(165, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(166, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(167, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(168, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(169, "And user enters \"BDDLifeSpaceProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeSpaceProj{time}"), el("gallery search")));
      await session.step(170, "And user picks \"Delete Project\" from the context menu of BDDLifeSpaceProj{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDLifeSpaceProj{time} gallery card"))));
      await session.step(171, "Then \"Are you sure?\" dialog should contain text 'Delete project \"BDDLifeSpaceProj{time}\"?'", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete project \"BDDLifeSpaceProj{time}\"?")));
      await session.step(172, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(173, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(174, "And 0 projects named \"BDDLifeSpaceProj{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDLifeSpaceProj{time}")));
      await session.step(175, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(176, "When user picks \"Delete Space\" from the context menu of BDDLifeSpaceRen{time} tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el(session.text("BDDLifeSpaceRen{time} tree node inside browse tree"))));
      await session.step(177, "Then \"Are you sure?\" dialog should contain text 'Delete space \"BDDLifeSpaceRen{time}\"?'", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete space \"BDDLifeSpaceRen{time}\"?")));
      await session.step(178, "And \"Are you sure?\" dialog should contain text \"This will delete space and its related data\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), "This will delete space and its related data"));
      await session.step(179, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(180, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(181, "And 0 spaces named \"BDDLifeSpaceRen{time}\" should be on the server", () => spacesOnServer(page, 0, session.text("BDDLifeSpaceRen{time}")));
    });
    run.finish();
  });
});
