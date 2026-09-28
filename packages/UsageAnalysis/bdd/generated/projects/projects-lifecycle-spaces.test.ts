/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-lifecycle-spaces.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, views.space, sharing.share-dialog, GROK-21025]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/minimized-viewers.js';
import '../../bindings/projects-copies.js';
import '../../bindings/projects-derived.js';
import '../../bindings/projects-regressions.js';
import '../../bindings/projects-sources.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, dragTo, enterInto, isExpanded, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldOffer, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, contextPanelShows, creationScriptHolds, dialogCloses, noProjectOnServer, noSpaceOnServer, pickSharingUser, projectsOnServer, reloadedByDataSync, savedWithDataSync, sharingPaneLists, signInAsSecond, signInAsSelf, spacesOnServer, urlShouldContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {fileOnServer, noDialogShows, openEndsIn} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project built from a file in a space, shared, and renamed with its space", () => {
  const session = feature(test, "features/projects/projects-lifecycle-spaces.feature", import.meta.url);
  test("A project built from a file in a space, shared, and renamed with its space", {tag: ["@journey", "@serial", "@realizes:views.projects", "@realizes:views.space", "@realizes:sharing.share-dialog", "@known-failure", "@realizes:GROK-21025"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 17, page);
    await session.step(34, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(36, "And no project named \"BDDLifeSpaceProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDLifeSpaceProj{time}")));
    await session.step(37, "And no project named \"BDDLifeSpaceProjRenamed{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDLifeSpaceProjRenamed{time}")));
    await session.step(38, "And no space named \"BDDLifeSpace{time}, BDDLifeSpaceRen{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDLifeSpace{time}, BDDLifeSpaceRen{time}")));
    await run.scenario("A space is created", async () => {
      await session.step(41, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(42, "Then Create Space dialog should be visible", () => shouldBe(page, el("Create Space dialog"), "visible"));
      await session.step(43, "When user enters \"BDDLifeSpace{time}\" into Name input in Create Space dialog", () => enterInto(page, session.text("BDDLifeSpace{time}"), el("Name input in Create Space dialog")));
      await session.step(44, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(45, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(46, "And 1 space named \"BDDLifeSpace{time}\" should be on the server", () => spacesOnServer(page, 1, session.text("BDDLifeSpace{time}")));
      await session.step(47, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(48, "Then BDDLifeSpace{time} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("BDDLifeSpace{time} tree node inside browse tree")), "visible"));
    });
    await run.scenario("The demo file is copied into the space", async () => {
      await session.step(51, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(52, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(53, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
      await session.step(54, "When user drags demog.csv link in gallery to BDDLifeSpace{time} tree node inside browse tree", () => dragTo(page, el("demog.csv link in gallery"), el(session.text("BDDLifeSpace{time} tree node inside browse tree"))));
      await session.step(55, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(56, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
      await session.step(57, "And choice input in Move entity dialog should offer \"Move, Link, Copy\"", () => shouldOffer(page, el("choice input in Move entity dialog"), "Move, Link, Copy"));
      await session.step(58, "When user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(59, "Then choice input in Move entity dialog should have value \"Copy\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Copy"));
      await session.step(60, "When user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(61, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(62, "And the file \"BDDLifeSpace{time}:Files/demog.csv\" should be on the server", () => fileOnServer(page, session.text("BDDLifeSpace{time}:Files/demog.csv")));
      await session.step(63, "And the file \"System:DemoFiles/demog.csv\" should be on the server", () => fileOnServer(page, "System:DemoFiles/demog.csv"));
      await session.step(64, "When user double-clicks on BDDLifeSpace{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDLifeSpace{time} tree node inside browse tree"))));
      await session.step(65, "Then the \"BDDLifeSpace{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDLifeSpace{time}")));
      await session.step(66, "And demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
      await session.step(67, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(68, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
      await session.step(69, "And demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
    });
    await run.scenario("The file opens from the space", async () => {
      await session.step(72, "When user double-clicks on BDDLifeSpace{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDLifeSpace{time} tree node inside browse tree"))));
      await session.step(73, "Then the \"BDDLifeSpace{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDLifeSpace{time}")));
      await session.step(74, "When user double-clicks on demog.csv link in gallery", () => doubleClickOn(page, el("demog.csv link in gallery")));
      await session.step(75, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(76, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(77, "And the page address should contain \"/file/BDDLifeSpace{time}.Files/demog.csv\"", () => urlShouldContain(page, session.text("/file/BDDLifeSpace{time}.Files/demog.csv")));
    });
    await run.scenario("The file from the space is saved as a project with Data sync", async () => {
      await session.step(80, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(81, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(82, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
      await session.step(83, "When user enters \"BDDLifeSpaceProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDLifeSpaceProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(84, "And user clicks on \"Creation script\" button in \"demog\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog")));
      await session.step(85, "Then \"demog\" project table in \"Save project\" dialog should contain text 'OpenFile(\"BDDLifeSpace{time}:Files/demog.csv\")'", () => shouldContainText(page, el("\"demog\" project table in \"Save project\" dialog"), session.text("OpenFile(\"BDDLifeSpace{time}:Files/demog.csv\")")));
      await session.step(86, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(87, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(88, "And an info balloon containing 'Project \"BDDLifeSpaceProj{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDLifeSpaceProj{time}\" uploaded")));
      await session.step(89, "And 1 project named \"BDDLifeSpaceProj{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDLifeSpaceProj{time}")));
      await session.step(90, "And the \"demog\" table of the \"BDDLifeSpaceProj{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog", session.text("BDDLifeSpaceProj{time}")));
      await session.step(91, "And the creation script of the \"demog\" table of the \"BDDLifeSpaceProj{time}\" project on the server should contain 'OpenFile(\"BDDLifeSpace{time}:Files/demog.csv\")'", () => creationScriptHolds(page, "demog", session.text("BDDLifeSpaceProj{time}"), session.text("OpenFile(\"BDDLifeSpace{time}:Files/demog.csv\")")));
      await session.step(92, "And \"Share BDDLifeSpaceProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDLifeSpaceProj{time}\" dialog")), "visible"));
      await session.step(93, "When user clicks on CANCEL button in \"Share BDDLifeSpaceProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDLifeSpaceProj{time}\" dialog"))));
      await session.step(94, "Then the \"Share BDDLifeSpaceProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDLifeSpaceProj{time}")));
      await session.step(95, "When user picks \"Close All\" from the context menu of left sidebar", () => pickFromContextMenu(page, "Close All", el("left sidebar")));
      await session.step(96, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Only the project is shared with the second account", async () => {
      await session.step(100, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(101, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(102, "And user enters \"BDDLifeSpaceProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeSpaceProj{time}"), el("gallery search")));
      await session.step(103, "And user picks \"Share...\" from the context menu of BDDLifeSpaceProj{time} gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("BDDLifeSpaceProj{time} gallery card"))));
      await session.step(104, "Then \"Share BDDLifeSpaceProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDLifeSpaceProj{time}\" dialog")), "visible"));
      await session.step(105, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(106, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDLifeSpaceProj{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDLifeSpaceProj{time}\" dialog"))));
      await session.step(107, "And user unchecks \"Send notifications\" input in \"Share BDDLifeSpaceProj{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDDLifeSpaceProj{time}\" dialog"))));
      await session.step(108, "And user clicks on OK button in \"Share BDDLifeSpaceProj{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDLifeSpaceProj{time}\" dialog"))));
      await session.step(109, "Then the \"Share BDDLifeSpaceProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDLifeSpaceProj{time}")));
      await session.step(110, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(111, "When user clicks on BDDLifeSpaceProj{time} gallery card", () => clickOn(page, el(session.text("BDDLifeSpaceProj{time} gallery card"))));
      await session.step(112, "Then the context panel should show \"BDDLifeSpaceProj{time}\"", () => contextPanelShows(page, session.text("BDDLifeSpaceProj{time}")));
      await session.step(113, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
    });
    await run.scenario("The second account opens the shared project (GROK-18345)", async () => {
      await session.step(116, "Given user signs in as the sharing user", () => signInAsSecond(page));
      await session.step(117, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(118, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(119, "And user enters \"BDDLifeSpaceProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeSpaceProj{time}"), el("gallery search")));
      await session.step(120, "And user clicks on Refresh icon in gallery toolbar", () => clickOn(page, el("Refresh icon in gallery toolbar")));
      await session.step(121, "And user double-clicks on BDDLifeSpaceProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDLifeSpaceProj{time} gallery card"))));
      await session.step(122, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(123, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(124, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(125, "And \"Data loading error\" dialog should be absent", () => shouldBe(page, el("\"Data loading error\" dialog"), "absent"));
      await session.step(126, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(127, "And no errors should have been logged", () => noErrors(page));
      await session.step(128, "When user picks \"Close All\" from the context menu of left sidebar", () => pickFromContextMenu(page, "Close All", el("left sidebar")));
      await session.step(129, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    });
    await run.scenario("The owner renames the space", async () => {
      await session.step(132, "Given user signs in as themselves again", () => signInAsSelf(page));
      await session.step(133, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(134, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(135, "When user picks \"Rename...\" from the context menu of BDDLifeSpace{time} tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("BDDLifeSpace{time} tree node inside browse tree"))));
      await session.step(136, "Then Rename project dialog should be visible", () => shouldBe(page, el("Rename project dialog"), "visible"));
      await session.step(137, "And Name input in Rename project dialog should have value \"BDDLifeSpace{time}\"", () => shouldHaveValue(page, el("Name input in Rename project dialog"), session.text("BDDLifeSpace{time}")));
      await session.step(138, "When user enters \"BDDLifeSpaceRen{time}\" into Name input in Rename project dialog", () => enterInto(page, session.text("BDDLifeSpaceRen{time}"), el("Name input in Rename project dialog")));
      await session.step(139, "And user clicks on OK button in Rename project dialog", () => clickOn(page, el("OK button in Rename project dialog")));
      await session.step(140, "Then the \"Rename project\" dialog should close", () => dialogCloses(page, "Rename project"));
      await session.step(141, "And 1 space named \"BDDLifeSpaceRen{time}\" should be on the server", () => spacesOnServer(page, 1, session.text("BDDLifeSpaceRen{time}")));
      await session.step(142, "And 0 spaces named \"BDDLifeSpace{time}\" should be on the server", () => spacesOnServer(page, 0, session.text("BDDLifeSpace{time}")));
    });
    await run.scenario("The owner opens the project after the space was renamed", async () => {
      await session.step(145, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(146, "And user enters \"BDDLifeSpaceProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeSpaceProj{time}"), el("gallery search")));
      await session.step(147, "And user clicks on Refresh icon in gallery toolbar", () => clickOn(page, el("Refresh icon in gallery toolbar")));
      await session.step(148, "And user double-clicks on BDDLifeSpaceProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDLifeSpaceProj{time} gallery card"))));
      await session.step(149, "Then the open should end in the \"demog\" view or the \"Data loading error\" dialog", () => openEndsIn(page, "demog", "Data loading error"));
    });
    await run.scenario("After the space rename, the owner's open does not fail on the old space name", async () => {
      await session.step(153, "Then no \"Data loading error\" dialog should show 'OpenFile(\"BDDLifeSpace{time}:Files/demog.csv\")'", () => noDialogShows(page, "Data loading error", session.text("OpenFile(\"BDDLifeSpace{time}:Files/demog.csv\")")));
      await session.step(154, "And the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(155, "And the table should have 5850 rows", () => rowCount(page, 5850));
    }, {knownFailure: true});
    await run.scenario("The second account opens the project after the space was renamed", async () => {
      await session.step(158, "Given user signs in as the sharing user", () => signInAsSecond(page));
      await session.step(159, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(160, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(161, "And user enters \"BDDLifeSpaceProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeSpaceProj{time}"), el("gallery search")));
      await session.step(162, "And user clicks on Refresh icon in gallery toolbar", () => clickOn(page, el("Refresh icon in gallery toolbar")));
      await session.step(163, "And user double-clicks on BDDLifeSpaceProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDLifeSpaceProj{time} gallery card"))));
      await session.step(164, "Then the open should end in the \"demog\" view or the \"Data loading error\" dialog", () => openEndsIn(page, "demog", "Data loading error"));
    });
    await run.scenario("After the space rename, the second account's open does not fail on the old space name", async () => {
      await session.step(168, "Then no \"Data loading error\" dialog should show 'OpenFile(\"BDDLifeSpace{time}:Files/demog.csv\")'", () => noDialogShows(page, "Data loading error", session.text("OpenFile(\"BDDLifeSpace{time}:Files/demog.csv\")")));
      await session.step(169, "And the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(170, "And the table should have 5850 rows", () => rowCount(page, 5850));
    }, {knownFailure: true});
    await run.scenario("The owner renames the project", async () => {
      await session.step(173, "Given user signs in as themselves again", () => signInAsSelf(page));
      await session.step(174, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(175, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(176, "And user enters \"BDDLifeSpaceProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeSpaceProj{time}"), el("gallery search")));
      await session.step(177, "And user picks \"Rename...\" from the context menu of BDDLifeSpaceProj{time} gallery card", () => pickFromContextMenu(page, "Rename...", el(session.text("BDDLifeSpaceProj{time} gallery card"))));
      await session.step(178, "Then Name input in Rename project dialog should have value \"BDDLifeSpaceProj{time}\"", () => shouldHaveValue(page, el("Name input in Rename project dialog"), session.text("BDDLifeSpaceProj{time}")));
      await session.step(179, "When user enters \"BDDLifeSpaceProjRenamed{time}\" into Name input in Rename project dialog", () => enterInto(page, session.text("BDDLifeSpaceProjRenamed{time}"), el("Name input in Rename project dialog")));
      await session.step(180, "And user clicks on OK button in Rename project dialog", () => clickOn(page, el("OK button in Rename project dialog")));
      await session.step(181, "Then the \"Rename project\" dialog should close", () => dialogCloses(page, "Rename project"));
      await session.step(182, "And 1 project named \"BDDLifeSpaceProjRenamed{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDLifeSpaceProjRenamed{time}")));
      await session.step(183, "And 0 projects named \"BDDLifeSpaceProj{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDLifeSpaceProj{time}")));
      await session.step(184, "When user enters \"BDDLifeSpaceProjRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeSpaceProjRenamed{time}"), el("gallery search")));
      await session.step(185, "And user clicks on Refresh icon in gallery toolbar", () => clickOn(page, el("Refresh icon in gallery toolbar")));
      await session.step(186, "Then BDDLifeSpaceProjRenamed{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDLifeSpaceProjRenamed{time} gallery card")), "visible"));
    });
    await run.scenario("The owner opens the renamed project", async () => {
      await session.step(189, "When user double-clicks on BDDLifeSpaceProjRenamed{time} gallery card", () => doubleClickOn(page, el(session.text("BDDLifeSpaceProjRenamed{time} gallery card"))));
      await session.step(190, "Then the open should end in the \"demog\" view or the \"Data loading error\" dialog", () => openEndsIn(page, "demog", "Data loading error"));
    });
    await run.scenario("The renamed project's open does not fail on the old space name for the owner", async () => {
      await session.step(194, "Then no \"Data loading error\" dialog should show 'OpenFile(\"BDDLifeSpace{time}:Files/demog.csv\")'", () => noDialogShows(page, "Data loading error", session.text("OpenFile(\"BDDLifeSpace{time}:Files/demog.csv\")")));
      await session.step(195, "And the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(196, "And the table should have 5850 rows", () => rowCount(page, 5850));
    }, {knownFailure: true});
    await run.scenario("The second account opens the renamed project", async () => {
      await session.step(199, "Given user signs in as the sharing user", () => signInAsSecond(page));
      await session.step(200, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(201, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(202, "And user enters \"BDDLifeSpaceProjRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeSpaceProjRenamed{time}"), el("gallery search")));
      await session.step(203, "And user clicks on Refresh icon in gallery toolbar", () => clickOn(page, el("Refresh icon in gallery toolbar")));
      await session.step(204, "And user double-clicks on BDDLifeSpaceProjRenamed{time} gallery card", () => doubleClickOn(page, el(session.text("BDDLifeSpaceProjRenamed{time} gallery card"))));
      await session.step(205, "Then the open should end in the \"demog\" view or the \"Data loading error\" dialog", () => openEndsIn(page, "demog", "Data loading error"));
    });
    await run.scenario("The renamed project's open does not fail on the old space name for the second account", async () => {
      await session.step(209, "Then no \"Data loading error\" dialog should show 'OpenFile(\"BDDLifeSpace{time}:Files/demog.csv\")'", () => noDialogShows(page, "Data loading error", session.text("OpenFile(\"BDDLifeSpace{time}:Files/demog.csv\")")));
      await session.step(210, "And the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(211, "And the table should have 5850 rows", () => rowCount(page, 5850));
    }, {knownFailure: true});
    await run.scenario("The owner deletes the project and the space", async () => {
      await session.step(214, "Given user signs in as themselves again", () => signInAsSelf(page));
      await session.step(215, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(216, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(217, "And user enters \"BDDLifeSpaceProjRenamed{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeSpaceProjRenamed{time}"), el("gallery search")));
      await session.step(218, "And user picks \"Delete Project\" from the context menu of BDDLifeSpaceProjRenamed{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDLifeSpaceProjRenamed{time} gallery card"))));
      await session.step(219, "Then \"Are you sure?\" dialog should contain text 'Delete project \"BDDLifeSpaceProjRenamed{time}\"?'", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete project \"BDDLifeSpaceProjRenamed{time}\"?")));
      await session.step(220, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(221, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(222, "And 0 projects named \"BDDLifeSpaceProjRenamed{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDLifeSpaceProjRenamed{time}")));
      await session.step(223, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(224, "When user picks \"Delete Space\" from the context menu of BDDLifeSpaceRen{time} tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el(session.text("BDDLifeSpaceRen{time} tree node inside browse tree"))));
      await session.step(225, "Then \"Are you sure?\" dialog should contain text 'Delete space \"BDDLifeSpaceRen{time}\"?'", () => shouldContainText(page, el("\"Are you sure?\" dialog"), session.text("Delete space \"BDDLifeSpaceRen{time}\"?")));
      await session.step(226, "And \"Are you sure?\" dialog should contain text \"This will delete space and its related data\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), "This will delete space and its related data"));
      await session.step(227, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(228, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(229, "And 0 spaces named \"BDDLifeSpaceRen{time}\" should be on the server", () => spacesOnServer(page, 0, session.text("BDDLifeSpaceRen{time}")));
    });
    run.finish();
  });
});
