/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-lifecycle-spaces.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, views.space, sharing.share-dialog, core#627]
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
import {clickOn, doubleClickOn, dragTo, enterInto, isExpanded, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldOffer, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, doubleClickOpensProject, noProjectOnServer, noSpaceOnServer, pickSharingUser, projectsOnServer, reloadedByDataSync, runningAccountSignedIn, savedWithDataSync, sharingPaneLists, sharingUserSignedIn, signInAsSelf, signInAsSharingUser, spaceOnServer, spacesOnServer, urlShouldContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project built from a file in a space, shared, and reopened after the space is renamed", () => {
  const session = feature(test, "features/projects/projects-lifecycle-spaces.feature", import.meta.url);
  test("A project built from a file in a space, shared, and reopened after the space is renamed", {tag: ["@journey", "@realizes:views.projects", "@realizes:views.space", "@realizes:sharing.share-dialog", "@known-failure:135", "@realizes:core#627"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(29, "And no project named \"BDDLifeSpaceProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDLifeSpaceProj{time}")));
    await session.step(30, "And no space named \"BDDLifeSpace{time}, BDDLifeSpaceRen{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDLifeSpace{time}, BDDLifeSpaceRen{time}")));
    await session.step(31, "And a space named \"BDDLifeSpace{time}\" is on the server", () => spaceOnServer(page, session.text("BDDLifeSpace{time}")));
    await run.scenario("The demo file is copied into the space", async () => {
      await session.step(34, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(35, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(36, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(37, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
      await session.step(38, "When user drags cars.csv link in gallery to BDDLifeSpace{time} tree node inside browse tree", () => dragTo(page, el("cars.csv link in gallery"), el(session.text("BDDLifeSpace{time} tree node inside browse tree"))));
      await session.step(39, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(40, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
      await session.step(41, "And choice input in Move entity dialog should offer \"Move, Link, Copy\"", () => shouldOffer(page, el("choice input in Move entity dialog"), "Move, Link, Copy"));
      await session.step(42, "When user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(43, "Then choice input in Move entity dialog should have value \"Copy\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Copy"));
      await session.step(44, "When user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(45, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(46, "When user double-clicks on BDDLifeSpace{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDLifeSpace{time} tree node inside browse tree"))));
      await session.step(47, "Then the \"BDDLifeSpace{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDLifeSpace{time}")));
      await session.step(48, "And cars.csv link in gallery should be visible", () => shouldBe(page, el("cars.csv link in gallery"), "visible"));
      await session.step(49, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(50, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
      await session.step(51, "And cars.csv link in gallery should be visible", () => shouldBe(page, el("cars.csv link in gallery"), "visible"));
    });
    await run.scenario("The file opens from the space", async () => {
      await session.step(54, "When user double-clicks on BDDLifeSpace{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDLifeSpace{time} tree node inside browse tree"))));
      await session.step(55, "Then the \"BDDLifeSpace{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDLifeSpace{time}")));
      await session.step(56, "When user double-clicks on cars.csv link in gallery", () => doubleClickOn(page, el("cars.csv link in gallery")));
      await session.step(57, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
      await session.step(58, "And the table should have 30 rows", () => rowCount(page, 30));
      await session.step(59, "And the page address should contain \"/file/BDDLifeSpace{time}.Files/cars.csv\"", () => urlShouldContain(page, session.text("/file/BDDLifeSpace{time}.Files/cars.csv")));
    });
    await run.scenario("The file from the space is saved as a project with Data sync", async () => {
      await session.step(62, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(63, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(64, "And \"Creation script\" button in \"cars\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"cars\" project table in \"Save project\" dialog"), "visible"));
      await session.step(65, "And Data sync switch in \"cars\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"cars\" project table in \"Save project\" dialog"), "checked"));
      await session.step(66, "When user enters \"BDDLifeSpaceProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDLifeSpaceProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(67, "And user clicks on \"Creation script\" button in \"cars\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"cars\" project table in \"Save project\" dialog")));
      await session.step(68, "Then \"cars\" project table in \"Save project\" dialog should contain text 'OpenFile(\"BDDLifeSpace{time}:Files/cars.csv\")'", () => shouldContainText(page, el("\"cars\" project table in \"Save project\" dialog"), session.text("OpenFile(\"BDDLifeSpace{time}:Files/cars.csv\")")));
      await session.step(69, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(70, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(71, "And an info balloon containing 'Project \"BDDLifeSpaceProj{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDLifeSpaceProj{time}\" uploaded")));
      await session.step(72, "And 1 project named \"BDDLifeSpaceProj{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDLifeSpaceProj{time}")));
      await session.step(73, "And the \"cars\" table of the \"BDDLifeSpaceProj{time}\" project should be saved with data sync", () => savedWithDataSync(page, "cars", session.text("BDDLifeSpaceProj{time}")));
      await session.step(74, "And \"Share BDDLifeSpaceProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDLifeSpaceProj{time}\" dialog")), "visible"));
      await session.step(75, "When user clicks on CANCEL button in \"Share BDDLifeSpaceProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDLifeSpaceProj{time}\" dialog"))));
      await session.step(76, "Then the \"Share BDDLifeSpaceProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDLifeSpaceProj{time}")));
      await session.step(77, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(78, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Only the project is shared with the second account", async () => {
      await session.step(82, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(83, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(84, "And user enters \"BDDLifeSpaceProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeSpaceProj{time}"), el("gallery search")));
      await session.step(85, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(86, "And user picks \"Share...\" from the context menu of BDDLifeSpaceProj{time} gallery card", () => pickFromContextMenu(page, "Share...", el(session.text("BDDLifeSpaceProj{time} gallery card"))));
      await session.step(87, "Then \"Share BDDLifeSpaceProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDLifeSpaceProj{time}\" dialog")), "visible"));
      await session.step(88, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(89, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDDLifeSpaceProj{time}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDDLifeSpaceProj{time}\" dialog"))));
      await session.step(90, "And user unchecks \"Send notifications\" input in \"Share BDDLifeSpaceProj{time}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDDLifeSpaceProj{time}\" dialog"))));
      await session.step(91, "And user clicks on OK button in \"Share BDDLifeSpaceProj{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDDLifeSpaceProj{time}\" dialog"))));
      await session.step(92, "Then the \"Share BDDLifeSpaceProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDLifeSpaceProj{time}")));
      await session.step(93, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(94, "When user clicks on BDDLifeSpaceProj{time} gallery card", () => clickOn(page, el(session.text("BDDLifeSpaceProj{time} gallery card"))));
      await session.step(95, "Then the context panel should show \"BDDLifeSpaceProj{time}\"", () => contextPanelShows(page, session.text("BDDLifeSpaceProj{time}")));
      await session.step(96, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
    });
    await run.scenario("The second account opens the shared project (GROK-18345)", async () => {
      await session.step(99, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(100, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(101, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(102, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(103, "And user enters \"BDDLifeSpaceProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeSpaceProj{time}"), el("gallery search")));
      await session.step(104, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(105, "And user double-clicks on BDDLifeSpaceProj{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDLifeSpaceProj{time} gallery card"))));
      await session.step(106, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
      await session.step(107, "And the table should have 30 rows", () => rowCount(page, 30));
      await session.step(108, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(109, "And \"Data loading error\" dialog should be absent", () => shouldBe(page, el("\"Data loading error\" dialog"), "absent"));
      await session.step(110, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(111, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(112, "And user signs in as themselves again", () => signInAsSelf(page));
      await session.step(113, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    });
    await run.scenario("The owner renames the space", async () => {
      await session.step(116, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(117, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(118, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(119, "When user picks \"Rename...\" from the context menu of BDDLifeSpace{time} tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("BDDLifeSpace{time} tree node inside browse tree"))));
      await session.step(120, "Then Rename project dialog should be visible", () => shouldBe(page, el("Rename project dialog"), "visible"));
      await session.step(121, "And Name input in Rename project dialog should have value \"BDDLifeSpace{time}\"", () => shouldHaveValue(page, el("Name input in Rename project dialog"), session.text("BDDLifeSpace{time}")));
      await session.step(122, "When user enters \"BDDLifeSpaceRen{time}\" into Name input in Rename project dialog", () => enterInto(page, session.text("BDDLifeSpaceRen{time}"), el("Name input in Rename project dialog")));
      await session.step(123, "And user clicks on OK button in Rename project dialog", () => clickOn(page, el("OK button in Rename project dialog")));
      await session.step(124, "Then the \"Rename project\" dialog should close", () => dialogCloses(page, "Rename project"));
      await session.step(125, "And 1 space named \"BDDLifeSpaceRen{time}\" should be on the server", () => spacesOnServer(page, 1, session.text("BDDLifeSpaceRen{time}")));
      await session.step(126, "And 0 spaces named \"BDDLifeSpace{time}\" should be on the server", () => spacesOnServer(page, 0, session.text("BDDLifeSpace{time}")));
      await session.step(127, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(128, "And user enters \"BDDLifeSpaceProj{time}\" into gallery search", () => enterInto(page, session.text("BDDLifeSpaceProj{time}"), el("gallery search")));
      await session.step(129, "And user clicks on Refresh icon in gallery toolbar", () => clickOn(page, el("Refresh icon in gallery toolbar")));
      await session.step(130, "Then BDDLifeSpaceProj{time} gallery card should be visible", () => shouldBe(page, el(session.text("BDDLifeSpaceProj{time} gallery card")), "visible"));
    });
    await run.scenario("After the space rename, the owner opens the project with its rows", async () => {
      await session.step(134, "When user double-clicks on BDDLifeSpaceProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDLifeSpaceProj{time} gallery card"))));
      await session.step(135, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
      await session.step(136, "And the table should have 30 rows", () => rowCount(page, 30));
      await session.step(137, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
    }, {knownFailure: true, line: 135});
    run.finish();
  });
});
