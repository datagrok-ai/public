/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-northwind-uploading-extra.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, file.menu.save.tables-as-project]
--- */
import {test} from '@playwright/test';
import '../../bindings/pkg-projects-copies.js';
import '../../bindings/pkg-projects-derived.js';
import '../../bindings/pkg-projects-sources.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {noErrorsButPreviewNoise} from '../../bindings/pkg-projects-regressions.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, doubleClickOn, enterInto, isExpanded, shouldBe, shouldContainText, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, dialogCloses, markOpenFrames, noProjectOnServer, projectsOnServer, savedAsSnapshot, savedWithDataSync, tableReloaded, tableSnapshot, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noTableLeft} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects uploaded from a Get Top 100 result of the NorthwindTest orders table", () => {
  const session = feature(test, "features/projects/projects-northwind-uploading-extra.feature", import.meta.url);
  test("A Get Top 100 result of the NorthwindTest orders is saved with Data sync ON [Sync=ON, Project=BDDNwUpXTopSync{time}, switch=checks, script=visible, saved=with data sync, how=reloaded by data sync]", {tag: ["@dev-only", "@serial", "@realizes:views.projects", "@realizes:file.menu.save.tables-as-project"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(32, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(36, "Given no project named \"BDDNwUpXTopSync{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDNwUpXTopSync{time}")));
    await session.step(37, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(38, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(39, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
    await session.step(40, "And Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree")));
    await session.step(41, "And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree")));
    await session.step(42, "When user picks \"Get Top 100\" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree", () => pickFromContextMenu(page, "Get Top 100", el("Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree")));
    await session.step(43, "Then the \"orders\" view should be current", () => viewIsCurrent(page, "orders"));
    await session.step(44, "And table \"orders\" should have 100 rows", () => tableRows(page, "orders", 100));
    await session.step(45, "When user marks the open tables as the frames in memory", () => markOpenFrames(page));
    await session.step(46, "And user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(47, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(48, "When user clicks on \"Creation script\" button in \"orders\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"orders\" project table in \"Save project\" dialog")));
    await session.step(49, "Then \"orders\" project table in \"Save project\" dialog should contain text \"limit = 100\"", () => shouldContainText(page, el("\"orders\" project table in \"Save project\" dialog"), "limit = 100"));
    await session.step(50, "When user enters \"BDDNwUpXTopSync{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDNwUpXTopSync{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(51, "And user checks Data sync switch in \"orders\" project table in \"Save project\" dialog", () => check(page, el("Data sync switch in \"orders\" project table in \"Save project\" dialog")));
    await session.step(52, "Then \"Creation script\" button in \"orders\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"orders\" project table in \"Save project\" dialog"), "visible"));
    await session.step(53, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(54, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(55, "And an info balloon containing 'Project \"BDDNwUpXTopSync{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDNwUpXTopSync{time}\" uploaded")));
    await session.step(56, "And 1 project named \"BDDNwUpXTopSync{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDNwUpXTopSync{time}")));
    await session.step(57, "And the \"orders\" table of the \"BDDNwUpXTopSync{time}\" project should be saved with data sync", () => savedWithDataSync(page, "orders", session.text("BDDNwUpXTopSync{time}")));
    await session.step(58, "And \"Share BDDNwUpXTopSync{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDNwUpXTopSync{time}\" dialog")), "visible"));
    await session.step(59, "When user clicks on CANCEL button in \"Share BDDNwUpXTopSync{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDNwUpXTopSync{time}\" dialog"))));
    await session.step(60, "Then the \"Share BDDNwUpXTopSync{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDNwUpXTopSync{time}")));
    await session.step(61, "And no errors but the project preview's should have been logged", () => noErrorsButPreviewNoise(page));
    await session.step(62, "When user picks \"Close All\" from the context menu of left sidebar", () => pickFromContextMenu(page, "Close All", el("left sidebar")));
    await session.step(63, "Then no table should be left in the workspace", () => noTableLeft(page));
    await session.step(64, "Given the browse panel is open", () => browsePanelOpen(page));
    await session.step(65, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(66, "And user enters \"BDDNwUpXTopSync{time}\" into gallery search", () => enterInto(page, session.text("BDDNwUpXTopSync{time}"), el("gallery search")));
    await session.step(67, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
    await session.step(68, "And user double-clicks on BDDNwUpXTopSync{time} gallery card", () => doubleClickOn(page, el(session.text("BDDNwUpXTopSync{time} gallery card"))));
    await session.step(69, "Then table \"orders\" should have been reloaded by data sync with 100 rows", () => tableReloaded(page, "orders", 100));
    await session.step(70, "And no errors should have been logged", () => noErrors(page));
    await session.step(71, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A Get Top 100 result of the NorthwindTest orders is saved with Data sync OFF [Sync=OFF, Project=BDDNwUpXTopNoSync{time}, switch=unchecks, script=hidden, saved=as a snapshot, how=loaded as a snapshot]", {tag: ["@dev-only", "@serial", "@realizes:views.projects", "@realizes:file.menu.save.tables-as-project"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(32, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(36, "Given no project named \"BDDNwUpXTopNoSync{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDNwUpXTopNoSync{time}")));
    await session.step(37, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(38, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(39, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
    await session.step(40, "And Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree")));
    await session.step(41, "And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree")));
    await session.step(42, "When user picks \"Get Top 100\" from the context menu of Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree", () => pickFromContextMenu(page, "Get Top 100", el("Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree")));
    await session.step(43, "Then the \"orders\" view should be current", () => viewIsCurrent(page, "orders"));
    await session.step(44, "And table \"orders\" should have 100 rows", () => tableRows(page, "orders", 100));
    await session.step(45, "When user marks the open tables as the frames in memory", () => markOpenFrames(page));
    await session.step(46, "And user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(47, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(48, "When user clicks on \"Creation script\" button in \"orders\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"orders\" project table in \"Save project\" dialog")));
    await session.step(49, "Then \"orders\" project table in \"Save project\" dialog should contain text \"limit = 100\"", () => shouldContainText(page, el("\"orders\" project table in \"Save project\" dialog"), "limit = 100"));
    await session.step(50, "When user enters \"BDDNwUpXTopNoSync{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDNwUpXTopNoSync{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(51, "And user unchecks Data sync switch in \"orders\" project table in \"Save project\" dialog", () => uncheck(page, el("Data sync switch in \"orders\" project table in \"Save project\" dialog")));
    await session.step(52, "Then \"Creation script\" button in \"orders\" project table in \"Save project\" dialog should be hidden", () => shouldBe(page, el("\"Creation script\" button in \"orders\" project table in \"Save project\" dialog"), "hidden"));
    await session.step(53, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(54, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(55, "And an info balloon containing 'Project \"BDDNwUpXTopNoSync{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDNwUpXTopNoSync{time}\" uploaded")));
    await session.step(56, "And 1 project named \"BDDNwUpXTopNoSync{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDNwUpXTopNoSync{time}")));
    await session.step(57, "And the \"orders\" table of the \"BDDNwUpXTopNoSync{time}\" project should be saved as a snapshot", () => savedAsSnapshot(page, "orders", session.text("BDDNwUpXTopNoSync{time}")));
    await session.step(58, "And \"Share BDDNwUpXTopNoSync{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDNwUpXTopNoSync{time}\" dialog")), "visible"));
    await session.step(59, "When user clicks on CANCEL button in \"Share BDDNwUpXTopNoSync{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDNwUpXTopNoSync{time}\" dialog"))));
    await session.step(60, "Then the \"Share BDDNwUpXTopNoSync{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDNwUpXTopNoSync{time}")));
    await session.step(61, "And no errors but the project preview's should have been logged", () => noErrorsButPreviewNoise(page));
    await session.step(62, "When user picks \"Close All\" from the context menu of left sidebar", () => pickFromContextMenu(page, "Close All", el("left sidebar")));
    await session.step(63, "Then no table should be left in the workspace", () => noTableLeft(page));
    await session.step(64, "Given the browse panel is open", () => browsePanelOpen(page));
    await session.step(65, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(66, "And user enters \"BDDNwUpXTopNoSync{time}\" into gallery search", () => enterInto(page, session.text("BDDNwUpXTopNoSync{time}"), el("gallery search")));
    await session.step(67, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
    await session.step(68, "And user double-clicks on BDDNwUpXTopNoSync{time} gallery card", () => doubleClickOn(page, el(session.text("BDDNwUpXTopNoSync{time} gallery card"))));
    await session.step(69, "Then table \"orders\" should have been loaded as a snapshot with 100 rows", () => tableSnapshot(page, "orders", 100));
    await session.step(70, "And no errors should have been logged", () => noErrors(page));
    await session.step(71, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
