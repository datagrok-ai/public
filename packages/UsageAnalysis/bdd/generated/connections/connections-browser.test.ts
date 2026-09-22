/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/connections/connections-browser.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {connectionHasChat, connectionHasNoChat, deleteChatOfConnection} from '../../bindings/connections.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded, pressKeyIn, rightClickOn, shouldBe, shouldContainText, typeInto, uncheck, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, connectionOnServer, contextPanelOpen, contextPanelShows, dialogCloses, pickSharingUser, sharingPaneLists, urlShouldContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuLists, noBalloons, noErrors, pickFromContextMenu, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A connection in the connections browser and its context panel", () => {
  const session = feature(test, "features/connections/connections-browser.feature", import.meta.url);
  test("A connection in the connections browser and its context panel", {tag: ["@connections", "@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(25, "And a \"Postgres\" connection named \"BDD-Conn-Browser-{run}\" is on the server", () => connectionOnServer(page, "Postgres", session.text("BDD-Conn-Browser-{run}")));
    await session.step(26, "And the context panel is open", () => contextPanelOpen(page));
    await run.scenario("The connections view searches by name and opens a card on the context panel", async () => {
      await session.step(29, "When user right-clicks on Databases---Postgres tree node inside browse tree", () => rightClickOn(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(30, "And user picks \"Browse connections\" from the open menu", () => pickFromOpenMenu(page, "Browse connections"));
      await session.step(31, "Then the \"Postgres\" view should be current", () => viewIsCurrent(page, "Postgres"));
      await session.step(32, "And the page address should contain \"/connections/Postgres\"", () => urlShouldContain(page, "/connections/Postgres"));
      await session.step(33, "When user types \"BDD-Conn-Browser-{run}\" into gallery search", () => typeInto(page, session.text("BDD-Conn-Browser-{run}"), el("gallery search")));
      await session.step(34, "Then \"BDD-Conn-Browser-{run}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDD-Conn-Browser-{run}\" link in gallery")), "visible"));
      await session.step(35, "And there should be 1 visible link in gallery", () => visibleCount(page, 1, el("link in gallery")));
      await session.step(36, "When user clicks on \"BDD-Conn-Browser-{run}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-Conn-Browser-{run}\" link in gallery"))));
      await session.step(37, "Then the context panel should show \"BDD-Conn-Browser-{run}\"", () => contextPanelShows(page, session.text("BDD-Conn-Browser-{run}")));
      await session.step(38, "And Details section in context panel should contain the text \"db.datagrok.ai\"", () => shouldContainText(page, el("Details section in context panel"), "db.datagrok.ai"));
      await session.step(39, "And Details section in context panel should contain the text \"northwind\"", () => shouldContainText(page, el("Details section in context panel"), "northwind"));
      await session.step(40, "And Details section in context panel should contain the text \"54322\"", () => shouldContainText(page, el("Details section in context panel"), "54322"));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The header arrow opens the connection's menu", async () => {
      await session.step(44, "When user clicks on \"context-arrow-down\" icon in context panel", () => clickOn(page, el("\"context-arrow-down\" icon in context panel")));
      await session.step(45, "Then the open menu should list \"Test connection\"", () => menuLists(page, "Test connection"));
      await session.step(46, "And the open menu should list \"Share...\"", () => menuLists(page, "Share..."));
      await session.step(47, "And the open menu should list \"Delete...\"", () => menuLists(page, "Delete..."));
      await session.step(48, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(49, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Sharing the connection shows on its Sharing pane", async () => {
      await session.step(52, "When user picks \"Share...\" from the context menu of \"BDD-Conn-Browser-{run}\" link in gallery", () => pickFromContextMenu(page, "Share...", el(session.text("\"BDD-Conn-Browser-{run}\" link in gallery"))));
      await session.step(53, "Then \"Share BDD-Conn-Browser-{run}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDD-Conn-Browser-{run}\" dialog")), "visible"));
      await session.step(54, "And \"Share BDD-Conn-Browser-{run}\" dialog should contain text \"Full access\"", () => shouldContainText(page, el(session.text("\"Share BDD-Conn-Browser-{run}\" dialog")), "Full access"));
      await session.step(55, "When user picks the sharing user in \"User, group, or email\" input in \"Share BDD-Conn-Browser-{run}\" dialog", () => pickSharingUser(page, el(session.text("\"User, group, or email\" input in \"Share BDD-Conn-Browser-{run}\" dialog"))));
      await session.step(56, "And user unchecks \"Send notifications\" input in \"Share BDD-Conn-Browser-{run}\" dialog", () => uncheck(page, el(session.text("\"Send notifications\" input in \"Share BDD-Conn-Browser-{run}\" dialog"))));
      await session.step(57, "And user clicks on OK button in \"Share BDD-Conn-Browser-{run}\" dialog", () => clickOn(page, el(session.text("OK button in \"Share BDD-Conn-Browser-{run}\" dialog"))));
      await session.step(58, "Then the \"Share BDD-Conn-Browser-{run}\" dialog should close", () => dialogCloses(page, session.text("Share BDD-Conn-Browser-{run}")));
      await session.step(59, "When user clicks on \"BDD-Conn-Browser-{run}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-Conn-Browser-{run}\" link in gallery"))));
      await session.step(60, "Then the sharing pane should list the sharing user", () => sharingPaneLists(page));
      await session.step(61, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The Activity pane records the connection", async () => {
      await session.step(64, "Then Activity section in context panel should be present", () => shouldBe(page, el("Activity section in context panel"), "present"));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A chat message is posted on the connection and removed with it", async () => {
      await session.step(68, "Given Chats section in context panel is expanded", () => isExpanded(page, el("Chats section in context panel")));
      await session.step(69, "When user types \"BDD chat {run}\" into chat post input in context panel", () => typeInto(page, session.text("BDD chat {run}"), el("chat post input in context panel")));
      await session.step(70, "And user presses Enter in chat post input in context panel", () => pressKeyIn(page, "Enter", el("chat post input in context panel")));
      await session.step(71, "Then Chats section in context panel should contain the text \"BDD chat {run}\"", () => shouldContainText(page, el("Chats section in context panel"), session.text("BDD chat {run}")));
      await session.step(72, "And the \"BDD-Conn-Browser-{run}\" connection should have a chat on the server", () => connectionHasChat(page, session.text("BDD-Conn-Browser-{run}")));
      await session.step(73, "When user deletes the chat of the \"BDD-Conn-Browser-{run}\" connection", () => deleteChatOfConnection(page, session.text("BDD-Conn-Browser-{run}")));
      await session.step(74, "Then the \"BDD-Conn-Browser-{run}\" connection should have no chat on the server", () => connectionHasNoChat(page, session.text("BDD-Conn-Browser-{run}")));
      await session.step(75, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
