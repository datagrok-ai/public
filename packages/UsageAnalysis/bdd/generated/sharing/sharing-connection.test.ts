/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sharing/sharing-connection.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sharing.share-dialog, core#977]
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
import {clickOn, isExpanded, shouldBe, shouldContainText, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, connectionOnServer, contextPanelOpen, contextPanelShows, dialogCloses, pickSharingUser, sharingPaneLists, sharingUserNoNotifications} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A connection shared from its Browse node, and what its Sharing pane says", () => {
  const session = feature(test, "features/sharing/sharing-connection.feature", import.meta.url);
  test("A connection shared from its Browse node, and what its Sharing pane says", {tag: ["@journey", "@realizes:sharing.share-dialog", "@known-failure:53", "@realizes:core#977"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the sharing user has no notifications, now and when the feature ends", () => sharingUserNoNotifications(page));
    await session.step(26, "And a \"Postgres\" connection named \"ShareConn1\" is on the server", () => connectionOnServer(page, "Postgres", "ShareConn1"));
    await session.step(27, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await run.scenario("The owner shares the connection from its Browse node", async () => {
      await session.step(31, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(32, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(33, "When user picks \"Share...\" from the context menu of Databases---Postgres---ShareConn1 tree node inside browse tree", () => pickFromContextMenu(page, "Share...", el("Databases---Postgres---ShareConn1 tree node inside browse tree")));
      await session.step(34, "Then \"Share ShareConn1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareConn1\" dialog"), "visible"));
      await session.step(35, "And share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
      await session.step(36, "When user picks the sharing user in \"User, group, or email\" input in \"Share ShareConn1\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share ShareConn1\" dialog")));
      await session.step(37, "And user unchecks \"Send notifications\" input in \"Share ShareConn1\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share ShareConn1\" dialog")));
      await session.step(38, "And user clicks on OK button in \"Share ShareConn1\" dialog", () => clickOn(page, el("OK button in \"Share ShareConn1\" dialog")));
      await session.step(39, "Then the \"Share ShareConn1\" dialog should close", () => dialogCloses(page, "Share ShareConn1"));
      await session.step(40, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(41, "When user clicks on Databases---Postgres---ShareConn1 tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---ShareConn1 tree node inside browse tree")));
      await session.step(42, "Then the context panel should show \"ShareConn1\"", () => contextPanelShows(page, "ShareConn1"));
      await session.step(43, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
      await session.step(44, "And Sharing section in context panel should contain text \"You are the owner\"", () => shouldContainText(page, el("Sharing section in context panel"), "You are the owner"));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
      await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The Sharing pane names the level of the connection's grant (core#977)", async () => {
      await session.step(50, "When user clicks on Databases---Postgres---ShareConn1 tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---ShareConn1 tree node inside browse tree")));
      await session.step(51, "Then the context panel should show \"ShareConn1\"", () => contextPanelShows(page, "ShareConn1"));
      await session.step(52, "And the sharing pane should list the sharing user", () => sharingPaneLists(page));
      await session.step(53, "And Sharing section in context panel should contain text \"can view and use\"", () => shouldContainText(page, el("Sharing section in context panel"), "can view and use"));
    }, {knownFailure: true, line: 53});
    run.finish();
  });
});
