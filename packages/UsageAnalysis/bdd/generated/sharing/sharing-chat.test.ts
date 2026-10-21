/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sharing/sharing-chat.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [sharing.share-dialog, views.queries]
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
import {clearField, clickOn, enterInto, isExpanded, pressKeyIn, shouldBe, shouldContainText, shouldHaveValue, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, expandShared, pickSharingUser, queryOnServer, refreshBrowse, runningAccountSignedIn, sharedShouldBe, sharingUserSignedIn, signInAsSelf, signInAsSharingUser, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Two people in the chat of a shared entity", () => {
  const session = feature(test, "features/sharing/sharing-chat.feature", import.meta.url);
  test("Two people in the chat of a shared entity", {tag: ["@journey", "@serial", "@realizes:sharing.share-dialog", "@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And a query \"ShareChat1\" on \"System:Datagrok\" reads \"select 1 as x\"", () => queryOnServer(page, "ShareChat1", "System:Datagrok", "select 1 as x"));
    await session.step(27, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(29, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(30, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
    await session.step(31, "When user clicks on Platform---Functions---Queries tree node inside browse tree", () => clickOn(page, el("Platform---Functions---Queries tree node inside browse tree")));
    await session.step(32, "Then the \"Queries\" view should be current", () => viewIsCurrent(page, "Queries"));
    await session.step(33, "When user enters \"ShareChat1\" into gallery search", () => enterInto(page, "ShareChat1", el("gallery search")));
    await session.step(34, "Then \"ShareChat1\" gallery card should be visible", () => shouldBe(page, el("\"ShareChat1\" gallery card"), "visible"));
    await run.scenario("The owner posts in the query's chat and shares the query", async () => {
      await session.step(37, "When user clicks on \"ShareChat1\" gallery card", () => clickOn(page, el("\"ShareChat1\" gallery card")));
      await session.step(38, "Then the context panel should show \"ShareChat1\"", () => contextPanelShows(page, "ShareChat1"));
      await session.step(39, "Given Chats section in context panel is expanded", () => isExpanded(page, el("Chats section in context panel")));
      await session.step(40, "When user types \"First comment from A {run}\" into chat post input in context panel", () => typeInto(page, session.text("First comment from A {run}"), el("chat post input in context panel")));
      await session.step(41, "And user presses Enter in chat post input in context panel", () => pressKeyIn(page, "Enter", el("chat post input in context panel")));
      await session.step(42, "Then Chats section in context panel should contain the text \"First comment from A {run}\"", () => shouldContainText(page, el("Chats section in context panel"), session.text("First comment from A {run}")));
      await session.step(43, "When user picks \"Share...\" from the context menu of \"ShareChat1\" gallery card", () => pickFromContextMenu(page, "Share...", el("\"ShareChat1\" gallery card")));
      await session.step(44, "Then \"Share ShareChat1\" dialog should be visible", () => shouldBe(page, el("\"Share ShareChat1\" dialog"), "visible"));
      await session.step(45, "When user picks the sharing user in \"User, group, or email\" input in \"Share ShareChat1\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share ShareChat1\" dialog")));
      await session.step(46, "And user unchecks \"Send notifications\" input in \"Share ShareChat1\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share ShareChat1\" dialog")));
      await session.step(47, "And user clicks on OK button in \"Share ShareChat1\" dialog", () => clickOn(page, el("OK button in \"Share ShareChat1\" dialog")));
      await session.step(48, "Then the \"Share ShareChat1\" dialog should close", () => dialogCloses(page, "Share ShareChat1"));
      await session.step(49, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
      await session.step(51, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The second account sees the post and replies", async () => {
      await session.step(54, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(55, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(56, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(57, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(58, "When user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(59, "Given \"My stuff\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff\" tree node inside browse tree")));
      await session.step(60, "And \"My stuff > Shared with me\" tree node inside browse tree is expanded", () => isExpanded(page, el("\"My stuff > Shared with me\" tree node inside browse tree")));
      await session.step(61, "When user expands \".\" shared by the running account", () => expandShared(page, "."));
      await session.step(62, "Then \"ShareChat1\" shared by the running account should be visible", () => sharedShouldBe(page, "ShareChat1", "visible"));
      await session.step(63, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
      await session.step(64, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
      await session.step(65, "When user clicks on Platform---Functions---Queries tree node inside browse tree", () => clickOn(page, el("Platform---Functions---Queries tree node inside browse tree")));
      await session.step(66, "Then the \"Queries\" view should be current", () => viewIsCurrent(page, "Queries"));
      await session.step(67, "When user enters \"ShareChat1\" into gallery search", () => enterInto(page, "ShareChat1", el("gallery search")));
      await session.step(68, "And user clicks on \"ShareChat1\" gallery card", () => clickOn(page, el("\"ShareChat1\" gallery card")));
      await session.step(69, "Then the context panel should show \"ShareChat1\"", () => contextPanelShows(page, "ShareChat1"));
      await session.step(70, "Given Chats section in context panel is expanded", () => isExpanded(page, el("Chats section in context panel")));
      await session.step(71, "Then Chats section in context panel should contain the text \"First comment from A {run}\"", () => shouldContainText(page, el("Chats section in context panel"), session.text("First comment from A {run}")));
      await session.step(72, "When user types \"Reply from B {run}\" into chat post input in context panel", () => typeInto(page, session.text("Reply from B {run}"), el("chat post input in context panel")));
      await session.step(73, "And user presses Enter in chat post input in context panel", () => pressKeyIn(page, "Enter", el("chat post input in context panel")));
      await session.step(74, "Then Chats section in context panel should contain the text \"Reply from B {run}\"", () => shouldContainText(page, el("Chats section in context panel"), session.text("Reply from B {run}")));
      await session.step(75, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(76, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
      await session.step(78, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The owner sees the reply", async () => {
      await session.step(81, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(82, "Then the running account should be signed in", () => runningAccountSignedIn(page));
      await session.step(83, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(84, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(85, "And Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
      await session.step(86, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
      await session.step(87, "When user clicks on Platform---Functions---Queries tree node inside browse tree", () => clickOn(page, el("Platform---Functions---Queries tree node inside browse tree")));
      await session.step(88, "Then the \"Queries\" view should be current", () => viewIsCurrent(page, "Queries"));
      await session.step(89, "When user enters \"ShareChat1\" into gallery search", () => enterInto(page, "ShareChat1", el("gallery search")));
      await session.step(90, "And user clicks on \"ShareChat1\" gallery card", () => clickOn(page, el("\"ShareChat1\" gallery card")));
      await session.step(91, "Then the context panel should show \"ShareChat1\"", () => contextPanelShows(page, "ShareChat1"));
      await session.step(92, "Given Chats section in context panel is expanded", () => isExpanded(page, el("Chats section in context panel")));
      await session.step(93, "Then Chats section in context panel should contain the text \"First comment from A {run}\"", () => shouldContainText(page, el("Chats section in context panel"), session.text("First comment from A {run}")));
      await session.step(94, "And Chats section in context panel should contain the text \"Reply from B {run}\"", () => shouldContainText(page, el("Chats section in context panel"), session.text("Reply from B {run}")));
      await session.step(95, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(96, "Then gallery search should have value \"\"", () => shouldHaveValue(page, el("gallery search"), ""));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
