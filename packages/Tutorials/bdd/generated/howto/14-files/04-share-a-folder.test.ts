/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/14-files/04-share-a-folder.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {noFileInHome, pickSharingUserShown} from '../../../bindings/howto-files.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, hoverOver, pressKey, shouldBe, shouldContainText, uncheck, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, connectionsOnServer, contextPanelClosed, dialogCloses, ownConnectionGone, refreshBrowse, sharingUserLevel, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {fixtureInUsersFiles} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {infoBalloonText, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Files", () => {
  const session = feature(test, "features/howto/14-files/04-share-a-folder.feature", import.meta.url);
  test("Share a folder", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And simple mode is off", () => simpleModeOff(page));
    await session.step(21, "And the user's own connection \"MyFilesAssays\" is removed now and at feature end", () => ownConnectionGone(page, "MyFilesAssays"));
    await session.step(22, "And no file \"Assays\" is in the user's files, now or when the feature ends", () => noFileInHome(page, "Assays"));
    await session.step(23, "And the \"features/howto/14-files/files/solubility-assay.csv\" file of the project is in the user's files as \"Assays/solubility-assay.csv\"", () => fixtureInUsersFiles(page, "features/howto/14-files/files/solubility-assay.csv", "Assays/solubility-assay.csv"));
    await session.step(24, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(25, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(26, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(27, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(28, "And user expands Files tree node inside browse tree", () => expand(page, el("Files tree node inside browse tree")));
    await session.step(29, "And user clicks on Files---My-files tree node inside browse tree", () => clickOn(page, el("Files---My-files tree node inside browse tree")));
    await session.step(30, "And user expands Files---My-files tree node inside browse tree", () => expand(page, el("Files---My-files tree node inside browse tree")));
    await session.step(31, "Then Files---My-files---Assays tree node inside browse tree should be visible", () => shouldBe(page, el("Files---My-files---Assays tree node inside browse tree"), "visible"));
    await session.step(32, "When the video says \"Files are shared by folder: right-click the folder\"", () => videoSays(page, "Files are shared by folder: right-click the folder"));
    await session.step(33, "And user picks \"Share folder...\" from the context menu of Files---My-files---Assays tree node inside browse tree", () => pickFromContextMenu(page, "Share folder...", el("Files---My-files---Assays tree node inside browse tree")));
    await session.step(34, "Then \"Share folder\" dialog should be visible", () => shouldBe(page, el("\"Share folder\" dialog"), "visible"));
    await session.step(35, "When the video says \"Pick a colleague or a group\"", () => videoSays(page, "Pick a colleague or a group"));
    await session.step(36, "And user picks the sharing user from the suggestions of \"User, group, or email\" input in \"Share folder\" dialog", () => pickSharingUserShown(page, el("\"User, group, or email\" input in \"Share folder\" dialog")));
    await session.step(37, "And the video says \"View and use: they can read the files, not change them\"", () => videoSays(page, "View and use: they can read the files, not change them"));
    await session.step(38, "And user clicks on share access selector", () => clickOn(page, el("share access selector")));
    await session.step(39, "Then \"View and use\" tree node inside privilege tree should be checked", () => shouldBe(page, el("\"View and use\" tree node inside privilege tree"), "checked"));
    await session.step(40, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(41, "Then share access selector should contain text \"View and use\"", () => shouldContainText(page, el("share access selector"), "View and use"));
    await session.step(43, "When user unchecks \"Send notifications\" input in \"Share folder\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share folder\" dialog")), undefined, "Notifications are optional; off here");
    await session.step(44, "And user clicks on OK button in \"Share folder\" dialog", () => clickOn(page, el("OK button in \"Share folder\" dialog")));
    await session.step(45, "Then the \"Share folder\" dialog should close", () => dialogCloses(page, "Share folder"));
    await session.step(46, "And an info balloon containing \"Shared\" should have been shown", () => infoBalloonText(page, "Shared"));
    await session.step(47, "And 1 connection named \"MyFilesAssays\" should be on the server", () => connectionsOnServer(page, 1, "MyFilesAssays"));
    await session.step(48, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(49, "Then \"My files/Assays\" tree node inside browse tree should be visible", () => shouldBe(page, el("\"My files/Assays\" tree node inside browse tree"), "visible"));
    await session.step(50, "When the video says \"The shared folder now shows under Files too\"", () => videoSays(page, "The shared folder now shows under Files too"));
    await session.step(51, "And user hovers over \"My files/Assays\" tree node inside browse tree", () => hoverOver(page, el("\"My files/Assays\" tree node inside browse tree")));
    await session.step(52, "And user expands Files---My-files tree node inside browse tree", () => expand(page, el("Files---My-files tree node inside browse tree")));
    await session.step(53, "And user picks \"Share folder...\" from the context menu of Files---My-files---Assays tree node inside browse tree", () => pickFromContextMenu(page, "Share folder...", el("Files---My-files---Assays tree node inside browse tree")));
    await session.step(54, "Then \"Share folder\" dialog should be visible", () => shouldBe(page, el("\"Share folder\" dialog"), "visible"));
    await session.step(55, "And the access level of the sharing user in \"Share folder\" dialog should be \"View and use\"", () => sharingUserLevel(page, el("\"Share folder\" dialog"), "View and use"));
    await session.step(56, "When the video says \"The folder's sharing lists who has access\"", () => videoSays(page, "The folder's sharing lists who has access"));
    await session.step(57, "And user clicks on CANCEL button in \"Share folder\" dialog", () => clickOn(page, el("CANCEL button in \"Share folder\" dialog")));
    await session.step(58, "Then the \"Share folder\" dialog should close", () => dialogCloses(page, "Share folder"));
    await session.step(59, "And no errors should have been logged", () => noErrors(page));
  });
});
