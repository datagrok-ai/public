/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/14-files/03-organize-your-files.feature
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
import {fileInHome, noFileInHome, noFileInHomeNow} from '../../../bindings/howto-files.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, dragTo, enterInto, expand, shouldBe, shouldHaveValue, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, dialogCloses, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {fixtureInUsersFiles} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Files", () => {
  const session = feature(test, "features/howto/14-files/03-organize-your-files.feature", import.meta.url);
  test("Organize your files", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And no file \"Solubility\" is in the user's files, now or when the feature ends", () => noFileInHome(page, "Solubility"));
    await session.step(16, "And the \"features/howto/14-files/files/solubility-assay.csv\" file of the project is in the user's files as \"solubility-assay.csv\"", () => fixtureInUsersFiles(page, "features/howto/14-files/files/solubility-assay.csv", "solubility-assay.csv"));
    await session.step(17, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(21, "And user expands Files tree node inside browse tree", () => expand(page, el("Files tree node inside browse tree")));
    await session.step(22, "And user clicks on Files---My-files tree node inside browse tree", () => clickOn(page, el("Files---My-files tree node inside browse tree")));
    await session.step(23, "Then the \"My files\" view should be current", () => viewIsCurrent(page, "My files"));
    await session.step(24, "And solubility-assay.csv link in gallery should be visible", () => shouldBe(page, el("solubility-assay.csv link in gallery"), "visible"));
    await session.step(25, "When the video says \"Browse > Files > My files: your own folder on the server\"", () => videoSays(page, "Browse > Files > My files: your own folder on the server"));
    await session.step(26, "And the video says \"Right-click the folder view to make a folder\"", () => videoSays(page, "Right-click the folder view to make a folder"));
    await session.step(27, "And user picks \"Create folder...\" from the context menu of gallery", () => pickFromContextMenu(page, "Create folder...", el("gallery")));
    await session.step(28, "Then \"Create folder\" dialog should be visible", () => shouldBe(page, el("\"Create folder\" dialog"), "visible"));
    await session.step(29, "When user enters \"Solubility\" into Name input in \"Create folder\" dialog", () => enterInto(page, "Solubility", el("Name input in \"Create folder\" dialog")));
    await session.step(30, "And user clicks on OK button in \"Create folder\" dialog", () => clickOn(page, el("OK button in \"Create folder\" dialog")));
    await session.step(31, "Then the \"Create folder\" dialog should close", () => dialogCloses(page, "Create folder"));
    await session.step(32, "And Solubility link in gallery should be visible", () => shouldBe(page, el("Solubility link in gallery"), "visible"));
    await session.step(33, "When the video says \"Drag a file onto a folder to move it there\"", () => videoSays(page, "Drag a file onto a folder to move it there"));
    await session.step(34, "And user expands Files---My-files tree node inside browse tree", () => expand(page, el("Files---My-files tree node inside browse tree")));
    await session.step(35, "And user drags Files---My-files---solubility-assay.csv tree node inside browse tree to Files---My-files---Solubility tree node inside browse tree", () => dragTo(page, el("Files---My-files---solubility-assay.csv tree node inside browse tree"), el("Files---My-files---Solubility tree node inside browse tree")));
    await session.step(36, "Then the file \"Solubility/solubility-assay.csv\" should be in the user's files on the server", () => fileInHome(page, "Solubility/solubility-assay.csv"));
    await session.step(37, "And the file \"solubility-assay.csv\" should not be in the user's files on the server", () => noFileInHomeNow(page, "solubility-assay.csv"));
    await session.step(38, "When user expands Files---My-files---Solubility tree node inside browse tree", () => expand(page, el("Files---My-files---Solubility tree node inside browse tree")));
    await session.step(39, "Then Files---My-files---Solubility---solubility-assay.csv tree node inside browse tree should be visible", () => shouldBe(page, el("Files---My-files---Solubility---solubility-assay.csv tree node inside browse tree"), "visible"));
    await session.step(40, "When the video says \"Right-click a file to rename it\"", () => videoSays(page, "Right-click a file to rename it"));
    await session.step(41, "And user picks \"Rename...\" from the context menu of Files---My-files---Solubility---solubility-assay.csv tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el("Files---My-files---Solubility---solubility-assay.csv tree node inside browse tree")));
    await session.step(42, "Then \"File name\" input in Rename dialog should have value \"solubility-assay.csv\"", () => shouldHaveValue(page, el("\"File name\" input in Rename dialog"), "solubility-assay.csv"));
    await session.step(43, "When user enters \"solubility-ph-study.csv\" into \"File name\" input in Rename dialog", () => enterInto(page, "solubility-ph-study.csv", el("\"File name\" input in Rename dialog")));
    await session.step(44, "And user clicks on OK button in Rename dialog", () => clickOn(page, el("OK button in Rename dialog")));
    await session.step(45, "Then the \"Rename\" dialog should close", () => dialogCloses(page, "Rename"));
    await session.step(46, "And the file \"Solubility/solubility-ph-study.csv\" should be in the user's files on the server", () => fileInHome(page, "Solubility/solubility-ph-study.csv"));
    await session.step(47, "And the file \"Solubility/solubility-assay.csv\" should not be in the user's files on the server", () => noFileInHomeNow(page, "Solubility/solubility-assay.csv"));
    await session.step(48, "And Files---My-files---Solubility---solubility-ph-study.csv tree node inside browse tree should be visible", () => shouldBe(page, el("Files---My-files---Solubility---solubility-ph-study.csv tree node inside browse tree"), "visible"));
    await session.step(49, "And no errors should have been logged", () => noErrors(page));
    await session.step(50, "When the video says \"The file has its new name, in its new folder\"", () => videoSays(page, "The file has its new name, in its new folder"));
  });
});
