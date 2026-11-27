/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/14-files/02-upload-a-file-to-your-files.feature
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
import {dropToUpload, fileInHome, noFileInHome} from '../../../bindings/howto-files.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, expand, hoverOver, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noTableLeft} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Files", () => {
  const session = feature(test, "features/howto/14-files/02-upload-a-file-to-your-files.feature", import.meta.url);
  test("Upload a file to your files", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And no file \"solubility-assay.csv\" is in the user's files, now or when the feature ends", () => noFileInHome(page, "solubility-assay.csv"));
    await session.step(17, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(21, "And user expands Files tree node inside browse tree", () => expand(page, el("Files tree node inside browse tree")));
    await session.step(22, "And user clicks on Files---My-files tree node inside browse tree", () => clickOn(page, el("Files---My-files tree node inside browse tree")));
    await session.step(23, "Then the \"My files\" view should be current", () => viewIsCurrent(page, "My files"));
    await session.step(24, "When the video says \"Browse > Files > My files: your own folder on the server\"", () => videoSays(page, "Browse > Files > My files: your own folder on the server"));
    await session.step(25, "And the video says \"Drag a file from your computer onto it\"", () => videoSays(page, "Drag a file from your computer onto it"));
    await session.step(26, "And user drops the \"features/howto/14-files/files/solubility-assay.csv\" file of the project onto gallery to upload it", () => dropToUpload(page, "features/howto/14-files/files/solubility-assay.csv", el("gallery")));
    await session.step(27, "Then solubility-assay.csv link in gallery should be visible", () => shouldBe(page, el("solubility-assay.csv link in gallery"), "visible"));
    await session.step(28, "And the file \"solubility-assay.csv\" should be in the user's files on the server", () => fileInHome(page, "solubility-assay.csv"));
    await session.step(29, "When the video says \"It is uploaded: in the folder and in the tree\"", () => videoSays(page, "It is uploaded: in the folder and in the tree"));
    await session.step(30, "And user hovers over solubility-assay.csv link in gallery", () => hoverOver(page, el("solubility-assay.csv link in gallery")));
    await session.step(31, "And user clicks on Files---My-files---solubility-assay.csv tree node inside browse tree", () => clickOn(page, el("Files---My-files---solubility-assay.csv tree node inside browse tree")));
    await session.step(32, "Then the \"solubility-assay\" view should be current", () => viewIsCurrent(page, "solubility-assay"));
    await session.step(33, "And no table should be left in the workspace", () => noTableLeft(page));
    await session.step(34, "When the video says \"One click previews it\"", () => videoSays(page, "One click previews it"));
    await session.step(35, "And user double-clicks on Files---My-files---solubility-assay.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---My-files---solubility-assay.csv tree node inside browse tree")));
    await session.step(36, "Then the table should have 16 rows", () => rowCount(page, 16));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
    await session.step(38, "When the video says \"A double-click opens it as a table\"", () => videoSays(page, "A double-click opens it as a table"));
  });
});
