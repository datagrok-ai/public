/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/14-files/07-write-a-text-file.feature
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
import {fileInHomeHolds, noFileInHome} from '../../../bindings/howto-files.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, expand, shouldBe, typeAtCaret, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, dialogCloses, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Files", () => {
  const session = feature(test, "features/howto/14-files/07-write-a-text-file.feature", import.meta.url);
  test("Write a text file", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And no file \"assay-notes.txt\" is in the user's files, now or when the feature ends", () => noFileInHome(page, "assay-notes.txt"));
    await session.step(15, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(16, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(17, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(18, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(19, "And user expands Files tree node inside browse tree", () => expand(page, el("Files tree node inside browse tree")));
    await session.step(20, "And user clicks on Files---My-files tree node inside browse tree", () => clickOn(page, el("Files---My-files tree node inside browse tree")));
    await session.step(21, "Then the \"My files\" view should be current", () => viewIsCurrent(page, "My files"));
    await session.step(22, "When the video says \"Right-click the folder view to make a new file\"", () => videoSays(page, "Right-click the folder view to make a new file"));
    await session.step(23, "And user picks \"Create file...\" from the context menu of gallery", () => pickFromContextMenu(page, "Create file...", el("gallery")));
    await session.step(24, "Then \"Create file\" dialog should be visible", () => shouldBe(page, el("\"Create file\" dialog"), "visible"));
    await session.step(25, "When user enters \"assay-notes.txt\" into Name input in \"Create file\" dialog", () => enterInto(page, "assay-notes.txt", el("Name input in \"Create file\" dialog")));
    await session.step(26, "And user clicks on OK button in \"Create file\" dialog", () => clickOn(page, el("OK button in \"Create file\" dialog")));
    await session.step(27, "Then the \"Create file\" dialog should close", () => dialogCloses(page, "Create file"));
    await session.step(28, "And assay-notes.txt link in gallery should be visible", () => shouldBe(page, el("assay-notes.txt link in gallery"), "visible"));
    await session.step(29, "When user expands Files---My-files tree node inside browse tree", () => expand(page, el("Files---My-files tree node inside browse tree")));
    await session.step(30, "And user clicks on Files---My-files---assay-notes.txt tree node inside browse tree", () => clickOn(page, el("Files---My-files---assay-notes.txt tree node inside browse tree")));
    await session.step(31, "Then the \"assay-notes.txt\" view should be current", () => viewIsCurrent(page, "assay-notes.txt"));
    await session.step(32, "When the video says \"A text file opens in an editor\"", () => videoSays(page, "A text file opens in an editor"));
    await session.step(33, "And user clicks on code editor", () => clickOn(page, el("code editor")));
    await session.step(34, "And user types \"Solubility assay, plate 3: re-run CMP-106 at pH 6.5\" at the caret", () => typeAtCaret(page, "Solubility assay, plate 3: re-run CMP-106 at pH 6.5"));
    await session.step(35, "Then SAVE button should be visible", () => shouldBe(page, el("SAVE button"), "visible"));
    await session.step(36, "When the video says \"SAVE shows once the text changes\"", () => videoSays(page, "SAVE shows once the text changes"));
    await session.step(37, "And user clicks on SAVE button", () => clickOn(page, el("SAVE button")));
    await session.step(38, "Then the file \"assay-notes.txt\" in the user's files should hold the text \"re-run CMP-106 at pH 6.5\"", () => fileInHomeHolds(page, "assay-notes.txt", "re-run CMP-106 at pH 6.5"));
    await session.step(39, "And SAVE button should be hidden", () => shouldBe(page, el("SAVE button"), "hidden"));
    await session.step(40, "And no errors should have been logged", () => noErrors(page));
    await session.step(41, "When the video says \"Saved to your files on the server\"", () => videoSays(page, "Saved to your files on the server"));
  });
});
