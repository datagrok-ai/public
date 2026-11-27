/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/14-files/05-search-files.feature
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
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, expand, hoverOver, shouldBe, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, packageInstalled, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Files", () => {
  const session = feature(test, "features/howto/14-files/05-search-files.feature", import.meta.url);
  test("Search files", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(17, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(21, "And user expands Files tree node inside browse tree", () => expand(page, el("Files tree node inside browse tree")));
    await session.step(22, "And user clicks on Files---Demo tree node inside browse tree", () => clickOn(page, el("Files---Demo tree node inside browse tree")));
    await session.step(23, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
    await session.step(24, "When the video says \"A folder has a search box: type part of a name\"", () => videoSays(page, "A folder has a search box: type part of a name"));
    await session.step(25, "And user types \"aspirin\" into gallery search", () => typeInto(page, "aspirin", el("gallery search")));
    await session.step(26, "Then aspirin-2d.sdf link in gallery should be visible", () => shouldBe(page, el("aspirin-2d.sdf link in gallery"), "visible"));
    await session.step(27, "And aspirin-3d.sdf link in gallery should be visible", () => shouldBe(page, el("aspirin-3d.sdf link in gallery"), "visible"));
    await session.step(28, "And aspirin.mol link in gallery should be visible", () => shouldBe(page, el("aspirin.mol link in gallery"), "visible"));
    await session.step(29, "When the video says \"It finds files in the folders inside it too\"", () => videoSays(page, "It finds files in the folders inside it too"));
    await session.step(30, "And user hovers over \"Search in subfolders\" icon inside gallery toolbar", () => hoverOver(page, el("\"Search in subfolders\" icon inside gallery toolbar")));
    await session.step(31, "And user clicks on \"Search in subfolders\" icon inside gallery toolbar", () => clickOn(page, el("\"Search in subfolders\" icon inside gallery toolbar")));
    await session.step(32, "Then aspirin-2d.sdf link in gallery should be absent", () => shouldBe(page, el("aspirin-2d.sdf link in gallery"), "absent"));
    await session.step(33, "When the video says \"Switched off, it looks in this folder only\"", () => videoSays(page, "Switched off, it looks in this folder only"));
    await session.step(34, "And user clicks on \"Search in subfolders\" icon inside gallery toolbar", () => clickOn(page, el("\"Search in subfolders\" icon inside gallery toolbar")));
    await session.step(35, "Then aspirin-2d.sdf link in gallery should be visible", () => shouldBe(page, el("aspirin-2d.sdf link in gallery"), "visible"));
    await session.step(36, "When the video says \"Double-click a file to open it\"", () => videoSays(page, "Double-click a file to open it"));
    await session.step(37, "And user double-clicks on aspirin-2d.sdf link in gallery", () => doubleClickOn(page, el("aspirin-2d.sdf link in gallery")));
    await session.step(38, "Then the table should have 1 row", () => rowCount(page, 1));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
});
