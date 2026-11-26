/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/03-more-viewers/10-tile-viewer.feature
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
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {dragTileToLane} from '../../../bindings/howto-more-viewers.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, selectIn, shouldBe, typeInto, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dockToViewEdge, dockedAtViewEdge} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · More viewers", () => {
  const session = feature(test, "features/howto/03-more-viewers/10-tile-viewer.feature", import.meta.url);
  test("Show rows as cards with the Tile Viewer and the Form", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(19, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(20, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "When the video says \"The Tile Viewer shows every row as a card\"", () => videoSays(page, "The Tile Viewer shows every row as a card"));
    await session.step(23, "And user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
    await session.step(24, "And user types \"tile\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "tile", el("viewer gallery search in \"Add Viewer\" dialog")));
    await session.step(25, "And user clicks on first \"Tile Viewer\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Tile Viewer\" card in \"Add Viewer\" dialog")));
    await session.step(26, "Then Tile Viewer viewer should be visible", () => shouldBe(page, el("Tile Viewer viewer"), "visible"));
    await session.step(28, "When user docks Tile Viewer viewer to the bottom edge of the view", () => dockToViewEdge(page, el("Tile Viewer viewer"), "bottom"), undefined, "Drag its title bar to the bottom, so the lanes have room");
    await session.step(29, "Then Tile Viewer viewer should be docked along the bottom edge of the view", () => dockedAtViewEdge(page, el("Tile Viewer viewer"), "bottom"));
    await session.step(31, "When user hovers over Tile Viewer viewer", () => hoverOver(page, el("Tile Viewer viewer")), undefined, "Split the cards into lanes by a column");
    await session.step(32, "And user clicks on settings icon of Tile Viewer viewer", () => clickOn(page, el("settings icon of Tile Viewer viewer")));
    await session.step(33, "And user selects \"SEVERITY\" in \"Lanes\" property in context panel", () => selectIn(page, "SEVERITY", el("\"Lanes\" property in context panel")));
    await session.step(34, "Then the \"lanes\" reading of Tile Viewer viewer should be 5", () => readingIs(page, "lanes", el("Tile Viewer viewer"), 5));
    await session.step(36, "When user drags the card of row 1 of Tile Viewer viewer to the \"High\" lane", () => dragTileToLane(page, 1, el("Tile Viewer viewer"), "High"), undefined, "Drag a card to another lane to change its value");
    await session.step(37, "Then the value of \"SEVERITY\" column in row 1 should be \"High\"", () => valueInRow(page, "SEVERITY", 1, "High"));
    await session.step(39, "When user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")), undefined, "Hide the context panel for more room");
    await session.step(40, "Then context panel should be hidden", () => shouldBe(page, el("context panel"), "hidden"));
    await session.step(41, "When the video says \"The Form viewer shows one row at a time\" at the top", () => videoSaysAt(page, "The Form viewer shows one row at a time", "top"));
    await session.step(42, "And user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
    await session.step(43, "And user types \"form\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "form", el("viewer gallery search in \"Add Viewer\" dialog")));
    await session.step(44, "And user clicks on first \"Form\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Form\" card in \"Add Viewer\" dialog")));
    await session.step(45, "Then form viewer should be visible", () => shouldBe(page, el("form viewer"), "visible"));
    await session.step(47, "When user clicks on the \"tile of row 2\" area of Tile Viewer viewer", () => clickArea(page, "tile of row 2", el("Tile Viewer viewer")), undefined, "Click a card: the form shows its row");
    await session.step(48, "Then the \"row\" reading of form viewer should be 2", () => readingIs(page, "row", el("form viewer"), 2));
    await session.step(50, "When user clicks on the \"next row\" area of form viewer", () => clickArea(page, "next row", el("form viewer")), undefined, "The arrows step through the rows");
    await session.step(51, "Then the \"row\" reading of form viewer should be 3", () => readingIs(page, "row", el("form viewer"), 3));
    await session.step(52, "And no errors should have been logged", () => noErrors(page));
  });
});
