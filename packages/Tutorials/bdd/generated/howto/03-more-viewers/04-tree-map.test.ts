/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/03-more-viewers/04-tree-map.feature
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
import {clickOn, hoverOver, selectIn, shouldBe, typeInto, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverArea, noErrors, readingIs, readingReads, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dockToViewEdge, dockedAtViewEdge, pickInColumnSelector, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · More viewers", () => {
  const session = feature(test, "features/howto/03-more-viewers/04-tree-map.feature", import.meta.url);
  test("Show a hierarchy with a tree map or a sunburst", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(21, "When the video says \"A tree map shows categories inside categories\"", () => videoSays(page, "A tree map shows categories inside categories"));
    await session.step(22, "And user clicks on tree map icon on toolbox", () => clickOn(page, el("tree map icon on toolbox")));
    await session.step(23, "Then tree map viewer should be added to the open tableview", () => viewerAdded(page, "tree map"));
    await session.step(25, "When user clicks on toolbox tab", () => clickOn(page, el("toolbox tab")), undefined, "Hide the toolbox to give the map room");
    await session.step(26, "Then toolbox should be hidden", () => shouldBe(page, el("toolbox"), "hidden"));
    await session.step(28, "When user docks tree map viewer to the right edge of the view", () => dockToViewEdge(page, el("tree map viewer"), "right"), undefined, "Drag its title bar to dock it on the right");
    await session.step(29, "Then tree map viewer should be docked along the right edge of the view", () => dockedAtViewEdge(page, el("tree map viewer"), "right"));
    await session.step(31, "When user selects \"SEX\" in first tree map split selector", () => selectIn(page, "SEX", el("first tree map split selector")), undefined, "Pick the top level here");
    await session.step(33, "And user selects \"DIS_POP\" in second tree map split selector", () => selectIn(page, "DIS_POP", el("second tree map split selector")), undefined, "and add the next level beside it");
    await session.step(34, "Then the \"levels\" reading of tree map viewer should be 2", () => readingIs(page, "levels", el("tree map viewer"), 2));
    await session.step(36, "When user picks \"WEIGHT\" in the \"color\" column selector of tree map viewer", () => pickInColumnSelector(page, "WEIGHT", "color", el("tree map viewer")), undefined, "Color the rectangles by the average of a column");
    await session.step(37, "Then the \"color column\" reading of tree map viewer should be \"WEIGHT\"", () => readingReads(page, "color column", el("tree map viewer"), "WEIGHT"));
    await session.step(39, "When user hovers over the \"leaf F | RA\" area of tree map viewer", () => hoverArea(page, "leaf F | RA", el("tree map viewer")), undefined, "Each rectangle is as large as its share of the rows");
    await session.step(41, "And user clicks on the \"leaf F | RA\" area of tree map viewer", () => clickArea(page, "leaf F | RA", el("tree map viewer")), undefined, "Click a rectangle to select its rows");
    await session.step(42, "Then 1941 rows should be selected", () => selectedRowCount(page, 1941));
    await session.step(43, "When the video says \"A sunburst shows the same levels as rings\" at the top", () => videoSaysAt(page, "A sunburst shows the same levels as rings", "top"));
    await session.step(44, "And user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
    await session.step(45, "And user types \"sunburst\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "sunburst", el("viewer gallery search in \"Add Viewer\" dialog")));
    await session.step(46, "And user clicks on first \"Sunburst\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Sunburst\" card in \"Add Viewer\" dialog")));
    await session.step(47, "Then sunburst viewer should be visible", () => shouldBe(page, el("sunburst viewer"), "visible"));
    await session.step(49, "When user hovers over sunburst viewer", () => hoverOver(page, el("sunburst viewer")), undefined, "Its levels are set in the settings");
    await session.step(50, "And user clicks on settings icon of sunburst viewer", () => clickOn(page, el("settings icon of sunburst viewer")));
    await session.step(51, "And user clicks on \"...\" button in \"Hierarchy\" property in context panel", () => clickOn(page, el("\"...\" button in \"Hierarchy\" property in context panel")));
    await session.step(52, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(53, "When user clicks on \"None\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"None\" link in \"Select columns...\" dialog")));
    await session.step(54, "And user toggles the \"DIS_POP\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "DIS_POP", el("\"Select columns...\" dialog")));
    await session.step(55, "And user toggles the \"SEX\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "SEX", el("\"Select columns...\" dialog")));
    await session.step(56, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(57, "Then the \"hierarchy columns\" reading of sunburst viewer should be \"SEX, DIS_POP\"", () => readingReads(page, "hierarchy columns", el("sunburst viewer"), "SEX, DIS_POP"));
    await session.step(59, "When user clicks on the \"segment M | Psoriasis\" area of sunburst viewer", () => clickArea(page, "segment M | Psoriasis", el("sunburst viewer")), undefined, "Click a segment to select its rows");
    await session.step(60, "Then 820 rows should be selected", () => selectedRowCount(page, 820));
    await session.step(61, "And no errors should have been logged", () => noErrors(page));
  });
});
