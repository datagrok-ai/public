/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/03-more-viewers/14-network-diagram.feature
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
import {clickOn, doubleClickOn, hoverOver, isExpanded, selectIn, shouldBe, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, browsePanelOpen, contextPanelClosed, contextPanelShows, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, propertyShouldBe, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInAreaSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · More viewers", () => {
  const session = feature(test, "features/howto/03-more-viewers/14-network-diagram.feature", import.meta.url);
  test("Show links between things with a network diagram", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(17, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(18, "And user double-clicks on Files---Demo---got-s1-edges.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---got-s1-edges.csv tree node inside browse tree")));
    await session.step(19, "And the \"got-s1-edges\" view should be current", () => viewIsCurrent(page, "got-s1-edges"));
    await session.step(20, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(23, "When the video says \"A network diagram: who is linked to whom\"", () => videoSays(page, "A network diagram: who is linked to whom"));
    await session.step(24, "And user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
    await session.step(25, "And user types \"network\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "network", el("viewer gallery search in \"Add Viewer\" dialog")));
    await session.step(26, "And user clicks on first \"Network diagram\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Network diagram\" card in \"Add Viewer\" dialog")));
    await session.step(27, "Then network diagram viewer should be visible", () => shouldBe(page, el("network diagram viewer"), "visible"));
    await session.step(29, "When user picks \"Source\" in the column selector at the \"node 1 selector\" area of network diagram viewer", () => pickInAreaSelector(page, "Source", "node 1 selector", el("network diagram viewer")), undefined, "Pick the two columns each link joins");
    await session.step(30, "And user picks \"Target\" in the column selector at the \"node 2 selector\" area of network diagram viewer", () => pickInAreaSelector(page, "Target", "node 2 selector", el("network diagram viewer")));
    await session.step(31, "Then the \"node 2 column\" reading of network diagram viewer should be \"Target\"", () => readingReads(page, "node 2 column", el("network diagram viewer"), "Target"));
    await session.step(33, "When user hovers over network diagram viewer", () => hoverOver(page, el("network diagram viewer")), undefined, "The rest is in the settings");
    await session.step(34, "And user clicks on settings icon of network diagram viewer", () => clickOn(page, el("settings icon of network diagram viewer")));
    await session.step(35, "Then the context panel should show \"Network diagram\"", () => contextPanelShows(page, "Network diagram"));
    await session.step(37, "When user selects \"Weight\" in \"Node1 Size\" property in context panel", () => selectIn(page, "Weight", el("\"Node1 Size\" property in context panel")), undefined, "Data: bigger nodes for characters with more links");
    await session.step(38, "Then \"node1SizeColumnName\" property of network diagram viewer should be \"Weight\"", () => propertyShouldBe(page, "node1SizeColumnName", el("network diagram viewer"), "Weight"));
    await session.step(40, "When user selects \"Weight\" in \"Edge Width\" property in context panel", () => selectIn(page, "Weight", el("\"Edge Width\" property in context panel")), undefined, "and links as wide as the sum of their Weight");
    await session.step(41, "And user selects \"sum\" in \"Edge Width Aggr Type\" property in context panel", () => selectIn(page, "sum", el("\"Edge Width Aggr Type\" property in context panel")));
    await session.step(42, "Then \"edgeWidthAggrType\" property of network diagram viewer should be \"sum\"", () => propertyShouldBe(page, "edgeWidthAggrType", el("network diagram viewer"), "sum"));
    await session.step(44, "When user clicks on the \"node \\\"Eddard Stark\\\"\" area of network diagram viewer", () => clickArea(page, "node \"Eddard Stark\"", el("network diagram viewer")), undefined, "Click a character to select the rows of their links");
    await session.step(45, "Then 33 rows should be selected", () => selectedRowCount(page, 33));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
  });
});
