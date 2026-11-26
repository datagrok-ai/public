/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/03-more-viewers/08-sankey.feature
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
import {clickOn, doubleClickOn, expand, hoverOver, isExpanded, selectIn, shouldBe, typeInto, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, browsePanelClosed, browsePanelOpen, contextPanelClosed, contextPanelShows, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, propertyShouldBe} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · More viewers", () => {
  const session = feature(test, "features/howto/03-more-viewers/08-sankey.feature", import.meta.url);
  test("Show flows with a Sankey diagram", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(17, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(18, "And user double-clicks on Files---Demo---energy_uk.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---energy_uk.csv tree node inside browse tree")));
    await session.step(19, "And the \"energy_uk\" view should be current", () => viewIsCurrent(page, "energy_uk"));
    await session.step(20, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(23, "When the video says \"A Sankey diagram shows amounts flowing between categories\"", () => videoSays(page, "A Sankey diagram shows amounts flowing between categories"));
    await session.step(24, "And user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")));
    await session.step(25, "And user types \"sankey\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "sankey", el("viewer gallery search in \"Add Viewer\" dialog")));
    await session.step(26, "And user clicks on first \"Sankey\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Sankey\" card in \"Add Viewer\" dialog")));
    await session.step(27, "Then sankey viewer should be visible", () => shouldBe(page, el("sankey viewer"), "visible"));
    await session.step(28, "And \"Source\" property of sankey viewer should be \"source\"", () => propertyShouldBe(page, "Source", el("sankey viewer"), "source"));
    await session.step(29, "When the video says \"It found the source, target and amount columns itself\" at the top", () => videoSaysAt(page, "It found the source, target and amount columns itself", "top"));
    await session.step(31, "And user clicks on the \"node Electricity grid\" area of sankey viewer", () => clickArea(page, "node Electricity grid", el("sankey viewer")), undefined, "Click a node to select the rows of its flows");
    await session.step(32, "Then 18 rows should be selected", () => selectedRowCount(page, 18));
    await session.step(34, "When user hovers over sankey viewer", () => hoverOver(page, el("sankey viewer")), undefined, "The rest is in the settings");
    await session.step(35, "And user clicks on settings icon of sankey viewer", () => clickOn(page, el("settings icon of sankey viewer")));
    await session.step(36, "Then the context panel should show \"Sankey\"", () => contextPanelShows(page, "Sankey"));
    await session.step(38, "When user expands \"Color\" category in context panel", () => expand(page, el("\"Color\" category in context panel")), undefined, "Color: color the links by where they come from");
    await session.step(39, "And user selects \"source\" in \"Link Color\" property in context panel", () => selectIn(page, "source", el("\"Link Color\" property in context panel")));
    await session.step(40, "Then \"Link Color Column Name\" property of sankey viewer should be \"source\"", () => propertyShouldBe(page, "Link Color Column Name", el("sankey viewer"), "source"));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
  });
});
