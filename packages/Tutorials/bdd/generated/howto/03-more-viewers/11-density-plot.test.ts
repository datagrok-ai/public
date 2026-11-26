/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/03-more-viewers/11-density-plot.feature
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
import {clickOn, enterInto, expand, hoverOver, selectIn, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, simpleModeOff, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noErrors, readingIs, readingReads, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · More viewers", () => {
  const session = feature(test, "features/howto/03-more-viewers/11-density-plot.feature", import.meta.url);
  test("See where points crowd with a density plot", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And user opens earthquakes dataset", () => openDataset(page, ds("earthquakes")));
    await session.step(15, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(16, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(17, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(18, "When the video says \"When points pile up, count them in bins\"", () => videoSays(page, "When points pile up, count them in bins"));
    await session.step(19, "And user clicks on density plot icon on toolbox", () => clickOn(page, el("density plot icon on toolbox")));
    await session.step(20, "Then density plot viewer should be added to the open tableview", () => viewerAdded(page, "density plot"));
    await session.step(22, "When user clicks on toolbox tab", () => clickOn(page, el("toolbox tab")), undefined, "Hide the toolbox to give the plot room");
    await session.step(23, "Then toolbox should be hidden", () => shouldBe(page, el("toolbox"), "hidden"));
    await session.step(25, "When user picks \"Longitude\" in the \"x\" column selector of density plot viewer", () => pickInColumnSelector(page, "Longitude", "x", el("density plot viewer")), undefined, "Pick the columns right on the plot");
    await session.step(26, "And user picks \"Latitude\" in the \"y\" column selector of density plot viewer", () => pickInColumnSelector(page, "Latitude", "y", el("density plot viewer")));
    await session.step(27, "Then the \"x column\" reading of density plot viewer should be \"Longitude\"", () => readingReads(page, "x column", el("density plot viewer"), "Longitude"));
    await session.step(29, "When user hovers over the \"densest bin\" area of density plot viewer", () => hoverArea(page, "densest bin", el("density plot viewer")), undefined, "Hover a bin for its range and count");
    await session.step(30, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
    await session.step(32, "When user hovers over density plot viewer", () => hoverOver(page, el("density plot viewer")), undefined, "The rest is in the settings");
    await session.step(33, "And user clicks on settings icon of density plot viewer", () => clickOn(page, el("settings icon of density plot viewer")));
    await session.step(34, "Then the context panel should show \"Density plot\"", () => contextPanelShows(page, "Density plot"));
    await session.step(36, "When user expands \"Misc\" category in context panel", () => expand(page, el("\"Misc\" category in context panel")), undefined, "Misc: more bins, finer detail");
    await session.step(37, "And user enters \"100\" into \"Bins\" property in context panel", () => enterInto(page, "100", el("\"Bins\" property in context panel")));
    await session.step(38, "Then the \"bins\" reading of density plot viewer should be 100", () => readingIs(page, "bins", el("density plot viewer"), 100));
    await session.step(40, "When user expands \"Color\" category in context panel", () => expand(page, el("\"Color\" category in context panel")), undefined, "Color: a log scale brings out the sparse areas too");
    await session.step(41, "And user selects \"logarithmic\" in \"Color Transform Type\" property in context panel", () => selectIn(page, "logarithmic", el("\"Color Transform Type\" property in context panel")));
    await session.step(42, "Then the \"color transform\" reading of density plot viewer should be \"logarithmic\"", () => readingReads(page, "color transform", el("density plot viewer"), "logarithmic"));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
  });
});
