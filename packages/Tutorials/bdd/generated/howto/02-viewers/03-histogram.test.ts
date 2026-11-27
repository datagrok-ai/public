/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/02-viewers/03-histogram.feature
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
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {viewerShare} from '../../../bindings/howto-viewers.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, dragSliderTo, expand, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPassesAll, filterPassesFewer} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, doubleClickArea, dragAreaBy, hasArea, hoverArea, noErrors, propertyShouldBe, readingBetween} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Viewers", () => {
  const session = feature(test, "features/howto/02-viewers/03-histogram.feature", import.meta.url);
  test("Configure a histogram", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(18, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(19, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(20, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(21, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","WEIGHT"]]), [["Value","WEIGHT"]]);
    await session.step(23, "And histogram viewer takes 65% of the space it shares with the grid", () => viewerShare(page, el("histogram viewer"), 65));
    await session.step(24, "When the video says \"A histogram shows how the values of one column are spread\"", () => videoSays(page, "A histogram shows how the values of one column are spread"));
    await session.step(26, "And user picks \"AGE\" in the \"value\" column selector of histogram viewer", () => pickInColumnSelector(page, "AGE", "value", el("histogram viewer")), undefined, "Pick the column under the plot");
    await session.step(27, "Then \"Value\" property of histogram viewer should be \"AGE\"", () => propertyShouldBe(page, "Value", el("histogram viewer"), "AGE"));
    await session.step(29, "When user drags the slider of histogram viewer to 12", () => dragSliderTo(page, el("histogram viewer"), 12), undefined, "The slider at the top sets the number of bins");
    await session.step(30, "Then the \"bins shown\" reading of histogram viewer should be between 10 and 14", () => readingBetween(page, "bins shown", el("histogram viewer"), 10, 14));
    await session.step(32, "When user hovers over the \"bin 6\" area of histogram viewer", () => hoverArea(page, "bin 6", el("histogram viewer")), undefined, "Hover a bar to see how many rows it holds");
    await session.step(34, "When user hovers over the \"view\" area of histogram viewer", () => hoverArea(page, "view", el("histogram viewer")), undefined, "The handles under the bars filter the table");
    await session.step(35, "And user drags the \"range max handle\" area of histogram viewer by 400 pixels to the left", () => dragAreaBy(page, "range max handle", el("histogram viewer"), 400, "left"));
    await session.step(36, "Then fewer than 5850 rows should pass the filter", () => filterPassesFewer(page, 5850));
    await session.step(38, "When user double-clicks on the \"range slider\" area of histogram viewer", () => doubleClickArea(page, "range slider", el("histogram viewer")), undefined, "Double-click the slider to show all rows again");
    await session.step(39, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(41, "When user picks \"SEX\" in the \"split\" column selector of histogram viewer", () => pickInColumnSelector(page, "SEX", "split", el("histogram viewer")), undefined, "Split by a category: one line for each");
    await session.step(42, "Then histogram viewer should have a \"line F\" area", () => hasArea(page, el("histogram viewer"), "line F"));
    await session.step(43, "And histogram viewer should have a \"line M\" area", () => hasArea(page, el("histogram viewer"), "line M"));
    await session.step(45, "When user hovers over the \"x axis\" area of histogram viewer", () => hoverArea(page, "x axis", el("histogram viewer")), undefined, "The gear opens every setting, grouped in sections");
    await session.step(46, "And user clicks on settings icon of histogram viewer", () => clickOn(page, el("settings icon of histogram viewer")));
    await session.step(47, "Then the context panel should show \"Histogram\"", () => contextPanelShows(page, "Histogram"));
    await session.step(49, "When user expands \"Category\" category in context panel", () => expand(page, el("\"Category\" category in context panel")), undefined, "Category: stack the categories in each bar");
    await session.step(50, "And user checks \"Split Stack\" property in context panel", () => check(page, el("\"Split Stack\" property in context panel")));
    await session.step(51, "Then histogram viewer should have a \"bin 8 | F\" area", () => hasArea(page, el("histogram viewer"), "bin 8 | F"));
    await session.step(53, "When user expands \"X Axis\" category in context panel", () => expand(page, el("\"X Axis\" category in context panel")), undefined, "X Axis: write the count over each bar");
    await session.step(54, "And user checks \"Show Values\" property in context panel", () => check(page, el("\"Show Values\" property in context panel")));
    await session.step(55, "Then histogram viewer should have a \"bin labels\" area", () => hasArea(page, el("histogram viewer"), "bin labels"));
    await session.step(56, "And no errors should have been logged", () => noErrors(page));
  });
});
