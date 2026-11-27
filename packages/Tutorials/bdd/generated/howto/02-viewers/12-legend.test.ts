/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/02-viewers/12-legend.feature
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
import {clickOn, clickOnHolding, hoverOver, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, legendItemColor, legendLists, noErrors, pickColorSwatch, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Viewers", () => {
  const session = feature(test, "features/howto/02-viewers/12-legend.feature", import.meta.url);
  test("Use the legend", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(19, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(20, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(21, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"],["Color","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["X","WEIGHT"],["Y","HEIGHT"],["Color","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await session.step(27, "And scatter plot viewer takes 65% of the space it shares with the grid", () => viewerShare(page, el("scatter plot viewer"), 65));
    await session.step(28, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
    await session.step(29, "When the video says \"The legend names the colors on a viewer\"", () => videoSays(page, "The legend names the colors on a viewer"));
    await session.step(31, "And user clicks on \"Other\" legend item in legend of scatter plot viewer", () => clickOn(page, el("\"Other\" legend item in legend of scatter plot viewer")), undefined, "Click a category to show only it on this plot");
    await session.step(32, "Then scatter plot viewer should show 59 rows", () => showsRows(page, el("scatter plot viewer"), 59));
    await session.step(33, "And \"Other\" legend item in legend of scatter plot viewer should be selected", () => shouldBe(page, el("\"Other\" legend item in legend of scatter plot viewer"), "selected"));
    await session.step(35, "And all rows should pass the filter", () => filterPassesAll(page), undefined, "The table keeps all its rows");
    await session.step(37, "When user clicks on \"Black\" legend item in legend of scatter plot viewer holding Control", () => clickOnHolding(page, el("\"Black\" legend item in legend of scatter plot viewer"), "Control"), undefined, "Ctrl+click adds another category");
    await session.step(38, "Then scatter plot viewer should show 84 rows", () => showsRows(page, el("scatter plot viewer"), 84));
    await session.step(40, "When user clicks on cross of \"Black\" legend item in legend of scatter plot viewer", () => clickOn(page, el("cross of \"Black\" legend item in legend of scatter plot viewer")), undefined, "The x takes one back; a click on the last one shows all");
    await session.step(41, "Then scatter plot viewer should show 59 rows", () => showsRows(page, el("scatter plot viewer"), 59));
    await session.step(42, "When user clicks on \"Other\" legend item in legend of scatter plot viewer", () => clickOn(page, el("\"Other\" legend item in legend of scatter plot viewer")));
    await session.step(43, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(45, "When user hovers over \"Caucasian\" legend item in legend of scatter plot viewer", () => hoverOver(page, el("\"Caucasian\" legend item in legend of scatter plot viewer")), undefined, "Point at an item and click the palette to recolor it");
    await session.step(46, "And user clicks on color picker icon", () => clickOn(page, el("color picker icon")));
    await session.step(47, "Then \"Caucasian\" dialog should be visible", () => shouldBe(page, el("\"Caucasian\" dialog"), "visible"));
    await session.step(48, "When user picks the color \"#9467BD\" in the color picker dialog", () => pickColorSwatch(page, "#9467BD"));
    await session.step(49, "And user clicks on OK button in \"Caucasian\" dialog", () => clickOn(page, el("OK button in \"Caucasian\" dialog")));
    await session.step(50, "Then the \"Caucasian\" item in the legend of scatter plot viewer should be colored \"#9467BD\"", () => legendItemColor(page, "Caucasian", el("scatter plot viewer"), "#9467BD"));
    await session.step(51, "And the \"view\" area of scatter plot viewer should contain the color \"#9467BD\"", () => areaColor(page, "view", el("scatter plot viewer"), "#9467BD"));
    await session.step(52, "And no errors should have been logged", () => noErrors(page));
  });
});
