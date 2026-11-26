/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/15-annotate-charts/04-color-by-conditions.feature
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
import {clickOn, enterInto, pressKey, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {colorCodedAs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, legendItemColor, legendLists, noErrors, pickColorSwatch, pickFromAreaContextMenu, readingsDiffer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Annotations", () => {
  const session = feature(test, "features/howto/15-annotate-charts/04-color-by-conditions.feature", import.meta.url);
  test("Color points by conditions", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(23, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(24, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"],["Color","AGE"]]), [["X","WEIGHT"],["Y","HEIGHT"],["Color","AGE"]]);
    await session.step(28, "And scatter plot viewer takes 60% of the space it shares with the grid", () => viewerShare(page, el("scatter plot viewer"), 60));
    await session.step(29, "When the video says \"Conditional colors: one color for each range of values\"", () => videoSays(page, "Conditional colors: one color for each range of values"));
    await session.step(31, "And user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")), undefined, "Right-click the column's header: Color Coding > Conditional");
    await session.step(32, "And user presses Escape", () => pressKey(page, "Escape"));
    await session.step(34, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"), undefined, "The plot takes the column's colors, one per range");
    await session.step(35, "And the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
    await session.step(37, "When user picks \"Color Coding > Edit...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header AGE", el("grid")), undefined, "Color Coding > Edit... sets the ranges");
    await session.step(38, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
    await session.step(40, "When user clicks on \"Bin...\" link in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"Bin...\" link in \"Color-coding: AGE\" dialog")), undefined, "Bin...: how many ranges");
    await session.step(41, "Then \"Edit\" dialog should be visible", () => shouldBe(page, el("\"Edit\" dialog"), "visible"));
    await session.step(42, "When user enters \"3\" into \"Number of bins\" input in \"Edit\" dialog", () => enterInto(page, "3", el("\"Number of bins\" input in \"Edit\" dialog")));
    await session.step(43, "And user clicks on OK button in \"Edit\" dialog", () => clickOn(page, el("OK button in \"Edit\" dialog")));
    await session.step(44, "Then \"Edit\" dialog should be hidden", () => shouldBe(page, el("\"Edit\" dialog"), "hidden"));
    await session.step(46, "When user clicks on \"Color 2\" button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("\"Color 2\" button in \"Color-coding: AGE\" dialog")), undefined, "Click a color to change it");
    await session.step(47, "Then \"41.67 - 65.33\" dialog should be visible", () => shouldBe(page, el("\"41.67 - 65.33\" dialog"), "visible"));
    await session.step(48, "When user picks the color \"#c7c7c7\" in the color picker dialog", () => pickColorSwatch(page, "#c7c7c7"));
    await session.step(49, "And user clicks on OK button in \"41.67 - 65.33\" dialog", () => clickOn(page, el("OK button in \"41.67 - 65.33\" dialog")));
    await session.step(51, "When user enters \"<40\" into \"Pattern 1\" text input in \"Color-coding: AGE\" dialog", () => enterInto(page, "<40", el("\"Pattern 1\" text input in \"Color-coding: AGE\" dialog")), undefined, "Type each range: under 40, 40 to 65, over 65");
    await session.step(52, "And user enters \"40-65\" into \"Pattern 2\" text input in \"Color-coding: AGE\" dialog", () => enterInto(page, "40-65", el("\"Pattern 2\" text input in \"Color-coding: AGE\" dialog")));
    await session.step(53, "And user enters \">65\" into \"Pattern 3\" text input in \"Color-coding: AGE\" dialog", () => enterInto(page, ">65", el("\"Pattern 3\" text input in \"Color-coding: AGE\" dialog")));
    await session.step(54, "And user clicks on CLOSE button in \"Color-coding: AGE\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: AGE\" dialog")));
    await session.step(56, "Then the legend of scatter plot viewer should list 3 items", () => legendLists(page, el("scatter plot viewer"), 3), undefined, "The grid and the plot follow the new ranges");
    await session.step(57, "And the \"40-65\" item in the legend of scatter plot viewer should be colored \"#c7c7c7\"", () => legendItemColor(page, "40-65", el("scatter plot viewer"), "#c7c7c7"));
    await session.step(58, "And \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
    await session.step(59, "And the \"color of cell 1 of AGE\" and \"color of cell 3 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of AGE", "color of cell 3 of AGE", el("grid")));
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
  });
});
