/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/scatter-plot-conditional-colors.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {colorCodedAs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, closeContextMenu, hasArea, legendItemColor, legendLists, pickColorSwatch, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Color a scatter plot by conditions, then change one condition's color", () => {
  const session = feature(test, "features/guides/scatter-plot-conditional-colors.feature", import.meta.url);
  test("Switch the color column to conditional coding and recolor one range from the legend", {tag: ["@guide", "@help:visualize/viewers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"],["Color","AGE"]]), [["X","WEIGHT"],["Y","HEIGHT"],["Color","AGE"]]);
    await session.step(19, "Then scatter plot viewer should have a \"color scale\" area", () => hasArea(page, el("scatter plot viewer"), "color scale"));
    await session.step(20, "When user picks \"Color Coding > Conditional\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header AGE", el("grid")));
    await session.step(21, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(22, "Then \"AGE\" column should be color-coded conditionally", () => colorCodedAs(page, "AGE", "conditionally"));
    await session.step(23, "And the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
    await session.step(24, "When user hovers over \"18 - 35.75\" legend item in legend of scatter plot viewer", () => hoverOver(page, el("\"18 - 35.75\" legend item in legend of scatter plot viewer")));
    await session.step(25, "And user clicks on color picker icon", () => clickOn(page, el("color picker icon")));
    await session.step(26, "Then \"18 - 35.75\" dialog should be visible", () => shouldBe(page, el("\"18 - 35.75\" dialog"), "visible"));
    await session.step(27, "When user picks the color \"#9467BD\" in the color picker dialog", () => pickColorSwatch(page, "#9467BD"));
    await session.step(28, "And user clicks on OK button in \"18 - 35.75\" dialog", () => clickOn(page, el("OK button in \"18 - 35.75\" dialog")));
    await session.step(29, "Then the \"18 - 35.75\" item in the legend of scatter plot viewer should be colored \"#9467BD\"", () => legendItemColor(page, "18 - 35.75", el("scatter plot viewer"), "#9467BD"));
    await session.step(30, "And the \"view\" area of scatter plot viewer should contain the color \"#9467BD\"", () => areaColor(page, "view", el("scatter plot viewer"), "#9467BD"));
  });
});
