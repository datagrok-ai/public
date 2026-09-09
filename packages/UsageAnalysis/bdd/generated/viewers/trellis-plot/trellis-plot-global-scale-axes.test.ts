/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-global-scale-axes.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {setInnerProperty} from '../../../bindings/trellis-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, closeContextMenu, dragAreaBy, hasArea, hasNoArea, hoverArea, noErrors, pickFromAreaContextMenu, readingAsRemembered, readingAtLeast, readingDiffers, readingIs, readingNotAsRemembered, readingReads, rememberReading, rightClickArea, setProperties, setProperty, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot global scale, axes and range sliders", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-global-scale-axes.feature", import.meta.url);
  test("Trellis plot global scale, axes and range sliders", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]));
    await session.step(22, "Then the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
    await session.step(23, "And the \"x axis sliders\" reading of trellis plot viewer should be 0", () => readingIs(page, "x axis sliders", el("trellis plot viewer"), 0));
    await run.scenario("Global Scale redraws every cell on every flip", async () => {
      await session.step(26, "Then trellis plot viewer should not have an \"x axis\" area", () => hasNoArea(page, el("trellis plot viewer"), "x axis"));
      await session.step(27, "When user sets \"Global Scale\" property of trellis plot viewer to \"true\"", () => setProperty(page, "Global Scale", el("trellis plot viewer"), "true"));
      await session.step(28, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should differ from before", () => readingDiffers(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(29, "And the \"cell signature M | Asian\" reading of trellis plot viewer should differ from before", () => readingDiffers(page, "cell signature M | Asian", el("trellis plot viewer")));
      await session.step(30, "And trellis plot viewer should have an \"x axis\" area", () => hasArea(page, el("trellis plot viewer"), "x axis"));
      await session.step(31, "And trellis plot viewer should have a \"y axis\" area", () => hasArea(page, el("trellis plot viewer"), "y axis"));
      await session.step(32, "When user sets \"Global Scale\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Global Scale", el("trellis plot viewer"), "false"));
      await session.step(33, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should differ from before", () => readingDiffers(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(34, "And the \"cell signature M | Asian\" reading of trellis plot viewer should differ from before", () => readingDiffers(page, "cell signature M | Asian", el("trellis plot viewer")));
      await session.step(35, "And trellis plot viewer should not have an \"x axis\" area", () => hasNoArea(page, el("trellis plot viewer"), "x axis"));
      await session.step(36, "When user sets \"Global Scale\" property of trellis plot viewer to \"true\"", () => setProperty(page, "Global Scale", el("trellis plot viewer"), "true"));
      await session.step(37, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should differ from before", () => readingDiffers(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(38, "And the \"cell signature M | Asian\" reading of trellis plot viewer should differ from before", () => readingDiffers(page, "cell signature M | Asian", el("trellis plot viewer")));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show X and Y Axes obey Always, Never and Auto", async () => {
      await session.step(42, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Show X Axes","Always"],["Show Y Axes","Always"]]));
      await session.step(45, "Then trellis plot viewer should have an \"x axis\" area", () => hasArea(page, el("trellis plot viewer"), "x axis"));
      await session.step(46, "And trellis plot viewer should have an \"x axis cell 1\" area", () => hasArea(page, el("trellis plot viewer"), "x axis cell 1"));
      await session.step(47, "And trellis plot viewer should have a \"y axis\" area", () => hasArea(page, el("trellis plot viewer"), "y axis"));
      await session.step(48, "And the \"x axis sliders\" reading of trellis plot viewer should be 3", () => readingIs(page, "x axis sliders", el("trellis plot viewer"), 3));
      await session.step(49, "And the \"y axis sliders\" reading of trellis plot viewer should be 5", () => readingIs(page, "y axis sliders", el("trellis plot viewer"), 5));
      await session.step(50, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Show X Axes","Never"],["Show Y Axes","Never"]]));
      await session.step(53, "Then trellis plot viewer should not have an \"x axis\" area", () => hasNoArea(page, el("trellis plot viewer"), "x axis"));
      await session.step(54, "And trellis plot viewer should not have a \"y axis\" area", () => hasNoArea(page, el("trellis plot viewer"), "y axis"));
      await session.step(55, "And the \"x axis sliders\" reading of trellis plot viewer should be 0", () => readingIs(page, "x axis sliders", el("trellis plot viewer"), 0));
      await session.step(56, "And the \"y axis sliders\" reading of trellis plot viewer should be 0", () => readingIs(page, "y axis sliders", el("trellis plot viewer"), 0));
      await session.step(57, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Show X Axes","Auto"],["Show Y Axes","Auto"]]));
      await session.step(60, "Then trellis plot viewer should have an \"x axis\" area", () => hasArea(page, el("trellis plot viewer"), "x axis"));
      await session.step(61, "And the \"x axis sliders\" reading of trellis plot viewer should be 3", () => readingIs(page, "x axis sliders", el("trellis plot viewer"), 3));
      await session.step(62, "And the \"y axis sliders\" reading of trellis plot viewer should be 5", () => readingIs(page, "y axis sliders", el("trellis plot viewer"), 5));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Range Sliders gates the sliders and keeps the strip", async () => {
      await session.step(66, "When user sets \"Show Range Sliders\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Show Range Sliders", el("trellis plot viewer"), "false"));
      await session.step(67, "Then the \"x axis sliders\" reading of trellis plot viewer should be 0", () => readingIs(page, "x axis sliders", el("trellis plot viewer"), 0));
      await session.step(68, "And the \"y axis sliders\" reading of trellis plot viewer should be 0", () => readingIs(page, "y axis sliders", el("trellis plot viewer"), 0));
      await session.step(69, "And trellis plot viewer should have an \"x axis\" area", () => hasArea(page, el("trellis plot viewer"), "x axis"));
      await session.step(70, "When user sets \"Show Range Sliders\" property of trellis plot viewer to \"true\"", () => setProperty(page, "Show Range Sliders", el("trellis plot viewer"), "true"));
      await session.step(71, "Then the \"x axis sliders\" reading of trellis plot viewer should be 3", () => readingIs(page, "x axis sliders", el("trellis plot viewer"), 3));
      await session.step(72, "And the \"y axis sliders\" reading of trellis plot viewer should be 5", () => readingIs(page, "y axis sliders", el("trellis plot viewer"), 5));
      await session.step(73, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Reset Inner Range Sliders is offered only while a slider exists", async () => {
      await session.step(76, "When user right-clicks on the \"view\" area of trellis plot viewer", () => rightClickArea(page, "view", el("trellis plot viewer")));
      await session.step(77, "Then \"Properties...\" menu item in context menu should be visible", () => shouldBe(page, el("\"Properties...\" menu item in context menu"), "visible"));
      await session.step(78, "And \"Reset Inner Range Sliders\" menu item in context menu should be visible", () => shouldBe(page, el("\"Reset Inner Range Sliders\" menu item in context menu"), "visible"));
      await session.step(79, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(80, "And user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Show X Axes","Never"],["Show Y Axes","Never"]]));
      await session.step(83, "And user right-clicks on the \"view\" area of trellis plot viewer", () => rightClickArea(page, "view", el("trellis plot viewer")));
      await session.step(84, "Then \"Properties...\" menu item in context menu should be visible", () => shouldBe(page, el("\"Properties...\" menu item in context menu"), "visible"));
      await session.step(85, "And \"Reset Inner Range Sliders\" menu item in context menu should be absent", () => shouldBe(page, el("\"Reset Inner Range Sliders\" menu item in context menu"), "absent"));
      await session.step(86, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(87, "And user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Show X Axes","Always"],["Show Y Axes","Always"]]));
      await session.step(90, "And user right-clicks on the \"view\" area of trellis plot viewer", () => rightClickArea(page, "view", el("trellis plot viewer")));
      await session.step(91, "Then \"Reset Inner Range Sliders\" menu item in context menu should be visible", () => shouldBe(page, el("\"Reset Inner Range Sliders\" menu item in context menu"), "visible"));
      await session.step(92, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(93, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The shared slider re-bounds every cell and Reset puts them back exactly", async () => {
      await session.step(96, "When user hovers over the \"cell body F | Caucasian\" area of trellis plot viewer", () => hoverArea(page, "cell body F | Caucasian", el("trellis plot viewer")));
      await session.step(97, "Then trellis plot viewer should have an \"x range slider\" area", () => hasArea(page, el("trellis plot viewer"), "x range slider"));
      await session.step(98, "And trellis plot viewer should have an \"x range slider max handle\" area", () => hasArea(page, el("trellis plot viewer"), "x range slider max handle"));
      await session.step(99, "When user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(100, "And user remembers the \"cell signature M | Asian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature M | Asian", el("trellis plot viewer")));
      await session.step(101, "And user drags the \"x range slider max handle\" area of trellis plot viewer by 60 pixels to the left", () => dragAreaBy(page, "x range slider max handle", el("trellis plot viewer"), 60, "left"));
      await session.step(102, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(103, "And the \"cell signature M | Asian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature M | Asian", el("trellis plot viewer")));
      await session.step(104, "When user picks \"Reset Inner Range Sliders\" from the context menu of the \"view\" area of trellis plot viewer", () => pickFromAreaContextMenu(page, "Reset Inner Range Sliders", "view", el("trellis plot viewer")));
      await session.step(105, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(106, "And the \"cell signature M | Asian\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature M | Asian", el("trellis plot viewer")));
      await session.step(107, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The wheel zooms a cell only once Allow Zoom says so", async () => {
      await session.step(110, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Global Scale","false"]]));
      await session.step(112, "And user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(113, "And user scrolls the mouse wheel down over the \"cell body F | Caucasian\" area of trellis plot viewer", () => wheelOverArea(page, "down", "cell body F | Caucasian", el("trellis plot viewer")));
      await session.step(114, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(115, "When user sets \"allowZoom\" inner property of trellis plot viewer to \"true\"", () => setInnerProperty(page, "allowZoom", el("trellis plot viewer"), "true"));
      await session.step(116, "And user scrolls the mouse wheel down over the \"cell body F | Caucasian\" area of trellis plot viewer", () => wheelOverArea(page, "down", "cell body F | Caucasian", el("trellis plot viewer")));
      await session.step(117, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(118, "When user sets \"allowZoom\" inner property of trellis plot viewer to \"false\"", () => setInnerProperty(page, "allowZoom", el("trellis plot viewer"), "false"));
      await session.step(119, "And user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(120, "And user scrolls the mouse wheel down over the \"cell body F | Caucasian\" area of trellis plot viewer", () => wheelOverArea(page, "down", "cell body F | Caucasian", el("trellis plot viewer")));
      await session.step(121, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(122, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The wheel leaves a bar chart cell alone too", async () => {
      await session.step(125, "When user sets \"Viewer Type\" property of trellis plot viewer to \"Bar chart\"", () => setProperty(page, "Viewer Type", el("trellis plot viewer"), "Bar chart"));
      await session.step(126, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Bar chart\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Bar chart"));
      await session.step(127, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 2", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 2));
      await session.step(128, "When user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(129, "And user scrolls the mouse wheel down over the \"cell body F | Caucasian\" area of trellis plot viewer", () => wheelOverArea(page, "down", "cell body F | Caucasian", el("trellis plot viewer")));
      await session.step(130, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(131, "When user sets \"Viewer Type\" property of trellis plot viewer to \"Scatter plot\"", () => setProperty(page, "Viewer Type", el("trellis plot viewer"), "Scatter plot"));
      await session.step(132, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Scatter plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Scatter plot"));
      await session.step(133, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
