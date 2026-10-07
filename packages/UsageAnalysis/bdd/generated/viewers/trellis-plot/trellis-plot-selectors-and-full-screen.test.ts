/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-selectors-and-full-screen.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {canvasColors} from '@datagrok-libraries/bdd/bindings/common/pixels';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, hasArea, hasNoArea, hoverArea, noErrors, propertyShouldBe, readingIs, resizeTo, restoreSize, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot selectors, control panel and the full-screen cell", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-selectors-and-full-screen.feature", import.meta.url);
  test("Trellis plot selectors, control panel and the full-screen cell", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(22, "Then the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
    await session.step(23, "And trellis plot viewer should have a \"x selectors\" area", () => hasArea(page, el("trellis plot viewer"), "x selectors"));
    await session.step(24, "And trellis plot viewer should have a \"y selectors\" area", () => hasArea(page, el("trellis plot viewer"), "y selectors"));
    await session.step(25, "And trellis plot viewer should have a \"control panel\" area", () => hasArea(page, el("trellis plot viewer"), "control panel"));
    await run.scenario("Each selector strip and the control panel leave the screen and come back", async () => {
      await session.step(28, "When user sets \"Show X Selectors\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Show X Selectors", el("trellis plot viewer"), "false"));
      await session.step(29, "Then trellis plot viewer should not have a \"x selectors\" area", () => hasNoArea(page, el("trellis plot viewer"), "x selectors"));
      await session.step(30, "And trellis plot viewer should have a \"y selectors\" area", () => hasArea(page, el("trellis plot viewer"), "y selectors"));
      await session.step(31, "When user sets \"Show Y Selectors\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Show Y Selectors", el("trellis plot viewer"), "false"));
      await session.step(32, "Then trellis plot viewer should not have a \"y selectors\" area", () => hasNoArea(page, el("trellis plot viewer"), "y selectors"));
      await session.step(33, "When user sets \"Show Control Panel\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Show Control Panel", el("trellis plot viewer"), "false"));
      await session.step(34, "Then trellis plot viewer should not have a \"control panel\" area", () => hasNoArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(35, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Show X Selectors","true"],["Show Y Selectors","true"],["Show Control Panel","true"]]), [["Show X Selectors","true"],["Show Y Selectors","true"],["Show Control Panel","true"]]);
      await session.step(39, "Then trellis plot viewer should have a \"x selectors\" area", () => hasArea(page, el("trellis plot viewer"), "x selectors"));
      await session.step(40, "And trellis plot viewer should have a \"y selectors\" area", () => hasArea(page, el("trellis plot viewer"), "y selectors"));
      await session.step(41, "And trellis plot viewer should have a \"control panel\" area", () => hasArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(42, "And the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An X selector strip switched off stays off through Auto Layout's shrink and restore", async () => {
      await session.step(46, "Then \"Auto Layout\" property of trellis plot viewer should be \"true\"", () => propertyShouldBe(page, "Auto Layout", el("trellis plot viewer"), "true"));
      await session.step(47, "When user sets \"Show X Selectors\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Show X Selectors", el("trellis plot viewer"), "false"));
      await session.step(48, "And user resizes trellis plot viewer to 240 by 240", () => resizeTo(page, el("trellis plot viewer"), 240, 240));
      await session.step(49, "Then trellis plot viewer should not have a \"control panel\" area", () => hasNoArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(50, "And trellis plot viewer should not have a \"y selectors\" area", () => hasNoArea(page, el("trellis plot viewer"), "y selectors"));
      await session.step(51, "When user restores the size of trellis plot viewer", () => restoreSize(page, el("trellis plot viewer")));
      await session.step(52, "Then trellis plot viewer should have a \"control panel\" area", () => hasArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(53, "And trellis plot viewer should have a \"y selectors\" area", () => hasArea(page, el("trellis plot viewer"), "y selectors"));
      await session.step(54, "And trellis plot viewer should not have a \"x selectors\" area", () => hasNoArea(page, el("trellis plot viewer"), "x selectors"));
      await session.step(55, "When user sets \"Show X Selectors\" property of trellis plot viewer to \"true\"", () => setProperty(page, "Show X Selectors", el("trellis plot viewer"), "true"));
      await session.step(56, "Then trellis plot viewer should have a \"x selectors\" area", () => hasArea(page, el("trellis plot viewer"), "x selectors"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The full-screen icon of a cell opens its viewer on its own", async () => {
      await session.step(60, "When user hovers over the \"cell body F | Caucasian\" area of trellis plot viewer", () => hoverArea(page, "cell body F | Caucasian", el("trellis plot viewer")));
      await session.step(61, "Then trellis plot viewer should have a \"full screen icon\" area", () => hasArea(page, el("trellis plot viewer"), "full screen icon"));
      await session.step(62, "When user clicks on the \"full screen icon\" area of trellis plot viewer", () => clickArea(page, "full screen icon", el("trellis plot viewer")));
      await session.step(63, "Then \"SEX: F, RACE: Caucasian\" dialog should be visible", () => shouldBe(page, el("\"SEX: F, RACE: Caucasian\" dialog"), "visible"));
      await session.step(64, "And the canvases of \"SEX: F, RACE: Caucasian\" dialog should be painted in at least 2 colors", () => canvasColors(page, el("\"SEX: F, RACE: Caucasian\" dialog"), 2));
      await session.step(65, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(66, "Then \"SEX: F, RACE: Caucasian\" dialog should be absent", () => shouldBe(page, el("\"SEX: F, RACE: Caucasian\" dialog"), "absent"));
      await session.step(67, "And the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
