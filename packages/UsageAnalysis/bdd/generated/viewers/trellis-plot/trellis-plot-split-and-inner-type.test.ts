/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-split-and-inner-type.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, eventFired, hasArea, hasNoArea, listenFor, noBalloons, noErrors, readingAsRemembered, readingAtLeast, readingDiffers, readingIs, readingNotAsRemembered, readingReads, readingsDiffer, rememberReading, setProperties, setProperty, showsRows, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {cellsWideTall, innerPropertyShouldBe, pickInnerViewer, setInnerProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot split columns and the inner viewer", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-split-and-inner-type.feature", import.meta.url);
  test("Trellis plot split columns and the inner viewer", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 22, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(24, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]));
    await session.step(28, "Then the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
    await run.scenario("The viewer is added clean", async () => {
      await session.step(31, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(32, "And the \"error\" reading of trellis plot viewer should be \"\"", () => readingReads(page, "error", el("trellis plot viewer"), ""));
      await session.step(33, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(34, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
      await session.step(35, "And trellis plot viewer should show 1000 rows", () => showsRows(page, el("trellis plot viewer"), 1000));
      await session.step(36, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Two split columns give eight cells", async () => {
      await session.step(40, "Then the \"x categories\" reading of trellis plot viewer should be 2", () => readingIs(page, "x categories", el("trellis plot viewer"), 2));
      await session.step(41, "And the \"y categories\" reading of trellis plot viewer should be 4", () => readingIs(page, "y categories", el("trellis plot viewer"), 4));
      await session.step(42, "And the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(43, "And the \"x categories shown\" reading of trellis plot viewer should be 2", () => readingIs(page, "x categories shown", el("trellis plot viewer"), 2));
      await session.step(44, "And the \"y categories shown\" reading of trellis plot viewer should be 4", () => readingIs(page, "y categories shown", el("trellis plot viewer"), 4));
      await session.step(45, "And trellis plot viewer should have a \"cell F | Caucasian\" area", () => hasArea(page, el("trellis plot viewer"), "cell F | Caucasian"));
      await session.step(46, "And trellis plot viewer should have a \"cell M | Asian\" area", () => hasArea(page, el("trellis plot viewer"), "cell M | Asian"));
      await session.step(47, "And trellis plot viewer should have a \"cell 1,1\" area", () => hasArea(page, el("trellis plot viewer"), "cell 1,1"));
      await session.step(48, "And trellis plot viewer should not have a \"cell F | Klingon\" area", () => hasNoArea(page, el("trellis plot viewer"), "cell F | Klingon"));
      await session.step(49, "And the \"distinct cell signatures\" reading of trellis plot viewer should be 8", () => readingIs(page, "distinct cell signatures", el("trellis plot viewer"), 8));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second X split column multiplies the cells and dropping it brings them back", async () => {
      await session.step(53, "When user sets \"X Column Names\" property of trellis plot viewer to \"SEX, CONTROL\"", () => setProperty(page, "X Column Names", el("trellis plot viewer"), "SEX, CONTROL"));
      await session.step(54, "Then the \"x categories\" reading of trellis plot viewer should be 4", () => readingIs(page, "x categories", el("trellis plot viewer"), 4));
      await session.step(55, "And the \"cells\" reading of trellis plot viewer should be 16", () => readingIs(page, "cells", el("trellis plot viewer"), 16));
      await session.step(56, "And the cells of trellis plot viewer should be 4 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 4, 4));
      await session.step(57, "And trellis plot viewer should have a \"cell F, true | Caucasian\" area", () => hasArea(page, el("trellis plot viewer"), "cell F, true | Caucasian"));
      await session.step(58, "When user sets \"X Column Names\" property of trellis plot viewer to \"SEX\"", () => setProperty(page, "X Column Names", el("trellis plot viewer"), "SEX"));
      await session.step(59, "Then the \"x categories\" reading of trellis plot viewer should be 2", () => readingIs(page, "x categories", el("trellis plot viewer"), 2));
      await session.step(60, "And the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(61, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bar chart picked in the control panel announces itself and draws a picture per cell [type=Bar chart, pictures=8]", async () => {
      await session.step(65, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(66, "When user picks \"Bar chart\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Bar chart", el("trellis plot viewer")));
      await session.step(67, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(68, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Bar chart\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Bar chart"));
      await session.step(69, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(70, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 8", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 8));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Histogram picked in the control panel announces itself and draws a picture per cell [type=Histogram, pictures=8]", async () => {
      await session.step(65, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(66, "When user picks \"Histogram\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Histogram", el("trellis plot viewer")));
      await session.step(67, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(68, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Histogram\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Histogram"));
      await session.step(69, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(70, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 8", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 8));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Line chart picked in the control panel announces itself and draws a picture per cell [type=Line chart, pictures=8]", async () => {
      await session.step(65, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(66, "When user picks \"Line chart\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Line chart", el("trellis plot viewer")));
      await session.step(67, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(68, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Line chart\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Line chart"));
      await session.step(69, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(70, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 8", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 8));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Pie chart picked in the control panel announces itself and draws a picture per cell [type=Pie chart, pictures=2]", async () => {
      await session.step(65, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(66, "When user picks \"Pie chart\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Pie chart", el("trellis plot viewer")));
      await session.step(67, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(68, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Pie chart\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Pie chart"));
      await session.step(69, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(70, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 2", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 2));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Box plot picked in the control panel announces itself and draws a picture per cell [type=Box plot, pictures=8]", async () => {
      await session.step(65, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(66, "When user picks \"Box plot\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Box plot", el("trellis plot viewer")));
      await session.step(67, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(68, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Box plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Box plot"));
      await session.step(69, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(70, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 8", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 8));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Density plot picked in the control panel announces itself and draws a picture per cell [type=Density plot, pictures=2]", async () => {
      await session.step(65, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(66, "When user picks \"Density plot\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Density plot", el("trellis plot viewer")));
      await session.step(67, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(68, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Density plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Density plot"));
      await session.step(69, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(70, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 2", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 2));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Summary picked in the control panel announces itself and draws a picture per cell [type=Summary, pictures=8]", async () => {
      await session.step(65, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(66, "When user picks \"Summary\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Summary", el("trellis plot viewer")));
      await session.step(67, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(68, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Summary\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Summary"));
      await session.step(69, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(70, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 8", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 8));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Sparklines picked in the control panel announces itself and draws a picture per cell [type=Sparklines, pictures=8]", async () => {
      await session.step(65, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(66, "When user picks \"Sparklines\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Sparklines", el("trellis plot viewer")));
      await session.step(67, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(68, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Sparklines\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Sparklines"));
      await session.step(69, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(70, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 8", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 8));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("PC Plot picked in the control panel announces itself and draws a picture per cell [type=PC Plot, pictures=8]", async () => {
      await session.step(65, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(66, "When user picks \"PC Plot\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "PC Plot", el("trellis plot viewer")));
      await session.step(67, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(68, "And the \"inner viewer type\" reading of trellis plot viewer should be \"PC Plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "PC Plot"));
      await session.step(69, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(70, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 8", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 8));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Heatmap picked in the control panel announces itself and draws a picture per cell [type=Heatmap, pictures=8]", async () => {
      await session.step(65, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(66, "When user picks \"Heatmap\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Heatmap", el("trellis plot viewer")));
      await session.step(67, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(68, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Heatmap\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Heatmap"));
      await session.step(69, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(70, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 8", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 8));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Scatter plot picked in the control panel announces itself and draws a picture per cell [type=Scatter plot, pictures=8]", async () => {
      await session.step(65, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(66, "When user picks \"Scatter plot\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Scatter plot", el("trellis plot viewer")));
      await session.step(67, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(68, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Scatter plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Scatter plot"));
      await session.step(69, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(70, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 8", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 8));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Summary inside takes its own visualization and redraws every cell [type=Summary, setting=visualization, value=circles]", async () => {
      await session.step(88, "When user picks \"Summary\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Summary", el("trellis plot viewer")));
      await session.step(89, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Summary\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Summary"));
      await session.step(90, "When user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(91, "And user sets \"visualization\" inner property of trellis plot viewer to \"circles\"", () => setInnerProperty(page, "visualization", el("trellis plot viewer"), "circles"));
      await session.step(92, "Then \"visualization\" inner property of trellis plot viewer should be \"circles\"", () => innerPropertyShouldBe(page, "visualization", el("trellis plot viewer"), "circles"));
      await session.step(93, "And the \"cell signature F | Caucasian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(94, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
      await session.step(96, "When user picks \"Scatter plot\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Scatter plot", el("trellis plot viewer")));
      await session.step(97, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Scatter plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Scatter plot"));
    });
    await run.scenario("The Sparklines inside takes its own sparklineType and redraws every cell [type=Sparklines, setting=sparklineType, value=Bar Chart]", async () => {
      await session.step(88, "When user picks \"Sparklines\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Sparklines", el("trellis plot viewer")));
      await session.step(89, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Sparklines\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Sparklines"));
      await session.step(90, "When user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(91, "And user sets \"sparklineType\" inner property of trellis plot viewer to \"Bar Chart\"", () => setInnerProperty(page, "sparklineType", el("trellis plot viewer"), "Bar Chart"));
      await session.step(92, "Then \"sparklineType\" inner property of trellis plot viewer should be \"Bar Chart\"", () => innerPropertyShouldBe(page, "sparklineType", el("trellis plot viewer"), "Bar Chart"));
      await session.step(93, "And the \"cell signature F | Caucasian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(94, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
      await session.step(96, "When user picks \"Scatter plot\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Scatter plot", el("trellis plot viewer")));
      await session.step(97, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Scatter plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Scatter plot"));
    });
    await run.scenario("The PC Plot inside takes its own colorColumnName and redraws every cell [type=PC Plot, setting=colorColumnName, value=SEX]", async () => {
      await session.step(88, "When user picks \"PC Plot\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "PC Plot", el("trellis plot viewer")));
      await session.step(89, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"PC Plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "PC Plot"));
      await session.step(90, "When user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(91, "And user sets \"colorColumnName\" inner property of trellis plot viewer to \"SEX\"", () => setInnerProperty(page, "colorColumnName", el("trellis plot viewer"), "SEX"));
      await session.step(92, "Then \"colorColumnName\" inner property of trellis plot viewer should be \"SEX\"", () => innerPropertyShouldBe(page, "colorColumnName", el("trellis plot viewer"), "SEX"));
      await session.step(93, "And the \"cell signature F | Caucasian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(94, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
      await session.step(96, "When user picks \"Scatter plot\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Scatter plot", el("trellis plot viewer")));
      await session.step(97, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Scatter plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Scatter plot"));
    });
    await run.scenario("The Points viewer inside takes its columns and every cell follows a cut to one", async () => {
      await session.step(106, "When user picks \"Heatmap\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Heatmap", el("trellis plot viewer")));
      await session.step(107, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Heatmap\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Heatmap"));
      await session.step(108, "When user sets \"columnNames\" inner property of trellis plot viewer to \"AGE,HEIGHT,WEIGHT\"", () => setInnerProperty(page, "columnNames", el("trellis plot viewer"), "AGE,HEIGHT,WEIGHT"));
      await session.step(109, "Then \"columnNames\" inner property of trellis plot viewer should be \"AGE,HEIGHT,WEIGHT\"", () => innerPropertyShouldBe(page, "columnNames", el("trellis plot viewer"), "AGE,HEIGHT,WEIGHT"));
      await session.step(110, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(111, "When user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(112, "And user sets \"columnNames\" inner property of trellis plot viewer to \"AGE\"", () => setInnerProperty(page, "columnNames", el("trellis plot viewer"), "AGE"));
      await session.step(113, "Then \"columnNames\" inner property of trellis plot viewer should be \"AGE\"", () => innerPropertyShouldBe(page, "columnNames", el("trellis plot viewer"), "AGE"));
      await session.step(114, "And the \"cell signature F | Caucasian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(115, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(116, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
      await session.step(117, "When user sets \"columnNames\" inner property of trellis plot viewer to \"AGE,HEIGHT,WEIGHT\"", () => setInnerProperty(page, "columnNames", el("trellis plot viewer"), "AGE,HEIGHT,WEIGHT"));
      await session.step(118, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(119, "And no errors should have been logged", () => noErrors(page));
      await session.step(120, "When user picks \"Scatter plot\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Scatter plot", el("trellis plot viewer")));
      await session.step(121, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Scatter plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Scatter plot"));
    });
    await run.scenario("A no-data excursion leaves no stale frame", async () => {
      await session.step(124, "When user sets \"Viewer Type\" property of trellis plot viewer to \"Pie chart\"", () => setProperty(page, "Viewer Type", el("trellis plot viewer"), "Pie chart"));
      await session.step(125, "Then the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(126, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 2", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 2));
      await session.step(127, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Y Column Names",""],["Viewer Type","Bar chart"]]));
      await session.step(130, "Then the \"y categories\" reading of trellis plot viewer should be 1", () => readingIs(page, "y categories", el("trellis plot viewer"), 1));
      await session.step(131, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Y Column Names","RACE"],["Viewer Type","Pie chart"]]));
      await session.step(134, "Then the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(135, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 2", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 2));
      await session.step(136, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
      await session.step(137, "And the \"cell signature F | Caucasian\" and \"cell signature M | Asian\" readings of trellis plot viewer should differ", () => readingsDiffer(page, "cell signature F | Caucasian", "cell signature M | Asian", el("trellis plot viewer")));
      await session.step(138, "And no errors should have been logged", () => noErrors(page));
      await session.step(139, "When user sets \"Viewer Type\" property of trellis plot viewer to \"Scatter plot\"", () => setProperty(page, "Viewer Type", el("trellis plot viewer"), "Scatter plot"));
      await session.step(140, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Scatter plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Scatter plot"));
    });
    await run.scenario("Hiding the control panel keeps the type", async () => {
      await session.step(143, "Given user sets \"Viewer Type\" property of trellis plot viewer to \"Scatter plot\"", () => setProperty(page, "Viewer Type", el("trellis plot viewer"), "Scatter plot"));
      await session.step(144, "And trellis plot viewer should have a \"control panel\" area", () => hasArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(145, "When user sets \"Show Control Panel\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Show Control Panel", el("trellis plot viewer"), "false"));
      await session.step(146, "Then trellis plot viewer should not have a \"control panel\" area", () => hasNoArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(147, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Scatter plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Scatter plot"));
      await session.step(148, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(149, "When user sets \"Show Control Panel\" property of trellis plot viewer to \"true\"", () => setProperty(page, "Show Control Panel", el("trellis plot viewer"), "true"));
      await session.step(150, "Then trellis plot viewer should have a \"control panel\" area", () => hasArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(151, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Gridlines follow the mode and the inner type", async () => {
      await session.step(154, "When user sets \"Show Gridlines\" property of trellis plot viewer to \"always\"", () => setProperty(page, "Show Gridlines", el("trellis plot viewer"), "always"));
      await session.step(155, "Then the \"gridlines\" reading of trellis plot viewer should be \"true\"", () => readingReads(page, "gridlines", el("trellis plot viewer"), "true"));
      await session.step(156, "When user sets \"Show Gridlines\" property of trellis plot viewer to \"never\"", () => setProperty(page, "Show Gridlines", el("trellis plot viewer"), "never"));
      await session.step(157, "Then the \"gridlines\" reading of trellis plot viewer should be \"false\"", () => readingReads(page, "gridlines", el("trellis plot viewer"), "false"));
      await session.step(158, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Show Gridlines","auto"],["Viewer Type","Scatter plot"]]));
      await session.step(161, "Then the \"gridlines\" reading of trellis plot viewer should be \"true\"", () => readingReads(page, "gridlines", el("trellis plot viewer"), "true"));
      await session.step(162, "When user sets \"Viewer Type\" property of trellis plot viewer to \"Bar chart\"", () => setProperty(page, "Viewer Type", el("trellis plot viewer"), "Bar chart"));
      await session.step(163, "Then the \"gridlines\" reading of trellis plot viewer should be \"false\"", () => readingReads(page, "gridlines", el("trellis plot viewer"), "false"));
      await session.step(164, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Viewer Type","Scatter plot"],["Show Gridlines","always"]]));
      await session.step(167, "Then the \"gridlines\" reading of trellis plot viewer should be \"true\"", () => readingReads(page, "gridlines", el("trellis plot viewer"), "true"));
      await session.step(168, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Changing an inner viewer setting redraws every cell", async () => {
      await session.step(171, "When user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(172, "And user sets \"xColumnName\" inner property of trellis plot viewer to \"AGE\"", () => setInnerProperty(page, "xColumnName", el("trellis plot viewer"), "AGE"));
      await session.step(173, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(174, "And the \"cell signature M | Caucasian\" reading of trellis plot viewer should differ from before", () => readingDiffers(page, "cell signature M | Caucasian", el("trellis plot viewer")));
      await session.step(175, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(176, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
      await session.step(177, "When user sets \"xColumnName\" inner property of trellis plot viewer to \"HEIGHT\"", () => setInnerProperty(page, "xColumnName", el("trellis plot viewer"), "HEIGHT"));
      await session.step(178, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(179, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
