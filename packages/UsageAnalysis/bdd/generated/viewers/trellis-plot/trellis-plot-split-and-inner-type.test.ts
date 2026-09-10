/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-split-and-inner-type.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, eventFired, hasArea, hasNoArea, listenFor, noBalloons, noErrors, readingAsRemembered, readingAtLeast, readingDiffers, readingIs, readingNotAsRemembered, readingReads, readingsDiffer, rememberReading, setProperties, setProperty, showsRows, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {cellsWideTall, pickInnerViewer, setInnerProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot split columns and the inner viewer", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-split-and-inner-type.feature", import.meta.url);
  test("Trellis plot split columns and the inner viewer", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 13, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]));
    await session.step(22, "Then the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
    await run.scenario("The viewer is added clean", async () => {
      await session.step(25, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(26, "And the \"error\" reading of trellis plot viewer should be \"\"", () => readingReads(page, "error", el("trellis plot viewer"), ""));
      await session.step(27, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(28, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
      await session.step(29, "And trellis plot viewer should show 1000 rows", () => showsRows(page, el("trellis plot viewer"), 1000));
      await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Two split columns give eight cells", async () => {
      await session.step(34, "Then the \"x categories\" reading of trellis plot viewer should be 2", () => readingIs(page, "x categories", el("trellis plot viewer"), 2));
      await session.step(35, "And the \"y categories\" reading of trellis plot viewer should be 4", () => readingIs(page, "y categories", el("trellis plot viewer"), 4));
      await session.step(36, "And the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(37, "And the \"x categories shown\" reading of trellis plot viewer should be 2", () => readingIs(page, "x categories shown", el("trellis plot viewer"), 2));
      await session.step(38, "And the \"y categories shown\" reading of trellis plot viewer should be 4", () => readingIs(page, "y categories shown", el("trellis plot viewer"), 4));
      await session.step(39, "And trellis plot viewer should have a \"cell F | Caucasian\" area", () => hasArea(page, el("trellis plot viewer"), "cell F | Caucasian"));
      await session.step(40, "And trellis plot viewer should have a \"cell M | Asian\" area", () => hasArea(page, el("trellis plot viewer"), "cell M | Asian"));
      await session.step(41, "And trellis plot viewer should have a \"cell 1,1\" area", () => hasArea(page, el("trellis plot viewer"), "cell 1,1"));
      await session.step(42, "And trellis plot viewer should not have a \"cell F | Klingon\" area", () => hasNoArea(page, el("trellis plot viewer"), "cell F | Klingon"));
      await session.step(43, "And the \"distinct cell signatures\" reading of trellis plot viewer should be 8", () => readingIs(page, "distinct cell signatures", el("trellis plot viewer"), 8));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bar chart picked in the control panel announces itself and draws a picture per cell [type=Bar chart, pictures=8]", async () => {
      await session.step(47, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(48, "When user picks \"Bar chart\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Bar chart", el("trellis plot viewer")));
      await session.step(49, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(50, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Bar chart\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Bar chart"));
      await session.step(51, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(52, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 8", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 8));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Histogram picked in the control panel announces itself and draws a picture per cell [type=Histogram, pictures=8]", async () => {
      await session.step(47, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(48, "When user picks \"Histogram\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Histogram", el("trellis plot viewer")));
      await session.step(49, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(50, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Histogram\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Histogram"));
      await session.step(51, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(52, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 8", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 8));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Line chart picked in the control panel announces itself and draws a picture per cell [type=Line chart, pictures=8]", async () => {
      await session.step(47, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(48, "When user picks \"Line chart\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Line chart", el("trellis plot viewer")));
      await session.step(49, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(50, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Line chart\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Line chart"));
      await session.step(51, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(52, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 8", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 8));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Pie chart picked in the control panel announces itself and draws a picture per cell [type=Pie chart, pictures=2]", async () => {
      await session.step(47, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(48, "When user picks \"Pie chart\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Pie chart", el("trellis plot viewer")));
      await session.step(49, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(50, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Pie chart\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Pie chart"));
      await session.step(51, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(52, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 2", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 2));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Box plot picked in the control panel announces itself and draws a picture per cell [type=Box plot, pictures=8]", async () => {
      await session.step(47, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(48, "When user picks \"Box plot\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Box plot", el("trellis plot viewer")));
      await session.step(49, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(50, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Box plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Box plot"));
      await session.step(51, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(52, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 8", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 8));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Density plot picked in the control panel announces itself and draws a picture per cell [type=Density plot, pictures=2]", async () => {
      await session.step(47, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(48, "When user picks \"Density plot\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Density plot", el("trellis plot viewer")));
      await session.step(49, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(50, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Density plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Density plot"));
      await session.step(51, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(52, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 2", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 2));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Scatter plot picked in the control panel announces itself and draws a picture per cell [type=Scatter plot, pictures=8]", async () => {
      await session.step(47, "Given user listens for \"d4-trellis-plot-viewer-type-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(48, "When user picks \"Scatter plot\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Scatter plot", el("trellis plot viewer")));
      await session.step(49, "Then \"d4-trellis-plot-viewer-type-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-viewer-type-changed", el("trellis plot viewer")));
      await session.step(50, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Scatter plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Scatter plot"));
      await session.step(51, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(52, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 8", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 8));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A no-data excursion leaves no stale frame", async () => {
      await session.step(66, "When user sets \"Viewer Type\" property of trellis plot viewer to \"Pie chart\"", () => setProperty(page, "Viewer Type", el("trellis plot viewer"), "Pie chart"));
      await session.step(67, "Then the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(68, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 2", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 2));
      await session.step(69, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Y Column Names",""],["Viewer Type","Bar chart"]]));
      await session.step(72, "Then the \"y categories\" reading of trellis plot viewer should be 1", () => readingIs(page, "y categories", el("trellis plot viewer"), 1));
      await session.step(73, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Y Column Names","RACE"],["Viewer Type","Pie chart"]]));
      await session.step(76, "Then the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(77, "And the \"distinct cell signatures\" reading of trellis plot viewer should be at least 2", () => readingAtLeast(page, "distinct cell signatures", el("trellis plot viewer"), 2));
      await session.step(78, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
      await session.step(79, "And the \"cell signature F | Caucasian\" and \"cell signature M | Asian\" readings of trellis plot viewer should differ", () => readingsDiffer(page, "cell signature F | Caucasian", "cell signature M | Asian", el("trellis plot viewer")));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
      await session.step(81, "When user sets \"Viewer Type\" property of trellis plot viewer to \"Scatter plot\"", () => setProperty(page, "Viewer Type", el("trellis plot viewer"), "Scatter plot"));
      await session.step(82, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Scatter plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Scatter plot"));
    });
    await run.scenario("Hiding the control panel keeps the type", async () => {
      await session.step(85, "Given user sets \"Viewer Type\" property of trellis plot viewer to \"Scatter plot\"", () => setProperty(page, "Viewer Type", el("trellis plot viewer"), "Scatter plot"));
      await session.step(86, "And trellis plot viewer should have a \"control panel\" area", () => hasArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(87, "When user sets \"Show Control Panel\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Show Control Panel", el("trellis plot viewer"), "false"));
      await session.step(88, "Then trellis plot viewer should not have a \"control panel\" area", () => hasNoArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(89, "And the \"inner viewer type\" reading of trellis plot viewer should be \"Scatter plot\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Scatter plot"));
      await session.step(90, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(91, "When user sets \"Show Control Panel\" property of trellis plot viewer to \"true\"", () => setProperty(page, "Show Control Panel", el("trellis plot viewer"), "true"));
      await session.step(92, "Then trellis plot viewer should have a \"control panel\" area", () => hasArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(93, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Gridlines follow the mode and the inner type", async () => {
      await session.step(96, "When user sets \"Show Gridlines\" property of trellis plot viewer to \"always\"", () => setProperty(page, "Show Gridlines", el("trellis plot viewer"), "always"));
      await session.step(97, "Then the \"gridlines\" reading of trellis plot viewer should be \"true\"", () => readingReads(page, "gridlines", el("trellis plot viewer"), "true"));
      await session.step(98, "When user sets \"Show Gridlines\" property of trellis plot viewer to \"never\"", () => setProperty(page, "Show Gridlines", el("trellis plot viewer"), "never"));
      await session.step(99, "Then the \"gridlines\" reading of trellis plot viewer should be \"false\"", () => readingReads(page, "gridlines", el("trellis plot viewer"), "false"));
      await session.step(100, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Show Gridlines","auto"],["Viewer Type","Scatter plot"]]));
      await session.step(103, "Then the \"gridlines\" reading of trellis plot viewer should be \"true\"", () => readingReads(page, "gridlines", el("trellis plot viewer"), "true"));
      await session.step(104, "When user sets \"Viewer Type\" property of trellis plot viewer to \"Bar chart\"", () => setProperty(page, "Viewer Type", el("trellis plot viewer"), "Bar chart"));
      await session.step(105, "Then the \"gridlines\" reading of trellis plot viewer should be \"false\"", () => readingReads(page, "gridlines", el("trellis plot viewer"), "false"));
      await session.step(106, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Viewer Type","Scatter plot"],["Show Gridlines","always"]]));
      await session.step(109, "Then the \"gridlines\" reading of trellis plot viewer should be \"true\"", () => readingReads(page, "gridlines", el("trellis plot viewer"), "true"));
      await session.step(110, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Changing an inner viewer setting redraws every cell", async () => {
      await session.step(113, "When user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(114, "And user sets \"xColumnName\" inner property of trellis plot viewer to \"AGE\"", () => setInnerProperty(page, "xColumnName", el("trellis plot viewer"), "AGE"));
      await session.step(115, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(116, "And the \"cell signature M | Caucasian\" reading of trellis plot viewer should differ from before", () => readingDiffers(page, "cell signature M | Caucasian", el("trellis plot viewer")));
      await session.step(117, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(118, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
      await session.step(119, "When user sets \"xColumnName\" inner property of trellis plot viewer to \"HEIGHT\"", () => setInnerProperty(page, "xColumnName", el("trellis plot viewer"), "HEIGHT"));
      await session.step(120, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(121, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
