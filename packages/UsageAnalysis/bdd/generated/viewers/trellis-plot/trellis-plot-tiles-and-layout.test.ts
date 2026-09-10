/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-tiles-and-layout.feature
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
import {shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, hasNoArea, hoverArea, noErrors, pointerAway, propertyShouldBe, readingIs, readingReads, resizeTo, restoreSize, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {cellsWideTall} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot tiled view and auto layout", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-tiles-and-layout.feature", import.meta.url);
  test("Trellis plot tiled view and auto layout", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","DIS_POP"],["Y Column Names",""],["Viewer Type","Scatter plot"],["Pack Categories","false"]]));
    await session.step(22, "Then the \"x categories\" reading of trellis plot viewer should be 6", () => readingIs(page, "x categories", el("trellis plot viewer"), 6));
    await session.step(23, "And \"Tiles\" property of trellis plot viewer should be \"true\"", () => propertyShouldBe(page, "Tiles", el("trellis plot viewer"), "true"));
    await run.scenario("Tiled view rebuilds the band into a padded rectangle", async () => {
      await session.step(26, "Then \"Tiles Per Row\" property of trellis plot viewer should be \"4\"", () => propertyShouldBe(page, "Tiles Per Row", el("trellis plot viewer"), "4"));
      await session.step(27, "And the cells of trellis plot viewer should be 4 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 4, 2));
      await session.step(28, "And the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(29, "And the \"cells drawn\" reading of trellis plot viewer should be 6", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 6));
      await session.step(30, "And the \"blank cells\" reading of trellis plot viewer should be 1", () => readingIs(page, "blank cells", el("trellis plot viewer"), 1));
      await session.step(31, "And trellis plot viewer should have a \"cell RA\" area", () => hasArea(page, el("trellis plot viewer"), "cell RA"));
      await session.step(32, "And trellis plot viewer should have a \"cell UC\" area", () => hasArea(page, el("trellis plot viewer"), "cell UC"));
      await session.step(33, "And the \"distinct cell signatures\" reading of trellis plot viewer should be 6", () => readingIs(page, "distinct cell signatures", el("trellis plot viewer"), 6));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Tiles Per Row 1 gives a single-column strip", async () => {
      await session.step(37, "When user sets \"Tiles Per Row\" property of trellis plot viewer to \"1\"", () => setProperty(page, "Tiles Per Row", el("trellis plot viewer"), "1"));
      await session.step(38, "Then the cells of trellis plot viewer should be 1 wide and 5 tall", () => cellsWideTall(page, el("trellis plot viewer"), 1, 5));
      await session.step(39, "And the \"cells\" reading of trellis plot viewer should be 5", () => readingIs(page, "cells", el("trellis plot viewer"), 5));
      await session.step(40, "And the \"cells drawn\" reading of trellis plot viewer should be 5", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 5));
      await session.step(41, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
      await session.step(42, "And trellis plot viewer should have a \"cell RA\" area", () => hasArea(page, el("trellis plot viewer"), "cell RA"));
      await session.step(43, "And trellis plot viewer should not have a \"cell UC\" area", () => hasNoArea(page, el("trellis plot viewer"), "cell UC"));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Tiles Per Row 3 gives three tiles a row", async () => {
      await session.step(47, "When user sets \"Tiles Per Row\" property of trellis plot viewer to \"3\"", () => setProperty(page, "Tiles Per Row", el("trellis plot viewer"), "3"));
      await session.step(48, "Then the cells of trellis plot viewer should be 3 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 3, 2));
      await session.step(49, "And the \"cells\" reading of trellis plot viewer should be 6", () => readingIs(page, "cells", el("trellis plot viewer"), 6));
      await session.step(50, "And the \"cells drawn\" reading of trellis plot viewer should be 6", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 6));
      await session.step(51, "And trellis plot viewer should have a \"cell UC\" area", () => hasArea(page, el("trellis plot viewer"), "cell UC"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Tiles off returns a single band", async () => {
      await session.step(55, "When user sets \"Tiles\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Tiles", el("trellis plot viewer"), "false"));
      await session.step(56, "Then the cells of trellis plot viewer should be 5 wide and 1 tall", () => cellsWideTall(page, el("trellis plot viewer"), 5, 1));
      await session.step(57, "And the \"cells\" reading of trellis plot viewer should be 5", () => readingIs(page, "cells", el("trellis plot viewer"), 5));
      await session.step(58, "And the \"cells drawn\" reading of trellis plot viewer should be 5", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 5));
      await session.step(59, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
      await session.step(60, "And trellis plot viewer should not have a \"cell UC\" area", () => hasNoArea(page, el("trellis plot viewer"), "cell UC"));
      await session.step(61, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Tiles","true"],["Tiles Per Row","4"]]));
      await session.step(64, "Then the cells of trellis plot viewer should be 4 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 4, 2));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Auto layout drops the control panel when the viewer gets short", async () => {
      await session.step(68, "When user resizes trellis plot viewer to 900 by 700", () => resizeTo(page, el("trellis plot viewer"), 900, 700));
      await session.step(69, "Then trellis plot viewer should have a \"control panel\" area", () => hasArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(70, "When user resizes trellis plot viewer to 900 by 280", () => resizeTo(page, el("trellis plot viewer"), 900, 280));
      await session.step(71, "Then trellis plot viewer should not have a \"control panel\" area", () => hasNoArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(72, "And \"Show Control Panel\" property of trellis plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Control Panel", el("trellis plot viewer"), "true"));
      await session.step(73, "When user restores the size of trellis plot viewer", () => restoreSize(page, el("trellis plot viewer")));
      await session.step(74, "Then trellis plot viewer should have a \"control panel\" area", () => hasArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(75, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With auto layout off the control panel stays at any size", async () => {
      await session.step(78, "When user sets \"Auto Layout\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Auto Layout", el("trellis plot viewer"), "false"));
      await session.step(79, "And user resizes trellis plot viewer to 900 by 280", () => resizeTo(page, el("trellis plot viewer"), 900, 280));
      await session.step(80, "Then trellis plot viewer should have a \"control panel\" area", () => hasArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(81, "When user sets \"Show Control Panel\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Show Control Panel", el("trellis plot viewer"), "false"));
      await session.step(82, "Then trellis plot viewer should not have a \"control panel\" area", () => hasNoArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(83, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Show Control Panel","true"],["Auto Layout","true"]]));
      await session.step(86, "And user restores the size of trellis plot viewer", () => restoreSize(page, el("trellis plot viewer")));
      await session.step(87, "Then trellis plot viewer should have a \"control panel\" area", () => hasArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The two label strips and the two selector strips drop on their own thresholds", async () => {
      await session.step(91, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","SEX"],["Y Column Names","RACE"]]));
      await session.step(94, "Then the \"x labels shown\" reading of trellis plot viewer should be 2", () => readingIs(page, "x labels shown", el("trellis plot viewer"), 2));
      await session.step(95, "And the \"y labels shown\" reading of trellis plot viewer should be 4", () => readingIs(page, "y labels shown", el("trellis plot viewer"), 4));
      await session.step(96, "And trellis plot viewer should have an \"x selectors\" area", () => hasArea(page, el("trellis plot viewer"), "x selectors"));
      await session.step(97, "And trellis plot viewer should have a \"y selectors\" area", () => hasArea(page, el("trellis plot viewer"), "y selectors"));
      await session.step(98, "When user resizes trellis plot viewer to 240 by 400", () => resizeTo(page, el("trellis plot viewer"), 240, 400));
      await session.step(99, "Then trellis plot viewer should not have an \"x selectors\" area", () => hasNoArea(page, el("trellis plot viewer"), "x selectors"));
      await session.step(100, "And trellis plot viewer should have a \"y selectors\" area", () => hasArea(page, el("trellis plot viewer"), "y selectors"));
      await session.step(101, "When user resizes trellis plot viewer to 900 by 200", () => resizeTo(page, el("trellis plot viewer"), 900, 200));
      await session.step(102, "Then the \"x labels shown\" reading of trellis plot viewer should be 2", () => readingIs(page, "x labels shown", el("trellis plot viewer"), 2));
      await session.step(103, "And the \"y labels shown\" reading of trellis plot viewer should be 0", () => readingIs(page, "y labels shown", el("trellis plot viewer"), 0));
      await session.step(104, "And trellis plot viewer should have an \"x selectors\" area", () => hasArea(page, el("trellis plot viewer"), "x selectors"));
      await session.step(105, "And trellis plot viewer should not have a \"y selectors\" area", () => hasNoArea(page, el("trellis plot viewer"), "y selectors"));
      await session.step(106, "When user resizes trellis plot viewer to 400 by 180", () => resizeTo(page, el("trellis plot viewer"), 400, 180));
      await session.step(107, "Then the \"x labels shown\" reading of trellis plot viewer should be 0", () => readingIs(page, "x labels shown", el("trellis plot viewer"), 0));
      await session.step(108, "And the \"y labels shown\" reading of trellis plot viewer should be 0", () => readingIs(page, "y labels shown", el("trellis plot viewer"), 0));
      await session.step(109, "And the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(110, "When user restores the size of trellis plot viewer", () => restoreSize(page, el("trellis plot viewer")));
      await session.step(111, "Then the \"x labels shown\" reading of trellis plot viewer should be 2", () => readingIs(page, "x labels shown", el("trellis plot viewer"), 2));
      await session.step(112, "And the \"y labels shown\" reading of trellis plot viewer should be 4", () => readingIs(page, "y labels shown", el("trellis plot viewer"), 4));
      await session.step(113, "And no errors should have been logged", () => noErrors(page));
      await session.step(114, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","DIS_POP"],["Y Column Names",""]]));
      await session.step(117, "Then the cells of trellis plot viewer should be 4 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 4, 2));
    });
    await run.scenario("Title and description", async () => {
      await session.step(120, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Show Title","true"],["Title","Demographics"]]));
      await session.step(123, "Then title of trellis plot viewer should have text \"Demographics\"", () => shouldHaveText(page, el("title of trellis plot viewer"), "Demographics"));
      await session.step(124, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Description","By race and sex"],["Description Visibility Mode","Always"],["Description Position","Top"]]));
      await session.step(128, "Then description of trellis plot viewer should have text \"By race and sex\"", () => shouldHaveText(page, el("description of trellis plot viewer"), "By race and sex"));
      await session.step(129, "And the \"description slot\" reading of trellis plot viewer should be \"top\"", () => readingReads(page, "description slot", el("trellis plot viewer"), "top"));
      await session.step(130, "When user sets \"Description Position\" property of trellis plot viewer to \"Bottom\"", () => setProperty(page, "Description Position", el("trellis plot viewer"), "Bottom"));
      await session.step(131, "Then the \"description slot\" reading of trellis plot viewer should be \"bottom\"", () => readingReads(page, "description slot", el("trellis plot viewer"), "bottom"));
      await session.step(132, "When user sets \"Description Position\" property of trellis plot viewer to \"Left\"", () => setProperty(page, "Description Position", el("trellis plot viewer"), "Left"));
      await session.step(133, "Then the \"description slot\" reading of trellis plot viewer should be \"left\"", () => readingReads(page, "description slot", el("trellis plot viewer"), "left"));
      await session.step(134, "When user sets \"Description Position\" property of trellis plot viewer to \"Right\"", () => setProperty(page, "Description Position", el("trellis plot viewer"), "Right"));
      await session.step(135, "Then the \"description slot\" reading of trellis plot viewer should be \"right\"", () => readingReads(page, "description slot", el("trellis plot viewer"), "right"));
      await session.step(136, "When user sets \"Description Visibility Mode\" property of trellis plot viewer to \"Never\"", () => setProperty(page, "Description Visibility Mode", el("trellis plot viewer"), "Never"));
      await session.step(137, "Then the \"description slot\" reading of trellis plot viewer should be \"\"", () => readingReads(page, "description slot", el("trellis plot viewer"), ""));
      await session.step(138, "And description of trellis plot viewer should be absent", () => shouldBe(page, el("description of trellis plot viewer"), "absent"));
      await session.step(139, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Title",""],["Description",""],["Description Visibility Mode","Auto"],["Description Position","Top"]]));
      await session.step(144, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The full screen icon lives in the cell the pointer is in", async () => {
      await session.step(147, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","SEX"],["Y Column Names","RACE"]]));
      await session.step(150, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(151, "When user moves the pointer away from trellis plot viewer", () => pointerAway(page, el("trellis plot viewer")));
      await session.step(152, "Then trellis plot viewer should not have a \"full screen icon\" area", () => hasNoArea(page, el("trellis plot viewer"), "full screen icon"));
      await session.step(153, "When user hovers over the \"cell body F | Caucasian\" area of trellis plot viewer", () => hoverArea(page, "cell body F | Caucasian", el("trellis plot viewer")));
      await session.step(154, "Then trellis plot viewer should have a \"full screen icon\" area", () => hasArea(page, el("trellis plot viewer"), "full screen icon"));
      await session.step(155, "When user moves the pointer away from trellis plot viewer", () => pointerAway(page, el("trellis plot viewer")));
      await session.step(156, "Then trellis plot viewer should not have a \"full screen icon\" area", () => hasNoArea(page, el("trellis plot viewer"), "full screen icon"));
      await session.step(157, "When user hovers over the \"cell body M | Asian\" area of trellis plot viewer", () => hoverArea(page, "cell body M | Asian", el("trellis plot viewer")));
      await session.step(158, "Then trellis plot viewer should have a \"full screen icon\" area", () => hasArea(page, el("trellis plot viewer"), "full screen icon"));
      await session.step(159, "When user sets \"Allow Viewer Full Screen\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Allow Viewer Full Screen", el("trellis plot viewer"), "false"));
      await session.step(160, "And user hovers over the \"cell body F | Caucasian\" area of trellis plot viewer", () => hoverArea(page, "cell body F | Caucasian", el("trellis plot viewer")));
      await session.step(161, "Then trellis plot viewer should not have a \"full screen icon\" area", () => hasNoArea(page, el("trellis plot viewer"), "full screen icon"));
      await session.step(162, "When user sets \"Allow Viewer Full Screen\" property of trellis plot viewer to \"true\"", () => setProperty(page, "Allow Viewer Full Screen", el("trellis plot viewer"), "true"));
      await session.step(163, "And user moves the pointer away from trellis plot viewer", () => pointerAway(page, el("trellis plot viewer")));
      await session.step(164, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
