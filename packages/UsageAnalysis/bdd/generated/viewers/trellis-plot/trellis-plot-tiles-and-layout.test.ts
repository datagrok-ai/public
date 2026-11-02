/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-tiles-and-layout.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot, core#1057, core#1200]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, isExpanded, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, hasNoArea, hoverArea, noErrors, pointerAway, propertyShouldBe, readingIs, readingReads, reportsNoError, resizeTo, restoreSize, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {cellsWideTall} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot tiled view and auto layout", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-tiles-and-layout.feature", import.meta.url);
  test("Trellis plot tiled view and auto layout", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot", "@edge-cases", "@realizes:core#1057", "@known-failure:195", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(24, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","DIS_POP"],["Y Column Names",""],["Viewer Type","Scatter plot"],["Pack Categories","false"]]), [["X Column Names","DIS_POP"],["Y Column Names",""],["Viewer Type","Scatter plot"],["Pack Categories","false"]]);
    await session.step(29, "Then the \"x categories\" reading of trellis plot viewer should be 6", () => readingIs(page, "x categories", el("trellis plot viewer"), 6));
    await session.step(30, "And \"Tiles\" property of trellis plot viewer should be \"true\"", () => propertyShouldBe(page, "Tiles", el("trellis plot viewer"), "true"));
    await run.scenario("Tiled view rebuilds the band into a padded rectangle", async () => {
      await session.step(33, "Then \"Tiles Per Row\" property of trellis plot viewer should be \"4\"", () => propertyShouldBe(page, "Tiles Per Row", el("trellis plot viewer"), "4"));
      await session.step(34, "And the cells of trellis plot viewer should be 4 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 4, 2));
      await session.step(35, "And the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(36, "And the \"cells drawn\" reading of trellis plot viewer should be 6", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 6));
      await session.step(37, "And the \"blank cells\" reading of trellis plot viewer should be 1", () => readingIs(page, "blank cells", el("trellis plot viewer"), 1));
      await session.step(38, "And trellis plot viewer should have a \"cell RA\" area", () => hasArea(page, el("trellis plot viewer"), "cell RA"));
      await session.step(39, "And trellis plot viewer should have a \"cell UC\" area", () => hasArea(page, el("trellis plot viewer"), "cell UC"));
      await session.step(40, "And the \"distinct cell signatures\" reading of trellis plot viewer should be 6", () => readingIs(page, "distinct cell signatures", el("trellis plot viewer"), 6));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Tiles Per Row 1 gives a single-column strip", async () => {
      await session.step(44, "When user sets \"Tiles Per Row\" property of trellis plot viewer to \"1\"", () => setProperty(page, "Tiles Per Row", el("trellis plot viewer"), "1"));
      await session.step(45, "Then the cells of trellis plot viewer should be 1 wide and 5 tall", () => cellsWideTall(page, el("trellis plot viewer"), 1, 5));
      await session.step(46, "And the \"cells\" reading of trellis plot viewer should be 5", () => readingIs(page, "cells", el("trellis plot viewer"), 5));
      await session.step(47, "And the \"cells drawn\" reading of trellis plot viewer should be 5", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 5));
      await session.step(48, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
      await session.step(49, "And trellis plot viewer should have a \"cell RA\" area", () => hasArea(page, el("trellis plot viewer"), "cell RA"));
      await session.step(50, "And trellis plot viewer should not have a \"cell UC\" area", () => hasNoArea(page, el("trellis plot viewer"), "cell UC"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Tiles Per Row 3 gives three tiles a row", async () => {
      await session.step(54, "When user sets \"Tiles Per Row\" property of trellis plot viewer to \"3\"", () => setProperty(page, "Tiles Per Row", el("trellis plot viewer"), "3"));
      await session.step(55, "Then the cells of trellis plot viewer should be 3 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 3, 2));
      await session.step(56, "And the \"cells\" reading of trellis plot viewer should be 6", () => readingIs(page, "cells", el("trellis plot viewer"), 6));
      await session.step(57, "And the \"cells drawn\" reading of trellis plot viewer should be 6", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 6));
      await session.step(58, "And trellis plot viewer should have a \"cell UC\" area", () => hasArea(page, el("trellis plot viewer"), "cell UC"));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Tiles off returns a single band", async () => {
      await session.step(62, "When user sets \"Tiles\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Tiles", el("trellis plot viewer"), "false"));
      await session.step(63, "Then the cells of trellis plot viewer should be 5 wide and 1 tall", () => cellsWideTall(page, el("trellis plot viewer"), 5, 1));
      await session.step(64, "And the \"cells\" reading of trellis plot viewer should be 5", () => readingIs(page, "cells", el("trellis plot viewer"), 5));
      await session.step(65, "And the \"cells drawn\" reading of trellis plot viewer should be 5", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 5));
      await session.step(66, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
      await session.step(67, "And trellis plot viewer should not have a \"cell UC\" area", () => hasNoArea(page, el("trellis plot viewer"), "cell UC"));
      await session.step(68, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Tiles","true"],["Tiles Per Row","4"]]), [["Tiles","true"],["Tiles Per Row","4"]]);
      await session.step(71, "Then the cells of trellis plot viewer should be 4 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 4, 2));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Auto layout drops the control panel when the viewer gets short", async () => {
      await session.step(75, "When user resizes trellis plot viewer to 900 by 700", () => resizeTo(page, el("trellis plot viewer"), 900, 700));
      await session.step(76, "Then trellis plot viewer should have a \"control panel\" area", () => hasArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(77, "When user resizes trellis plot viewer to 900 by 280", () => resizeTo(page, el("trellis plot viewer"), 900, 280));
      await session.step(78, "Then trellis plot viewer should not have a \"control panel\" area", () => hasNoArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(79, "And \"Show Control Panel\" property of trellis plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Control Panel", el("trellis plot viewer"), "true"));
      await session.step(80, "When user restores the size of trellis plot viewer", () => restoreSize(page, el("trellis plot viewer")));
      await session.step(81, "Then trellis plot viewer should have a \"control panel\" area", () => hasArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With auto layout off the control panel stays at any size", async () => {
      await session.step(85, "When user sets \"Auto Layout\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Auto Layout", el("trellis plot viewer"), "false"));
      await session.step(86, "And user resizes trellis plot viewer to 900 by 280", () => resizeTo(page, el("trellis plot viewer"), 900, 280));
      await session.step(87, "Then trellis plot viewer should have a \"control panel\" area", () => hasArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(88, "When user sets \"Show Control Panel\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Show Control Panel", el("trellis plot viewer"), "false"));
      await session.step(89, "Then trellis plot viewer should not have a \"control panel\" area", () => hasNoArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(90, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Show Control Panel","true"],["Auto Layout","true"]]), [["Show Control Panel","true"],["Auto Layout","true"]]);
      await session.step(93, "And user restores the size of trellis plot viewer", () => restoreSize(page, el("trellis plot viewer")));
      await session.step(94, "Then trellis plot viewer should have a \"control panel\" area", () => hasArea(page, el("trellis plot viewer"), "control panel"));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The two label strips and the two selector strips drop on their own thresholds", async () => {
      await session.step(98, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","SEX"],["Y Column Names","RACE"]]), [["X Column Names","SEX"],["Y Column Names","RACE"]]);
      await session.step(101, "Then the \"x labels shown\" reading of trellis plot viewer should be 2", () => readingIs(page, "x labels shown", el("trellis plot viewer"), 2));
      await session.step(102, "And the \"y labels shown\" reading of trellis plot viewer should be 4", () => readingIs(page, "y labels shown", el("trellis plot viewer"), 4));
      await session.step(103, "And trellis plot viewer should have an \"x selectors\" area", () => hasArea(page, el("trellis plot viewer"), "x selectors"));
      await session.step(104, "And trellis plot viewer should have a \"y selectors\" area", () => hasArea(page, el("trellis plot viewer"), "y selectors"));
      await session.step(105, "When user resizes trellis plot viewer to 240 by 400", () => resizeTo(page, el("trellis plot viewer"), 240, 400));
      await session.step(106, "Then trellis plot viewer should not have an \"x selectors\" area", () => hasNoArea(page, el("trellis plot viewer"), "x selectors"));
      await session.step(107, "And trellis plot viewer should have a \"y selectors\" area", () => hasArea(page, el("trellis plot viewer"), "y selectors"));
      await session.step(108, "When user resizes trellis plot viewer to 900 by 200", () => resizeTo(page, el("trellis plot viewer"), 900, 200));
      await session.step(109, "Then the \"x labels shown\" reading of trellis plot viewer should be 2", () => readingIs(page, "x labels shown", el("trellis plot viewer"), 2));
      await session.step(110, "And the \"y labels shown\" reading of trellis plot viewer should be 0", () => readingIs(page, "y labels shown", el("trellis plot viewer"), 0));
      await session.step(111, "And trellis plot viewer should have an \"x selectors\" area", () => hasArea(page, el("trellis plot viewer"), "x selectors"));
      await session.step(112, "And trellis plot viewer should not have a \"y selectors\" area", () => hasNoArea(page, el("trellis plot viewer"), "y selectors"));
      await session.step(113, "When user resizes trellis plot viewer to 400 by 180", () => resizeTo(page, el("trellis plot viewer"), 400, 180));
      await session.step(114, "Then the \"x labels shown\" reading of trellis plot viewer should be 0", () => readingIs(page, "x labels shown", el("trellis plot viewer"), 0));
      await session.step(115, "And the \"y labels shown\" reading of trellis plot viewer should be 0", () => readingIs(page, "y labels shown", el("trellis plot viewer"), 0));
      await session.step(116, "And the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(117, "When user restores the size of trellis plot viewer", () => restoreSize(page, el("trellis plot viewer")));
      await session.step(118, "Then the \"x labels shown\" reading of trellis plot viewer should be 2", () => readingIs(page, "x labels shown", el("trellis plot viewer"), 2));
      await session.step(119, "And the \"y labels shown\" reading of trellis plot viewer should be 4", () => readingIs(page, "y labels shown", el("trellis plot viewer"), 4));
      await session.step(120, "And no errors should have been logged", () => noErrors(page));
      await session.step(121, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","DIS_POP"],["Y Column Names",""]]), [["X Column Names","DIS_POP"],["Y Column Names",""]]);
      await session.step(124, "Then the cells of trellis plot viewer should be 4 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 4, 2));
    });
    await run.scenario("Title and description", async () => {
      await session.step(127, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Show Title","true"],["Title","Demographics"]]), [["Show Title","true"],["Title","Demographics"]]);
      await session.step(130, "Then title of trellis plot viewer should have text \"Demographics\"", () => shouldHaveText(page, el("title of trellis plot viewer"), "Demographics"));
      await session.step(131, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Description","By race and sex"],["Description Visibility Mode","Always"],["Description Position","Top"]]), [["Description","By race and sex"],["Description Visibility Mode","Always"],["Description Position","Top"]]);
      await session.step(135, "Then description of trellis plot viewer should have text \"By race and sex\"", () => shouldHaveText(page, el("description of trellis plot viewer"), "By race and sex"));
      await session.step(136, "And the \"description slot\" reading of trellis plot viewer should be \"top\"", () => readingReads(page, "description slot", el("trellis plot viewer"), "top"));
      await session.step(137, "When user sets \"Description Position\" property of trellis plot viewer to \"Bottom\"", () => setProperty(page, "Description Position", el("trellis plot viewer"), "Bottom"));
      await session.step(138, "Then the \"description slot\" reading of trellis plot viewer should be \"bottom\"", () => readingReads(page, "description slot", el("trellis plot viewer"), "bottom"));
      await session.step(139, "When user sets \"Description Position\" property of trellis plot viewer to \"Left\"", () => setProperty(page, "Description Position", el("trellis plot viewer"), "Left"));
      await session.step(140, "Then the \"description slot\" reading of trellis plot viewer should be \"left\"", () => readingReads(page, "description slot", el("trellis plot viewer"), "left"));
      await session.step(141, "When user sets \"Description Position\" property of trellis plot viewer to \"Right\"", () => setProperty(page, "Description Position", el("trellis plot viewer"), "Right"));
      await session.step(142, "Then the \"description slot\" reading of trellis plot viewer should be \"right\"", () => readingReads(page, "description slot", el("trellis plot viewer"), "right"));
      await session.step(143, "When user sets \"Description Visibility Mode\" property of trellis plot viewer to \"Never\"", () => setProperty(page, "Description Visibility Mode", el("trellis plot viewer"), "Never"));
      await session.step(144, "Then the \"description slot\" reading of trellis plot viewer should be \"\"", () => readingReads(page, "description slot", el("trellis plot viewer"), ""));
      await session.step(145, "And description of trellis plot viewer should be absent", () => shouldBe(page, el("description of trellis plot viewer"), "absent"));
      await session.step(146, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Title",""],["Description",""],["Description Visibility Mode","Auto"],["Description Position","Top"]]), [["Title",""],["Description",""],["Description Visibility Mode","Auto"],["Description Position","Top"]]);
      await session.step(151, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The full screen icon lives in the cell the pointer is in", async () => {
      await session.step(154, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","SEX"],["Y Column Names","RACE"]]), [["X Column Names","SEX"],["Y Column Names","RACE"]]);
      await session.step(157, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
      await session.step(158, "When user moves the pointer away from trellis plot viewer", () => pointerAway(page, el("trellis plot viewer")));
      await session.step(159, "Then trellis plot viewer should not have a \"full screen icon\" area", () => hasNoArea(page, el("trellis plot viewer"), "full screen icon"));
      await session.step(160, "When user hovers over the \"cell body F | Caucasian\" area of trellis plot viewer", () => hoverArea(page, "cell body F | Caucasian", el("trellis plot viewer")));
      await session.step(161, "Then trellis plot viewer should have a \"full screen icon\" area", () => hasArea(page, el("trellis plot viewer"), "full screen icon"));
      await session.step(162, "When user moves the pointer away from trellis plot viewer", () => pointerAway(page, el("trellis plot viewer")));
      await session.step(163, "Then trellis plot viewer should not have a \"full screen icon\" area", () => hasNoArea(page, el("trellis plot viewer"), "full screen icon"));
      await session.step(164, "When user hovers over the \"cell body M | Asian\" area of trellis plot viewer", () => hoverArea(page, "cell body M | Asian", el("trellis plot viewer")));
      await session.step(165, "Then trellis plot viewer should have a \"full screen icon\" area", () => hasArea(page, el("trellis plot viewer"), "full screen icon"));
      await session.step(166, "When user sets \"Allow Viewer Full Screen\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Allow Viewer Full Screen", el("trellis plot viewer"), "false"));
      await session.step(167, "And user hovers over the \"cell body F | Caucasian\" area of trellis plot viewer", () => hoverArea(page, "cell body F | Caucasian", el("trellis plot viewer")));
      await session.step(168, "Then trellis plot viewer should not have a \"full screen icon\" area", () => hasNoArea(page, el("trellis plot viewer"), "full screen icon"));
      await session.step(169, "When user sets \"Allow Viewer Full Screen\" property of trellis plot viewer to \"true\"", () => setProperty(page, "Allow Viewer Full Screen", el("trellis plot viewer"), "true"));
      await session.step(170, "And user moves the pointer away from trellis plot viewer", () => pointerAway(page, el("trellis plot viewer")));
      await session.step(171, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Tiles Per Row 0 on a one-column trellis lays the tiles out one a row and logs no error (core#1057)", async () => {
      await session.step(175, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","DIS_POP"],["Y Column Names",""],["Tiles Per Row","4"]]), [["X Column Names","DIS_POP"],["Y Column Names",""],["Tiles Per Row","4"]]);
      await session.step(179, "Then the \"y categories\" reading of trellis plot viewer should be 1", () => readingIs(page, "y categories", el("trellis plot viewer"), 1));
      await session.step(180, "And the cells of trellis plot viewer should be 4 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 4, 2));
      await session.step(181, "When user sets \"Tiles Per Row\" property of trellis plot viewer to \"0\"", () => setProperty(page, "Tiles Per Row", el("trellis plot viewer"), "0"));
      await session.step(182, "Then the cells of trellis plot viewer should be 1 wide and 5 tall", () => cellsWideTall(page, el("trellis plot viewer"), 1, 5));
      await session.step(183, "And trellis plot viewer should report no error", () => reportsNoError(page, el("trellis plot viewer")));
      await session.step(184, "And no errors should have been logged", () => noErrors(page));
      await session.step(185, "When user sets \"Tiles Per Row\" property of trellis plot viewer to \"4\"", () => setProperty(page, "Tiles Per Row", el("trellis plot viewer"), "4"));
      await session.step(186, "Then the cells of trellis plot viewer should be 4 wide and 2 tall", () => cellsWideTall(page, el("trellis plot viewer"), 4, 2));
    });
    await run.scenario("A Tiles Per Row of 11 typed in the context panel outside its declared bounds is not taken (core#1200)", async () => {
      await session.step(190, "Given user sets \"Tiles Per Row\" property of trellis plot viewer to \"4\"", () => setProperty(page, "Tiles Per Row", el("trellis plot viewer"), "4"));
      await session.step(191, "When user clicks on settings icon of trellis plot viewer", () => clickOn(page, el("settings icon of trellis plot viewer")));
      await session.step(192, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(193, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
      await session.step(194, "When user enters \"11\" into \"Tiles Per Row\" property in context panel", () => enterInto(page, "11", el("\"Tiles Per Row\" property in context panel")));
      await session.step(195, "Then \"Tiles Per Row\" property of trellis plot viewer should be \"4\"", () => propertyShouldBe(page, "Tiles Per Row", el("trellis plot viewer"), "4"));
    }, {knownFailure: true, line: 195});
    run.finish();
  });
});
