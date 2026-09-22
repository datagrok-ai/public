/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-persistence.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid]
--- */
import {test} from '@playwright/test';
import '../../../bindings/grid.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag, makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {colorCodedCategorically, colorConditional, colorLinkedTo} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaAtLeastTall, clickArea, doubleClickArea, dragAreaBy, dragAreaToArea, hasNoArea, loadLayout, noErrors, pickFromAreaContextMenu, propertyShouldBe, readingAsRemembered, readingDiffers, readingHigher, readingIs, readingNotAsRemembered, readingReads, readingsDiffer, readingsEqual, rememberReading, saveLayoutToServer, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingNotContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid appearance and geometry across a layout and a project", () => {
  const session = feature(test, "features/viewers/grid/grid-persistence.feature", import.meta.url);
  test("Grid appearance and geometry across a layout and a project", {tag: ["@journey", "@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(47, "Given user is logged in", () => loggedIn(page));
    await session.step(48, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(49, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(50, "And \"Frozen Columns\" property of grid should be \"1\"", () => propertyShouldBe(page, "Frozen Columns", el("grid"), "1"));
    await run.scenario("Four colour codings, the row height and the missing-value colour are set", async () => {
      await session.step(53, "When user picks \"Color Coding > Linear\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linear", "header AGE", el("grid")));
      await session.step(54, "And user picks \"Color Coding > Conditional\" from the context menu of the \"header HEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Conditional", "header HEIGHT", el("grid")));
      await session.step(55, "And user colors \"HEIGHT\" column conditionally:", () => colorConditional(page, "HEIGHT", [["<160","#0000FF"],[">180","#FF0000"]]), [["<160","#0000FF"],[">180","#FF0000"]]);
      await session.step(58, "And user picks \"Color Coding > Categorical\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header SEX", el("grid")));
      await session.step(59, "Then \"AGE\" column should have tag \".color-coding-type\" equal to \"Linear\"", () => columnTag(page, "AGE", ".color-coding-type", "Linear"));
      await session.step(60, "And \"HEIGHT\" column should have tag \".color-coding-type\" equal to \"Conditional\"", () => columnTag(page, "HEIGHT", ".color-coding-type", "Conditional"));
      await session.step(61, "And \"SEX\" column should be color-coded categorically", () => colorCodedCategorically(page, "SEX"));
      await session.step(62, "When user picks \"Color Coding > Linked\" from the context menu of the \"header WEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Linked", "header WEIGHT", el("grid")));
      await session.step(63, "And user picks \"Color Coding > Edit...\" from the context menu of the \"header WEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Edit...", "header WEIGHT", el("grid")));
      await session.step(64, "Then \"Color-coding: WEIGHT\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: WEIGHT\" dialog"), "visible"));
      await session.step(65, "When user selects \"SEX\" in \"Source column\" input in \"Color-coding: WEIGHT\" dialog", () => selectIn(page, "SEX", el("\"Source column\" input in \"Color-coding: WEIGHT\" dialog")));
      await session.step(66, "And user clicks on CLOSE button in \"Color-coding: WEIGHT\" dialog", () => clickOn(page, el("CLOSE button in \"Color-coding: WEIGHT\" dialog")));
      await session.step(67, "Then the coloring of \"WEIGHT\" column should be linked to \"SEX\" column", () => colorLinkedTo(page, "WEIGHT", "SEX"));
      await session.step(68, "And the \"color of cell 1 of WEIGHT\" and \"color of cell 1 of SEX\" readings of grid should be the same", () => readingsEqual(page, "color of cell 1 of WEIGHT", "color of cell 1 of SEX", el("grid")));
      await session.step(69, "And the \"color of cell 4 of WEIGHT\" and \"color of cell 4 of SEX\" readings of grid should be the same", () => readingsEqual(page, "color of cell 4 of WEIGHT", "color of cell 4 of SEX", el("grid")));
      await session.step(70, "And the \"color of cell 1 of WEIGHT\" and \"color of cell 4 of WEIGHT\" readings of grid should differ", () => readingsDiffer(page, "color of cell 1 of WEIGHT", "color of cell 4 of WEIGHT", el("grid")));
      await session.step(71, "When user sets \"Row Height\" property of grid to \"40\"", () => setProperty(page, "Row Height", el("grid"), "40"));
      await session.step(72, "And user sets \"Missing Value Color\" property of grid to \"#FFAAAA\"", () => setProperty(page, "Missing Value Color", el("grid"), "#FFAAAA"));
      await session.step(73, "Then the \"row height\" reading of grid should be 40", () => readingIs(page, "row height", el("grid"), 40));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The min and max stats rows leave the columns as they were", async () => {
      await session.step(77, "When user remembers the \"column order\" reading of grid", () => rememberReading(page, "column order", el("grid")));
      await session.step(78, "And user picks \"Add > Column Stats > min\" from the context menu of the \"cell 2 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Add > Column Stats > min", "cell 2 of AGE", el("grid")));
      await session.step(79, "And user picks \"Add > Column Stats > max\" from the context menu of the \"cell 2 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Add > Column Stats > max", "cell 2 of AGE", el("grid")));
      await session.step(80, "Then the \"column order\" reading of grid should be as remembered", () => readingAsRemembered(page, "column order", el("grid")));
      await session.step(81, "And grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Moving, hiding, widening, pinning and sorting set the geometry", async () => {
      await session.step(85, "When user drags the \"header HEIGHT\" area of grid to the \"header DEMOG\" area", () => dragAreaToArea(page, "header HEIGHT", el("grid"), "header DEMOG"));
      await session.step(86, "Then the \"column order\" reading of grid should differ from before", () => readingDiffers(page, "column order", el("grid")));
      await session.step(87, "When user picks \"Order or Hide Columns...\" from the context menu of the \"cell 2 of AGE\" area of grid", () => pickFromAreaContextMenu(page, "Order or Hide Columns...", "cell 2 of AGE", el("grid")));
      await session.step(88, "Then the \"text of cell 4 of __name\" reading of Grid viewer in Order or Hide Columns dialog should be \"RACE\"", () => readingReads(page, "text of cell 4 of __name", el("Grid viewer in Order or Hide Columns dialog"), "RACE"));
      await session.step(89, "And the \"text of cell 4 of x\" reading of Grid viewer in Order or Hide Columns dialog should be \"true\"", () => readingReads(page, "text of cell 4 of x", el("Grid viewer in Order or Hide Columns dialog"), "true"));
      await session.step(90, "When user clicks on the \"cell 4 of x\" area of Grid viewer in Order or Hide Columns dialog", () => clickArea(page, "cell 4 of x", el("Grid viewer in Order or Hide Columns dialog")));
      await session.step(91, "Then the \"text of cell 4 of x\" reading of Grid viewer in Order or Hide Columns dialog should be \"false\"", () => readingReads(page, "text of cell 4 of x", el("Grid viewer in Order or Hide Columns dialog"), "false"));
      await session.step(92, "When user clicks on CLOSE button in Order or Hide Columns dialog", () => clickOn(page, el("CLOSE button in Order or Hide Columns dialog")));
      await session.step(93, "Then grid should not have a \"header RACE\" area", () => hasNoArea(page, el("grid"), "header RACE"));
      await session.step(94, "And the \"column order\" reading of grid should not contain \"RACE\"", () => readingNotContains(page, "column order", el("grid"), "RACE"));
      await session.step(95, "When user drags the \"column resizer AGE\" area of grid by 60 pixels to the right", () => dragAreaBy(page, "column resizer AGE", el("grid"), 60, "right"));
      await session.step(96, "Then the \"column width of AGE\" reading of grid should be higher than before", () => readingHigher(page, "column width of AGE", el("grid")));
      await session.step(97, "When user picks \"Pin > Pin Column\" from the context menu of the \"header SEX\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Column", "header SEX", el("grid")));
      await session.step(98, "Then \"Frozen Columns\" property of grid should be \"2\"", () => propertyShouldBe(page, "Frozen Columns", el("grid"), "2"));
      await session.step(99, "When user picks \"Pin > Pin Row\" from the context menu of the \"cell 2 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Row", "cell 2 of USUBJID", el("grid")));
      await session.step(100, "And user picks \"Pin > Pin Row\" from the context menu of the \"cell 4 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Row", "cell 4 of USUBJID", el("grid")));
      await session.step(101, "And user makes row 298 current", () => makeRowCurrent(page, 298));
      await session.step(102, "And user makes row 300 current", () => makeRowCurrent(page, 300));
      await session.step(103, "And user picks \"Pin > Pin Row\" from the context menu of the \"cell 298 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Row", "cell 298 of USUBJID", el("grid")));
      await session.step(104, "Then the \"pinned rows\" reading of grid should be 3", () => readingIs(page, "pinned rows", el("grid"), 3));
      await session.step(105, "When user double-clicks on the \"header AGE\" area of grid", () => doubleClickArea(page, "header AGE", el("grid")));
      await session.step(106, "And user double-clicks on the \"header AGE\" area of grid", () => doubleClickArea(page, "header AGE", el("grid")));
      await session.step(107, "Then the \"sort column\" reading of grid should be \"AGE\"", () => readingReads(page, "sort column", el("grid"), "AGE"));
      await session.step(108, "And the \"sort direction\" reading of grid should be \"ascending\"", () => readingReads(page, "sort direction", el("grid"), "ascending"));
      await session.step(109, "And the \"color of cell 298 of HEIGHT\" reading of grid should be \"#ffaaaa\"", () => readingReads(page, "color of cell 298 of HEIGHT", el("grid"), "#ffaaaa"));
      await session.step(110, "And the \"color of cell 4 of HEIGHT\" reading of grid should be \"#ff0000\"", () => readingReads(page, "color of cell 4 of HEIGHT", el("grid"), "#ff0000"));
      await session.step(111, "And the \"color of cell 2 of HEIGHT\" reading of grid should be \"#0000ff\"", () => readingReads(page, "color of cell 2 of HEIGHT", el("grid"), "#0000ff"));
      await session.step(112, "And the \"color of cell 4 of AGE\" and \"color of cell 298 of AGE\" readings of grid should differ", () => readingsDiffer(page, "color of cell 4 of AGE", "color of cell 298 of AGE", el("grid")));
      await session.step(113, "And the \"color of cell 4 of SEX\" and \"color of cell 298 of SEX\" readings of grid should differ", () => readingsDiffer(page, "color of cell 4 of SEX", "color of cell 298 of SEX", el("grid")));
      await session.step(114, "And the \"color of cell 298 of WEIGHT\" and \"color of cell 298 of SEX\" readings of grid should be the same", () => readingsEqual(page, "color of cell 298 of WEIGHT", "color of cell 298 of SEX", el("grid")));
      await session.step(115, "When user remembers the \"column order\" reading of grid", () => rememberReading(page, "column order", el("grid")));
      await session.step(116, "And user remembers the \"column width of AGE\" reading of grid", () => rememberReading(page, "column width of AGE", el("grid")));
      await session.step(117, "And user remembers the \"row height\" reading of grid", () => rememberReading(page, "row height", el("grid")));
      await session.step(118, "And user remembers the \"color of cell 4 of AGE\" reading of grid", () => rememberReading(page, "color of cell 4 of AGE", el("grid")));
      await session.step(119, "And user remembers the \"color of cell 298 of AGE\" reading of grid", () => rememberReading(page, "color of cell 298 of AGE", el("grid")));
      await session.step(120, "And user remembers the \"color of cell 4 of SEX\" reading of grid", () => rememberReading(page, "color of cell 4 of SEX", el("grid")));
      await session.step(121, "And user remembers the \"color of cell 298 of SEX\" reading of grid", () => rememberReading(page, "color of cell 298 of SEX", el("grid")));
      await session.step(122, "And user remembers the \"color of cell 4 of WEIGHT\" reading of grid", () => rememberReading(page, "color of cell 4 of WEIGHT", el("grid")));
      await session.step(123, "And user remembers the \"color of cell 298 of WEIGHT\" reading of grid", () => rememberReading(page, "color of cell 298 of WEIGHT", el("grid")));
      await session.step(124, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A layout saved to the server restores everything over a fresh view of the table", async () => {
      await session.step(127, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(128, "And user closes all views", () => closeAllViews(page));
      await session.step(129, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(130, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
      await session.step(131, "Then the \"column order\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "column order", el("grid")));
      await session.step(132, "And the \"column width of AGE\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "column width of AGE", el("grid")));
      await session.step(133, "And the \"row height\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "row height", el("grid")));
      await session.step(134, "And the \"color of cell 4 of AGE\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "color of cell 4 of AGE", el("grid")));
      await session.step(135, "And \"Frozen Columns\" property of grid should be \"1\"", () => propertyShouldBe(page, "Frozen Columns", el("grid"), "1"));
      await session.step(136, "And the \"pinned rows\" reading of grid should be 0", () => readingIs(page, "pinned rows", el("grid"), 0));
      await session.step(137, "And the \"sort column\" reading of grid should be \"\"", () => readingReads(page, "sort column", el("grid"), ""));
      await session.step(138, "When user loads the saved layout", () => loadLayout(page));
      await session.step(139, "Then scatter plot viewer should be absent", () => shouldBe(page, el("scatter plot viewer"), "absent"));
      await session.step(140, "And the \"column order\" reading of grid should be as remembered", () => readingAsRemembered(page, "column order", el("grid")));
      await session.step(141, "And the \"column width of AGE\" reading of grid should be as remembered", () => readingAsRemembered(page, "column width of AGE", el("grid")));
      await session.step(142, "And the \"row height\" reading of grid should be as remembered", () => readingAsRemembered(page, "row height", el("grid")));
      await session.step(143, "And the \"cell 4 of AGE\" area of grid should be at least 36 pixels tall", () => areaAtLeastTall(page, "cell 4 of AGE", el("grid"), 36));
      await session.step(144, "And \"Frozen Columns\" property of grid should be \"2\"", () => propertyShouldBe(page, "Frozen Columns", el("grid"), "2"));
      await session.step(145, "And the \"pinned rows\" reading of grid should be 3", () => readingIs(page, "pinned rows", el("grid"), 3));
      await session.step(146, "And the \"sort column\" reading of grid should be \"AGE\"", () => readingReads(page, "sort column", el("grid"), "AGE"));
      await session.step(147, "And the \"sort direction\" reading of grid should be \"ascending\"", () => readingReads(page, "sort direction", el("grid"), "ascending"));
      await session.step(148, "And the \"color of cell 4 of AGE\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 4 of AGE", el("grid")));
      await session.step(149, "And the \"color of cell 298 of AGE\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 298 of AGE", el("grid")));
      await session.step(150, "And the \"color of cell 4 of HEIGHT\" reading of grid should be \"#ff0000\"", () => readingReads(page, "color of cell 4 of HEIGHT", el("grid"), "#ff0000"));
      await session.step(151, "And the \"color of cell 298 of HEIGHT\" reading of grid should be \"#ffaaaa\"", () => readingReads(page, "color of cell 298 of HEIGHT", el("grid"), "#ffaaaa"));
      await session.step(152, "And the \"color of cell 2 of HEIGHT\" reading of grid should be \"#0000ff\"", () => readingReads(page, "color of cell 2 of HEIGHT", el("grid"), "#0000ff"));
      await session.step(153, "And the \"color of cell 4 of SEX\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 4 of SEX", el("grid")));
      await session.step(154, "And the \"color of cell 298 of SEX\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 298 of SEX", el("grid")));
      await session.step(155, "And the \"color of cell 4 of WEIGHT\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 4 of WEIGHT", el("grid")));
      await session.step(156, "And the \"color of cell 298 of WEIGHT\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 298 of WEIGHT", el("grid")));
      await session.step(157, "And the coloring of \"WEIGHT\" column should be linked to \"SEX\" column", () => colorLinkedTo(page, "WEIGHT", "SEX"));
      await session.step(158, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A project round-trip restores everything", async () => {
      await session.step(161, "When user saves the current view as project \"zz-grid-persistence\"", () => saveAsProject(page, "zz-grid-persistence"));
      await session.step(162, "And user closes all views", () => closeAllViews(page));
      await session.step(163, "And user opens the \"zz-grid-persistence\" project", () => openProject(page, "zz-grid-persistence"));
      await session.step(164, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
      await session.step(165, "And the \"column order\" reading of grid should be as remembered", () => readingAsRemembered(page, "column order", el("grid")));
      await session.step(166, "And the \"column width of AGE\" reading of grid should be as remembered", () => readingAsRemembered(page, "column width of AGE", el("grid")));
      await session.step(167, "And the \"row height\" reading of grid should be as remembered", () => readingAsRemembered(page, "row height", el("grid")));
      await session.step(168, "And the \"cell 4 of AGE\" area of grid should be at least 36 pixels tall", () => areaAtLeastTall(page, "cell 4 of AGE", el("grid"), 36));
      await session.step(169, "And \"Frozen Columns\" property of grid should be \"2\"", () => propertyShouldBe(page, "Frozen Columns", el("grid"), "2"));
      await session.step(170, "And the \"pinned rows\" reading of grid should be 3", () => readingIs(page, "pinned rows", el("grid"), 3));
      await session.step(171, "And the \"sort column\" reading of grid should be \"AGE\"", () => readingReads(page, "sort column", el("grid"), "AGE"));
      await session.step(172, "And the \"sort direction\" reading of grid should be \"ascending\"", () => readingReads(page, "sort direction", el("grid"), "ascending"));
      await session.step(173, "And the \"color of cell 4 of AGE\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 4 of AGE", el("grid")));
      await session.step(174, "And the \"color of cell 298 of AGE\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 298 of AGE", el("grid")));
      await session.step(175, "And the \"color of cell 4 of HEIGHT\" reading of grid should be \"#ff0000\"", () => readingReads(page, "color of cell 4 of HEIGHT", el("grid"), "#ff0000"));
      await session.step(176, "And the \"color of cell 298 of HEIGHT\" reading of grid should be \"#ffaaaa\"", () => readingReads(page, "color of cell 298 of HEIGHT", el("grid"), "#ffaaaa"));
      await session.step(177, "And the \"color of cell 2 of HEIGHT\" reading of grid should be \"#0000ff\"", () => readingReads(page, "color of cell 2 of HEIGHT", el("grid"), "#0000ff"));
      await session.step(178, "And the \"color of cell 4 of SEX\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 4 of SEX", el("grid")));
      await session.step(179, "And the \"color of cell 298 of SEX\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 298 of SEX", el("grid")));
      await session.step(180, "And the \"color of cell 4 of WEIGHT\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 4 of WEIGHT", el("grid")));
      await session.step(181, "And the \"color of cell 298 of WEIGHT\" reading of grid should be as remembered", () => readingAsRemembered(page, "color of cell 298 of WEIGHT", el("grid")));
      await session.step(182, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
