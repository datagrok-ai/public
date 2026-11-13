/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-layout.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart, public#2350]
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
import {clickOn, rightClickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {categoricalColorIs, clearSelection, colorOff, noColorCoding, selectWhereIs, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, dragSelectionOverArea, legendItemColor, lessHighlight, loadLayout, moreHighlight, noErrors, painted, pickColorSwatch, propertiesShouldBe, propertyShouldBe, readingIs, readingReads, reportsNoError, saveLayout, saveLayoutToServer, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart through a saved layout and a saved project", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-layout.feature", import.meta.url);
  test("Line chart through a saved layout and a saved project", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart", "@realizes:public#2350"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(23, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","HEIGHT"]]), [["xColumnName","AGE"],["yColumnNames","HEIGHT"]]);
    await session.step(26, "Then line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await run.scenario("The configured chart comes back from a layout applied after the viewer was closed", async () => {
      await session.step(29, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xColumnName","STARTED"],["yColumnNames","AGE, HEIGHT"],["splitColumnNames","SEX"],["multiAxis","true"],["lineWidth","3"],["interpolation","Spline"]]), [["xColumnName","STARTED"],["yColumnNames","AGE, HEIGHT"],["splitColumnNames","SEX"],["multiAxis","true"],["lineWidth","3"],["interpolation","Spline"]]);
      await session.step(36, "Then the \"x column\" reading of line chart viewer should be \"STARTED\"", () => readingReads(page, "x column", el("line chart viewer"), "STARTED"));
      await session.step(37, "And the \"y columns\" reading of line chart viewer should be \"AGE, HEIGHT\"", () => readingReads(page, "y columns", el("line chart viewer"), "AGE, HEIGHT"));
      await session.step(38, "And the \"split columns\" reading of line chart viewer should be 1", () => readingIs(page, "split columns", el("line chart viewer"), 1));
      await session.step(39, "And the \"multi axis\" reading of line chart viewer should be \"true\"", () => readingReads(page, "multi axis", el("line chart viewer"), "true"));
      await session.step(40, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(41, "And user clicks on close icon of line chart viewer", () => clickOn(page, el("close icon of line chart viewer")));
      await session.step(42, "Then line chart viewer should be absent", () => shouldBe(page, el("line chart viewer"), "absent"));
      await session.step(43, "When user loads the saved layout", () => loadLayout(page));
      await session.step(44, "Then line chart viewer should be visible", () => shouldBe(page, el("line chart viewer"), "visible"));
      await session.step(45, "And properties of line chart viewer should be:", () => propertiesShouldBe(page, el("line chart viewer"), [["xColumnName","STARTED"],["yColumnNames","AGE, HEIGHT"],["splitColumnNames","SEX"],["multiAxis","true"],["lineWidth","3"],["interpolation","Spline"]]), [["xColumnName","STARTED"],["yColumnNames","AGE, HEIGHT"],["splitColumnNames","SEX"],["multiAxis","true"],["lineWidth","3"],["interpolation","Spline"]]);
      await session.step(52, "And the \"x column\" reading of line chart viewer should be \"STARTED\"", () => readingReads(page, "x column", el("line chart viewer"), "STARTED"));
      await session.step(53, "And the \"y columns\" reading of line chart viewer should be \"AGE, HEIGHT\"", () => readingReads(page, "y columns", el("line chart viewer"), "AGE, HEIGHT"));
      await session.step(54, "And the \"split columns\" reading of line chart viewer should be 1", () => readingIs(page, "split columns", el("line chart viewer"), 1));
      await session.step(55, "And the \"categories\" reading of line chart viewer should be 2", () => readingIs(page, "categories", el("line chart viewer"), 2));
      await session.step(56, "And the \"multi axis\" reading of line chart viewer should be \"true\"", () => readingReads(page, "multi axis", el("line chart viewer"), "true"));
      await session.step(57, "And the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
      await session.step(58, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(59, "And line chart viewer should be painted", () => painted(page, el("line chart viewer")));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Selection checkboxes keep their states through a layout", async () => {
      await session.step(63, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xColumnName","AGE"],["yColumnNames","HEIGHT"],["splitColumnNames",""],["multiAxis","false"]]), [["xColumnName","AGE"],["yColumnNames","HEIGHT"],["splitColumnNames",""],["multiAxis","false"]]);
      await session.step(68, "And user drags a selection box over the \"plot\" area of line chart viewer", () => dragSelectionOverArea(page, "plot", el("line chart viewer")));
      await session.step(69, "Then some rows should be selected", () => someSelected(page));
      await session.step(70, "When user sets \"showSelectedRows\" property of line chart viewer to \"false\"", () => setProperty(page, "showSelectedRows", el("line chart viewer"), "false"));
      await session.step(71, "Then line chart viewer should show less selection highlight than before", () => lessHighlight(page, el("line chart viewer")));
      await session.step(72, "When user sets \"showSelectedRows\" property of line chart viewer to \"true\"", () => setProperty(page, "showSelectedRows", el("line chart viewer"), "true"));
      await session.step(73, "Then line chart viewer should show more selection highlight than before", () => moreHighlight(page, el("line chart viewer")));
      await session.step(74, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["showSelectedRows","false"],["showCurrentRowLine","true"],["showMouseOverCategory","false"],["showMouseOverRowLine","false"]]), [["showSelectedRows","false"],["showCurrentRowLine","true"],["showMouseOverCategory","false"],["showMouseOverRowLine","false"]]);
      await session.step(79, "And user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(80, "And user clicks on close icon of line chart viewer", () => clickOn(page, el("close icon of line chart viewer")));
      await session.step(81, "Then line chart viewer should be absent", () => shouldBe(page, el("line chart viewer"), "absent"));
      await session.step(82, "When user loads the saved layout", () => loadLayout(page));
      await session.step(83, "Then line chart viewer should be visible", () => shouldBe(page, el("line chart viewer"), "visible"));
      await session.step(84, "And properties of line chart viewer should be:", () => propertiesShouldBe(page, el("line chart viewer"), [["showSelectedRows","false"],["showCurrentRowLine","true"],["showMouseOverCategory","false"],["showMouseOverRowLine","false"]]), [["showSelectedRows","false"],["showCurrentRowLine","true"],["showMouseOverCategory","false"],["showMouseOverRowLine","false"]]);
      await session.step(89, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["showSelectedRows","true"],["showCurrentRowLine","false"],["showMouseOverCategory","true"],["showMouseOverRowLine","true"]]), [["showSelectedRows","true"],["showCurrentRowLine","false"],["showMouseOverCategory","true"],["showMouseOverRowLine","true"]]);
      await session.step(94, "And user clears the row selection", () => clearSelection(page));
      await session.step(95, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Pack Categories off and Multi Axis on come back from a layout", async () => {
      await session.step(98, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["yColumnNames","AGE, HEIGHT"],["packCategories","false"],["multiAxis","true"]]), [["yColumnNames","AGE, HEIGHT"],["packCategories","false"],["multiAxis","true"]]);
      await session.step(102, "Then the \"multi axis\" reading of line chart viewer should be \"true\"", () => readingReads(page, "multi axis", el("line chart viewer"), "true"));
      await session.step(103, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(104, "And user clicks on close icon of line chart viewer", () => clickOn(page, el("close icon of line chart viewer")));
      await session.step(105, "Then line chart viewer should be absent", () => shouldBe(page, el("line chart viewer"), "absent"));
      await session.step(106, "When user loads the saved layout", () => loadLayout(page));
      await session.step(107, "Then line chart viewer should be visible", () => shouldBe(page, el("line chart viewer"), "visible"));
      await session.step(108, "And properties of line chart viewer should be:", () => propertiesShouldBe(page, el("line chart viewer"), [["packCategories","false"],["multiAxis","true"]]), [["packCategories","false"],["multiAxis","true"]]);
      await session.step(111, "And the \"multi axis\" reading of line chart viewer should be \"true\"", () => readingReads(page, "multi axis", el("line chart viewer"), "true"));
      await session.step(112, "And the \"charts\" reading of line chart viewer should be 1", () => readingIs(page, "charts", el("line chart viewer"), 1));
      await session.step(113, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["packCategories","true"],["multiAxis","false"],["yColumnNames","HEIGHT"]]), [["packCategories","true"],["multiAxis","false"],["yColumnNames","HEIGHT"]]);
      await session.step(117, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A layout with a categorical X and Row Source Selected loads with nothing selected (public#2350)", async () => {
      await session.step(121, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xColumnName","RACE"],["yColumnNames","AGE"],["rowSource","Selected"]]), [["xColumnName","RACE"],["yColumnNames","AGE"],["rowSource","Selected"]]);
      await session.step(125, "And user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
      await session.step(126, "Then the \"rows shown\" reading of line chart viewer should be 15", () => readingIs(page, "rows shown", el("line chart viewer"), 15));
      await session.step(127, "When user saves the layout of the current table view", () => saveLayout(page));
      await session.step(128, "And user clears the row selection", () => clearSelection(page));
      await session.step(129, "Then the \"rows shown\" reading of line chart viewer should be 0", () => readingIs(page, "rows shown", el("line chart viewer"), 0));
      await session.step(130, "When user loads the saved layout", () => loadLayout(page));
      await session.step(131, "Then line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
      await session.step(132, "And \"rowSource\" property of line chart viewer should be \"Selected\"", () => propertyShouldBe(page, "rowSource", el("line chart viewer"), "Selected"));
      await session.step(133, "And the \"rows shown\" reading of line chart viewer should be 0", () => readingIs(page, "rows shown", el("line chart viewer"), 0));
      await session.step(134, "When user selects rows where \"RACE\" is \"Black\"", () => selectWhereIs(page, "RACE", "Black"));
      await session.step(135, "Then the \"rows shown\" reading of line chart viewer should be 27", () => readingIs(page, "rows shown", el("line chart viewer"), 27));
      await session.step(136, "And the \"x categories\" reading of line chart viewer should be 1", () => readingIs(page, "x categories", el("line chart viewer"), 1));
      await session.step(137, "When user clears the row selection", () => clearSelection(page));
      await session.step(138, "And user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xColumnName","AGE"],["yColumnNames","HEIGHT"],["rowSource","Filtered"]]), [["xColumnName","AGE"],["yColumnNames","HEIGHT"],["rowSource","Filtered"]]);
      await session.step(142, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A colour picked in the legend comes back from a layout and from a project", async () => {
      await session.step(145, "When user sets \"splitColumnNames\" property of line chart viewer to \"SEX\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "SEX"));
      await session.step(146, "Then legend of line chart viewer should be visible", () => shouldBe(page, el("legend of line chart viewer"), "visible"));
      await session.step(147, "When user right-clicks on \"F\" legend item in legend of line chart viewer", () => rightClickOn(page, el("\"F\" legend item in legend of line chart viewer")));
      await session.step(148, "Then \"F\" dialog should be visible", () => shouldBe(page, el("\"F\" dialog"), "visible"));
      await session.step(149, "When user picks the color \"#9467BD\" in the color picker dialog", () => pickColorSwatch(page, "#9467BD"));
      await session.step(150, "And user clicks on OK button in \"F\" dialog", () => clickOn(page, el("OK button in \"F\" dialog")));
      await session.step(151, "Then the categorical color of \"F\" in \"SEX\" column should be \"#9467BD\"", () => categoricalColorIs(page, "F", "SEX", "#9467BD"));
      await session.step(152, "And the \"F\" item in the legend of line chart viewer should be colored \"#9467BD\"", () => legendItemColor(page, "F", el("line chart viewer"), "#9467BD"));
      await session.step(153, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(154, "And user clicks on close icon of line chart viewer", () => clickOn(page, el("close icon of line chart viewer")));
      await session.step(155, "Then line chart viewer should be absent", () => shouldBe(page, el("line chart viewer"), "absent"));
      await session.step(156, "When user removes the coloring of \"SEX\" column", () => colorOff(page, "SEX"));
      await session.step(157, "Then \"SEX\" column should have no color coding", () => noColorCoding(page, "SEX"));
      await session.step(158, "When user loads the saved layout", () => loadLayout(page));
      await session.step(159, "Then line chart viewer should be visible", () => shouldBe(page, el("line chart viewer"), "visible"));
      await session.step(160, "And the categorical color of \"F\" in \"SEX\" column should be \"#9467BD\"", () => categoricalColorIs(page, "F", "SEX", "#9467BD"));
      await session.step(161, "And the \"F\" item in the legend of line chart viewer should be colored \"#9467BD\"", () => legendItemColor(page, "F", el("line chart viewer"), "#9467BD"));
      await session.step(162, "When user saves the current view as project \"bdd-line-chart-layout\"", () => saveAsProject(page, "bdd-line-chart-layout"));
      await session.step(163, "And user closes all views", () => closeAllViews(page));
      await session.step(164, "And user opens the \"bdd-line-chart-layout\" project", () => openProject(page, "bdd-line-chart-layout"));
      await session.step(165, "Then line chart viewer should be visible", () => shouldBe(page, el("line chart viewer"), "visible"));
      await session.step(166, "And legend of line chart viewer should be visible", () => shouldBe(page, el("legend of line chart viewer"), "visible"));
      await session.step(167, "And the categorical color of \"F\" in \"SEX\" column should be \"#9467BD\"", () => categoricalColorIs(page, "F", "SEX", "#9467BD"));
      await session.step(168, "And the \"F\" item in the legend of line chart viewer should be colored \"#9467BD\"", () => legendItemColor(page, "F", el("line chart viewer"), "#9467BD"));
      await session.step(169, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
