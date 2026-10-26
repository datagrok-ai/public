/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-look-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot]
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
import {colorLinearOver, colorOff, filterPassesAll, filterPassesFewer, noColorCoding, passingAsRemembered, rememberPassing} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, lessInk, loadLayout, moreInk, noErrors, pickFromContextMenu, readingAsRemembered, readingDoesNotRead, readingHigher, readingIs, rememberReading, saveLayout, scaleNarrower, scaleWider, setProperties, setProperty, showsRows, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragRangeHandle} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot regressions — markers, colors, layouts and the project", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-look-regressions.feature", import.meta.url);
  test("Scatter plot regressions — markers, colors, layouts and the project", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(20, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await run.scenario("Marker Default Size changes the markers (public#2285)", async () => {
      await session.step(24, "When user sets \"markerDefaultSize\" property of scatter plot viewer to \"12\"", () => setProperty(page, "markerDefaultSize", el("scatter plot viewer"), "12"));
      await session.step(25, "Then scatter plot viewer should have more ink than before", () => moreInk(page, el("scatter plot viewer")));
      await session.step(26, "When user sets \"markerDefaultSize\" property of scatter plot viewer to \"2\"", () => setProperty(page, "markerDefaultSize", el("scatter plot viewer"), "2"));
      await session.step(27, "Then scatter plot viewer should have less ink than before", () => lessInk(page, el("scatter plot viewer")));
      await session.step(28, "When user sets \"markerDefaultSize\" property of scatter plot viewer to \"\"", () => setProperty(page, "markerDefaultSize", el("scatter plot viewer"), ""));
      await session.step(29, "Then scatter plot viewer should have more ink than before", () => moreInk(page, el("scatter plot viewer")));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Over 20 color categories Regression Per Category draws one common line (public#3695)", async () => {
      await session.step(35, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Show Regression Line","true"],["Regression Per Category","true"],["Color","RACE"]]), [["Show Regression Line","true"],["Regression Per Category","true"],["Color","RACE"]]);
      await session.step(39, "Then the \"regression lines\" reading of scatter plot viewer should be 4", () => readingIs(page, "regression lines", el("scatter plot viewer"), 4));
      await session.step(40, "When user sets \"Color\" property of scatter plot viewer to \"USUBJID\"", () => setProperty(page, "Color", el("scatter plot viewer"), "USUBJID"));
      await session.step(41, "Then the \"regression lines\" reading of scatter plot viewer should be 1", () => readingIs(page, "regression lines", el("scatter plot viewer"), 1));
      await session.step(42, "And scatter plot viewer should have a \"regression line\" area", () => hasArea(page, el("scatter plot viewer"), "regression line"));
      await session.step(43, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Color",""],["Show Regression Line","false"]]), [["Color",""],["Show Regression Line","false"]]);
      await session.step(46, "Then the \"regression lines\" reading of scatter plot viewer should be 0", () => readingIs(page, "regression lines", el("scatter plot viewer"), 0));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Color Min and Color Max narrow the color scale (public#2658)", async () => {
      await session.step(51, "When user sets \"Color\" property of scatter plot viewer to \"AGE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "AGE"));
      await session.step(52, "Then the \"color scale min\" reading of scatter plot viewer should be 18", () => readingIs(page, "color scale min", el("scatter plot viewer"), 18));
      await session.step(53, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["colorMin","30"],["colorMax","50"]]), [["colorMin","30"],["colorMax","50"]]);
      await session.step(56, "Then the color scale of scatter plot viewer should cover a narrower range than before", () => scaleNarrower(page, el("scatter plot viewer")));
      await session.step(57, "And the \"color scale min\" reading of scatter plot viewer should be 30", () => readingIs(page, "color scale min", el("scatter plot viewer"), 30));
      await session.step(58, "And the \"color scale max\" reading of scatter plot viewer should be 50", () => readingIs(page, "color scale max", el("scatter plot viewer"), 50));
      await session.step(59, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["colorMin",""],["colorMax",""]]), [["colorMin",""],["colorMax",""]]);
      await session.step(62, "Then the \"color scale min\" reading of scatter plot viewer should be 18", () => readingIs(page, "color scale min", el("scatter plot viewer"), 18));
      await session.step(63, "When user sets \"Color\" property of scatter plot viewer to \"\"", () => setProperty(page, "Color", el("scatter plot viewer"), ""));
      await session.step(64, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Removing the grid's linear coloring gives the color scale its own range back (public#1434)", async () => {
      await session.step(68, "When user colors \"AGE\" column linearly from \"#0000FF\" to \"#FF0000\" over 30 to 40", () => colorLinearOver(page, "AGE", "#0000FF", "#FF0000", 30, 40));
      await session.step(69, "And user sets \"Color\" property of scatter plot viewer to \"AGE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "AGE"));
      await session.step(70, "Then the \"color scale min\" reading of scatter plot viewer should be 30", () => readingIs(page, "color scale min", el("scatter plot viewer"), 30));
      await session.step(71, "When user removes the coloring of \"AGE\" column", () => colorOff(page, "AGE"));
      await session.step(72, "Then \"AGE\" column should have no color coding", () => noColorCoding(page, "AGE"));
      await session.step(73, "And the color scale of scatter plot viewer should cover a wider range than before", () => scaleWider(page, el("scatter plot viewer")));
      await session.step(74, "And the \"color scale min\" reading of scatter plot viewer should be 18", () => readingIs(page, "color scale min", el("scatter plot viewer"), 18));
      await session.step(75, "When user sets \"Color\" property of scatter plot viewer to \"\"", () => setProperty(page, "Color", el("scatter plot viewer"), ""));
      await session.step(76, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An explicit X window comes back from a layout (public#2402)", async () => {
      await session.step(80, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xMin","100"],["xMax","200"]]), [["xMin","100"],["xMax","200"]]);
      await session.step(83, "Then the \"x axis min\" reading of scatter plot viewer should be 100", () => readingIs(page, "x axis min", el("scatter plot viewer"), 100));
      await session.step(84, "When user saves the layout of the current table view", () => saveLayout(page));
      await session.step(85, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xMin",""],["xMax",""]]), [["xMin",""],["xMax",""]]);
      await session.step(88, "Then the \"x axis min\" reading of scatter plot viewer should not be \"100\"", () => readingDoesNotRead(page, "x axis min", el("scatter plot viewer"), "100"));
      await session.step(89, "When user loads the saved layout", () => loadLayout(page));
      await session.step(90, "Then the \"x axis min\" reading of scatter plot viewer should be 100", () => readingIs(page, "x axis min", el("scatter plot viewer"), 100));
      await session.step(91, "And the \"x axis max\" reading of scatter plot viewer should be 200", () => readingIs(page, "x axis max", el("scatter plot viewer"), 200));
      await session.step(92, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xMin",""],["xMax",""]]), [["xMin",""],["xMax",""]]);
      await session.step(95, "Then the \"x axis min\" reading of scatter plot viewer should not be \"100\"", () => readingDoesNotRead(page, "x axis min", el("scatter plot viewer"), "100"));
      await session.step(96, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The window a zoom left comes back with the project (core#375)", async () => {
      await session.step(101, "When user drags the min handle of the \"x\" range slider of scatter plot viewer by 40 pixels", () => dragRangeHandle(page, "min", "x", el("scatter plot viewer"), 40));
      await session.step(102, "Then the \"x axis min\" reading of scatter plot viewer should be higher than before", () => readingHigher(page, "x axis min", el("scatter plot viewer")));
      await session.step(103, "And fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(104, "When user remembers the \"x axis min\" reading of scatter plot viewer", () => rememberReading(page, "x axis min", el("scatter plot viewer")));
      await session.step(105, "And user remembers how many rows pass the filter", () => rememberPassing(page));
      await session.step(106, "And user saves the current view as project \"bdd-scatter-zoom-{time}\"", () => saveAsProject(page, session.text("bdd-scatter-zoom-{time}")));
      await session.step(107, "And user closes all views", () => closeAllViews(page));
      await session.step(108, "And user opens the \"bdd-scatter-zoom-{time}\" project", () => openProject(page, session.text("bdd-scatter-zoom-{time}")));
      await session.step(109, "Then the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(110, "And the \"x axis min\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "x axis min", el("scatter plot viewer")));
      await session.step(111, "And as many rows as remembered should pass the filter", () => passingAsRemembered(page));
      await session.step(112, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(113, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(114, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
