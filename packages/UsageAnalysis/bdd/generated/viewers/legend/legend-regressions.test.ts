/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/legend/legend-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.legend, public#1209, public#2505, public#2965, public#3289]
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
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {categoricalColorIs, colorOff, noColorCoding} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, legendDocked, legendItemColor, legendLists, loadLayout, noErrors, pickColorSwatch, saveLayout, setProperties, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Legend regressions on the scatter plot", () => {
  const session = feature(test, "features/viewers/legend/legend-regressions.feature", import.meta.url);
  test("Legend regressions on the scatter plot", {tag: ["@journey", "@viewers", "@realizes:viewers.legend", "@realizes:public#1209", "@realizes:public#2505", "@realizes:public#2965", "@realizes:public#3289"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(19, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await run.scenario("A viewer Filter that empties the plot and lets rows back brings the legend back (public#1209)", async () => {
      await session.step(24, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Color","RACE"],["Legend Visibility","Auto"]]), [["Color","RACE"],["Legend Visibility","Auto"]]);
      await session.step(27, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(28, "When user sets \"Filter\" property of scatter plot viewer to '${RACE} == \"none\"'", () => setProperty(page, "Filter", el("scatter plot viewer"), "${RACE} == \"none\""));
      await session.step(29, "Then scatter plot viewer should show 0 rows", () => showsRows(page, el("scatter plot viewer"), 0));
      await session.step(30, "When user sets \"Filter\" property of scatter plot viewer to '${RACE} in [\"Asian\", \"Black\"]'", () => setProperty(page, "Filter", el("scatter plot viewer"), "${RACE} in [\"Asian\", \"Black\"]"));
      await session.step(31, "Then scatter plot viewer should show 38 rows", () => showsRows(page, el("scatter plot viewer"), 38));
      await session.step(32, "And the legend of scatter plot viewer should list 2 items", () => legendLists(page, el("scatter plot viewer"), 2));
      await session.step(33, "And \"Asian\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"Asian\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(34, "And \"Caucasian\" legend item in legend of scatter plot viewer should be absent", () => shouldBe(page, el("\"Caucasian\" legend item in legend of scatter plot viewer"), "absent"));
      await session.step(35, "When user saves the layout of the current table view", () => saveLayout(page));
      await session.step(36, "And user loads the saved layout", () => loadLayout(page));
      await session.step(37, "Then scatter plot viewer should show 38 rows", () => showsRows(page, el("scatter plot viewer"), 38));
      await session.step(38, "And the legend of scatter plot viewer should list 2 items", () => legendLists(page, el("scatter plot viewer"), 2));
      await session.step(39, "And \"Black\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(40, "And \"Caucasian\" legend item in legend of scatter plot viewer should be absent", () => shouldBe(page, el("\"Caucasian\" legend item in legend of scatter plot viewer"), "absent"));
      await session.step(41, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Filter",""],["Color",""]]), [["Filter",""],["Color",""]]);
      await session.step(44, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A viewer Filter set after Markers drops the filtered-out categories from the marker legend (public#2505)", async () => {
      await session.step(50, "When user sets \"Markers\" property of scatter plot viewer to \"RACE\"", () => setProperty(page, "Markers", el("scatter plot viewer"), "RACE"));
      await session.step(51, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(52, "When user sets \"Filter\" property of scatter plot viewer to '${RACE} in [\"Asian\", \"Black\"]'", () => setProperty(page, "Filter", el("scatter plot viewer"), "${RACE} in [\"Asian\", \"Black\"]"));
      await session.step(53, "Then scatter plot viewer should show 38 rows", () => showsRows(page, el("scatter plot viewer"), 38));
      await session.step(54, "And the legend of scatter plot viewer should list 2 items", () => legendLists(page, el("scatter plot viewer"), 2));
      await session.step(55, "And \"Caucasian\" legend item in legend of scatter plot viewer should be absent", () => shouldBe(page, el("\"Caucasian\" legend item in legend of scatter plot viewer"), "absent"));
      await session.step(56, "And \"Other\" legend item in legend of scatter plot viewer should be absent", () => shouldBe(page, el("\"Other\" legend item in legend of scatter plot viewer"), "absent"));
      await session.step(57, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Filter",""],["Markers",""]]), [["Filter",""],["Markers",""]]);
      await session.step(60, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Color and Markers on the same column a color is picked from the legend (public#2965)", async () => {
      await session.step(66, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Color","RACE"],["Markers","RACE"],["Legend Visibility","Always"]]), [["Color","RACE"],["Markers","RACE"],["Legend Visibility","Always"]]);
      await session.step(70, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(71, "When user hovers over \"Asian\" legend item in legend of scatter plot viewer", () => hoverOver(page, el("\"Asian\" legend item in legend of scatter plot viewer")));
      await session.step(72, "And user clicks on color picker icon", () => clickOn(page, el("color picker icon")));
      await session.step(73, "Then \"Asian\" dialog should be visible", () => shouldBe(page, el("\"Asian\" dialog"), "visible"));
      await session.step(74, "When user picks the color \"#9467BD\" in the color picker dialog", () => pickColorSwatch(page, "#9467BD"));
      await session.step(75, "And user clicks on OK button in \"Asian\" dialog", () => clickOn(page, el("OK button in \"Asian\" dialog")));
      await session.step(76, "Then \"Asian\" dialog should be absent", () => shouldBe(page, el("\"Asian\" dialog"), "absent"));
      await session.step(77, "And the categorical color of \"Asian\" in \"RACE\" column should be \"#9467BD\"", () => categoricalColorIs(page, "Asian", "RACE", "#9467BD"));
      await session.step(78, "And the \"Asian\" item in the legend of scatter plot viewer should be colored \"#9467BD\"", () => legendItemColor(page, "Asian", el("scatter plot viewer"), "#9467BD"));
      await session.step(79, "When user removes the coloring of \"RACE\" column", () => colorOff(page, "RACE"));
      await session.step(80, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Color",""],["Markers",""],["Legend Visibility","Auto"]]), [["Color",""],["Markers",""],["Legend Visibility","Auto"]]);
      await session.step(84, "Then \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Color column of a thousand categories leaves the Markers section on the screen (public#3289)", async () => {
      await session.step(90, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Color","USUBJID"],["Markers","SEX"],["Legend Visibility","Always"]]), [["Color","USUBJID"],["Markers","SEX"],["Legend Visibility","Always"]]);
      await session.step(94, "Then the legend of scatter plot viewer should be docked", () => legendDocked(page, el("scatter plot viewer")));
      await session.step(95, "And \"F\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"F\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(96, "And \"M\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"M\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(97, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Color",""],["Markers",""],["Legend Visibility","Auto"]]), [["Color",""],["Markers",""],["Legend Visibility","Auto"]]);
      await session.step(101, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(102, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
