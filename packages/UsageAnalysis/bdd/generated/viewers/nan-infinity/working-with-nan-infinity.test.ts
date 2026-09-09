/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/nan-infinity/working-with-nan-infinity.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {columnIncomplete} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {setCell} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, loadLayout, noErrors, painted, paintedInColors, propertyShouldBe, readingAsRemembered, readingAtLeast, readingFinite, rememberReading, repainted, saveLayout, setProperties, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A NaN and an Infinity in the plotted columns", () => {
  const session = feature(test, "features/viewers/nan-infinity/working-with-nan-infinity.feature", import.meta.url);
  test("A NaN and an Infinity in the plotted columns", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","HEIGHT"],["Y","WEIGHT"],["Size","AGE"],["Color","RACE"],["Show Regression Line","true"]]));
    await session.step(19, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(20, "And \"HEIGHT\" column should have missing values", () => columnIncomplete(page, "HEIGHT"));
    await run.scenario("A NaN in the X column drops its row and leaves the axes finite", async () => {
      await session.step(23, "When user sets \"HEIGHT\" column in row 1 to \"NaN\"", () => setCell(page, "HEIGHT", 1, "NaN"));
      await session.step(24, "Then scatter plot viewer should show 871 rows", () => showsRows(page, el("scatter plot viewer"), 871));
      await session.step(25, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(26, "And the \"x axis min\" reading of scatter plot viewer should be a finite number", () => readingFinite(page, "x axis min", el("scatter plot viewer")));
      await session.step(27, "And the \"x axis max\" reading of scatter plot viewer should be a finite number", () => readingFinite(page, "x axis max", el("scatter plot viewer")));
      await session.step(28, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An Infinity in the Y column keeps the Y axis on a finite range", async () => {
      await session.step(31, "When user sets \"WEIGHT\" column in row 2 to \"Infinity\"", () => setCell(page, "WEIGHT", 2, "Infinity"));
      await session.step(32, "Then scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(33, "And the \"y axis min\" reading of scatter plot viewer should be a finite number", () => readingFinite(page, "y axis min", el("scatter plot viewer")));
      await session.step(34, "And the \"y axis max\" reading of scatter plot viewer should be a finite number", () => readingFinite(page, "y axis max", el("scatter plot viewer")));
      await session.step(35, "And the \"regression lines\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "regression lines", el("scatter plot viewer"), 1));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The square marker and the regression line survive both values", async () => {
      await session.step(39, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Marker Type","square"]]));
      await session.step(41, "Then \"Marker Type\" property of scatter plot viewer should be \"square\"", () => propertyShouldBe(page, "Marker Type", el("scatter plot viewer"), "square"));
      await session.step(42, "And scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(43, "And scatter plot viewer should be painted in at least 2 colors", () => paintedInColors(page, el("scatter plot viewer"), 2));
      await session.step(44, "And the \"regression lines\" reading of scatter plot viewer should be at least 1", () => readingAtLeast(page, "regression lines", el("scatter plot viewer"), 1));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A layout round-trip brings the plot back over the same data", async () => {
      await session.step(48, "When user remembers the \"rows shown\" reading of scatter plot viewer", () => rememberReading(page, "rows shown", el("scatter plot viewer")));
      await session.step(49, "And user saves the layout of the current table view", () => saveLayout(page));
      await session.step(50, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Marker Type","circle"],["X","AGE"]]));
      await session.step(53, "Then \"X\" property of scatter plot viewer should be \"AGE\"", () => propertyShouldBe(page, "X", el("scatter plot viewer"), "AGE"));
      await session.step(54, "When user loads the saved layout", () => loadLayout(page));
      await session.step(55, "Then \"X\" property of scatter plot viewer should be \"HEIGHT\"", () => propertyShouldBe(page, "X", el("scatter plot viewer"), "HEIGHT"));
      await session.step(56, "And \"Marker Type\" property of scatter plot viewer should be \"square\"", () => propertyShouldBe(page, "Marker Type", el("scatter plot viewer"), "square"));
      await session.step(57, "And the \"rows shown\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "rows shown", el("scatter plot viewer")));
      await session.step(58, "And the \"x axis max\" reading of scatter plot viewer should be a finite number", () => readingFinite(page, "x axis max", el("scatter plot viewer")));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The two values put back restore every row", async () => {
      await session.step(62, "When user sets \"HEIGHT\" column in row 1 to \"174.705\"", () => setCell(page, "HEIGHT", 1, "174.705"));
      await session.step(63, "And user sets \"WEIGHT\" column in row 2 to \"64\"", () => setCell(page, "WEIGHT", 2, "64"));
      await session.step(64, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(65, "And the \"y axis max\" reading of scatter plot viewer should be a finite number", () => readingFinite(page, "y axis max", el("scatter plot viewer")));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
