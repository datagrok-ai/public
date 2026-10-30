/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot-3d/scatter-plot-3d-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot3d, core#1027]
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
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterBetween, filterPasses, removeColumn, renameColumn, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, boundTable, noBalloons, noErrors, propertyShouldBe, readingAsRemembered, readingIs, readingNotAsRemembered, rememberReading, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("3D scatter plot edge cases: no rows, columns removed or renamed, a table swap", () => {
  const session = feature(test, "features/viewers/scatter-plot-3d/scatter-plot-3d-edge-cases.feature", import.meta.url);
  test("A filter that removes every row leaves the plot empty without an error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot3d"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a 3d scatter plot viewer with:", () => addViewerWith(page, "3d scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["zColumnName","HEIGHT"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["zColumnName","HEIGHT"]]);
    await session.step(19, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 872", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 872));
    await session.step(22, "When user remembers the \"scene signature\" reading of 3d scatter plot viewer", () => rememberReading(page, "scene signature", el("3d scatter plot viewer")));
    await session.step(23, "And user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(24, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(25, "And the \"rows shown\" reading of 3d scatter plot viewer should be 0", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 0));
    await session.step(26, "And the \"scene signature\" reading of 3d scatter plot viewer should not be as remembered", () => readingNotAsRemembered(page, "scene signature", el("3d scatter plot viewer")));
    await session.step(27, "When user resets the filter", () => resetFilter(page));
    await session.step(28, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 872", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 872));
    await session.step(29, "And the \"scene signature\" reading of 3d scatter plot viewer should be as remembered", () => readingAsRemembered(page, "scene signature", el("3d scatter plot viewer")));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
    await session.step(31, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Renaming the axis columns is followed by the plot's settings", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot3d"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a 3d scatter plot viewer with:", () => addViewerWith(page, "3d scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["zColumnName","HEIGHT"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["zColumnName","HEIGHT"]]);
    await session.step(19, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 872", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 872));
    await session.step(34, "When user renames \"AGE\" column to \"AGE_YEARS\"", () => renameColumn(page, "AGE", "AGE_YEARS"));
    await session.step(35, "And user renames \"HEIGHT\" column to \"HEIGHT_CM\"", () => renameColumn(page, "HEIGHT", "HEIGHT_CM"));
    await session.step(36, "Then \"xColumnName\" property of 3d scatter plot viewer should be \"AGE_YEARS\"", () => propertyShouldBe(page, "xColumnName", el("3d scatter plot viewer"), "AGE_YEARS"));
    await session.step(37, "And \"zColumnName\" property of 3d scatter plot viewer should be \"HEIGHT_CM\"", () => propertyShouldBe(page, "zColumnName", el("3d scatter plot viewer"), "HEIGHT_CM"));
    await session.step(38, "And the \"rows shown\" reading of 3d scatter plot viewer should be 872", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 872));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("Removing the Size column keeps the plot working without an error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot3d", "@realizes:core#1027"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a 3d scatter plot viewer with:", () => addViewerWith(page, "3d scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["zColumnName","HEIGHT"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["zColumnName","HEIGHT"]]);
    await session.step(19, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 872", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 872));
    await session.step(45, "When user sets \"Size\" property of 3d scatter plot viewer to \"WEIGHT\"", () => setProperty(page, "Size", el("3d scatter plot viewer"), "WEIGHT"));
    await session.step(46, "And user removes \"WEIGHT\" column", () => removeColumn(page, "WEIGHT"));
    await session.step(47, "Then no errors should have been logged", () => noErrors(page));
    await session.step(48, "When user remembers the \"scene signature\" reading of 3d scatter plot viewer", () => rememberReading(page, "scene signature", el("3d scatter plot viewer")));
    await session.step(49, "And user sets \"Marker Type\" property of 3d scatter plot viewer to \"box\"", () => setProperty(page, "Marker Type", el("3d scatter plot viewer"), "box"));
    await session.step(50, "Then the \"scene signature\" reading of 3d scatter plot viewer should not be as remembered", () => readingNotAsRemembered(page, "scene signature", el("3d scatter plot viewer")));
  });
  test("Removing the Color column keeps the chosen X, Y and Z", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot3d", "@realizes:core#1027"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a 3d scatter plot viewer with:", () => addViewerWith(page, "3d scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["zColumnName","HEIGHT"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["zColumnName","HEIGHT"]]);
    await session.step(19, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 872", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 872));
    await session.step(56, "When user sets \"Color\" property of 3d scatter plot viewer to \"RACE\"", () => setProperty(page, "Color", el("3d scatter plot viewer"), "RACE"));
    await session.step(57, "And user removes \"RACE\" column", () => removeColumn(page, "RACE"));
    await session.step(58, "Then \"xColumnName\" property of 3d scatter plot viewer should be \"AGE\"", () => propertyShouldBe(page, "xColumnName", el("3d scatter plot viewer"), "AGE"));
    await session.step(59, "And \"yColumnName\" property of 3d scatter plot viewer should be \"WEIGHT\"", () => propertyShouldBe(page, "yColumnName", el("3d scatter plot viewer"), "WEIGHT"));
    await session.step(60, "And \"zColumnName\" property of 3d scatter plot viewer should be \"HEIGHT\"", () => propertyShouldBe(page, "zColumnName", el("3d scatter plot viewer"), "HEIGHT"));
    await session.step(61, "And no errors should have been logged", () => noErrors(page));
  });
  test("A table swap and back keeps the on-viewer column selectors", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot3d", "@realizes:core#1027"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a 3d scatter plot viewer with:", () => addViewerWith(page, "3d scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["zColumnName","HEIGHT"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["zColumnName","HEIGHT"]]);
    await session.step(19, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 872", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 872));
    await session.step(67, "Given user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(68, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
    await session.step(69, "When user sets \"table\" property of 3d scatter plot viewer to \"cars\"", () => setProperty(page, "table", el("3d scatter plot viewer"), "cars"));
    await session.step(70, "Then 3d scatter plot viewer should be bound to table \"cars\"", () => boundTable(page, el("3d scatter plot viewer"), "cars"));
    await session.step(71, "When user sets \"table\" property of 3d scatter plot viewer to \"demog-1000\"", () => setProperty(page, "table", el("3d scatter plot viewer"), "demog-1000"));
    await session.step(72, "Then 3d scatter plot viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("3d scatter plot viewer"), "demog-1000"));
    await session.step(73, "And X column input in 3d scatter plot viewer should be visible", () => shouldBe(page, el("X column input in 3d scatter plot viewer"), "visible"));
    await session.step(74, "And Y column input in 3d scatter plot viewer should be visible", () => shouldBe(page, el("Y column input in 3d scatter plot viewer"), "visible"));
    await session.step(75, "And Z column input in 3d scatter plot viewer should be visible", () => shouldBe(page, el("Z column input in 3d scatter plot viewer"), "visible"));
    await session.step(76, "And Color column input in 3d scatter plot viewer should be visible", () => shouldBe(page, el("Color column input in 3d scatter plot viewer"), "visible"));
  });
});
