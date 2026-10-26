/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart]
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
import {addCalculated, selectWhereIs, selectWhereOneOf} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, legendLists, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart legend regressions", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-regressions.feature", import.meta.url);
  test("The split legend lists no Y / split pair that has no data (core#825)", {tag: ["@viewers", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "When user adds a calculated column \"H_F\" with formula \"if(${SEX} == \\\"F\\\", ${HEIGHT}, null)\"", () => addCalculated(page, "H_F", "if(${SEX} == \"F\", ${HEIGHT}, null)"));
    await session.step(17, "And user adds a calculated column \"W_M\" with formula \"if(${SEX} == \\\"M\\\", ${WEIGHT}, null)\"", () => addCalculated(page, "W_M", "if(${SEX} == \"M\", ${WEIGHT}, null)"));
    await session.step(18, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","H_F, W_M"],["splitColumnNames","SEX"],["multiAxis","true"],["yGlobalScale","true"]]), [["xColumnName","AGE"],["yColumnNames","H_F, W_M"],["splitColumnNames","SEX"],["multiAxis","true"],["yGlobalScale","true"]]);
    await session.step(24, "Then the legend of line chart viewer should list 2 items", () => legendLists(page, el("line chart viewer"), 2));
    await session.step(25, "And \"H_F / F\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"H_F / F\" legend item in legend of line chart viewer"), "visible"));
    await session.step(26, "And \"W_M / M\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"W_M / M\" legend item in legend of line chart viewer"), "visible"));
    await session.step(27, "And \"H_F / M\" legend item in legend of line chart viewer should be absent", () => shouldBe(page, el("\"H_F / M\" legend item in legend of line chart viewer"), "absent"));
    await session.step(28, "And \"W_M / F\" legend item in legend of line chart viewer should be absent", () => shouldBe(page, el("\"W_M / F\" legend item in legend of line chart viewer"), "absent"));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
  });
  test("Under Row Source = Selected the split legend lists the categories of the selected rows (public#1651)", {tag: ["@viewers", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(33, "When user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["splitColumnNames","RACE"],["Row Source","Selected"]]), [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["splitColumnNames","RACE"],["Row Source","Selected"]]);
    await session.step(38, "And user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
    await session.step(39, "Then the legend of line chart viewer should list 1 item", () => legendLists(page, el("line chart viewer"), 1));
    await session.step(40, "And \"Asian\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"Asian\" legend item in legend of line chart viewer"), "visible"));
    await session.step(41, "When user selects rows where \"RACE\" is one of \"Asian, Black\"", () => selectWhereOneOf(page, "RACE", "Asian, Black"));
    await session.step(42, "Then the legend of line chart viewer should list 2 items", () => legendLists(page, el("line chart viewer"), 2));
    await session.step(43, "And \"Black\" legend item in legend of line chart viewer should be visible", () => shouldBe(page, el("\"Black\" legend item in legend of line chart viewer"), "visible"));
    await session.step(44, "And \"Caucasian\" legend item in legend of line chart viewer should be absent", () => shouldBe(page, el("\"Caucasian\" legend item in legend of line chart viewer"), "absent"));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
  });
});
