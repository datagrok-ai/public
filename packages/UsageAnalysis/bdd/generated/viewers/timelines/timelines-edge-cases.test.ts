/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/timelines/timelines-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.timelines, public#617, public#1973]
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
import {columnType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterBetween, filterPasses, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, legendLists, noErrors, painted, readingIs, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Timelines edge cases: a cleared color, long lane names, no rows", () => {
  const session = feature(test, "features/viewers/timelines/timelines-edge-cases.feature", import.meta.url);
  test("Clearing Color takes the legend away without an error", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.timelines", "@realizes:public#617"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(11, "And user opens a table \"events\" with:", () => openTableOf(page, "events", [["subject","start","end","grade","dose"],["s1","2020-01-01","2020-01-10","mild","1"],["s2","2020-01-05","","mild","2"],["s3","","2020-01-20","high","3"],["s4","2020-01-03","2020-01-08","high","4"]]), [["subject","start","end","grade","dose"],["s1","2020-01-01","2020-01-10","mild","1"],["s2","2020-01-05","","mild","2"],["s3","","2020-01-20","high","3"],["s4","2020-01-03","2020-01-08","high","4"]]);
    await session.step(17, "And user adds a timelines viewer with:", () => addViewerWith(page, "timelines", [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]), [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]);
    await session.step(21, "Then \"start\" column should have type \"datetime\"", () => columnType(page, "start", "datetime"));
    await session.step(22, "And timelines viewer should be painted", () => painted(page, el("timelines viewer")));
    await session.step(26, "When user sets \"colorColumnName\" property of timelines viewer to \"grade\"", () => setProperty(page, "colorColumnName", el("timelines viewer"), "grade"));
    await session.step(27, "Then the legend of timelines viewer should list 2 items", () => legendLists(page, el("timelines viewer"), 2));
    await session.step(28, "When user sets \"colorColumnName\" property of timelines viewer to \"\"", () => setProperty(page, "colorColumnName", el("timelines viewer"), ""));
    await session.step(29, "Then legend of timelines viewer should be hidden", () => shouldBe(page, el("legend of timelines viewer"), "hidden"));
    await session.step(30, "And timelines viewer should be painted", () => painted(page, el("timelines viewer")));
    await session.step(31, "And no errors should have been logged", () => noErrors(page));
  });
  test("A filter that removes every row leaves no interval and no error", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.timelines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(11, "And user opens a table \"events\" with:", () => openTableOf(page, "events", [["subject","start","end","grade","dose"],["s1","2020-01-01","2020-01-10","mild","1"],["s2","2020-01-05","","mild","2"],["s3","","2020-01-20","high","3"],["s4","2020-01-03","2020-01-08","high","4"]]), [["subject","start","end","grade","dose"],["s1","2020-01-01","2020-01-10","mild","1"],["s2","2020-01-05","","mild","2"],["s3","","2020-01-20","high","3"],["s4","2020-01-03","2020-01-08","high","4"]]);
    await session.step(17, "And user adds a timelines viewer with:", () => addViewerWith(page, "timelines", [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]), [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]);
    await session.step(21, "Then \"start\" column should have type \"datetime\"", () => columnType(page, "start", "datetime"));
    await session.step(22, "And timelines viewer should be painted", () => painted(page, el("timelines viewer")));
    await session.step(34, "When user filters rows where \"dose\" is between 50 and 60", () => filterBetween(page, "dose", 50, 60));
    await session.step(35, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(36, "And the \"rows shown\" reading of timelines viewer should be 0", () => readingIs(page, "rows shown", el("timelines viewer"), 0));
    await session.step(37, "And the \"intervals\" reading of timelines viewer should be 0", () => readingIs(page, "intervals", el("timelines viewer"), 0));
    await session.step(38, "When user resets the filter", () => resetFilter(page));
    await session.step(39, "Then the \"rows shown\" reading of timelines viewer should be 4", () => readingIs(page, "rows shown", el("timelines viewer"), 4));
    await session.step(40, "And no errors should have been logged", () => noErrors(page));
  });
  test("Lanes with 200-character names are drawn (public#1973)", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.timelines", "@realizes:public#1973"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(11, "And user opens a table \"events\" with:", () => openTableOf(page, "events", [["subject","start","end","grade","dose"],["s1","2020-01-01","2020-01-10","mild","1"],["s2","2020-01-05","","mild","2"],["s3","","2020-01-20","high","3"],["s4","2020-01-03","2020-01-08","high","4"]]), [["subject","start","end","grade","dose"],["s1","2020-01-01","2020-01-10","mild","1"],["s2","2020-01-05","","mild","2"],["s3","","2020-01-20","high","3"],["s4","2020-01-03","2020-01-08","high","4"]]);
    await session.step(17, "And user adds a timelines viewer with:", () => addViewerWith(page, "timelines", [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]), [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]);
    await session.step(21, "Then \"start\" column should have type \"datetime\"", () => columnType(page, "start", "datetime"));
    await session.step(22, "And timelines viewer should be painted", () => painted(page, el("timelines viewer")));
    await session.step(44, "Given user opens a table \"long-lanes\" with:", () => openTableOf(page, "long-lanes", [["subject","start","end"],["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa1","2020-01-01","2020-01-10"],["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa2","2020-01-05","2020-01-12"]]), [["subject","start","end"],["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa1","2020-01-01","2020-01-10"],["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa2","2020-01-05","2020-01-12"]]);
    await session.step(48, "And user adds a timelines viewer with:", () => addViewerWith(page, "timelines", [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]), [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]);
    await session.step(52, "Then the \"lanes\" reading of timelines viewer should be 2", () => readingIs(page, "lanes", el("timelines viewer"), 2));
    await session.step(53, "And the \"intervals\" reading of timelines viewer should be 2", () => readingIs(page, "intervals", el("timelines viewer"), 2));
    await session.step(54, "And timelines viewer should be painted", () => painted(page, el("timelines viewer")));
    await session.step(55, "And no errors should have been logged", () => noErrors(page));
  });
});
