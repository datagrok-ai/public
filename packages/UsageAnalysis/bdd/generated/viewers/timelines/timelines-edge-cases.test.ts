/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/timelines/timelines-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.timelines]
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

test.describe("Timelines edge cases: missing dates, a cleared color, long lane names, no rows", () => {
  const session = feature(test, "features/viewers/timelines/timelines-edge-cases.feature", import.meta.url);
  test("Events with a missing date are drawn only with Show Open Intervals on", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.timelines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And user opens a table \"events\" with:", () => openTableOf(page, "events", [["subject","start","end","grade","dose"],["s1","2020-01-01","2020-01-10","mild","1"],["s2","2020-01-05","","mild","2"],["s3","","2020-01-20","high","3"],["s4","2020-01-03","2020-01-08","high","4"]]), [["subject","start","end","grade","dose"],["s1","2020-01-01","2020-01-10","mild","1"],["s2","2020-01-05","","mild","2"],["s3","","2020-01-20","high","3"],["s4","2020-01-03","2020-01-08","high","4"]]);
    await session.step(20, "And user adds a timelines viewer with:", () => addViewerWith(page, "timelines", [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]), [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]);
    await session.step(24, "Then \"start\" column should have type \"datetime\"", () => columnType(page, "start", "datetime"));
    await session.step(25, "And timelines viewer should be painted", () => painted(page, el("timelines viewer")));
    await session.step(28, "When user sets \"showOpenIntervals\" property of timelines viewer to \"false\"", () => setProperty(page, "showOpenIntervals", el("timelines viewer"), "false"));
    await session.step(29, "Then the \"intervals\" reading of timelines viewer should be 2", () => readingIs(page, "intervals", el("timelines viewer"), 2));
    await session.step(30, "When user sets \"showOpenIntervals\" property of timelines viewer to \"true\"", () => setProperty(page, "showOpenIntervals", el("timelines viewer"), "true"));
    await session.step(31, "Then the \"intervals\" reading of timelines viewer should be 4", () => readingIs(page, "intervals", el("timelines viewer"), 4));
    await session.step(32, "When user sets \"showOpenIntervals\" property of timelines viewer to \"false\"", () => setProperty(page, "showOpenIntervals", el("timelines viewer"), "false"));
    await session.step(33, "Then the \"intervals\" reading of timelines viewer should be 2", () => readingIs(page, "intervals", el("timelines viewer"), 2));
    await session.step(34, "And no errors should have been logged", () => noErrors(page));
  });
  test("Clearing Color takes the legend away without an error", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.timelines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And user opens a table \"events\" with:", () => openTableOf(page, "events", [["subject","start","end","grade","dose"],["s1","2020-01-01","2020-01-10","mild","1"],["s2","2020-01-05","","mild","2"],["s3","","2020-01-20","high","3"],["s4","2020-01-03","2020-01-08","high","4"]]), [["subject","start","end","grade","dose"],["s1","2020-01-01","2020-01-10","mild","1"],["s2","2020-01-05","","mild","2"],["s3","","2020-01-20","high","3"],["s4","2020-01-03","2020-01-08","high","4"]]);
    await session.step(20, "And user adds a timelines viewer with:", () => addViewerWith(page, "timelines", [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]), [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]);
    await session.step(24, "Then \"start\" column should have type \"datetime\"", () => columnType(page, "start", "datetime"));
    await session.step(25, "And timelines viewer should be painted", () => painted(page, el("timelines viewer")));
    await session.step(37, "When user sets \"colorColumnName\" property of timelines viewer to \"grade\"", () => setProperty(page, "colorColumnName", el("timelines viewer"), "grade"));
    await session.step(38, "Then the legend of timelines viewer should list 2 items", () => legendLists(page, el("timelines viewer"), 2));
    await session.step(39, "When user sets \"colorColumnName\" property of timelines viewer to \"\"", () => setProperty(page, "colorColumnName", el("timelines viewer"), ""));
    await session.step(40, "Then legend of timelines viewer should be hidden", () => shouldBe(page, el("legend of timelines viewer"), "hidden"));
    await session.step(41, "And timelines viewer should be painted", () => painted(page, el("timelines viewer")));
    await session.step(42, "And no errors should have been logged", () => noErrors(page));
  });
  test("A filter that removes every row leaves no interval and no error", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.timelines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And user opens a table \"events\" with:", () => openTableOf(page, "events", [["subject","start","end","grade","dose"],["s1","2020-01-01","2020-01-10","mild","1"],["s2","2020-01-05","","mild","2"],["s3","","2020-01-20","high","3"],["s4","2020-01-03","2020-01-08","high","4"]]), [["subject","start","end","grade","dose"],["s1","2020-01-01","2020-01-10","mild","1"],["s2","2020-01-05","","mild","2"],["s3","","2020-01-20","high","3"],["s4","2020-01-03","2020-01-08","high","4"]]);
    await session.step(20, "And user adds a timelines viewer with:", () => addViewerWith(page, "timelines", [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]), [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]);
    await session.step(24, "Then \"start\" column should have type \"datetime\"", () => columnType(page, "start", "datetime"));
    await session.step(25, "And timelines viewer should be painted", () => painted(page, el("timelines viewer")));
    await session.step(45, "When user filters rows where \"dose\" is between 50 and 60", () => filterBetween(page, "dose", 50, 60));
    await session.step(46, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(47, "And the \"rows shown\" reading of timelines viewer should be 0", () => readingIs(page, "rows shown", el("timelines viewer"), 0));
    await session.step(48, "And the \"intervals\" reading of timelines viewer should be 0", () => readingIs(page, "intervals", el("timelines viewer"), 0));
    await session.step(49, "When user resets the filter", () => resetFilter(page));
    await session.step(50, "Then the \"rows shown\" reading of timelines viewer should be 4", () => readingIs(page, "rows shown", el("timelines viewer"), 4));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
  test("Lanes with 200-character names are drawn (public#1973)", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.timelines"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And user opens a table \"events\" with:", () => openTableOf(page, "events", [["subject","start","end","grade","dose"],["s1","2020-01-01","2020-01-10","mild","1"],["s2","2020-01-05","","mild","2"],["s3","","2020-01-20","high","3"],["s4","2020-01-03","2020-01-08","high","4"]]), [["subject","start","end","grade","dose"],["s1","2020-01-01","2020-01-10","mild","1"],["s2","2020-01-05","","mild","2"],["s3","","2020-01-20","high","3"],["s4","2020-01-03","2020-01-08","high","4"]]);
    await session.step(20, "And user adds a timelines viewer with:", () => addViewerWith(page, "timelines", [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]), [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]);
    await session.step(24, "Then \"start\" column should have type \"datetime\"", () => columnType(page, "start", "datetime"));
    await session.step(25, "And timelines viewer should be painted", () => painted(page, el("timelines viewer")));
    await session.step(54, "Given user opens a table \"long-lanes\" with:", () => openTableOf(page, "long-lanes", [["subject","start","end"],["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa1","2020-01-01","2020-01-10"],["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa2","2020-01-05","2020-01-12"]]), [["subject","start","end"],["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa1","2020-01-01","2020-01-10"],["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa2","2020-01-05","2020-01-12"]]);
    await session.step(58, "And user adds a timelines viewer with:", () => addViewerWith(page, "timelines", [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]), [["splitByColumnName","subject"],["startColumnName","start"],["endColumnName","end"]]);
    await session.step(62, "Then the \"lanes\" reading of timelines viewer should be 2", () => readingIs(page, "lanes", el("timelines viewer"), 2));
    await session.step(63, "And the \"intervals\" reading of timelines viewer should be 2", () => readingIs(page, "intervals", el("timelines viewer"), 2));
    await session.step(64, "And timelines viewer should be painted", () => painted(page, el("timelines viewer")));
    await session.step(65, "And no errors should have been logged", () => noErrors(page));
  });
});
