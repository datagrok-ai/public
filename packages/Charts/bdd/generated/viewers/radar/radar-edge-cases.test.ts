/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/radar/radar-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.radar]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {deleteSelected, filterBetween, filterPasses, removeColumn, rowCount, selectAllRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, noErrors, painted, propertyShouldBe, repainted} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Radar edge cases — Values columns removed from the table, an empty row set, every row deleted", () => {
  const session = feature(test, "features/viewers/radar/radar-edge-cases.feature", import.meta.url);
  test("Values columns removed from the table drop out of the axes, and with none left the Radar says why", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.radar"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a radar viewer", () => addViewer(page, "radar"));
    await session.step(17, "Then \"Values\" property of radar viewer should be \"AGE, HEIGHT, WEIGHT\"", () => propertyShouldBe(page, "Values", el("radar viewer"), "AGE, HEIGHT, WEIGHT"));
    await session.step(18, "And radar viewer should be painted", () => painted(page, el("radar viewer")));
    await session.step(21, "When user removes \"HEIGHT\" column", () => removeColumn(page, "HEIGHT"));
    await session.step(22, "Then \"Values\" property of radar viewer should be \"AGE, WEIGHT\"", () => propertyShouldBe(page, "Values", el("radar viewer"), "AGE, WEIGHT"));
    await session.step(23, "And \"The Radar viewer requires a minimum of 1 numerical column.\" text in radar viewer should be absent", () => shouldBe(page, el("\"The Radar viewer requires a minimum of 1 numerical column.\" text in radar viewer"), "absent"));
    await session.step(24, "And radar viewer should be painted", () => painted(page, el("radar viewer")));
    await session.step(25, "And no errors should have been logged", () => noErrors(page));
    await session.step(26, "When user removes \"AGE\" column", () => removeColumn(page, "AGE"));
    await session.step(27, "And user removes \"WEIGHT\" column", () => removeColumn(page, "WEIGHT"));
    await session.step(28, "Then \"The Radar viewer requires a minimum of 1 numerical column.\" text in radar viewer should be visible", () => shouldBe(page, el("\"The Radar viewer requires a minimum of 1 numerical column.\" text in radar viewer"), "visible"));
    await session.step(29, "And \"Values\" property of radar viewer should be \"\"", () => propertyShouldBe(page, "Values", el("radar viewer"), ""));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
  });
  test("With every row filtered out the Radar redraws without the column message and logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.radar"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a radar viewer", () => addViewer(page, "radar"));
    await session.step(17, "Then \"Values\" property of radar viewer should be \"AGE, HEIGHT, WEIGHT\"", () => propertyShouldBe(page, "Values", el("radar viewer"), "AGE, HEIGHT, WEIGHT"));
    await session.step(18, "And radar viewer should be painted", () => painted(page, el("radar viewer")));
    await session.step(33, "When user filters rows where \"AGE\" is between 200 and 300", () => filterBetween(page, "AGE", 200, 300));
    await session.step(34, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(35, "And radar viewer should have repainted", () => repainted(page, el("radar viewer")));
    await session.step(36, "And \"The Radar viewer requires a minimum of 1 numerical column.\" text in radar viewer should be absent", () => shouldBe(page, el("\"The Radar viewer requires a minimum of 1 numerical column.\" text in radar viewer"), "absent"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
  test("Once every row of the table is deleted the Radar redraws without the column message and logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.radar"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a radar viewer", () => addViewer(page, "radar"));
    await session.step(17, "Then \"Values\" property of radar viewer should be \"AGE, HEIGHT, WEIGHT\"", () => propertyShouldBe(page, "Values", el("radar viewer"), "AGE, HEIGHT, WEIGHT"));
    await session.step(18, "And radar viewer should be painted", () => painted(page, el("radar viewer")));
    await session.step(40, "When user selects all rows", () => selectAllRows(page));
    await session.step(41, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(42, "Then the table should have 0 rows", () => rowCount(page, 0));
    await session.step(43, "And radar viewer should have repainted", () => repainted(page, el("radar viewer")));
    await session.step(44, "And \"The Radar viewer requires a minimum of 1 numerical column.\" text in radar viewer should be absent", () => shouldBe(page, el("\"The Radar viewer requires a minimum of 1 numerical column.\" text in radar viewer"), "absent"));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
  });
  test("A row with NaN or Infinity under Column normalization leaves the rest of the Radar drawn and logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:charts.viewer.radar"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a radar viewer", () => addViewer(page, "radar"));
    await session.step(17, "Then \"Values\" property of radar viewer should be \"AGE, HEIGHT, WEIGHT\"", () => propertyShouldBe(page, "Values", el("radar viewer"), "AGE, HEIGHT, WEIGHT"));
    await session.step(18, "And radar viewer should be painted", () => painted(page, el("radar viewer")));
    await session.step(52, "Given user opens a table \"non-finite\" with:", () => openTableOf(page, "non-finite", [["name","a","b","c"],["one","1","4","7"],["two","NaN","5","8"],["three","Infinity","6","9"],["four","-Infinity","5","7"],["five","3","4","8"]]), [["name","a","b","c"],["one","1","4","7"],["two","NaN","5","8"],["three","Infinity","6","9"],["four","-Infinity","5","7"],["five","3","4","8"]]);
    await session.step(59, "When user adds a radar viewer", () => addViewer(page, "radar"));
    await session.step(60, "Then \"Normalization\" property of radar viewer should be \"Column\"", () => propertyShouldBe(page, "Normalization", el("radar viewer"), "Column"));
    await session.step(61, "And \"Values\" property of radar viewer should be \"a, b, c\"", () => propertyShouldBe(page, "Values", el("radar viewer"), "a, b, c"));
    await session.step(62, "And \"The Radar viewer requires a minimum of 1 numerical column.\" text in radar viewer should be absent", () => shouldBe(page, el("\"The Radar viewer requires a minimum of 1 numerical column.\" text in radar viewer"), "absent"));
    await session.step(63, "And radar viewer should be painted", () => painted(page, el("radar viewer")));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
  });
});
