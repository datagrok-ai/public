/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tile-viewer/tile-viewer-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tile-viewer]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {addTileViewerWith} from '../../../bindings/tile-viewer.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterBetween, filterPasses, removeColumn, renameColumn, resetFilter, selectNoRows, selectWhereIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hasArea, hasNoArea, noErrors, propertyShouldBe, readingIs, readingReads, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Tile viewer edge cases: no rows, one row, the lanes column removed or renamed", () => {
  const session = feature(test, "features/viewers/tile-viewer/tile-viewer-edge-cases.feature", import.meta.url);
  test("A filter that removes every row leaves no card and no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tile-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","RACE"]]), [["Lanes Column Name","RACE"]]);
    await session.step(17, "Then tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
    await session.step(18, "And the \"lane names\" reading of tile viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "lane names", el("tile viewer"), "Asian, Black, Caucasian, Other"));
    await session.step(21, "When user filters rows where \"AGE\" is between 500 and 600", () => filterBetween(page, "AGE", 500, 600));
    await session.step(22, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(23, "And tile viewer should show 0 rows", () => showsRows(page, el("tile viewer"), 0));
    await session.step(24, "And the \"tiles\" reading of tile viewer should be 0", () => readingIs(page, "tiles", el("tile viewer"), 0));
    await session.step(25, "And tile viewer should not have a \"tile of row 1\" area", () => hasNoArea(page, el("tile viewer"), "tile of row 1"));
    await session.step(26, "When user resets the filter", () => resetFilter(page));
    await session.step(27, "Then tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
    await session.step(28, "And tile viewer should have a \"tile of row 1\" area", () => hasArea(page, el("tile viewer"), "tile of row 1"));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
  });
  test("Row Source Selected with nothing selected shows no card and no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tile-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","RACE"]]), [["Lanes Column Name","RACE"]]);
    await session.step(17, "Then tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
    await session.step(18, "And the \"lane names\" reading of tile viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "lane names", el("tile viewer"), "Asian, Black, Caucasian, Other"));
    await session.step(32, "Given user selects no rows", () => selectNoRows(page));
    await session.step(33, "When user sets \"Row Source\" property of tile viewer to \"Selected\"", () => setProperty(page, "Row Source", el("tile viewer"), "Selected"));
    await session.step(34, "Then tile viewer should show 0 rows", () => showsRows(page, el("tile viewer"), 0));
    await session.step(35, "And the \"tiles\" reading of tile viewer should be 0", () => readingIs(page, "tiles", el("tile viewer"), 0));
    await session.step(36, "When user selects rows where \"RACE\" is \"Asian\"", () => selectWhereIs(page, "RACE", "Asian"));
    await session.step(37, "Then tile viewer should show 15 rows", () => showsRows(page, el("tile viewer"), 15));
    await session.step(38, "And the \"lane of row 10\" reading of tile viewer should be \"Asian\"", () => readingReads(page, "lane of row 10", el("tile viewer"), "Asian"));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("A one-row table shows its one card", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tile-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","RACE"]]), [["Lanes Column Name","RACE"]]);
    await session.step(17, "Then tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
    await session.step(18, "And the \"lane names\" reading of tile viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "lane names", el("tile viewer"), "Asian, Black, Caucasian, Other"));
    await session.step(42, "Given user opens a table \"one-row\" with:", () => openTableOf(page, "one-row", [["name","group","x"],["a","g1","42"]]), [["name","group","x"],["a","g1","42"]]);
    await session.step(45, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","group"]]), [["Lanes Column Name","group"]]);
    await session.step(47, "Then tile viewer should show 1 rows", () => showsRows(page, el("tile viewer"), 1));
    await session.step(48, "And the \"tiles\" reading of tile viewer should be 1", () => readingIs(page, "tiles", el("tile viewer"), 1));
    await session.step(49, "And the \"lane of row 1\" reading of tile viewer should be \"g1\"", () => readingReads(page, "lane of row 1", el("tile viewer"), "g1"));
    await session.step(50, "And no errors should have been logged", () => noErrors(page));
  });
  test("The lanes column renamed is followed and the lanes stay", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tile-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","RACE"]]), [["Lanes Column Name","RACE"]]);
    await session.step(17, "Then tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
    await session.step(18, "And the \"lane names\" reading of tile viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "lane names", el("tile viewer"), "Asian, Black, Caucasian, Other"));
    await session.step(53, "When user renames \"RACE\" column to \"ORIGIN\"", () => renameColumn(page, "RACE", "ORIGIN"));
    await session.step(54, "Then \"Lanes Column Name\" property of tile viewer should be \"ORIGIN\"", () => propertyShouldBe(page, "Lanes Column Name", el("tile viewer"), "ORIGIN"));
    await session.step(55, "And the \"lane names\" reading of tile viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "lane names", el("tile viewer"), "Asian, Black, Caucasian, Other"));
    await session.step(56, "And the \"lane of row 10\" reading of tile viewer should be \"Asian\"", () => readingReads(page, "lane of row 10", el("tile viewer"), "Asian"));
    await session.step(57, "And no errors should have been logged", () => noErrors(page));
  });
  test("The lanes column removed while the lanes are shown logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.tile-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","RACE"]]), [["Lanes Column Name","RACE"]]);
    await session.step(17, "Then tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
    await session.step(18, "And the \"lane names\" reading of tile viewer should be \"Asian, Black, Caucasian, Other\"", () => readingReads(page, "lane names", el("tile viewer"), "Asian, Black, Caucasian, Other"));
    await session.step(60, "When user removes \"RACE\" column", () => removeColumn(page, "RACE"));
    await session.step(61, "Then the table should not have a column \"RACE\"", () => hasNoColumn(page, "RACE"));
    await session.step(62, "And tile viewer should show 1000 rows", () => showsRows(page, el("tile viewer"), 1000));
    await session.step(63, "And no errors should have been logged", () => noErrors(page));
  });
});
