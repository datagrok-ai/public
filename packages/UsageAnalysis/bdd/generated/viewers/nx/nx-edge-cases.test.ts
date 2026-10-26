/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/nx/nx-edge-cases.feature
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
import {columnComplete, columnIncomplete, hasNoColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, pickFromAreaContextMenu, reportsNoError, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("NX edge cases — a column a calculated column reads is removed while a viewer shows the result", () => {
  const session = feature(test, "features/viewers/nx/nx-edge-cases.feature", import.meta.url);
  test("Removing the column a calculated column reads empties the result and names it in a balloon", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "When user adds a calculated column \"AGE2\" with formula \"${AGE} * 2\"", () => addCalculated(page, "AGE2", "${AGE} * 2"));
    await session.step(16, "Then the value of \"AGE2\" column in row 1 should be \"52\"", () => valueInRow(page, "AGE2", 1, "52"));
    await session.step(17, "And \"AGE2\" column should have no missing values", () => columnComplete(page, "AGE2"));
    await session.step(18, "When user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE2"],["yColumnName","WEIGHT"]]), [["xColumnName","AGE2"],["yColumnName","WEIGHT"]]);
    await session.step(21, "And user picks \"Remove\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header AGE", el("grid")));
    await session.step(22, "Then the table should not have a column \"AGE\"", () => hasNoColumn(page, "AGE"));
    await session.step(23, "And a warning balloon containing \"AGE2 uses the removed column and lost its values\" should have been shown", () => warningBalloonText(page, "AGE2 uses the removed column and lost its values"));
    await session.step(24, "And the value of \"AGE2\" column in row 1 should be \"\"", () => valueInRow(page, "AGE2", 1, ""));
    await session.step(25, "And \"AGE2\" column should have missing values", () => columnIncomplete(page, "AGE2"));
    await session.step(26, "And scatter plot viewer should report no error", () => reportsNoError(page, el("scatter plot viewer")));
    await session.step(27, "And no errors should have been logged", () => noErrors(page));
  });
});
