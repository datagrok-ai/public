/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/general/table-manager-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {noTableLeft, noTableViewOpen} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Table Manager edge cases: no table open", () => {
  const session = feature(test, "features/general/table-manager-edge-cases.feature", import.meta.url);
  test("With no table open the manager docks with an empty list and no error", {tag: ["@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And no table view is open", () => noTableViewOpen(page));
    await session.step(10, "Then no table should be left in the workspace", () => noTableLeft(page));
    await session.step(13, "Then \"Tables\" dock panel should be absent", () => shouldBe(page, el("\"Tables\" dock panel"), "absent"));
    await session.step(14, "When user presses Alt+T", () => pressKey(page, "Alt+T"));
    await session.step(15, "Then \"Tables\" dock panel should be visible", () => shouldBe(page, el("\"Tables\" dock panel"), "visible"));
    await session.step(16, "And the \"rows\" reading of Grid viewer in \"Tables\" dock panel should be 0", () => readingIs(page, "rows", el("Grid viewer in \"Tables\" dock panel"), 0));
    await session.step(17, "And no errors should have been logged", () => noErrors(page));
    await session.step(18, "When user presses Alt+T", () => pressKey(page, "Alt+T"));
    await session.step(19, "Then \"Tables\" dock panel should be absent", () => shouldBe(page, el("\"Tables\" dock panel"), "absent"));
    await session.step(20, "And no errors should have been logged", () => noErrors(page));
    await session.step(21, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
