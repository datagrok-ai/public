/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/general/table-manager.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Table Manager is docked and closed with Alt+T", () => {
  const session = feature(test, "features/general/table-manager.feature", import.meta.url);
  test("Alt+T docks the Table Manager and closes it again", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(21, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await session.step(22, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(23, "Then the open table views should be exactly \"demog, smiles, spgi-100\"", () => openTableViewsExactly(page, "demog, smiles, spgi-100"));
    await session.step(26, "Then \"Tables\" dock panel should be absent", () => shouldBe(page, el("\"Tables\" dock panel"), "absent"));
    await session.step(27, "When user presses Alt+T", () => pressKey(page, "Alt+T"));
    await session.step(28, "Then \"Tables\" dock panel should be visible", () => shouldBe(page, el("\"Tables\" dock panel"), "visible"));
    await session.step(29, "And Grid viewer in \"Tables\" dock panel should be visible", () => shouldBe(page, el("Grid viewer in \"Tables\" dock panel"), "visible"));
    await session.step(30, "When user presses Alt+T", () => pressKey(page, "Alt+T"));
    await session.step(31, "Then \"Tables\" dock panel should be absent", () => shouldBe(page, el("\"Tables\" dock panel"), "absent"));
    await session.step(32, "And no errors should have been logged", () => noErrors(page));
  });
});
