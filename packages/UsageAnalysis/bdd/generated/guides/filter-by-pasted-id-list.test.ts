/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/filter-by-pasted-id-list.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {pasteIntoCardSearch} from '../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterIsExactlyAnyOf, filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter a table down to a list of IDs", () => {
  const session = feature(test, "features/guides/filter-by-pasted-id-list.feature", import.meta.url);
  test("Paste a list of IDs into the search of the Id filter card", {tag: ["@guide", "@help:visualize/viewers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(15, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
    await session.step(16, "And user adds a card for \"Id\" to the filter panel", () => addCardFor(page, "Id"));
    await session.step(17, "And user hovers over \"Id\" filter card", () => hoverOver(page, el("\"Id\" filter card")));
    await session.step(18, "And user clicks on search icon of \"Id\" filter card", () => clickOn(page, el("search icon of \"Id\" filter card")));
    await session.step(19, "And user pastes \"CAST-634783\\nCAST-634790\\nCAST-634812\\nCAST-634851\\nCAST-634880\" into the search of the \"Id\" filter card", () => pasteIntoCardSearch(page, "CAST-634783\\nCAST-634790\\nCAST-634812\\nCAST-634851\\nCAST-634880", "Id"));
    await session.step(20, "Then 5 rows should pass the filter", () => filterPasses(page, 5));
    await session.step(21, "And the filter should pass exactly the rows where \"Id\" is one of \"CAST-634783, CAST-634790, CAST-634812, CAST-634851, CAST-634880\"", () => filterIsExactlyAnyOf(page, "Id", "CAST-634783, CAST-634790, CAST-634812, CAST-634851, CAST-634880"));
    await session.step(22, "And grid should show 5 rows", () => showsRows(page, el("grid"), 5));
  });
});
