/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/zz-scratch/probe.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {addTileViewer} from '../../bindings/tile-viewer.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, pointerAway} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("scratch probe", () => {
  const session = feature(test, "features/zz-scratch/probe.feature", import.meta.url);
  test("tile label tooltip", {tag: ["@viewers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(4, "Given user is logged in", () => loggedIn(page));
    await session.step(5, "And user opens a table \"long-values\" with:", () => openTableOf(page, "long-values", [["id","note","measured_concentration_of_compound_in_plasma"],["1","The first value of this column is long enough that no field of a card can show it whole, so its end is cut off here.","5"],["2","short","6"]]), [["id","note","measured_concentration_of_compound_in_plasma"],["1","The first value of this column is long enough that no field of a card can show it whole, so its end is cut off here.","5"],["2","short","6"]]);
    await session.step(9, "And user adds a tile viewer", () => addTileViewer(page));
    await session.step(10, "When user hovers over the \"label measured_concentration_of_compound_in_plasma of row 1\" area of tile viewer", () => hoverArea(page, "label measured_concentration_of_compound_in_plasma of row 1", el("tile viewer")));
    await session.step(11, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
    await session.step(12, "And tooltip should contain text \"measured_concentration_of_compound_in_plasma\"", () => shouldContainText(page, el("tooltip"), "measured_concentration_of_compound_in_plasma"));
    await session.step(13, "When user moves the pointer away from tile viewer", () => pointerAway(page, el("tile viewer")));
    await session.step(14, "And user hovers over the \"label note of row 1\" area of tile viewer", () => hoverArea(page, "label note of row 1", el("tile viewer")));
    await session.step(15, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
    await session.step(16, "When user hovers over the \"field note of row 1\" area of tile viewer", () => hoverArea(page, "field note of row 1", el("tile viewer")));
    await session.step(17, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
  });
});
