/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/forms/forms-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.forms]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnType, makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, noErrors, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Forms viewer — regressions of closed issues", () => {
  const session = feature(test, "features/viewers/forms/forms-regressions.feature", import.meta.url);
  test("An empty integer cell shows an empty field on the card (public#3561)", {tag: ["@viewers", "@realizes:viewers.forms"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens a table \"ints\" with:", () => openTableOf(page, "ints", [["id","n"],["1","5"],["2",""],["3","7"]]), [["id","n"],["1","5"],["2",""],["3","7"]]);
    await session.step(18, "Then \"n\" column should have type \"int\"", () => columnType(page, "n", "int"));
    await session.step(19, "When user adds a forms viewer", () => addViewer(page, "forms"));
    await session.step(20, "Then forms viewer should be visible", () => shouldBe(page, el("forms viewer"), "visible"));
    await session.step(21, "When user makes row 1 current", () => makeRowCurrent(page, 1));
    await session.step(22, "Then the \"n of current card\" reading of forms viewer should be \"5\"", () => readingReads(page, "n of current card", el("forms viewer"), "5"));
    await session.step(23, "When user makes row 2 current", () => makeRowCurrent(page, 2));
    await session.step(24, "Then the \"current record\" reading of forms viewer should be 2", () => readingIs(page, "current record", el("forms viewer"), 2));
    await session.step(25, "And the \"id of current card\" reading of forms viewer should be \"2\"", () => readingReads(page, "id of current card", el("forms viewer"), "2"));
    await session.step(26, "And the \"n of current card\" reading of forms viewer should be \"\"", () => readingReads(page, "n of current card", el("forms viewer"), ""));
    await session.step(27, "When user makes row 3 current", () => makeRowCurrent(page, 3));
    await session.step(28, "Then the \"n of current card\" reading of forms viewer should be \"7\"", () => readingReads(page, "n of current card", el("forms viewer"), "7"));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
  });
});
