import {Page, expect} from '@playwright/test';

// Balloon observation. The product publishes no balloon event on the JS API side — neither
// AppEvents nor grok.events carries one — so the DOM is the only channel available. A
// MutationObserver is still a subscription: it fires on the mutation itself, never on a timer,
// which is what keeps these helpers on the top rung of the wait ladder.
//
// Classes: the product builds `div('d4-balloon,$type')` (d4/.../widgets/balloon/balloon.dart:28),
// so a balloon carries TWO separate classes and the compound `.d4-balloon.error` is what matches.
// The hyphenated `.d4-balloon-error` and `.grok-balloon-error` exist nowhere in the product; a
// spec reading either counted zero under every product behaviour.
//
// Lifetime differs by who raised the balloon, which decides whether an absence reading can be
// trusted at all. Raised from a plugin through grok.shell.error/warning it auto-hides after ~5 s
// (xamgle/.../interop/grok_api.dart:42-44 passes `autoHide: options['autoHide'] ?? true`); raised
// inside the Dart client through Balloon.error/Balloon.warning it is STICKY (balloon.dart:128,131
// default `autoHide: false`). So an instantaneous count is honest only about the moment it is
// taken — for a claim about a whole window, arm the recorder and read what it collected.

export const BALLOON = '.d4-balloon';
export const ERROR_BALLOON = '.d4-balloon.error';
export const WARNING_BALLOON = '.d4-balloon.warning';

// Every probe balloon this module raises carries this prefix, so a reader can tell its own probe
// from product output without the caller having to remember to exclude it.
const PROBE_MARKER = 'channel-probe ';

export interface RecordedBalloon {
  cls: string;
  text: string;
}

// Records every balloon that enters the document from this point on. Re-arming discards what the
// previous window collected, so each call scopes the next reading to what follows it.
export async function armBalloonRecorder(page: Page): Promise<void> {
  await page.evaluate(() => {
    const w = window as any;
    if (w.__balloonObs) w.__balloonObs.disconnect();
    w.__balloonSeen = [];
    // The identity ledger deliberately SURVIVES re-arming, while the window list above does not.
    // Balloon.show re-appends the container to body on every show (balloon.dart:26), so every show
    // re-delivers the container — and with it every balloon still standing inside it. Clearing the
    // ledger on re-arm made those older balloons recordable again, and the next window attributed
    // them to whatever action it was judging (measured 2026-08-23: armed with one sticky balloon
    // standing, raised one, read back two). A sticky balloon is the realistic case, since Dart-side
    // Balloon.error/warning default to autoHide: false.
    w.__balloonNodes = w.__balloonNodes ?? [];
    const record = (el: Element): void => {
      // querySelectorAll, not querySelector: a re-appended container holds every standing balloon,
      // and taking only the first would drop a genuinely new one arriving alongside an old.
      const found = el.classList?.contains('d4-balloon') ? [el] : Array.from(el.querySelectorAll?.('.d4-balloon') ?? []);
      for (const b of found) {
        if (w.__balloonNodes.indexOf(b) >= 0) continue;
        w.__balloonNodes.push(b);
        // textContent, not innerText: a balloon that already left the DOM has no layout.
        w.__balloonSeen.push({cls: String(b.className), text: String(b.textContent || '').slice(0, 300)});
      }
    };
    w.__balloonObs = new MutationObserver((muts: MutationRecord[]) => {
      for (const m of muts)
        for (const n of Array.from(m.addedNodes))
          if (n.nodeType === 1) record(n as Element);
    });
    w.__balloonObs.observe(document.body, {childList: true, subtree: true});
  });
}

export async function readRecordedBalloons(page: Page): Promise<RecordedBalloon[]> {
  return page.evaluate(() => ((window as any).__balloonSeen ?? []) as RecordedBalloon[]);
}

// Raises a balloon carrying a token unique to this call, waits for THAT balloon (not merely some
// balloon) to appear, then dismisses it by clicking and waits for its removal. Both waits are
// subscriptions raced against a failure cap, so a green run pays the mechanism's real cost —
// measured 2026-08-23 at 161-256 ms — instead of the 5 s auto-hide the previous form waited out.
//
// Two hazards this guards, both measured on 2026-08-23:
//  - the token match. `.first()` on a shared selector matches whatever was appended first, so an
//    ambient balloon (dev raises `Debugging packages: "..."` at login) could satisfy the probe
//    while the probe's own balloon never rendered.
//  - the parked cursor. A balloon inserted under a stationary pointer receives mouseEnter, and
//    that CANCELS its auto-hide timer (balloon.dart:70-73) — a probe balloon left standing >12 s
//    was observed this way. The pointer is moved off the balloon corner after the click.
// The two caps are deliberately different. A JS-API balloon auto-hides at exactly 5 000 ms
// (grok_api.dart:44 passes `timeout ?? 5`), so a removal cap of 5 s would race the auto-hide: a
// click that MISSED would still see the node leave, and the dismissal assertion would pass for the
// wrong reason. The removal cap stays far under that, which keeps it a real discriminator — the
// mechanism itself was measured at 142-160 ms (the product's own 100 ms closeOnMouseUp delay plus
// observer turnaround), so 1.5 s is roughly ten times the observed cost.
const DISMISS_CAP_MS = 1_500;

export async function proveBalloonChannel(page: Page, label: string, cap = 5_000): Promise<void> {
  const token = `${PROBE_MARKER}${label} ${Date.now()}`;
  const appeared = await page.evaluate(async ({tok, ms, dismissMs}) => {
    const w = window as any;
    const hit = (): Element | undefined => Array.from(document.querySelectorAll('.d4-balloon.warning'))
      .find((b) => (b.textContent || '').indexOf(tok) >= 0);
    // Armed BEFORE the raise, so a balloon that renders synchronously cannot slip past the observer.
    const seen = new Promise<Element | null>((resolve) => {
      const already = hit();
      if (already) {
        resolve(already);
        return;
      }
      const obs = new MutationObserver(() => {
        const el = hit();
        if (el) {
          obs.disconnect();
          resolve(el);
        }
      });
      obs.observe(document.body, {childList: true, subtree: true});
      setTimeout(() => {
        obs.disconnect();
        resolve(null);
      }, ms);
    });
    w.grok.shell.warning(tok);
    const node = await seen;
    if (node === null) return false;
    w.__probeGone = new Promise<boolean>((resolve) => {
      if (!node.isConnected) {
        resolve(true);
        return;
      }
      const obs = new MutationObserver(() => {
        if (!node.isConnected) {
          obs.disconnect();
          resolve(true);
        }
      });
      obs.observe(document.body, {childList: true, subtree: true});
      setTimeout(() => {
        obs.disconnect();
        resolve(false);
      }, dismissMs);
    });
    return true;
  }, {tok: token, ms: cap, dismissMs: DISMISS_CAP_MS});
  expect(appeared,
    `the balloon channel must render "${token}"; if it does not, every "no balloon" reading in this spec `
    + 'is a reading of a channel nothing has been shown to arrive on').toBe(true);

  // Clicked, not waited out: closeOnMouseUp removes the node ~100 ms after mouseUp with no fade
  // (balloon.dart:83-94). x=8 keeps clear of the copy and close icons, measured at x=245 and x=264
  // in a 281 px balloon; a warning carries no Report icon at all (it is built only for type
  // 'error', balloon.dart:37-43).
  await page.locator(WARNING_BALLOON, {hasText: token}).first().click({position: {x: 8, y: 8}});
  await page.mouse.move(0, 0);

  const dismissed = await page.evaluate(() => (window as any).__probeGone as Promise<boolean>);
  expect(dismissed,
    `the probe balloon "${token}" must leave the DOM when clicked; still standing means the click missed `
    + 'its body, and the cancelled auto-hide timer will leave it there for the rest of the run').toBe(true);
}

// Instantaneous — never a retrying matcher. `expect(locator).toHaveCount(0)` polls for up to the
// configured expect timeout (15 s, playwright.config.ts) and goes green the moment a balloon
// auto-hides, so a plugin-raised error standing at read time would be waited out and reported as
// absence. Ambient balloons are excluded BY NAME, never by shortening the window.
export async function expectNoErrorBalloon(page: Page, what: string, ignore: RegExp[] = []): Promise<void> {
  const found = await page.evaluate((sel) => Array.from(document.querySelectorAll(sel))
    .map((b) => String(b.textContent || '').slice(0, 300)), ERROR_BALLOON);
  const kept = found.filter((t) => !ignore.some((re) => re.test(t)));
  expect(kept, `${what}: no error balloon may be standing; found ${JSON.stringify(found)}`).toEqual([]);
}

// The window claim rather than the instant claim: what the recorder collected since it was armed.
// Use this when the action being judged runs long enough that a plugin-raised balloon could have
// auto-hidden before an instantaneous read.
export async function expectNoBalloonSinceArmed(
  page: Page, what: string, cls: RegExp = /error/, ignore: RegExp[] = [],
): Promise<void> {
  const all = await readRecordedBalloons(page);
  const kept = all.filter((b) => cls.test(b.cls)
    && b.text.indexOf(PROBE_MARKER) !== 0
    && !ignore.some((re) => re.test(b.text)));
  expect(kept, `${what}: no balloon matching ${cls} may have been raised since the recorder was armed; `
    + `recorded ${JSON.stringify(all)}`).toEqual([]);
}

// Arms the recorder and then proves, INSIDE the window it just opened, that this recorder instance
// is live — by raising a probe and requiring the recorder to have caught it.
//
// The distinction from calling `armBalloonRecorder` and `proveBalloonChannel` separately is the
// whole point. `armBalloonRecorder` disconnects the previous observer before installing a new one,
// so an install that fails leaves no observer at all and every later read returns an empty list
// under any product behaviour. A probe raised BEFORE the arm proves the channel but says nothing
// about the observer the reading will actually come from; a neighbouring step's successful reading
// proves a DIFFERENT observer, installed by a different call. Only a probe inside this window
// makes a later empty reading mean "nothing was raised" rather than "nothing was watching".
//
// The probe's own balloon stays in the recorded list and is excluded from absence assertions by its
// marker, so callers need not remember to filter it.
export async function armBalloonRecorderProved(page: Page, label: string): Promise<void> {
  await armBalloonRecorder(page);
  await proveBalloonChannel(page, label);
  const recorded = await readRecordedBalloons(page);
  const probes = recorded.filter((b) => b.text.indexOf(PROBE_MARKER) === 0);
  expect(probes.length, `the balloon recorder armed for "${label}" did not record its own probe, so an `
    + `empty reading from it would prove nothing; recorded ${JSON.stringify(recorded)}`).toBeGreaterThan(0);
}
