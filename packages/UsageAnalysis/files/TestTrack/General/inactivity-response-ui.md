---
feature: general
target_layer: manual-only
coverage_type: regression
produced_from: split
original_path: public/packages/UsageAnalysis/files/TestTrack/General/inactivity-response.md
split_date: 2026-06-16
related_bugs: []
manual_only_reason: |
  The premise is a real 20-minute idle period before resuming interaction. A
  literal 20-minute wait is impractical in CI, and faking the server-side
  session/idle timeout (token refresh, websocket reconnect, auth expiry) cannot
  be reliably reproduced from the client without exercising the exact wall-clock
  path the scenario is meant to validate. The post-idle interaction steps
  (navigate, open file, run query/script, watch console) are themselves
  automatable, but are meaningless without the genuine idle precondition. Keep
  as manual (or a long nightly soak), not a standard Playwright spec.
---

# Platform behavior after 20 minutes of inactivity

Verifies that the platform stays stable and responsive after being left idle
for 20 minutes, and that resuming interaction (navigating, opening a file,
running a query and a script) does not produce console errors.

1. Open platform. (User is logged in). Begin a timer to track 20 minutes of inactivity on the platform.
2. Do not interact. Leave the platform idle without any user interaction for the duration of the 20-minute timer.
3. Post-onactivity interaction. After 20 minutes, interact with the platform by performing the following actions:
  * Navigate to any section (e.g., Dashboards, Databases, or Plugins). Review console outputs for errors that might have occurred during the period of inactivity or during the first interaction afterward.
  * Open a data file or project.
  * Run a query and a script.
4. Observe if the platform is responsive. Verify that there are no unexpected delays or errors. 
5. Expected result:
* The platform remains stable and responsive after 20 minutes of inactivity.
* No errors or unexpected behavior should occur when resuming interaction with the platform.
* Console outputs should not show any errors related to the inactivity period.

---
{
  "order": 3
}