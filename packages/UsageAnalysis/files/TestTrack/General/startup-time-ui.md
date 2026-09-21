---
feature: general
realizes_atlas: []
priority: p2
target_layer: manual-only
coverage_type: perf
manual_only_reason: |
  This is a performance threshold check, not a functional test. The 2-4s budget
  is environment-dependent (cold vs warm cache, CI worker contention, network to
  the API), and "ready for user interaction" is not a crisp, assertable signal —
  encoding a hard 2-4s bound in a Playwright spec would be flaky and would fail
  for reasons unrelated to a real regression. Belongs to perf monitoring / a
  dedicated startup-time metric, not the functional regression suite.
related_bugs: []
---

# Platform startup time

Verifies that the platform starts and becomes ready for user interaction
within an acceptable time frame (2-4 seconds) after a fresh, cold-cache
launch.

- Close all running platform instances. 
- Start the platform from a fresh launch (cleared cache, new browser tab).
- Measure the time from the launch action to the point where the platform is ready for user interaction.
- **Expected Result**: The platform should load and become interactive within **2–4 seconds**.

Notes:
- Ignore minor loading of non-blocking background processes.

---
{
  "order": 9,
  "datasets": []
}
