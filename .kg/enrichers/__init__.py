"""LLM-driven enrichers run AFTER deterministic extractors.

Each enricher has a 3-step lifecycle:
  1. prepare()  — gather context per unit (e.g. per package), write input JSON
  2. (external) — agent/LLM consumes input, writes output JSON
  3. apply()    — validate output, emit JSONL nodes/edges into the canonical store

Step 2 is intentionally outside Python — it can be Claude Code agents,
the Anthropic SDK, or batch jobs. The contract is the input/output JSON shape.
"""
