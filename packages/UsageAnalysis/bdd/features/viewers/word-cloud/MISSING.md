# Word cloud — what the features cannot say yet

## Won't do

- Rotation Step 0 (public#4424): `rotationStep` is declared `min: 1` (core#1200's clamp, ApiTests `AI: Viewers: JsViewer`)
  and the render clamps it too.
- A removed word column (core#1136): Charts `Word cloud` "Removed word column"; a renamed one, no rows, empty strings,
  500 and 501 values, values over 50 characters: ApiTests `AI: Viewers: Shared rules`, Charts `Viewers: no rows`, low
  value.
