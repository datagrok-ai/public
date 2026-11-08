# Scripts — what the features leave out

## Still open

- core#456 — closing the script editor with an unsaved edit asks nothing; claim when done: the confirmation on
  close (`scripts-edit`).

## Won't do

- Interrupting a running script: JavaScript scripts end at once, and server-language scripts are not run in features.
- The run dialog's input checks (empty, out of range, wrong type, a validator): ApiTests `Widgets: InputForm validation
  and events`, `Widgets: InputForm skipLogic`; an empty required Table input keeps OK disabled (core#1155): xamgle
  `FuncParamEditor | DataFrameParam` (a feature cannot know which tables earlier features left in the page).
- Header parsing (extra whitespace, a malformed `input:` line): parser tests, not features.
