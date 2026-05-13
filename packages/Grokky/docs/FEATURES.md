## Smart Assistant / NLP-driven UI

- Natural-language query interface: Allow users to ask “Which compounds have high toxicity risk and low solubility?” or
  “Show me molecules with IC50 < 10 nM and Lipophilicity < 3” and translate into filters, visualizations, models.
- Search integration (help, etc)
    - Classification of search string (id, entity, question, request, etc)
    - Using LLMs to match the search-integrated query templates
- Querying databases (sort of already implemented but can be improved)
    - Automatic collection of “useful” words
    - Collecting statistics on the db
    - Mapping query parameters to the database columns
- Invoking functions and chains of functions
    - Filtering out available functions
- Contextual suggestions: Based on the semantic types detected (molecule, gene sequence, plate ID, etc), provide
  suggestions: “Would you like to run an ADMET prediction model?” or “Try clustering these compounds on similarity”
  automatically.
- Auto-narrative insights: For dashboards or result sets, auto-generate a short narrative (“In the hit list you see two
  structural clusters; cluster A shows significantly more low-toxicity compounds than cluster B”) to aid non-technical
  users.
- AI for chemists for managing molecules (“strip salts”), etc
- Working with dataframes (average height by race)
- Voice-driven chats

## Self-guided database access

There are multiple ways to access data in Datagrok, and multiple ways where AI can be useful. A lot depends on how
much power and uncertainty you want to give users, and how much are you willing to invest to annotate your existing
databases with metadata to help guide the LLMs.

- Manually written parameterized SQL, annotated with description and search patterns to help LLMs associate it with the
  user query and pass parameters, potentially using transformation functions to the parameters (like resolving compound
  ids and drug names to SMILES, and then passing to the pre-defined SQL query)
- Visual editor - optimized UI for creating joins, filters, aggregations, pivots etc visually, with the UI suggesting
  sensible options based on the knowledge of the schema (column names and types, foreign keys, etc). In the end, a visual
  query produces a parameterized SQL.
- AI-driven query building and ad-hoc querying. You point to a database schema and write a query in English, and LLM
  produces the SQL query, along with the visualization suitable to best interpret the results.

All of the approaches outlined above integrate into the AI-driven search engine that understands the nature of the
question, and routes it accordingly.

## Help
- Answer questions using our wiki (just like DeepWiki). Fully integrated into the search system - the AI recognizes that
  it’s a generic question about the platform, and answers it.
- Automatic troubleshooting
  Get an error -> suggest how to fix (using our knowledge base, etc)

## Vibe-coded apps (on the fly, by users)
- Produce apps based on the user input. Example: “I want a sketcher on the left, and predicted ADME properties appear on
  the right as I sketch the molecule. When I click the SUBMIT button, the molecule should get registered in MolTrack”.
  Amazingly, this is almost working out of the box.

## Agents

Long-performing, collaborative tasks with memory

## MCP

- Connect to external MCP servers, so that Datagrok can make use of them when executing user tasks. Availability of the
  MCP servers is subject to standard Datagrok privileges, so you can easily specify who gets access to what.
- Expose the Datagrok MCP server(s) for others to use

## Ideas

### Different models

- Super-fast old-school, regex-based NLP for obvious things to keep UI snappy
- Client-side, Chrome-integrated (fast and does not make network calls - great for security) Gemini Nano for quick
  triaging/classifying of requests
- ChatGPT for more complex reasoning, function chaining, etc
- Enterprise customers just bring their OpenAI token - your data stays within your infrastructure
