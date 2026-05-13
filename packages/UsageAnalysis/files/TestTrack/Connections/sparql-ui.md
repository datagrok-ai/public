# SPARQL — manual UI checks

This is the **manual companion** to `sparql.md`. The autotest
(`sparql.test.ts`) covers most of the scenario end-to-end via the global
**New Connection** button (which surfaces a Data Source dropdown that lists
every provider, including Sparql). The carve-out is the bit the autotest
deliberately skips — the manual scenario calls for clicking the `...` ("Show
more") footer under Browse > Databases to surface the Sparql provider as a
right-clickable tree node, and then choosing **Add new connection** from its
context menu. On dev that footer behavior is brittle and the dialog the
button surfaces is identical to the New Connection one, so the autotest
takes the New-Connection path. The Test button is also absent on the Sparql
dialog on dev — the autotest goes straight to OK.

The autotest covers:

- Open the New Connection dialog, switch Data Source to **Sparql**,
  fill **Name** = `test_sparql`, **Endpoint**, **Requires Server** = true,
  leave Prefixes empty
- Click **OK** → connection saved with `dataSource = "Sparql"`
- Right-click `test_sparql` in the Browse tree → **Delete...** → DELETE
- Server-side `dapi.connections.filter(...)` returns nothing afterwards

## Pre-conditions

- Logged into Datagrok dev/release; no `test_sparql` connection from a previous run

## Steps — `Show more` discoverability

1. Open **Browse > Databases**
2. Click the `...` footer at the bottom of the Databases group
3. Verify the list expands to show every provider, including **Sparql**, even
   though it has no saved connection
4. Right-click **Sparql** → **Add new connection** opens the same dialog the
   global New Connection button does
5. (The rest of the flow — fill / OK / right-click delete — is what the
   autotest exercises end-to-end.)

## Steps — TEST button on Sparql dialog (release/cloud)

If the Sparql dialog has a **TEST** button on your env (it doesn't on dev):

1. Fill in the connection fields per `sparql.md`
2. Click **TEST**
3. Wait for a balloon — green for success, red with the SPARQL endpoint error otherwise

## What to look for

- `...` footer click reliably reveals every provider with no error balloons
- Sparql provider's right-click menu has **Add new connection / Add connection...**
- Sparql Add-Connection dialog carries Endpoint / Requires Server / Prefixes inputs
- After OK, the new connection's tree node appears under Sparql immediately
