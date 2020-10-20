<!-- TITLE: Tests: Elastic search -->
<!-- SUBTITLE: -->

# Tests: Elastic search

Elastic search allows search for functions and other main menu items, as well as looking for search query inside help pages.

## Testing scenario

1. Open elastic search field (last item of main menu or use "Alt + Q")
   * Search field expanded
   
1. Enter *"Open"* to search field
   * In drop-down list the main menu actions is offered (Open Local File, Open Web File, Open Test Dataset, Open Layout Gallery)
   * Drop-down list contains help pages, including test scenarios
   * For menu actions that have icon, same icon is displayed in drop-down list
   * Help pages marked withs "?"
   
1. Click on "Open Local File" from drop-down list (Or go to it using ↓ and click "Enter")
   * An explorer opens to select a local file
   * Search field cleared

1. Close explorer ("Open Local File" dialog)   

1. Open elastic search field (last item of main menu or use "Alt + Q")
   * Previous search keys entered in search field (*"Open"*)

1. Clear the search field

1. Enter *"Join"* to search field
   * In  drop-down list there is main menu item "Join Tables", which is not active for pressing (not enough tables open)
  
1. Click on "Join Tables" with icon "?" from drop-down list (Or go to it using ↓ and click "Enter")
   * "Join Tables" is open in context help

1. Clear the search field

1. Enter *"eq"* to search field
   * In  drop-down list there is actions in name of which is "eq"
   * Actions in list are marked with "lightning" icon
  
1. Click on "eq" from drop-down list (Or go to it using ↓ and click "Enter")
   * "eq" action is open on [Property Panel](../overview/navigation.md#properties)
