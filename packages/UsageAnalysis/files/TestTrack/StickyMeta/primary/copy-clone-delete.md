
### Copy / clone / move objects with metadata	
1. Open SPGI table (with some sticky meta added by you).
2. Verify metadata presence for:
   - Cloned table
   - New view opened
   - Saved as a project and reopened
   - Moved project to some Space and opened from there
   - Imported/exported
- Expected Result: 
   - No metadata is lost
   - Copy/cloned object shows same metadata; moving between Spaces preserves metadata	

### Delete metadata and verify removal	
1. On a cell with metadata, open Context Panel → Sticky Meta
2. Delete fields rating and notes. Save
3. Refresh page and/or reopen table	
- Expected Result: 
   - Metadata values are deleted; no blue marker appears if there is no more sticky meta for the cell	
   - After reload, metadata is absent	

### Persistency after refresh, relogin, server restart	
1. Add metadata to some objects. 
2. Perform:
• Refresh browser tab
• Logout & login again
3. Check metadata on same objects. 	
- Expected Result: 
   - Metadata remains intact in all scenarios.	
   - No metadata is lost.


---
{
  "order": 3
}
