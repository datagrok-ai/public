1. Browse navigation tree. Navigate to the file Files/Demo/bio/peptides.csv and double click it in browse tree. 
-  Make sure this element is enabled : #rootDiv > div.layout-workarea-wrapper > div.layout-workarea > div.layout-status-bar > div.d4-flex-row.d4-global-status-panel > div > div:nth-child(1). You can on/off it by clocking on <div class="windows-manager-toggle active"><i class="grok-icon fal fa-window-maximize"></i></div> 
Expected : data should open in view. this means that inside the element '<div class="tab-handle-list-container" style=""><i class="grok-icon far fa-search grok-smart-bar-icon" name="icon-search"></i><input type="text" class="ui-input-editor grok-smart-bar-input" placeholder="Alt+Q to search commands" style="visibility: hidden;"><i class="grok-icon fal fa-chevron-down" name="icon-chevron-down" style="display: none;"></i><div class="tab-handle disable-selection" name="view-handle: Home" data-home="true"><div class="tab-handle-text"><i class="grok-icon fal fa-home d4-view-icon" name="icon-home"></i></div><div class="tab-handle-close-button"><i class="fal fa-times"></i></div></div><div class="tab-handle disable-selection tab-handle-selected" name="view-handle: peptides"><div class="tab-handle-text"><span>peptides</span></div><div class="tab-handle-close-button"><i class="fal fa-times"></i></div></div></div>' should be present tab-handle-text = peptides. 
2. Make sure that each amino acid is rendered with a different color. 
This means than letters visible for column AlignedSequence are visually colored in different colors. 
3. Make sure context panel is opened. Click on the AlignedSequence column title.
4. Check that in **Context Panel** present these elements:
- #elementContent > div > div > div > div:nth-child(11) > div.d4-accordion-pane-header.expanded
- <div class="d4-accordion-pane-header expanded" name="div-section--Details">Details</div>
5. On the **Context Panel**, expand each tab d4-accordion-pane-header should become expanded for tab you are checking. If a tab contains some not expanded parts - expand them all where possible. 
6. Make sure there are no Error messages inside expanded tabs in context panel.
