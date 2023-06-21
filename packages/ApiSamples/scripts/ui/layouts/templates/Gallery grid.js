//Layout: Gallery grid

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showContextPanel = false; //*hide context panel

let gallery = ui.divH([],'grok-gallery-grid');
for (let i =0; i<40; i++){
  gallery.append(ui.card(ui.p('item-'+i.toString())))
}

let view = grok.shell.newView('Gallery grid', [
  gallery
],'layout');

view.box = true;

$('.layout div').css('border','1px dashed var(--grey-2)');