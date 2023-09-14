//Datagrok has pre-build CSS helper classes to manipulate most common parameters such as:
//margins, paddings, text size, text weight, etc.
//To call them, use ui.css. 

let div1 = ui.div('Child 1');

div1.classList.add(
  ui.css.padding.large,    // Padding
  ui.css.border.all,       // Border
  ui.css.background.white, // Background color
  ui.css.margin.medium,    // Margin
  ui.css.textSize.medium,  // Text size
  ui.css.textWeight.bold,  // Text weight
  ui.css.shadow.small      // Shadow
);

ui.div([
  'Parent', 
   div1
], ui.css.background.light); //add light grey background color