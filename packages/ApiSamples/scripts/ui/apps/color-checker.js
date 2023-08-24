//name: color-checker
//description: color contrast and blindness check
//language: javascript

let view = grok.shell.newView('HexToRGB');
view.box = true;

let contrastBox = ui.div();
contrastBox.style.padding = '5px';
contrastBox.style.maxWidth = '150px';
contrastBox.style.marginBottom = '10px';

let hexActive = ui.boolInput('', true, v=>{
  if (v)
    colorHex.enabled = true;
  else
    colorHex.enabled = false;
  if ((v)&&(rgbActive.value)){
    rgbActive.value = false;
  }
});

let rgbActive = ui.boolInput('', false, v=>{
  if (v)
    colorRgb.enabled = true;
  else
    colorRgb.enabled = false;
  if ((v)&&(hexActive.value)){
    hexActive.value = false;
  } 
});

let colorBox = ui.div();
colorBox.style.marginBottom = '10px';
colorBox.style.padding = '5px';
colorBox.style.maxWidth = '150px';

let colorBlidnessBox = ui.div();
colorBlidnessBox.style.marginBottom = '10px';
colorBlidnessBox.style.maxWidth = '150px';

let inputColor = ui.element('input');
inputColor.type = 'color';
inputColor.style.width = '40px';
inputColor.style.height = '20px';
inputColor.style.marginLeft = '10px';

inputColor.onchange = (e) =>{
  if (colorHex.enabled){
    colorHex.value = inputColor.value.substring(1);
  }else{
    let rgb = hexToRGB(inputColor.value.substring(1));
    colorRgb.value = rgb.R+','+rgb.G+','+rgb.B;
  }
  
}


let colorHex = ui.stringInput('HEX:', '', v=>{
  let result = colorHex.value.match(/[\d\w]{6}/);
  colorBox.innerHTML = '';
  contrastBox.innerHTML = '';
  saveBox.innerHTML = '';
  contrastBox.style.backgroundColor = 'white';
  colorBox.style.backgroundColor = 'white';
  
  if (result!=null){
    let rgb = hexToRGB(colorHex.value);
    if (colorHex.enabled){
      colorRgb.value = rgb.R+','+rgb.G+','+rgb.B;
      calcRatio(rgb);
    }
  }
  if ((colorHex.value =='')&&(colorRgb.value.length>0))
    colorRgb.value = '';
});

let colorRgb = ui.stringInput('RGB:', '', v=>{
  let result = colorRgb.value.match(/([\d]{1,3}\,[\d]{1,3}\,[\d]{1,3}){1,12}/);
  colorBox.innerHTML = '';
  contrastBox.innerHTML = '';
  saveBox.innerHTML = '';
  contrastBox.style.backgroundColor = 'white';
  colorBox.style.backgroundColor = 'white';

  if (result!=null){
	let hex = rgbToHex(result[0]);
  	let rgbresult = colorRgb.value.split(',');
    let rgb = {
      'R': rgbresult[0],
      'G': rgbresult[1],
      'B': rgbresult[2]
    };
    if (colorRgb.enabled){
      colorHex.value = hex;
      calcRatio(rgb);
    }
  }
  if (colorRgb.value == '')
    colorHex.value = '';
});
colorRgb.enabled = false;
colorRgb.input.setAttribute('maxlength','11');

let copyRGB = ui.button(ui.iconFA('copy'), ()=>{
  copyToClipboard(colorRgb.value);
  grok.shell.info(colorRgb.value);
})

let copyHEX = ui.button(ui.iconFA('copy'), ()=>{
  copyToClipboard(colorHex.value);
  grok.shell.info(colorHex.value);
})

view.box = true;
let palette = ui.panel([
  ui.h2('Palette')
]);
let saveBox = ui.div();
view.append(ui.splitH([
  ui.panel([
    ui.divH([
      hexActive,
      colorHex,
      copyHEX
    ]),
    ui.divH([
      rgbActive,
      colorRgb,
      copyRGB
    ]),
    ui.divH([
      ui.div([],{style:{width:'35px'}}),
      ui.label('Picker'),
      inputColor
    ],{style:{alignItems:'center',margin:'10px 0'}}),
    ui.div([
      ui.label('Color:'),
      colorBox,
      ui.label('Contrast:'),
      contrastBox,
      ui.label('Color blindness:'),
      colorBlidnessBox,
      saveBox
    ],{style:{margin:'10px 0 0 35px'}})
  ]),
  palette
]));


//calculation process

function calcRatio(rgb){
  
  let textColor = colorBalance(rgb);
  let flag = '';
  colorBox.innerHTML = 'Text color';
  colorBlidnessBox.innerHTML = '';
  colorBox.style.color = textColor;
  colorBox.style.backgroundColor = '#'+colorHex.value;

  if (textColor == 'black'){
    textColor = [0,0,0];
  }else{
    textColor = [255,255,255];
  }  

  let ratio = contrast(rgb, textColor).toFixed(2);
  let ratioResult = ui.h1(ratio);
  ratioResult.style.margin = '0px';
  if(ratio>7){
    flag = 'Good';
    contrastBox.style.color = 'var(--green-2)'; 
  	contrastBox.style.backgroundColor = 'var(--green-1)';
  }else if ((ratio>4.5)&&(ratio<7)){
    flag = 'Avarage';
    contrastBox.style.color = 'var(--orange-2)'; 
  	contrastBox.style.backgroundColor = 'var(--orange-1)';
  }else if(ratio<4.5){
    flag = 'Bad';
    contrastBox.style.color = 'var(--red-2)'; 
  	contrastBox.style.backgroundColor = 'var(--red-1)';
  }
  
  contrastBox.append(
    ui.div([
      ratioResult,
      ui.divText(flag)
    ]));
    
  let protanopia = ColorMatrix({R:rgb.R,G:rgb.G,B:rgb.B,A:1},Blind('Protanopia'));
  let protanomaly = ColorMatrix({R:rgb.R,G:rgb.G,B:rgb.B,A:1},Blind('Protanomaly'));
  let deuteranopia = ColorMatrix({R:rgb.R,G:rgb.G,B:rgb.B,A:1},Blind('Deuteranopia'));
  let deuteranomaly = ColorMatrix({R:rgb.R,G:rgb.G,B:rgb.B,A:1},Blind('Deuteranomaly'));
  let tritanopia = ColorMatrix({R:rgb.R,G:rgb.G,B:rgb.B,A:1},Blind('Tritanopia'));
  let tritanomaly = ColorMatrix({R:rgb.R,G:rgb.G,B:rgb.B,A:1},Blind('Tritanomaly'));
  let achromatopsia = ColorMatrix({R:rgb.R,G:rgb.G,B:rgb.B,A:1},Blind('Achromatopsia'));
  let achromatomaly = ColorMatrix({R:rgb.R,G:rgb.G,B:rgb.B,A:1},Blind('Achromatomaly'));
  
  colorBlidnessBox.append(
    ui.divH([
      ui.div(['Protanopia ('+contrast(protanopia, textColor).toFixed(2)+')'],{style:{margin:'0 10px'}}),
    ],{style:{minWidth:'200px',color:'rgb('+textColor+')',alignItems:'center',padding:'5px 0', backgroundColor:'rgb('+protanopia.R+','+protanopia.G+','+protanopia.B+')'}}),
    ui.divH([
      ui.div(['Protanomaly ('+contrast(protanomaly, textColor).toFixed(2)+')'],{style:{margin:'0 10px'}}),
    ],{style:{minWidth:'200px',color:'rgb('+textColor+')',alignItems:'center',padding:'5px 0', backgroundColor:'rgb('+protanomaly.R+','+protanomaly.G+','+protanomaly.B+')'}}),
    ui.divH([
      ui.div(['Deuteranopia ('+contrast(deuteranopia, textColor).toFixed(2)+')'],{style:{margin:'0 10px'}}),
    ],{style:{minWidth:'200px',color:'rgb('+textColor+')',alignItems:'center',padding:'5px 0', backgroundColor:'rgb('+deuteranopia.R+','+deuteranopia.G+','+deuteranopia.B+')'}}),
    ui.divH([
      ui.div(['Deuteranomaly ('+contrast(deuteranomaly, textColor).toFixed(2)+')'],{style:{margin:'0 10px'}}),
    ],{style:{minWidth:'200px',color:'rgb('+textColor+')',alignItems:'center',padding:'5px 0', backgroundColor:'rgb('+deuteranomaly.R+','+deuteranomaly.G+','+deuteranomaly.B+')'}}),
    ui.divH([
      ui.div(['Achromatopsia ('+contrast(achromatopsia, textColor).toFixed(2)+')'],{style:{margin:'0 10px'}}),
    ],{style:{minWidth:'200px',color:'rgb('+textColor+')',alignItems:'center',padding:'5px 0', backgroundColor:'rgb('+achromatopsia.R+','+achromatopsia.G+','+achromatopsia.B+')'}}),
    ui.divH([
      ui.div(['Achromatomaly ('+contrast(achromatomaly, textColor).toFixed(2)+')'],{style:{margin:'0 10px'}}),
    ],{style:{minWidth:'200px',color:'rgb('+textColor+')',alignItems:'center',padding:'5px 0', backgroundColor:'rgb('+achromatomaly.R+','+achromatomaly.G+','+achromatomaly.B+')'}}),
  );
  
  saveBox.append(ui.bigButton('Save color',()=>{
    let color = ui.divH([
      ui.div(['#'+colorHex.value],{style:{marginLeft:'10px'}}),
      ui.div([ratio+' - '+flag],{style:{margin:'0 20px'}}),
    ],{style:{alignItems:'center', minWidth:'200px', padding:'5px 0', backgroundColor:'#'+colorHex.value}});
    $(color).attr('data-id',colorHex.value);
    $(color).css('color','rgb('+textColor+')');
    
    palette.append(ui.divH([
      color,
      ui.button(ui.iconFA('copy'),()=>{
          copyToClipboard($(color).children().first().text().substring(1));
  		  grok.shell.info($(color).children().first().text().substring(1));
      }),
      ui.button(ui.iconFA('trash-alt'),()=>{
        $(color).parent().detach();
      })
    ]))
  }))
}

//color converting
function hexToRGB(h) {  
  let aRgbHex = [];
  if (h.length==3){
  	aRgbHex = h.match(/.{1,1}/g);
    aRgbHex[0] = aRgbHex[0].toString().repeat(2);
    aRgbHex[1] = aRgbHex[1].toString().repeat(2);
    aRgbHex[2] = aRgbHex[2].toString().repeat(2);
  }
  else if (h.length==6)
    aRgbHex = h.match(/.{1,2}/g);
  var aRgb = {
    'R': parseInt(aRgbHex[0], 16),
    'G': parseInt(aRgbHex[1], 16),
    'B': parseInt(aRgbHex[2], 16)
  };
  return aRgb
}

function rgbToHex(color) {  
  let rgb = color.split(',')
  let hex = parseInt(rgb[0]).toString(16) + parseInt(rgb[1]).toString(16) + parseInt(rgb[2]).toString(16);
  return hex
}

//colorbalance
function colorBalance(Color){
   let luma = ((0.299 * Color.R) + (0.587 * Color.G) + (0.114 * Color.B)) / 255;
   return luma > 0.5 ? 'black' : 'white';
}

//copy to clipboard
function copyToClipboard(text) {
    var dummy = document.createElement("textarea");
    document.body.appendChild(dummy);
    dummy.value = text;
    dummy.select();
    document.execCommand("copy");
    document.body.removeChild(dummy);
}

//Contrast ratio
function luminance(r, g, b) {
    var a = [r, g, b].map(function (v) {
        v /= 255;
        return v <= 0.03928
            ? v / 12.92
            : Math.pow( (v + 0.055) / 1.055, 2.4 );
    });
    return a[0] * 0.2126 + a[1] * 0.7152 + a[2] * 0.0722;
}
function contrast(rgbBg, rgbText) {
    var lum1 = luminance(rgbBg.R, rgbBg.G, rgbBg.B);
    var lum2 = luminance(rgbText[0], rgbText[1], rgbText[2]);
    var brightest = Math.max(lum1, lum2);
    var darkest = Math.min(lum1, lum2);
    return (brightest + 0.05)
         / (darkest + 0.05);
}

//*****************
function ColorMatrix(o,m) { // takes: RGBA object, Matrix array

    function fu(n) { return(n<0?0:(n<255?n:255)); }

    var r=((o.R*m[0])+(o.G*m[1])+(o.B*m[2])+(o.A*m[3])+m[4]);
    var g=((o.R*m[5])+(o.G*m[6])+(o.B*m[7])+(o.A*m[8])+m[9]);
    var b=((o.R*m[10])+(o.G*m[11])+(o.B*m[12])+(o.A*m[13])+m[14]);
    var a=((o.R*m[15])+(o.G*m[16])+(o.B*m[17])+(o.A*m[18])+m[19]);
    
    return({'R':fu(r), 'G':fu(g), 'B':fu(b), 'A':fu(a)});
    
};

function Blind (v) { // this function just returns the Matrix

    return({'Normal':[1,0,0,0,0, 0,1,0,0,0, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,1],
            'Protanopia':[0.567,0.433,0,0,0, 0.558,0.442,0,0,0, 0,0.242,0.758,0,0, 0,0,0,1,0, 0,0,0,0,1],
            'Protanomaly':[0.817,0.183,0,0,0, 0.333,0.667,0,0,0, 0,0.125,0.875,0,0, 0,0,0,1,0, 0,0,0,0,1],
            'Deuteranopia':[0.625,0.375,0,0,0, 0.7,0.3,0,0,0, 0,0.3,0.7,0,0, 0,0,0,1,0, 0,0,0,0,1],
            'Deuteranomaly':[0.8,0.2,0,0,0, 0.258,0.742,0,0,0, 0,0.142,0.858,0,0, 0,0,0,1,0, 0,0,0,0,1],
            'Tritanopia':[0.95,0.05,0,0,0, 0,0.433,0.567,0,0, 0,0.475,0.525,0,0, 0,0,0,1,0, 0,0,0,0,1],
            'Tritanomaly':[0.967,0.033,0,0,0, 0,0.733,0.267,0,0, 0,0.183,0.817,0,0, 0,0,0,1,0, 0,0,0,0,1],
            'Achromatopsia':[0.299,0.587,0.114,0,0, 0.299,0.587,0.114,0,0, 0.299,0.587,0.114,0,0, 0,0,0,1,0, 0,0,0,0,1],
            'Achromatomaly':[0.618,0.320,0.062,0,0, 0.163,0.775,0.062,0,0, 0.163,0.320,0.516,0,0,0,0,0,1,0,0,0,0,0]}[v]);

};