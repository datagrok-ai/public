/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';

export class IconTool extends DG.ViewBase {
  constructor(name) {
    super(name);

    this.name = name;

    const fileUrl = ui.input.string('', {value: ''});
    fileUrl.readOnly = true;
    fileUrl.enabled = false;
    //@ts-ignore
    fileUrl.input.placeholder = '';

    const browseBtn = ui.button(ui.iconFA('folder-open'), ()=> $(browseFile).trigger('click'));
    const clearImage = ui.button(ui.iconFA('trash-alt'), ()=> deleteImageFile());
    const fileButtons = ui.buttonsInput([
      //@ts-ignore
      ui.divH([
        //@ts-ignore
        fileUrl.input,
        browseBtn,
        clearImage,
      ]),
    ]);

    $(fileButtons).css('margin', '0');
    $(fileButtons).children('label').text('Add image');

    const blendMode = ui.input.choice('Blend mode', {value: 'normal', items: ['normal', 'multiply', 'darken', 'lighter', 'screen', 'overlay', 'color-dodge', 'color-burn', 'hard-light', 'soft-light', 'difference', 'exclusion', 'hue', 'color', 'luminosity'], onValueChanged: (value) => setBlendMode(value)});
    const text = ui.input.string('Add text', {value: '', onValueChanged: (value) => setText(value, fontColorPicker.value, fontSize.value, fontWeight.value)});
    const fontSize = ui.input.int('Text size', {value: 60, onValueChanged: (value) => setText(text.value, fontColorPicker.value, value, fontWeight.value)});
    $(fontSize.input).attr('type', 'number');
    const fontWeight = ui.input.choice('Font weight', {value: 'normal', items: ['lighter', 'normal', 'bold'], onValueChanged: (value) => setText(text.value, fontColorPicker.value, fontSize.value, value)});
    fontWeight.input.style.width='100%';

    const textX = ui.element('input', 'textX');
    $(textX).attr('type', 'range');
    $(textX).attr('min', '-100');
    $(textX).attr('max', '100');
    $(textX).attr('value', '00');
    $(textX).css('margin-left', '5px');

    const textY = ui.element('input', 'textY');
    $(textY).attr('type', 'range');
    $(textY).attr('min', '-100');
    $(textY).attr('max', '100');
    $(textY).attr('value', '00');
    $(textY).css('margin-left', '5px');
    //@ts-ignore
    const imageMargin = ui.input.int('Margins', {value: 0,
      onValueChanged: (value) => setImage(image, imageRange.value, value)});
    $(imageMargin.input).attr('type', 'number');
    const imageRange = ui.element('input', 'imageRange');
    $(imageRange).attr('type', 'range');
    $(imageRange).attr('min', '0');
    $(imageRange).attr('max', '100');
    $(imageRange).attr('value', '100');

    const removeColorPicker = ui.element('input', 'colorPicker');
    const bgColorPicker = ui.element('input', 'colorPicker');
    const fontColorPicker = ui.element('input', 'colorPicker');
    $(bgColorPicker).attr('type', 'color');
    $(removeColorPicker).attr('type', 'color');
    $(removeColorPicker).css('visibility', 'hidden');
    $(removeColorPicker).css('margin-left', '-20px');
    $(fontColorPicker).attr('type', 'color');
    $(bgColorPicker).val('#efefef');

    const infoIcon = ui.iconFA('info-circle');
    infoIcon.style.color='var(--blue-1)';
    infoIcon.style.marginLeft='10px';

    const removeColorBtn = ui.button(ui.iconFA('eye-dropper'), ()=>$(removeColorPicker).trigger('click'));
    removeColorBtn.style.backgroundColor='var(--steel-1)';
    removeColorBtn.style.border='1px solid var(--steel-2)';
    removeColorBtn.style.width='28px';
    removeColorBtn.style.borderRadius='100%';
    $(removeColorBtn).children('i').removeClass('fal');
    $(removeColorBtn).children('i').addClass('fas');

    const removeColor = ui.buttonsInput([
      removeColorBtn,
      //@ts-ignore
      ui.tooltip.bind(infoIcon, 'Set color or use the eyedropper to remove the color from the image'),
      //@ts-ignore
      removeColorPicker,
    ]);
    $(removeColor).children('label').text('Remove color');
    $(removeColor).children('label').css('padding-top', '4px');
    $(removeColor).css('margin', '0');

    const setColors = ui.buttonsInput([
      //@ts-ignore
      ui.divH([
        ui.divText('Font: ', {style: {color: 'var(--grey-4)', margin: '0 5px 0 0'}}),
        fontColorPicker,
        ui.divText('BG: ', {style: {color: 'var(--grey-4)', margin: '0 5px 0 10px'}}),
        bgColorPicker,
      ], {style: {alignItems: 'center'}}),
    ]);
    $(setColors).children('label').text('Colors');

    const setTransperency = ui.buttonsInput([
      //@ts-ignore
      imageRange,
    ]);
    $(setTransperency).children('label').text('Transperency');
    $(setTransperency).css('margin', '0');

    const textPosition = ui.buttonsInput([
      //@ts-ignore
      ui.divH([
        ui.divText('x: ', {style: {color: 'var(--grey-4)'}}),
        textX,
      ]),
      //@ts-ignore
      ui.divH([
        ui.divText('y: ', {style: {color: 'var(--grey-4)'}}),
        textY,
      ]),
    ]);
    $(textPosition).children('label').text('Text position');


    removeColorPicker.onchange = (e) =>{
      //@ts-ignore
      const rgbColor = convertHex(removeColorPicker.value);
      //@ts-ignore
      removeImageBG(image, rgbColor, imageRange.value, imageMargin.value);
      //@ts-ignore
      removeColorPicker.value = null;
    };

    bgColorPicker.onchange = (e) =>{
      //@ts-ignore
      drawIcon(bgColorPicker.value);
    };

    fontColorPicker.onchange = (e) =>{
      //@ts-ignore
      setText(text.value, fontColorPicker.value, fontSize.value, fontWeight.value);
    };

    imageRange.onchange = (e)=>{
      //@ts-ignore
      setImage(image, imageRange.value, imageMargin.value);
    };

    textX.onchange = (e)=>{
      //@ts-ignore
      fX=textX.value;
      //@ts-ignore
      setText(text.value, fontColorPicker.value, fontSize.value, fontWeight.value);
    };

    textY.onchange = (e)=>{
      //@ts-ignore
      fY=textY.value;
      //@ts-ignore
      setText(text.value, fontColorPicker.value, fontSize.value, fontWeight.value);
    };

    const dwnl = ui.bigButton('Download', ()=> saveOutputIcon());

    const image = ui.element('img', 'image');
    $(image).hide();

    const image1 = ui.element('img', 'image');
    $(image1).hide();
    const image2 = ui.element('img', 'image');
    $(image2).hide();
    const image3 = ui.element('img', 'image');
    $(image3).hide();

    const browseFile = ui.element('input', 'browseFile');
    $(browseFile).attr('type', 'file');
    $(browseFile).hide();


    browseFile.onchange = (e) =>{
      //@ts-ignore
      fileUrl.value = $(browseFile).val().replace(/^.*[\\\/]/, '');
      importImageFile(e);
    };


    const canvaBG = ui.element('canvas', 'canva');
    $(canvaBG).attr('width', '300');
    $(canvaBG).attr('height', '300');
    $(canvaBG).css('position', 'absolute');
    $(canvaBG).css('z-index', '1');

    const canvaImage = ui.element('canvas', 'canva');
    $(canvaImage).attr('width', '300');
    $(canvaImage).attr('height', '300');
    $(canvaImage).css('position', 'absolute');
    $(canvaImage).css('z-index', '2');

    const canvaText = ui.element('canvas', 'canva');
    $(canvaText).attr('width', '300');
    $(canvaText).attr('height', '300');
    $(canvaText).css('position', 'absolute');
    $(canvaText).css('z-index', '3');

    const canvaOutPut = ui.element('canvas', 'canva');
    $(canvaOutPut).attr('width', '300');
    $(canvaOutPut).attr('height', '300');
    $(canvaOutPut).hide();

    const canva = ui.div([
      canvaBG,
      canvaImage,
      canvaText,
    ]);

    $(canva).css('padding', '10px');
    $(canva).css('border', '1px dashed var(--grey-3)');
    $(canva).css('position', 'relative');
    $(canva).css('width', '320px');
    $(canva).css('height', '320px');

    const radius = 60;
    const x = 0;
    const y = 0;
    const width = 300;
    const height = 300;

    let fX = 0;
    let fY = 0;

    drawIcon('#efefef');

    this.root.appendChild(ui.divH([
      ui.divV([
        canva,
        canvaOutPut,
      ]),
      ui.divV([
        browseFile,
        ui.div([
          fileButtons,
          imageMargin,
          setTransperency,
          removeColor,
          text,
          fontSize,
          textPosition,
          fontWeight,
          setColors,
          ui.buttonsInput([dwnl]),
        ], {classes: 'ui-form'}),
        image,
        image1,
        image2,
        image3,
      ]),
    ]));
    function setText(string, color, size, fontWeight) {
      //@ts-ignore
      const ctx1 = canvaText.getContext('2d');
      //@ts-ignore
      const ctxOutput = canvaOutPut.getContext('2d');
      //@ts-ignore
      ctx1.clearRect(x, y, canvaText.width, canvaText.height);
      iconMask(ctx1);
      ctx1.clip();
      ctx1.fillStyle = color;
      ctx1.textAlign = 'center';
      ctx1.font = fontWeight+' '+size+'px Roboto';
      //@ts-ignore
      ctx1.fillText(string, 150-fX, (canvaText.height/2)+(size/3)-fY);
      //@ts-ignore
      image3.src = canvaText.toDataURL();
      ctxOutput.drawImage(image3, x, y, width, height);
    }

    function setImage(image, transperency, margin) {
      transperency = transperency/100;
      //@ts-ignore
      const ctx1 = canvaImage.getContext('2d');
      //@ts-ignore
      const ctxOutput = canvaOutPut.getContext('2d');
      //@ts-ignore
      ctx1.clearRect(x, y, canvaImage.width, canvaImage.height);
      iconMask(ctx1);
      ctx1.clip();
      ctx1.globalAlpha = transperency;
      ctx1.drawImage(image, x+(margin/2), y+(margin/2), width-margin, height-margin);
      //@ts-ignore
      image2.src = canvaImage.toDataURL();
      ctxOutput.drawImage(image2, x+(margin/2), y+(margin/2), width-margin, height-margin);
    }
    function setBlendMode(value) {
      //@ts-ignore
      const ctx1 = canvaImage.getContext('2d');
      //@ts-ignore
      const ctxOutput = canvaOutPut.getContext('2d');
      mergeCanvas();
      ctxOutput.globalCompositeOperation = value;
      //@ts-ignore
      const imgdata = ctxOutput.getImageData(x+(margin/2), y+(margin/2), width-margin, height-margin);
      //@ts-ignore
      ctx1.clearRect(x, y, canvaImage.width, canvaImage.height);
      iconMask(ctx1);
      ctx1.clip();
      //@ts-ignore
      ctx1.putImageData(imgdata, x+(margin/2), y+(margin/2));
      //@ts-ignore
      image2.src = canvaImage.toDataURL();
      //@ts-ignore
      ctxOutput.drawImage(image2, x+(margin/2), y+(margin/2), width-margin, height-margin);
    }

    function removeImageBG(image, color, transperency, margin) {
      //@ts-ignore
      const ctx1 = canvaImage.getContext('2d');
      //@ts-ignore
      const ctxOutput = canvaOutPut.getContext('2d');
      const imgdata = ctx1.getImageData(x+(margin/2), y+(margin/2), width-margin, height-margin);
      const data = imgdata.data;
      //@ts-ignore
      ctx1.clearRect(x, y, canvaImage.width, canvaImage.height);
      iconMask(ctx1);
      ctx1.clip();

      for (let i = 0, n = data.length; i <n; i += 4) {
        const r = data[i];
        const g = data[i+1];
        const b = data[i+2];
        if ((r >= color.r-20 && r<=color.r+20) && (g >= color.g-20 && g<=color.g+20) && (b >= color.b-20 && b<=color.b+20))
          data[i+3] = '0.0';

        if ((r >= color.r-40 && r<=color.r+40) && (g >= color.g-40 && g<=color.g+40) && (b >= color.b-40 && b<=color.b+40))
          data[i+3] = '0.1';

        if ((r >= color.r-60 && r<=color.r+60) && (g >= color.g-60 && g<=color.g+60) && (b >= color.b-60 && b<=color.b+60))
          data[i+3] = '0.2';

        if ((r >= color.r-80 && r<=color.r+80) && (g >= color.g-80 && g<=color.g+80) && (b >= color.b-80 && b<=color.b+80))
          data[i+3] = '0.3';

        if ((r >= color.r-100 && r<=color.r+100) && (g >= color.g-100 && g<=color.g+100) && (b >= color.b-100 && b<=color.b+100))
          data[i+3] = '0.4';
      }
      ctx1.putImageData(imgdata, x+(margin/2), y+(margin/2));
      //@ts-ignore
      image.src = canvaImage.toDataURL();
      //@ts-ignore
      image2.src = canvaImage.toDataURL();
      ctxOutput.drawImage(image2, x+(margin/2), y+(margin/2), width-margin, height-margin);
    }

    function drawIcon(color) {
      if (color === undefined)
        color='#efefef';

      //@ts-ignore
      const ctx = canvaBG.getContext('2d');
      //@ts-ignore
      const ctxOutput = canvaOutPut.getContext('2d');
      //@ts-ignore
      ctx.clearRect(x, y, canvaImage.width, canvaImage.height);
      ctx.fillStyle = color;
      ctx.beginPath();
      ctx.moveTo(x + radius, y);
      ctx.lineTo(x + width - radius, y);
      ctx.quadraticCurveTo(x + width, y, x + width, y + radius);
      ctx.lineTo(x + width, y + height - radius);
      ctx.quadraticCurveTo(x + width, y + height, x + width - radius, y + height);
      ctx.lineTo(x + radius, y + height);
      ctx.quadraticCurveTo(x, y + height, x, y + height - radius);
      ctx.lineTo(x, y + radius);
      ctx.quadraticCurveTo(x, y, x + radius, y);
      ctx.closePath();
      ctx.fill();
      //@ts-ignore
      image1.src = canvaBG.toDataURL();
      ctxOutput.drawImage(image1, x, y, width, height);
    }

    function iconMask(ctx) {
      ctx.beginPath();
      ctx.moveTo(x + radius, y);
      ctx.lineTo(x + width - radius, y);
      ctx.quadraticCurveTo(x + width, y, x + width, y + radius);
      ctx.lineTo(x + width, y + height - radius);
      ctx.quadraticCurveTo(x + width, y + height, x + width - radius, y + height);
      ctx.lineTo(x + radius, y + height);
      ctx.quadraticCurveTo(x, y + height, x, y + height - radius);
      ctx.lineTo(x, y + radius);
      ctx.quadraticCurveTo(x, y, x + radius, y);
      ctx.closePath();
    }

    //Upload, delete and save image

    function importImageFile(e) {
      const file = e.target.files[0];
      if (!file)
        return;

      const reader = new FileReader();
      reader.onload = function(e) {
        const contents = e.target.result;
        //@ts-ignore
        image.src=contents;
      };
      reader.onloadend = function(e) {
        //@ts-ignore
        setImage(image, imageRange.value, imageMargin.value);
      };
      reader.readAsDataURL(file);
    }

    function mergeCanvas() {
      //@ts-ignore
      const ctxOutput = canvaOutPut.getContext('2d');
      ctxOutput.clearRect(x, y, width, height);
      iconMask(ctxOutput);
      ctxOutput.clip();

      ctxOutput.drawImage(image1, x, y, width, height);
      ctxOutput.drawImage(image2, x, y, width, height);
      ctxOutput.drawImage(image3, x, y, width, height);
    }

    function saveOutputIcon() {
      mergeCanvas();

      const outputLink = ui.element('a');
      outputLink.innerHTML = 'download image';
      //@ts-ignore
      outputLink.href = canvaOutPut.toDataURL();
      //@ts-ignore
      outputLink.download = 'package.png';
      $(outputLink).trigger('click');
    }

    function deleteImageFile() {
      //@ts-ignore
      const ctx = canvaImage.getContext('2d');
      ctx.clearRect(x, y, width, height);
      //@ts-ignore
      image.src='';
      //@ts-ignore
      browseFile.value='';
      fileUrl.value='';
      //@ts-ignore
      imageRange.value = 100;
      //@ts-ignore
      setImage(image, imageRange.value, imageMargin.value);
    }

    function convertHex(hex) {
      hex = hex.replace('#', '');
      const r = parseInt(hex.substring(0, 2), 16);
      const g = parseInt(hex.substring(2, 4), 16);
      const b = parseInt(hex.substring(4, 6), 16);
      const result = {r, g, b};
      return result;
    }
  }
}
