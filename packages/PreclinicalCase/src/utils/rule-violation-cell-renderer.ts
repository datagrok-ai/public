import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

@grok.decorators.cellRenderer({
  name: 'sdiscRuleViolationRenderer',
  cellType: 'sdisc-rule-violation',
})
export class SdiscRuleViolationCellRenderer extends DG.GridCellRenderer {
  get name(): string {
    return 'sdiscRuleViolationRenderer';
  }

  get cellType(): string {
    return 'sdisc-rule-violation';
  }

  get defaultWidth(): number {
    return 30;
  }

  get defaultHeight(): number {
    return 30;
  }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle,
  ): void {
    let isTooltip = false;
    if (x === 0 && y === 0 && w === this.defaultWidth && h === this.defaultHeight)
      isTooltip = true;
    const value = gridCell.cell.value;
    if (!value)
      return;

    let errors;
    try {
      errors = typeof value === 'string' ? JSON.parse(value) : value;
    }
    catch (e) {
      g.font = '12px sans-serif';
      g.fillText(String(value), x + 10, y + 15);
      return;
    }

    if (!Array.isArray(errors) || errors.length === 0)
      return;

    if (!isTooltip) {
      g.font = '12px sans-serif';
      g.fillText(String(errors.length), x + 10, y + 15);
      return;
    }

    const tooltipCanvasWidth = 300;
    const lineHeight = 15;
    const separatorHeight = 1;
    const separatorMargin = 5;
    const padding = 5;

    const mandatoryKeys = ['ruleID', 'core_id', 'value', 'values', 'variables', 'message'];
    const array: {[key: string]: any}[] = [];
    let tooltipHeight = 0;
    for (const error of errors) {
      const rowsByKey: {[key: string]: any} = {};
      for (const key of Object.keys(error)) {
        tooltipHeight += lineHeight;
        rowsByKey[key] = [];
        if (mandatoryKeys.includes(key)) {
          g.font = 'bold 12px sans-serif';
          const keyWithPadding = g.measureText(key).width + padding * 2;
          const valueText = String(error[key] || '');
          const words = valueText.split(key === 'message' ? ' ' : ',');
          let line = '';

          g.font = '12px sans-serif';
          for (const word of words) {
            const testLine = line + word + ' ';
            const metrics = g.measureText(testLine);
            if (metrics.width > tooltipCanvasWidth - keyWithPadding && line !== '') {
              rowsByKey[key].push(line.trim());
              tooltipHeight += lineHeight;
              line = word + ' ';
            } else
              line = testLine;
          }
          if (line.trim()) {
            rowsByKey[key].push(line.trim());
            tooltipHeight += lineHeight;
          }
        }
      }
      array.push(rowsByKey);
    }
    if (errors.length > 1)
      tooltipHeight += (errors.length - 1) * (separatorHeight + separatorMargin);

    g.save();
    g.canvas.width = tooltipCanvasWidth;
    g.canvas.style.width = `${tooltipCanvasWidth}px`;
    g.canvas.height = tooltipHeight;
    g.canvas.style.height = `${tooltipHeight}px`;

    let currentY = 0;

    const labelColor = cellStyle.textColorHtml || '#000000';
    const textColor = 'grey';
    const separatorColor = '#E0E0E0';

    errors.forEach((error, errorIndex) => {
      if (errorIndex > 0) {
        currentY += separatorMargin;
        g.strokeStyle = separatorColor;
        g.lineWidth = separatorHeight;
        g.beginPath();
        g.moveTo(x + padding, currentY);
        g.lineTo(x + tooltipCanvasWidth - padding, currentY);
        g.stroke();
        currentY += separatorMargin + separatorHeight;
      }

      const rowsByKey = array[errorIndex];

      Object.keys(error).forEach((key) => {
        if (!key && !error[key] || !mandatoryKeys.includes(key))
          return;

        currentY += lineHeight;

        g.fillStyle = labelColor;
        g.font = 'bold 12px sans-serif';
        const keyText = key ? `${key}: ` : '';
        const keyY = currentY;
        g.fillText(keyText, x + padding, keyY);

        const keyWidth = keyText ? g.measureText(keyText).width : 0;
        const valueX = x + padding + keyWidth;

        g.fillStyle = textColor;
        g.font = '12px sans-serif';
        let valueY = keyY;

        if (rowsByKey[key] && rowsByKey[key].length > 0) {
          rowsByKey[key].forEach((row: string) => {
            g.fillText(row, valueX, valueY);
            valueY += lineHeight;
          });
        } else {
          const valueText = error[key].toString() || 'N/A';
          g.fillText(valueText, valueX, valueY);
        }

        currentY = valueY;
      });
    });

    g.restore();
  }
}
