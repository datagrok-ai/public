import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

/**
 * Custom cell renderer for sdisc-rule-violation semantic type.
 * Displays validation rule violations in the format:
 * Rule Id: <id>
 * Message: <message>
 * Value: <value>
 */
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
    } catch (e) {
      // If parsing fails, try to display as-is
      g.font = '12px sans-serif';
      g.fillText(String(value), x + 10, y + 15);
      return;
    }

    if (!Array.isArray(errors) || errors.length === 0)
      return;

    //in case we are in grid - only show number of violated rules
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

    //pre-calculate tootip height
    const mandatoryKeys = ['ruleID', 'core_id', 'value', 'values', 'variables', 'message'];
    const array = [];
    let tooltipHeight = 0;
    for (const error of errors) {
      const rowsByKey = {};
      for (const key of Object.keys(error)) {
        //advance to next line
        tooltipHeight += lineHeight;
        rowsByKey[key] = [];
        if (mandatoryKeys.includes(key)) {
          g.font = 'bold 12px sans-serif';
          const keyWithPadding = g.measureText(key).width + padding * 2;
          //split value by words
          // Wrap long values and track the final Y position
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
          // Push the last line
          if (line.trim()) {
            rowsByKey[key].push(line.trim());
            tooltipHeight += lineHeight;
          }
        }
      }
      array.push(rowsByKey);
    }
    //consider distance between violated rules in case more than one
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
      // Draw horizontal separator before each error (except the first)
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

      // Get pre-calculated rows for this error
      const rowsByKey = array[errorIndex];

      Object.keys(error).forEach((key) => {
        if (!key && !error[key] || !mandatoryKeys.includes(key))
          return;

        // Advance to next line before drawing
        currentY += lineHeight;

        // Draw key (label)
        g.fillStyle = labelColor;
        g.font = 'bold 12px sans-serif';
        const keyText = key ? `${key}: ` : '';
        const keyY = currentY;
        g.fillText(keyText, x + padding, keyY);

        // Calculate position for value (after key)
        const keyWidth = keyText ? g.measureText(keyText).width : 0;
        const valueX = x + padding + keyWidth;

        // Draw value using pre-calculated split rows
        g.fillStyle = textColor;
        g.font = '12px sans-serif';
        let valueY = keyY; // Start value at same Y as key

        // Use pre-calculated split rows if available
        if (rowsByKey[key] && rowsByKey[key].length > 0) {
          // Draw each pre-calculated row
          rowsByKey[key].forEach((row: string) => {
            g.fillText(row, valueX, valueY);
            valueY += lineHeight;
          });
        } else {
          // Fallback: use original value if not pre-calculated
          const valueText = error[key].toString() || 'N/A';
          g.fillText(valueText, valueX, valueY);
        }

        // Always update currentY to the bottom of the value
        currentY = valueY;
      });
    });

    g.restore();
  }
}
