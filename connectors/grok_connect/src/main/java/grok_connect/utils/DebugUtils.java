package grok_connect.utils;

import org.apache.commons.io.FileUtils;

import java.io.File;
import java.math.BigDecimal;
import java.math.MathContext;

public class DebugUtils {
    public static byte[] substituteBlob(String path) {
        try {
            return FileUtils.readFileToByteArray(new File(path));
        }
        catch (Exception e){
            System.out.println(e.toString());
        }
        return null;
    }

    public static class NumericColumnStats {
        public BigDecimal min;
        public BigDecimal max;
        public BigDecimal mean;
        public BigDecimal valuesCounter = new BigDecimal(0);

        public void updateStats(Object value) {
            try {
                BigDecimal bigDecimalValue = new BigDecimal(value.toString());
                max = getMaxForColumn(max, bigDecimalValue);
                min = getMinForColumn(min, bigDecimalValue);
                mean = getMeanOnStep(bigDecimalValue);
            }
            catch (NumberFormatException e)  { }
        }

        public BigDecimal getMeanOnStep (BigDecimal value) {
            valuesCounter = valuesCounter.add(new BigDecimal(1));

            if (mean == null)
                mean = new BigDecimal(0);

            return mean = mean.add(value.subtract(mean).divide(valuesCounter, MathContext.DECIMAL32));
        }

        public static BigDecimal getMaxForColumn(BigDecimal max, BigDecimal value) {
            return max == null ? value : value.compareTo(max) == 1 ? value : max;
        }

        public static BigDecimal getMinForColumn(BigDecimal min, BigDecimal value) {
            return min == null ? value : value.compareTo(min) == -1 ? value : min;
        }
    }


}
