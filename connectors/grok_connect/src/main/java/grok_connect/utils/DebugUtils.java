package grok_connect.utils;

import org.apache.commons.io.FileUtils;

import java.io.File;

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
}
