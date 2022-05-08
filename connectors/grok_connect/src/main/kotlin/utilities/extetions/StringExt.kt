package utilities.extetions

fun String?.substringOrEmpty(startIndex: Int, endIndex: Int): String {
    return try {
        this?.substring(startIndex, endIndex).orEmpty()
    } catch (e: Exception) {
        ""
    }
}