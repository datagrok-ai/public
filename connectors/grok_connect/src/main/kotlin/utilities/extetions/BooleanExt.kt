package utilities.extetions

fun Boolean.takeIfTrue() = this.takeIf { it }
fun Boolean.takeIfFalse() = this.takeIf { !it }