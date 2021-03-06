# from https://github.com/r-lib/usethis/blob/master/R/style.R
# as of Dec 31, 2017
bullet <- function(lines, bullet) {
  lines <- paste0(bullet, " ", lines)
  cat_line(lines)
}

todo_bullet <- function() crayon::red(clisymbols::symbol$bullet)

todo <- function(...) {
  bullet(paste0(...), bullet = todo_bullet())
}

done <- function(...) {
  bullet(paste0(...), bullet = crayon::green(clisymbols::symbol$tick))
}

failure <- function(...) {
  bullet(paste0(...), bullet = crayon::red(clisymbols::symbol$cross))
}

stop_failure <- function(...) {
  stop(paste0(crayon::red(clisymbols::symbol$cross), " ", ..., "\n"), call. = FALSE)
}

code_block <- function(..., copy = interactive()) {
  block <- paste0("  ", c(...), collapse = "\n")
  if (copy && clipr::clipr_available()) {
    clipr::write_clip(paste0(c(...), collapse = "\n"))
    message("Copying code to clipboard:")
  }
  cat_line(crayon::make_style("darkgrey")(block))
}

cat_line <- function(...) {
  cat(..., "\n", sep = "")
}

field <- function(...) {
  x <- paste0(...)
  crayon::green(x)
}
value <- function(...) {
  x <- paste0(...)
  crayon::blue(encodeString(x, quote = "'"))
}

code <- function(...) {
  x <- paste0(...)
  crayon::make_style("darkgrey")(encodeString(x, quote = "`"))
}

collapse <- function(x, sep = ", ") {
  paste0(x, collapse = sep)
}
