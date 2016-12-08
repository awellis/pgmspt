library("stringr")

fix_inline_citations <- function(text) {
    re_inline_year <- "[(]\\d{4}[)]"
    re_author <- "[[:alpha:]- ]+"
    re_author_year <- paste(re_author, re_inline_year)
    re_ampersand <- " & "
    re_ampersand_author_year <- sprintf("%s(?=%s)", re_ampersand, re_author_year)
    str_replace_all(text, re_ampersand_author_year, " and ")
}

do_fix <- c(
    "Jones & Name (2005) found...",
    "Jones & Hyphen-Name (2005) found...",
    "Jones & Space Name (2005) found...",
    "Marge, Maggie, & Lisa (2005) found...")


do_not_fix <- c(
    "...have been found (Jones & Name, 2005)",
    "...have been found (Jones & Hyphen-Name, 2005)",
    "...have been found (Jones & Space Name, 2005)",
    "...have been found (Marge, Maggie, & Lisa, 2005)")
