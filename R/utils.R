#' Read SnapGene files and convert to R list structure
#' Required packages: xml2, httr, stringr, Biostrings, IRanges, bit64

#' Convert raw bytes to integer (big endian)
#' @param bytes Raw vector of bytes
#' @return Integer value
raw_to_int <- function(bytes) {
  if (length(bytes) == 4) {
    # 4-byte integer
    sum(256^(3:0) * as.integer(bytes))
  } else if (length(bytes) == 2) {
    # 2-byte short
    sum(256^(1:0) * as.integer(bytes))
  } else {
    stop("Unsupported number of bytes")
  }
}

#' Parse HTML content
#' @param val String containing HTML content
#' @return Cleaned text string
parse_html <- function(val) {
  if (is.character(val)) {
    # Remove HTML tags, normalize whitespace, replace quotes
    cleaned <- val %>%
      stringr::str_replace_all("<[^>]+>", " ") %>%
      stringr::str_replace_all("\n", " ") %>%
      stringr::str_replace_all('"', "'") %>%
      stringr::str_trim() %>%
      stringr::str_squish()
    return(cleaned)
  }
  return(val)
}

#' Parse dictionary/list recursively
#' @param obj List object potentially containing HTML content
#' @return List with parsed HTML content
parse_dict <- function(obj) {
  if (is.list(obj)) {
    for (key in names(obj)) {
      if (is.character(obj[[key]])) {
        obj[[key]] <- parse_html(obj[[key]])
      } else if (is.list(obj[[key]])) {
        obj[[key]] <- parse_dict(obj[[key]])
      }
    }
  }
  return(obj)
}

#' Extract attribute value from XML node attributes
#' @param node XML node
#' @param attr_name Name of attribute to extract
#' @param default Default value if attribute not found
#' @return Attribute value or default
get_attr <- function(node, attr_name, default = NULL) {
  if (!is.null(node) && !is.null(node$.attrs) && attr_name %in% names(node$.attrs)) {
    return(node$.attrs[[attr_name]])
  }
  return(default)
}
