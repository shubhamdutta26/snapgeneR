#' Convert SnapGene file to R list
#' @param filepath Path to SnapGene .dna file
#' @return List containing parsed SnapGene data
read_snapgene <- function(filepath) {
  tryCatch({
    con <- file(filepath, "rb")

    # Check file format
    first_byte <- readBin(con, "raw", n=1)
    if (!identical(first_byte, as.raw(0x09))) {
      stop("Wrong format for a SnapGene file! First byte mismatch")
    }

    # Read document properties
    length_bytes <- readBin(con, "raw", n=4)
    length <- raw_to_int(length_bytes)

    title_raw <- readBin(con, "raw", n=8)
    title <- rawToChar(title_raw)

    if (length != 14 || title != "SnapGene") {
      stop(sprintf("Wrong format for a SnapGene file! Length: %d, Title: %s", length, title))
    }

    # Initialize data structure
    isDNA_bytes <- readBin(con, "raw", n=2)
    exportVersion_bytes <- readBin(con, "raw", n=2)
    importVersion_bytes <- readBin(con, "raw", n=2)

    data <- list(
      isDNA = raw_to_int(isDNA_bytes),
      exportVersion = raw_to_int(exportVersion_bytes),
      importVersion = raw_to_int(importVersion_bytes),
      features = list()
    )

    # Read blocks until end of file
    while (length(raw_data <- readBin(con, "raw", n=1)) > 0) {
      block_type <- as.integer(raw_data)
      block_size_bytes <- readBin(con, "raw", n=4)
      block_size <- raw_to_int(block_size_bytes)

      if (block_type == 0) {
        # DNA sequence block
        props_raw <- readBin(con, "raw", n=1)
        props <- as.integer(props_raw)

        data$dna <- list(
          topology = if (bitwAnd(props, 0x01) > 0) "circular" else "linear",
          strandedness = if (bitwAnd(props, 0x02) > 0) "double" else "single",
          damMethylated = bitwAnd(props, 0x04) > 0,
          dcmMethylated = bitwAnd(props, 0x08) > 0,
          ecoKIMethylated = bitwAnd(props, 0x10) > 0,
          length = block_size - 1
        )

        seq_raw <- readBin(con, "raw", n=block_size-1)
        data$seq <- rawToChar(seq_raw)

      } else if (block_type == 6) {
        # Notes block
        block_raw <- readBin(con, "raw", n=block_size)
        block_content <- rawToChar(block_raw)
        note_data <- xml2::read_xml(block_content) %>%
          xml2::as_list() %>%
          parse_dict()
        data$notes <- note_data$Notes

      } else if (block_type == 10) {
        # Features block
        block_raw <- readBin(con, "raw", n=block_size)
        block_content <- rawToChar(block_raw)
        features_data <- xml2::read_xml(block_content) %>%
          xml2::as_list()

        features <- features_data$Features$Feature
        if (!is.null(features) && !is.list(features[[1]])) {
          features <- list(features)
        }

        if (!is.null(features)) {
          for (feature in features) {
            segments <- feature$Segment
            if (!is.null(segments) && !is.list(segments[[1]])) {
              segments <- list(segments)
            }

            # Process segments
            segments_ranges <- lapply(segments, function(seg) {
              range_str <- get_attr(seg, "range")
              if (!is.null(range_str)) {
                range <- as.numeric(strsplit(range_str, "-")[[1]])
                return(sort(range))
              }
              return(NULL)
            })

            segments_ranges <- segments_ranges[!sapply(segments_ranges, is.null)]

            # Process qualifiers
            qualifiers <- feature$Q
            if (!is.null(qualifiers) && !is.list(qualifiers[[1]])) {
              qualifiers <- list(qualifiers)
            }

            parsed_qualifiers <- list()
            if (!is.null(qualifiers)) {
              for (qualifier in qualifiers) {
                q_name <- get_attr(qualifier, "name")
                if (is.null(q_name)) next

                q_value <- qualifier$V

                if (is.null(q_value)) {
                  next
                } else if (is.list(q_value)) {
                  if (length(names(q_value[[1]])) == 1) {
                    parsed_qualifiers[[q_name]] <- lapply(q_value, function(v) {
                      parse_html(v[[1]])
                    })
                  } else {
                    parsed_qualifiers[[q_name]] <- lapply(q_value, function(v) {
                      setNames(parse_html(v[[1]]), v[[2]])
                    })
                  }
                } else {
                  parsed_qualifiers[[q_name]] <- parse_html(q_value)
                }
              }
            }

            # Add label if missing
            if (is.null(parsed_qualifiers$label)) {
              parsed_qualifiers$label <- get_attr(feature, "name", "")
            }

            # Process notes
            if (is.null(parsed_qualifiers$note)) {
              parsed_qualifiers$note <- list()
            } else if (!is.list(parsed_qualifiers$note)) {
              parsed_qualifiers$note <- list(parsed_qualifiers$note)
            }

            # Add color note
            first_segment_color <- get_attr(segments[[1]], "color", "#000000")
            parsed_qualifiers$note <- c(
              parsed_qualifiers$note,
              paste("color:", first_segment_color)
            )

            # Create feature entry
            if (length(segments_ranges) > 0) {
              data$features[[length(data$features) + 1]] <- list(
                start = min(sapply(segments_ranges, `[`, 1)) - 1,
                end = max(sapply(segments_ranges, `[`, 2)),
                strand = switch(get_attr(feature, "directionality", "0"),
                                "0" = ".",
                                "1" = "+",
                                "2" = "-",
                                "3" = "=",
                                "."),
                type = get_attr(feature, "type", "misc_feature"),
                name = get_attr(feature, "name", ""),
                color = first_segment_color,
                textColor = "black",
                segments = segments,
                row = 0,
                isOrf = FALSE,
                qualifiers = parsed_qualifiers
              )
            }
          }
        }
      } else {
        # Skip unknown blocks
        readBin(con, "raw", n=block_size)
      }
    }

    return(data)

  }, finally = {
    if (exists("con")) close(con)
  })
}
