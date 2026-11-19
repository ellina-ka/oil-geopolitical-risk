# Data access and validation.
# Raw FRED files are written once to data/raw/ and never modified; every
# download is logged (date, rows, md5) in data/raw/MANIFEST.csv.

fred_url <- function(id) {
  paste0("https://fred.stlouisfed.org/graph/fredgraph.csv?id=", id)
}

#' Download a FRED series (or reuse the cached raw file) and return a clean
#' data.frame(date, value, series).
fetch_fred <- function(id, raw_dir = "data/raw", refresh = FALSE) {
  stopifnot(is.character(id), length(id) == 1, nzchar(id))
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(raw_dir, paste0(id, ".csv"))

  if (refresh || !file.exists(path)) {
    tmp <- tempfile(fileext = ".csv")
    ok <- tryCatch({
      utils::download.file(fred_url(id), tmp, quiet = TRUE, mode = "wb")
      TRUE
    }, error = function(e) {
      warning("Download failed for ", id, ": ", conditionMessage(e))
      FALSE
    })
    if (ok) {
      read_fred_csv(tmp, id) # validate before replacing the cached copy
      file.copy(tmp, path, overwrite = TRUE)
      log_download(id, path, raw_dir)
    } else if (!file.exists(path)) {
      stop("No data for ", id, ": download failed and no cached copy at ", path)
    }
  }
  read_fred_csv(path, id)
}

log_download <- function(id, path, raw_dir) {
  manifest <- file.path(raw_dir, "MANIFEST.csv")
  row <- data.frame(
    series = id,
    downloaded_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    rows = nrow(utils::read.csv(path)),
    md5 = unname(tools::md5sum(path)),
    url = fred_url(id)
  )
  utils::write.table(row, manifest, sep = ",", row.names = FALSE,
                     col.names = !file.exists(manifest),
                     append = file.exists(manifest))
}

#' Parse a FRED csv. FRED marks missing observations with "." or blanks;
#' these are dropped (holidays), so returns computed later span the gap.
read_fred_csv <- function(path, id) {
  x <- utils::read.csv(path, stringsAsFactors = FALSE,
                       na.strings = c(".", ""), check.names = FALSE)
  if (ncol(x) != 2) stop(path, ": expected 2 columns, found ", ncol(x))
  names(x) <- c("date", "value")
  x$date <- as.Date(x$date)
  if (anyNA(x$date)) stop(path, ": unparseable dates")
  x$value <- suppressWarnings(as.numeric(x$value))
  x <- x[!is.na(x$value), , drop = FALSE]
  if (nrow(x) == 0) stop(path, ": no numeric observations")
  if (anyDuplicated(x$date)) stop(path, ": duplicated dates")
  if (is.unsorted(x$date)) stop(path, ": dates not sorted")
  x$series <- id
  rownames(x) <- NULL
  x
}

#' Read and check the hand-coded event list.
read_events <- function(path = "data/events.csv") {
  ev <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("event_id", "event_date", "description", "direction", "source_url")
  missing_cols <- setdiff(required, names(ev))
  if (length(missing_cols)) stop("events file lacks: ", paste(missing_cols, collapse = ", "))
  ev$event_date <- as.Date(ev$event_date)
  if (anyNA(ev$event_date)) stop("events file: unparseable event_date")
  if (anyDuplicated(ev$event_id)) stop("events file: duplicated event_id")
  bad_dir <- setdiff(unique(ev$direction), c("escalation", "de-escalation"))
  if (length(bad_dir)) stop("unknown direction value(s): ", paste(bad_dir, collapse = ", "))
  no_src <- ev$event_id[is.na(ev$source_url) | !nzchar(ev$source_url)]
  if (length(no_src)) warning("events without a source: ", paste(no_src, collapse = ", "))
  ev[order(ev$event_date), ]
}
