version 1.1

workflow csv_overlap {
  input {
    String query_csv
    String haystack_csv
  }

  ####################################################################
  # STEP 1 — Split CSVs
  ####################################################################
  Array[String] q_raw = split(query_csv, ",")
  Array[String] h_raw = split(haystack_csv, ",")

  ####################################################################
  # STEP 2 — Trim whitespace (case-sensitive)
  # NOTE: This uses the ONLY legal comprehension form.
  ####################################################################
  Array[String] q_trim = [for x in q_raw: trim(x)]
  Array[String] h_trim = [for x in h_raw: trim(x)]

  ####################################################################
  # STEP 3 — Replace empty tokens with a sentinel string
  ####################################################################
  Array[String] q_norm = [
    for x in q_trim:
      if x == "" then "___NULL___" else x
  ]

  Array[String] h_norm = [
    for x in h_trim:
      if x == "" then "___NULL___" else x
  ]

  ####################################################################
  # STEP 4 — Compute dirty overlap by mapping mismatches to ""
  ####################################################################
  Array[String] overlap_dirty = [
    for x in q_norm:
      if (x != "___NULL___" && x in h_norm) then x else ""
  ]

  ####################################################################
  # STEP 5 — Clean empty values (this is valid syntax)
  ####################################################################
  Array[String] overlap_list = [
    for x in overlap_dirty:
      if x != "" then x else ""
  ]

  ####################################################################
  # STEP 6 — Final formatted outputs
  ####################################################################
  String  overlap_csv = sep(",", overlap_list)
  Boolean has_overlap = length(overlap_list) > 0

  output {
    Array[String] overlap_list = overlap_list
    String        overlap_csv  = overlap_csv
    Boolean       has_overlap  = has_overlap
  }
}
