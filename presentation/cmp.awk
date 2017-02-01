BEGIN {
  FS="\t"
}
NR != 1 { it[$2, $1] += $3; anz[$2, $1]++ }
END {
  for (combined in it) {
    split(combined, separate, SUBSEP)
    it[combined] /= anz[combined]
    print separate[2], it[combined] | "sort -g > it_" separate[1] ".csv"
  }
}
