## function to retrieve raw count data for specific contigs and return it in
## tall, tidy form

extract_and_tidy <- function(contig, raw_data, exp_des) {
  require(reshape2)
  the_row <- which(rownames(raw_data) %in% contig)
  jDat <- cbind(exp_des, t(raw_data[the_row, ]))
  names(jDat) <- gsub("WPW_Inoculation_Trinity_C500_", "", names(jDat))
  jDat <- melt(jDat,
               id.vars = names(exp_des),
               variable.name = 'contig', value.name = 'count')
  return(jDat)
}

