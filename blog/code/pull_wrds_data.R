# --- Header -------------------------------------------------------------------
# See LICENSE file for details
#
# This code pulls data from WRDS
# ------------------------------------------------------------------------------

suppressMessages({
  library(RPostgres)
  library(DBI)
})

if (!exists("cfg")) source("code/read_config.R")

save_wrds_data <- function(df, fname) {
  if(file.exists(fname)) {
    file.rename(
      fname,
      paste0(
        substr(fname, 1, nchar(fname) - 4),
        "_",
        format(file.info(fname)$mtime, "%Y-%m-%d_%H_%M_%S"),".rds")
    )
  }
  saveRDS(df, fname)
}

# --- Connect to WRDS ----------------------------------------------------------

wrds <- dbConnect(
  Postgres(),
  host = 'wrds-pgdata.wharton.upenn.edu',
  port = 9737,
  user = cfg$wrds_user,
  password = cfg$wrds_pwd,
  sslmode = 'require',
  dbname = 'wrds'
)

message("Logged on to WRDS ...")

# --- Specify filters and variables --------------------------------------------

stat_vars <- c(
  "gvkey", "prirow", "loc", "fic", "sic", "spcindcd", "ipodate", "dldte"
)
stat_var_str <- paste(stat_vars, collapse = ", ")

dyn_vars <- c(
  "gvkey",  "conm", "fyear", "datadate", "iid", "indfmt", "sich",
  "consol", "popsrc", "datafmt", "curcd", "exchg", "fyr",
  "act", "ap", "aqc", "at", "ceq", "che", "cogs",
  "dlc", "dp", "dpc",  "ib", "ibc", "lct", "nicon", "oancf",
  "ivncf", "fincf", "oiadp", "pi", "sale"
)

dyn_var_str <- paste(dyn_vars, collapse = ", ")
dyn_filter <- paste(
  "consol='C' and (indfmt='INDL') and datafmt='HIST_STD'",
  "and popsrc='I'"
)

# --- Pull Global data --------------------------------------------------

message("Pulling global static data ... ", appendLF = FALSE)
res <- dbSendQuery(wrds, paste("select ", stat_var_str, "from COMP.G_COMPANY"))
cstat_global_static <- dbFetch(res,n=-1)
dbClearResult(res)
message("done!")

message("Pulling global annual financial data ... ", appendLF = FALSE)
res <- dbSendQuery(wrds, paste(
  "select", dyn_var_str, "from COMP.G_FUNDA", "where", dyn_filter
))
cstat_global_dynamic <- dbFetch(res, n=-1)
dbClearResult(res)
message("done!")

message("Pulling currency exchange data ... ", appendLF = FALSE)
res <- dbSendQuery(wrds, paste(
  "SELECT * FROM comp.g_exrt_dly"
))
cstat_exrate_daily <- dbFetch(res, n=-1)
dbClearResult(res)
message("done!")

message(
  "Pulling global securities data for fiscal year ends ... ",
  appendLF = FALSE
)
res <- dbSendQuery(
  wrds,
  paste(
    "select link.gvkey, secd2.isin, secd2.curcdd, secd2.iid,",
    "fyedate, secd2.datadate, secd2.cshoc, secd2.prccd, secd2.trfd,",
    "secd2.ajexdi from",
    "(select funda1.gvkey, max(secd1.datadate) as datadate,",
    "funda1.datadate as fyedate",
    "from comp.g_secd as secd1 left join comp.g_funda as funda1 on",
    "secd1.gvkey = funda1.gvkey and funda1.datadate - secd1.datadate >= 0",
    "and funda1.datadate - secd1.datadate < 7",
    "group by funda1.gvkey, funda1.datadate) as link",
    "inner join comp.g_secd as secd2 on",
    "link.gvkey = secd2.gvkey and link.datadate = secd2.datadate"
  )
)

cstat_global_sec <- dbFetch(res,n=-1)
dbClearResult(res)
message("done!")

dbDisconnect(wrds)
message("Disconnected from WRDS")

# --- Store data ------------------------------------------------------

cstat_fund <- merge(cstat_global_static, cstat_global_dynamic, by="gvkey")
save_wrds_data(cstat_exrate_daily, "data/pulled/cstat_exrate_daily.rds")
save_wrds_data(cstat_fund, "data/pulled/cstat_global_fund.rds")
save_wrds_data(cstat_global_sec, "data/pulled/cstat_global_sec.rds")

