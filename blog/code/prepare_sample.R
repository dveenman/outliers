# --- Header -------------------------------------------------------------------
# See LICENSE file for details 
#
# Creates an earnings/cash flow/accrual/fiscal year returns sample based on
# Compustat Global sample
# ------------------------------------------------------------------------------
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# Implemented sample types are "World" and "Eurozone"
sample_type <- "Eurozone"
message(sprintf("Preparing %s sample... ", sample_type), appendLF = FALSE)

mleadlag <- function(x, n, ts_id) {
  pos <- match(as.numeric(ts_id) + n, as.numeric(ts_id))
  x[pos]
}

if (sample_type == "World") {
  fx <- readRDS("data/pulled/cstat_exrate_daily.rds")
  fx_usd <-  fx %>%
    filter(tocurd == "USD") %>%
    select(datadate, gbp_to_usd = exratd) %>%
    left_join(fx, by = "datadate") %>%
    mutate(
      fx_to_usd = gbp_to_usd/exratd 
    ) %>%
    select(datadate, curcd = tocurd, fx_to_usd) %>%
    arrange(curcd, datadate)
  
  sec <- readRDS("data/pulled/cstat_global_sec.rds") %>% 
    filter(!is.na(prccd)) %>%
    left_join(
      fx_usd %>% rename(curcdd = curcd), by = c("curcdd", "datadate")
    ) %>%
    mutate(prccd = prccd * fx_to_usd) %>%
    select(-fx_to_usd)
  
  fund <- readRDS("data/pulled/cstat_global_fund.rds") %>% 
    arrange(gvkey, fyr, fyear) %>%
    filter(at > 0, sale > 0) %>%
    left_join(fx_usd, by = c("curcd", "datadate")) %>%
    mutate(across(act:sale, ~.x * fx_to_usd)) %>%
    select(-fx_to_usd)
} else if (sample_type == "Eurozone") {
  sec <- readRDS("data/pulled/cstat_global_sec.rds") %>% 
    filter(!is.na(prccd)) %>%
    filter(curcdd == "EUR")
  
  locs_euro <- read_csv(
    "data/external/eea_ctries_2021-01-01.csv", col_types = "cccclll"
  ) %>% select(loc = "iso3c", euro) %>% filter(euro)

  fund <- readRDS("data/pulled/cstat_global_fund.rds") %>% 
    arrange(gvkey, fyr, fyear) %>%
    filter(at > 0, sale > 0)  %>%
    filter(curcd == "EUR") %>%
    left_join(locs_euro, by = "loc") %>% 
    filter(euro) %>% select(-euro)
} else stop("Unknown sample type.")
    
fund %>%
  group_by(gvkey, datadate) %>%
  filter(n() > 1) -> dup_cases  

if (nrow(dup_cases) > 0) stop(
  "Compustat Global has duplicate financial data at the gvkey, datadate", 
  "level. Duplicates are stored in 'dup_cases'"
)

fund %>%
  group_by(gvkey, fyear, fyr) %>%
  filter(n() > 1) -> dup_cases  

if (nrow(dup_cases) > 0) stop(
  "Compustat Global has duplicate financial data at the gvkey, fyear, fyr level.",
  "Duplicates are stored in 'dup_cases'"
)

sec <- sec %>% 
  rename(
    prcdate = datadate,
    datadate = fyedate
  ) %>%
  arrange(gvkey, datadate) %>%
  distinct() 

sec %>%
  group_by(gvkey, datadate, iid) %>%
  filter(n() > 1) -> dup_cases  

if (nrow(dup_cases) > 0) stop(
  "Compustat Global has duplicate securities data at the gvkey,",
  "datadate, iid level. Duplicates are stored in 'dup_cases'"
)


# Close to BLZ (JAR, 2016)

smp_raw <- fund %>%
  left_join(sec, by = c("gvkey", "datadate", "iid")) %>%
  arrange(gvkey, fyr, fyear) %>%
  group_by(gvkey, fyr) %>%
  mutate(
    avg_at = 0.5*(at + mleadlag(at, -1, fyear)),
    lag_at = mleadlag(at, -1, fyear),
    lag_sale = mleadlag(sale, -1, fyear),
    sales_gr = log(sale)/log(lag_sale),
    ear = ib/avg_at,
    tacc_cf = (ib - oancf)/avg_at,
    tacc_bs = (act - mleadlag(act, -1, fyear) - 
                 (che - mleadlag(che, -1, fyear)) - 
                 ((lct - mleadlag(lct, -1, fyear)) - 
                    (dlc - mleadlag(dlc, -1, fyear))) - dp)/avg_at,
    tacc = ifelse(!is.na(tacc_cf), tacc_cf, tacc_bs),
    cfo_cf = oancf/avg_at,
    cfo_bs = ear - tacc_bs,
    cfo = ifelse(!is.na(cfo_cf), cfo_cf, cfo_bs),
    sales = sale/avg_at,
    ret = ((prccd/ajexdi) * trfd)/mleadlag((prccd/ajexdi) * trfd, -1, fyear) - 1 
  ) %>%
  select(
    gvkey, fyear, fyr, conm, loc, sic, avg_at, sales, sales_gr, ear, cfo, tacc, 
    ret, cfo_cf, cfo_bs, tacc_cf, tacc_bs, 
    at, lag_at, sale, lag_sale, ib, oancf
  ) 

# Deal with 120 cases where we have duplicate gvkey, fyear observations 
# Choose the series with more firm observations, grouped by fye

smp_raw_distinct <- smp_raw %>%
  filter(!is.na(avg_at)) %>%
  group_by(gvkey, fyr) %>%
  mutate(nobs = n()) %>%
  group_by(gvkey, fyear) %>%
  filter(n() ==  1 | nobs == max(nobs)) %>%
  filter(row_number() == 1) %>%
  ungroup()

smp_raw_distinct %>%
  group_by(gvkey, fyear) %>%
  filter(n() > 1) -> dup_cases

if (nrow(dup_cases) > 0) stop(
  "Distinct sample has duplicate observations data at the gvkey,",
  "fyear level. This should not happen."
)

ff12 <- read_csv(
  "data/external/fama_french_12_industries.csv", col_types = "cc"
)

smp <- smp_raw_distinct %>%
  left_join(ff12, by = "sic") %>%
  mutate(country = countrycode::countrycode(loc, "iso3c", "country.name")) %>% 
  filter_at(vars(ear, cfo, ret, sales, sales_gr), all_vars(!is.na(.))) %>%
  mutate(ln_ta = log(avg_at)) %>%
  select(
    loc, country, gvkey, conm, sic, ff12_ind, fyear, 
    ln_ta, sales, sales_gr, ear, cfo, tacc, ret,
    cfo_cf, cfo_bs, tacc_cf, tacc_bs, 
    at, lag_at, sale, lag_sale, ib, oancf
  ) %>%
  arrange(loc, gvkey, fyear)

smp %>%
  group_by(gvkey, fyear) %>%
  filter(n() > 1) -> dup_cases

if (nrow(dup_cases) > 0) stop(
  "Sample has duplicate observations data at the gvkey,",
  "fear level. Duplicates are stored in 'dup_cases'"
)

if (sample_type == "World") {
  smp <- smp %>% 
    filter(fyear < 2020, !is.na(ff12_ind), ff12_ind != "Finance")
} else if (sample_type == "Eurozone") {
  smp <- smp %>% 
    filter(
      fyear >= 2005 & fyear < 2020, !is.na(ff12_ind), ff12_ind != "Finance"
    )
}

saveRDS(smp, "data/generated/sample.rds")
foreign::write.dta(smp,"data/generated/sample.dta")
message("done")

