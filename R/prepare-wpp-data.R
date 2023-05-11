library(readxl)

wpp_col_name = c("Index",   "Variant", "Location", "Notes", "Loc_Code", "ISO_Alpha_3", "ISO_Alpha_2", "SDMX_Code", "Type", "Parent", "Year", 0:100)
wpp_col_type = c("numeric", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", rep("text", 101))

wpp_pop_xlsx = "data/WPP2022_POP_F01_2_POPULATION_SINGLE_AGE_MALE.xlsx"
wpp_pop_data = dplyr::bind_rows(data.frame(read_excel(wpp_pop_xlsx, sheet="Estimates",      skip=16, col_types=wpp_col_type)),
                                data.frame(read_excel(wpp_pop_xlsx, sheet="Medium variant", skip=16, col_types=wpp_col_type)))
colnames(wpp_pop_data) = wpp_col_name
wpp_pop_data = wpp_pop_data[!is.na(wpp_pop_data$ISO_Alpha_3),] # restrict to national estimates
wpp_pop_long = reshape2::melt(wpp_pop_data, id.vars=c("ISO_Alpha_3", "Year"), measure.vars=sprintf("%d", 0:100), variable.name="Age", value.name="Value")
wpp_pop_long$Value = 1000 * as.numeric(wpp_pop_long$Value)
wpp_pop_long$Age = as.numeric(as.character(wpp_pop_long$Age))

mc_pop_data = wpp_pop_long
usethis::use_data(mc_pop_data, overwrite=TRUE)
