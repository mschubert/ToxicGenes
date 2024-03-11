library(dplyr)
library(ggplot2)

# set SYNAPSE_USERNAME, SYNAPSE_PASSWORD in env vars
# request access to the 'genie nsclc'/'genie crc' on synapse
nsclc = genieBPC::pull_data_synapse("NSCLC", version="v2.0-public")[[1]]
#crc = genieBPC::pull_data_synapse("CRC", version="v2.0-public")

# no rbm14 samples
cna = nsclc$cna |>
    filter(Hugo_Symbol %in% c("CCND1", "RBM14")) |>
    tidyr::pivot_longer(-Hugo_Symbol, names_to="sample", values_to="copy")

clin = nsclc$pt_char |> # prob get this from ca_dx_index [stage_dx, dob_ca_dx_days/mos/yrs]
    transmute(sex = factor(naaccr_sex_code),
              age = age_last_fu_yrs,
              status = as.integer(!is.na(age_death_yrs)),
              time = ??)
