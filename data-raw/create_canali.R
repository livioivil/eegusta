elec_position <- eegUtils::import_chans("data-raw/canali.ced")

usethis::use_data(elec_position, overwrite = TRUE, compress = "xz")
