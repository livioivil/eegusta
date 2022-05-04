#' @title Converts eeg_data or eeg_epochs (eegUtils) to eeg_lst (eeguana) object
#' @description Builds an eeg_lst object composed of two data.table objects and a tibble one. All three are linked by a unique identifier .id. Amplitude values and timestamps appear in the signal table. Triggers, blinks, artifact rejection markings, and other events logged by the EEG recording software appear in the events table. Segment information and recording IDs appear in the segments tibble.
#' @usage eegUtils2eeguana(data, reference, unit, resolution)
#' @param data eegUtils data
#' @param reference reference
#' @param unit unit
#' @param resolution resolution
#' @author Angela Andreella
#' @return Returns an eeg_lst (eeguana package)
#' @export
#' @importFrom eeguana sample_int
#' @importFrom data.table setnames
#' @importFrom data.table setcolorder
#' @importFrom data.table setattr
#' @importFrom data.table setkey
#' @importFrom data.table data.table
#' @importFrom dplyr tibble 
#' @importFrom eeguana eeg_lst
#' @importFrom eeguana channel_dbl
#' @importFrom data.table as.data.table
#' 
eegUtils2eeguana <- function(data, reference = "", unit = "microvolt", resolution = 1){
  
  .sample <- sample_int(data$timings$sample,sampling_rate = data$srate)
  
  .id <- as.integer(data$timings$epoch)
  
  nch <- length(data$chan_info$electrode)
  ch_n <- sapply(c(1:nch), function(x) paste0("Ch",x))
  
  
  set_chnl<-function(i){
    channel_dbl(
      values = as.numeric(unlist(data$signals[,i])),
      .x = data$chan_info$cart_x[i],
      .y = data$chan_info$cart_y[i],
      .z = data$chan_info$cart_z[i],
      reference = reference,
      resolution = resolution,
      number = ch_n[i],
      unit = unit,
      radius = ifelse(!is.null(data$chan_info$radius[i]),data$chan_info$radius[i], NA),
      theta = ifelse(!is.null(data$chan_info$theta[i]),data$chan_info$radius[i], NA),
      phi = ifelse(!is.null(data$chan_info$phi[i]),data$chan_info$radius[i], NA)
    )}
  
  
  .signal=lapply(1:length(data$signals),set_chnl)
  names(.signal)=names(data$signals)
  
  .signal<- as.data.table(.signal)
  .signal <- cbind(.signal,.id,.sample)
  
  # ####signal tbl
  
  # The signal table is organised into columns representing timestamps
  # .sample and individual electrodes. Each .sample corresponds to
  # 1 sample in the original recording, i.e. if the sampling rate of the EEG
  # recording is 500 Hz, then each .sample corresponds to 2 milliseconds.
  # These timestamps correspond to .initial in the events table, which
  # displays only the timestamps where logged events began. 
  #TODO!!!!!
  #.signal <- data.table(.id,.sample,sapply(nam, function(x) assign(x,get(x))))
  # .signal
  .signal[, .id := .id][, .sample := .sample]
  setnames(.signal, make_names(colnames(.signal)))
  setcolorder(.signal, c(".id", ".sample"))
  setattr(.signal, "class", c("signal_tbl", class(.signal)))
  setkey(.signal, .id, .sample)
  
  ####.events tbl
  ###.id Integers indicating to which group the row of the signal matrix belongs.
  ## .initial = dati$events$event_onset
  ## .final =  as_sample_int(dati$events$event_onset + dati$events$event_time, sampling_rate =  dati$srate, unit = "s") - 1L,
  ## .channel = dati$chan_info$electrode
  ## descriptions_dt = data.table::data.table()
  
  .events <- new_events_tbl (
    .id = data$events$epoch,
    #.initial = dati$events$event_onset %>%
    #  as_sample_int(sampling_rate = dati$srate, unit = "s"),
    .initial = data$events$event_onset,
    #.initial =  init_events,
    #.final =  init_events,
    #.final =  round(dati$events$event_time* dati$srate) %>% as.integer() +init_events,
    .final = data$events$event_onset + data$srate - 1L,
    #.final = as_sample_int(dati$events$event_onset, 
    #                      sampling_rate =  dati$srate, unit = "s"),
    #.final = rep(1,length(dati$events$epoch)) %>%
    #  as_sample_int(sampling_rate = dati$srate, unit = "s"), 
    .channel = NA_character_,
    descriptions_dt = data.table(
      # .recording = dati$epochs$recording, 
      .type = rep("Stimulus",length(data$events$event_type)),
      .description = paste0(data$events$event_type)))
  #                                          .subj = dati$events$subj))
  
  
  setattr(.events$.initial,"sampling_rate",data$srate)
  setattr(.events$.initial,"class","sample_int")
  setattr(.events$.final,"sampling_rate",data$srate)
  setattr(.events$.final,"class","sample_int")
  setkey(.events, .id)
  #segments_tbl <- dplyr::tibble(.id =  dati$epochs$epoch, .recording = dati$epochs$recording)
  segments_tbl <- tibble(data$epochs)

  setnames(segments_tbl,gsub("^epoch$",".id",names(segments_tbl)))
  setnames(segments_tbl,gsub("^recording$",".recording",names(segments_tbl)))
  segments_tbl$segment=1:nrow(segments_tbl)


  #type = rep("Stimulus", length(dati$epochs$recording)))
  
  #segments_tbl <- validate_segments(segments_tbl)
  
  data <-   eeg_lst(
    signal_tbl = .signal,
    events_tbl = .events,
    segments_tbl = segments_tbl)
  
  return(data)
}