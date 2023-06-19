#' Reads aplitudes from a text delimited file 
#' @param file file's name
#' @param time_column name of the column that contains the time 
#' @param srate sampling rate. 256 is the default
#' @param chan_info a \code{tibble} with channel locations. It has the following columns \code{electrode, radius, theta, phi, cart_x, cart_y, cart_z, x, y}(\code{NULL} by default)
#' @param .unit if "milliseconds" it divides the time by 1000
#' @export
#' @param file file's name
#' @return returns a tibble object with columns \code{.channel,radius,theta,phi,.x,.y,.z}.
#' @author Livio Finos
##################
read_delim2eeg_data <- function(
    file,
    time_column="time",
    srate=256,
    chan_info=NULL,
    .unit = "milliseconds"){
  
  d=read_tsv(file)
  
  # file=strsplit(file,"/")[[1]]
  # file=file[length(file)]
  # 
  # file=gsub("\\.txt$","",file)
  # 
  # info=str_split(file,"_")[[1]]
  # info=list(info[1],info[2])
  # info[[1]]=as.numeric(info[[1]])
  # 
  timings <- tibble::tibble(sample = 1:nrow(d))
  timings$time <- d[,time_column]
  if(.unit == "milliseconds") timings$time <- timings$time/ 1000
  d[time_column] <- NULL
  sigs <- tibble::as_tibble(d)
  rm(d)
  
  event_table <- tibble::tibble(event_onset = which(timings$time==0), 
                                event_time = 0, 
                                event_type = file)
  
  epochs <- tibble::new_tibble(list(epoch = 1, participant_id = file, 
                                    recording = file), nrow = 1, class = "epoch_info")
  if(is.null(chan_info))
    chan_info=names(sigs)
  D <- eegUtils:::eeg_data(data = sigs, srate = srate, events = event_table, 
                           timings = timings,chan_info = chan_info, epochs = epochs)
  
  
  D <-  epoch_data( D ,
                    events = file,
                    epoch_labels = file,
                    time_lim = range(D$timings$time),
                    baseline = c(min(D$timings$time),0)
  )
  
  D
}
