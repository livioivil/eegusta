#' Reads channel locations in elp format
#' @export
#' @param file file's name
#' @return returns a tibble object with columns \code{.channel,radius,theta,phi,.x,.y,.z}.
#' 
# There are standards but the electrodes are not always in the standard 10-20 position. I created a dataset with those locations for 32 electrodes, layout_32_1020,  (which was what I needed so far). In the intro vignette, they weren't inside the vhdr (sometimes they are) so I had to load them:
# https://bnicenboim.github.io/eeguana/articles/intro.html
# 
# channels_tbl(faces) <- select(channels_tbl(faces), .channel) %>%
# left_join(layout_32_1020)

read_elp <- function(file){
  loc <- readr::read_csv(file,comment = "//")
  loc=as.data.frame(loc)
  
  loc=loc[-(1:(grep("%S",loc[,1])[1]-1)),]
  loc=matrix(loc,byrow = TRUE,ncol = 3)
#  str(loc)
  
  loc=loc[,2:3]
  
  loc[,1]=gsub("%N","",loc[,1])
  loc[,1]=gsub("\\t","",loc[,1])
  
  temp=matrix(as.numeric(unlist(strsplit(loc[,2],"\\t"))),byrow = TRUE,ncol=3)
  temp=data.frame(temp)
  temp=temp[,c(2,1,3)]
  loc=data.frame(loc[,1],as.numeric(NA),as.numeric(NA),as.numeric(NA),temp)
  names(loc)=c(".channel","radius","theta","phi",".x",".y",".z")
  loc=tibble::as_tibble(loc)
  return(loc)
}
# .channel radius theta   phi    .x    .y    .z
# <chr>     <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
