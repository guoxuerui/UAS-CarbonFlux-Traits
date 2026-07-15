#' Resample speclib object to Micasense RedEdge M bands, including transformation to data.frame
#' 
#' @param spec Speclib-object with full radiometric resolution
#' 
#' @import hsdar
#' 
#' @author Xuerui

.resample_mica_fact <- function() {
  
  mica<- data.frame(read.table("Mica_center.csv", col.names = "center"))
  mica$center<- c(475,560,668,840,717)
  mica$fwhm <- c(20,20,10,40,10)
  
  
  function(spec) {
    # spectra for this batch of parameters pars
    spec<-simulate_spectra
    spec_df <-spectralResampling(spec, sensor = mica)@spectra@spectra_ma %>%
      as.data.frame()
    # name the columns = spectral bands
    names(spec_df) <- paste0("B", 1:nrow(mica))
    
    spec_df 
  }
}

resample_mica<- .resample_mica_fact()
