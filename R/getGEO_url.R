#' Replace getGEO function
#'
#' @param GSE GSE number of GEO database
#' @return GSE dataset download address
#' @keywords Replace getGEO function
#' @export


getGEO_url = function(GSE){
  gse=GSE   # 需要更改GSE
  mirror='tencent'

  gse=toupper(gse)

  down=ifelse(as.numeric(gsub('GSE','',gse))<1000,
              paste0('/GEOmirror/GSEnnn/',gse,
                     '_eSet.Rdata'),
              paste0('/GEOmirror/',
                     gsub('[0-9][0-9][0-9]$','nnn',gse),'/',gse,
                     '_eSet.Rdata'))
  if(mirror=='tencent'){
    up='http://49.235.27.111'
  }


  tpf=paste0(gse, '_eSet.Rdata')

  result = paste0(up,down)
  return(result)
}

