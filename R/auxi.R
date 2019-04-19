#Ler parametros
definir_parametros = function(char,data){
  termos = char %>%
    stringr::str_replace_all(" ","") %>%
    stringr::str_replace_all(stringr::fixed(")"),"+") %>%
    stringr::str_replace_all(stringr::fixed("("),"+") %>%
    stringr::str_split("[(*),(+),(\\-),(^),(\\/)]") %>%
    unlist
  excluir=c("log","exp","sqrt","sen","cos","")
  covar = termos[termos %in% names(data)] %>% unique
  parametros = termos[!termos %in% c(names(data),excluir)]
  parametros=suppressWarnings(parametros[is.na(as.numeric(parametros))] %>% unique)
  return(list(covar=covar,parametros=parametros))
}

estruturar = function(char,lista,asArray=F){
  c0=list()
  c0$ori=c("log","exp","sqrt","sen","cos") %>% {.[order(nchar(.),decreasing = T)]}
  c0$int=paste0("|",1:length(c0$ori))
  c1=list()
  c1$ori=lista[[1]] %>% {.[order(nchar(.),decreasing = T)]}
  c1$int=paste0("~",1:length(c1$ori))
  c1$final=paste0("data$",c1$ori)
  c2=list()
  if(asArray==T) indices=lista[[2]] %>% {order(nchar(.),decreasing = T)}
  c2$ori=lista[[2]] %>% {.[order(nchar(.),decreasing = T)]}
  c2$int=paste0("¨",1:length(c2$ori))
  c2$final=paste0("theta$",c2$ori)
  if(asArray==T) c2$final=paste0("theta[",indices,"]")
  k_max=max(length(c0$ori),length(c1$ori),length(c2$ori))

  for(k in 1:k_max){
    if(k<=length(c0$ori))
      char=stringr::str_replace_all(char,c0$ori[k],c0$int[k])
    if(k<=length(c1$ori))
      char=stringr::str_replace_all(char,c1$ori[k],c1$int[k])
    if(k<=length(c2$ori))
      char=stringr::str_replace_all(char,c2$ori[k],c2$int[k])}

  for(k in k_max:1){
    if(k<=length(c1$ori))
      char=stringr::str_replace_all(char,stringr::fixed(c1$int[k]),c1$final[k])
    if(k<=length(c2$ori))
      char=stringr::str_replace_all(char,stringr::fixed(c2$int[k]),c2$final[k])
    if(k<=length(c0$ori))
      char=stringr::str_replace_all(char,stringr::fixed(c0$int[k]),c0$ori[k])}

  return(paste0("(",char,")"))  }

#' Pipe operator
#'
#' See \code{\link[magrittr]{\%>\%}} for more details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL
