#' @title  Normal regression model with general parametrization
#'
#' @description Determine the estimates of the parameters of a
#' general parametrization model with a Newthon Raphson method
#'
#'
#' @param formula a nonlinear model formula including variables and parameters
#' @param formula_var a nonlinear model formula for the diagonal of covariance matrix
#' @param data A data frame in which to evaluate the variables in \code{formula} and \code{formula_var}.
#' Can also be a list or an environment, but not a matrix
#' @param theta_ini A named list of starting estimates. When \code{start} is missing, a very cheap guess for \code{start} is tried and parameters names are automatically identified
#' @param control A list of control parameters. See 'Details'
#' @param bias_correction Logical. Should a bias correction be estimated for the parameters?
#'
#'
#' @details The control argument is a list that can supply any of the following components:
#'  \itemize{
#'  \item max_it Maximum number of iterations
#'  \item parada Maximum difference between iterations to consider convergence
#'}
#'
#'  @examples
#' \dontrun{
#'
#' fit <- reg_geral(y~alfa+X1^(gama)+beta*X2,~100+sigma*X3,data=data,
#' theta_ini=list(alfa=100,beta=1,gama=0.5,sigma=0.1),bias_correction = T)
#'
#' }
#'
#' @export


reg_geral=function(formula=NULL,
                   formula_var=NULL,
                   data,
                   theta_ini=NULL,
                   max_it=500,
                   method="NR",
                   parada=0.0001,
                   bias_correction=T,
                   prop=1,
                   verbose=1,
                   simular=F){

  if(is.null(formula)){
    stop('arguments "formula" is missing, with no default',call. = F)}
  if(is.null(formula_var)){
    stop('arguments "formula_var" is missing, with no default',call. = F)}
  if(!inherits(formula,"formula"))
      stop("formula has a indefinite format",call. = F)
  if(!inherits(formula_var,"formula"))
      stop("formula_var has a indefinite format",call. = F)
  if(!is.null(theta_ini)){
    if(!is.list(theta_ini))
      stop("theta_ini must be a list")
    theta=names(theta_ini)}

    resposta=data[,as.character(formula)[2]]
    mu=as.character(formula)[3]
    S=tail(as.character(formula_var),1)

  cl <- match.call()
  n=nrow(data)
  q=1
  par=definir_parametros(c(mu,S),data)
  covar=par[[1]]
  if(is.null(theta_ini)){
    theta=par[[2]]
  }else{
    par[[2]]=theta
  }

  if(verbose!=0) cat("Parametros do modelo:",paste0(theta,collapse=", "),"\n")
  if(verbose!=0) cat("Covariaveis do modelo:",paste0(covar,collapse=","),"\n")

  l=paste0("-0.5*log(",S,") -0.5*(",mu,"-resposta)^2/(",S,")")

  Di=array()
  Vi=array()
  Pd=array()
  ld=array()
  a=matrix(NA,nrow=length(theta),ncol=length(theta))
  C=matrix(NA,nrow=length(theta),ncol=length(theta))
  for(k in 1:length(theta)){
    Di[k]=as.character(Deriv::Deriv(mu,theta[k]))
    Vi[k]=as.character(Deriv::Deriv(S,theta[k]))
    ld[k]=as.character(Deriv::Deriv(l,theta[k]))
    Pd[k]=as.character(Deriv::Deriv(paste0("((",S,")^0.5)"),theta[k]))
    for(j in 1:length(theta)){
      a[k,j]=as.character(Deriv::Deriv(Di[k],theta[j]))
      C[k,j]=as.character(Deriv::Deriv(Vi[k],theta[j]))}
  }

  Fi=rbind(Di,Vi)
  aC=rbind(c(a),c(C))

  eval(parse(text=paste0(
    "gerar_l = function(theta, data,resposta){
    N=nrow(data)
    -0.5*sum(log(rep(1,N)*",estruturar(S,par,asArray=T),") +
    ((rep(1,N)*",estruturar(mu,par,asArray=T),")-resposta)^2/(rep(1,N)*",estruturar(S,par,asArray=T),"))}"
    )))

  eval(parse(text=paste0(
    "log_like=function(theta, data){
    sum(dnorm(c(",paste(resposta,collapse=","),"),",
    paste("rep(1,",length(resposta),")*",estruturar(mu,par,asArray=T),collapse = ","),
    ",sqrt(",paste("rep(1,",length(resposta),")*",estruturar(S,par,asArray=T),collapse = ","),"),log = T))
    }"
  )))


  eval(parse(text=paste0(
    "gerar_mu = function(theta, data){
    N=nrow(data)
    c(",
    paste("rep(1,N)*",estruturar(mu,par),collapse = ",")
    ,")}"
  )))

  eval(parse(text=paste0(
    "gerar_sigma = function(theta, data){
    N=nrow(data)
    Matrix::.sparseDiagonal(x=c(",
    paste("rep(1,N)*",estruturar(S,par),collapse = ","),
    "),n=N)}"
  )))

  eval(parse(text=paste0(
    "gerar_D = function(theta, data){
    N=nrow(data)
    Matrix::Matrix(c(",
    paste("rep(1,N)*",estruturar(t(Di),par),collapse = ","),
    "), ncol=N,byrow=T) %>% as.vector %>% Matrix::Matrix(ncol=",
    length(theta),", sparse=T,byrow=T)}"
  )))

  eval(parse(text=paste0(
    "gerar_V = function(theta, data){
    N=nrow(data)
    Matrix::Matrix(c(",
    paste("rep(1,N)*",estruturar(t(Vi),par),collapse = ","),
    "), ncol=N,byrow=T) %>% as.vector %>% Matrix::Matrix(ncol=",
    length(theta),", sparse=T,byrow=T)}"
  )))

  gerar_J = function(D){
    dt=dplyr::tibble(valor=as.vector(D),
                     id=rep(1:nrow(D),ncol(D)),
                     par=rep(1:ncol(D),each=nrow(D)))
    dt2=dplyr::left_join(dt,dt,by="id") %>%
      dplyr::mutate(total=2*valor.x*valor.y,
                    aux=1) %>%
      dplyr::select(-valor.x,-valor.y)
    dplyr::bind_rows(dt2,dt2 %>% dplyr::mutate(aux=0, total=0)) %>%
      dplyr::arrange(par.x,par.y,id,aux) %>%
      .$total %>%
      Matrix::Matrix(nrow=2*nrow(D),sparse=T)
  }

  eval(parse(text=paste0(
    "gerar_G = function(theta, data){
    N=nrow(data)
    Matrix::Matrix(c(",
    paste("rep(1,N)*",estruturar(t(aC),par),collapse = ","),
    "), ncol=N,byrow=T) %>% as.vector %>% Matrix::Matrix(ncol=",
    ncol(aC),", sparse=T,byrow=T)}"
  )))

  eval(parse(text=paste0(
    "gerar_F = function(theta, data){
    N=nrow(data)
    Matrix::Matrix(c(",
    paste("rep(1,N)*",estruturar(Matrix::t(Fi),par),collapse = ","),
    "), ncol=N,byrow=T) %>% as.vector %>% Matrix::Matrix(ncol=",
    ncol(Fi),", sparse=T,byrow=T)}"
  )))

  correcao_vies=function(){

    D=gerar_D(theta_val,data)
    G=gerar_G(theta_val,data)
    J=gerar_J(D)

    phi=-0.5*(G+J)

    K=solve(Matrix::t(Fn)%*%Hn%*%Fn,tol=1e-2000)
    e=phi%*%as.vector(K)

    K%*%Matrix::t(Fn)%*%Hn%*%e}

  gerar_H = function(sigma){
    Matrix::diag(sigma) %>%
    {rbind(1/.,0.5*(.^-2))} %>%
      as.vector() %>%
      {Matrix::.sparseDiagonal(x=.,n=length(.))}
  }

  # gerar_T = function(theta,data){
  #   z=resposta-gerar_mu(theta,data)
  #   sigma=gerar_sigma(theta,data)
  #   D=gerar_D(theta,data)
  #   V=gerar_V(theta,data)
  #   D+z*V/sigma
  # }


  gerar_s=function(media,sigma,theta,p=prop,all=T){
    res=resposta-media
    s2=-(Matrix::diag(sigma)-res^2)
    s=c(rbind(res,s2))
    if(all==F) return(s)
    S=Fn %*% unlist(theta) + s*p
    return(S)}

  #skovgaard:
  # eval(parse(text=paste0(
  #   "gerar_Pd = function(theta, data){
  #   N=nrow(data)
  #   Matrix::Matrix(c(",
  #   paste("rep(1,N)*",estruturar(Pd,par),collapse = ",")
  #   ,"), nrow=N, sparse=T)}"
  # )))

  eval(parse(text=paste0(
    "gerar_C = function(theta, data){
    N=nrow(data)
    Matrix::Matrix(c(",
    paste("rep(1,N)*",estruturar(c(C),par),collapse = ",")
    ,"), nrow=N, sparse=T)}"
  )))

  eval(parse(text=paste0(
    "gerar_D2 = function(theta, data){
    N=nrow(data)
    Matrix::Matrix(c(",
    paste("rep(1,N)*",estruturar(c(a),par),collapse = ",")
    ,"), nrow=N, sparse=T)}"
  )))

  gerar_M2=function(media,sigma){
    z=resposta-media
    c(-Matrix::diag(sigma),
           -2*z*Matrix::diag(sigma),
           -2*z*Matrix::diag(sigma),
           -4*(z^2)*Matrix::diag(sigma)+2*(Matrix::diag(sigma))^2) %>%
      matrix(ncol=length(z),byrow = TRUE) %>%
      apply(2,function(x) Matrix::Matrix(x,ncol=2,byrow=TRUE)) %>%
      Matrix::bdiag()
  }


  interacao=list()
  if(!is.null(theta_ini))
    theta_ini= unlist(theta_ini,use.names = F)
  if(is.null(theta_ini))
    theta_ini= rep(1,length(theta))
  theta_val = data.frame(matrix(theta_ini,nrow=1,dimnames = list(NULL,theta)))

  if(method=="optim"){
    logit=function(p){
      p=0.5*(p+1)
      p=log(p/(1-p))
      ifelse(is.infinite(p),NA,p)}

    opt = pso::psoptim(rep(NA,length(theta)),
                        fn=function(par){
                          aux=-log_like(logit(par),data)
                          if(is.na(aux)) aux=Inf
                          return(aux)},
                        control=list(maxit=max_it))

    k=opt$counts[2]
    theta_val_novo = logit(opt$par)
    names(theta_val_novo)=theta
    theta_val=data.frame(valor=logit(opt$par),nome=theta) %>% tidyr::spread(nome,valor)
  }

  if(method=="NR"){
    k=1
    dif=1
  while(dif>parada & k <= max_it){

    if(verbose>1) cat(paste0("Iteração ",k,"\n"))

    sigma=gerar_sigma(theta_val,data)
    media=gerar_mu(theta_val,data)
    Fn=gerar_F(theta_val,data)
    Hn=gerar_H(sigma)
    Sn=gerar_s(media,sigma,theta_val)

    theta_val_novo = (solve(Matrix::t(Fn)%*%Hn%*%Fn,tol=1e-2000) %*% Matrix::t(Fn)%*%Hn%*%Sn)[,1]
    names(theta_val_novo) = theta

    interacao[[k]]=theta_val_novo
    k=k+1
    dif=max(abs(theta_val-theta_val_novo)/abs(theta_val))
    theta_val[1,]=theta_val_novo

    if(is.na(dif) | dif>10000000) stop("Não convergiu",call. = F)

  }
  if(k >= max_it | is.na(dif)) stop("Não convergiu após o número máximo de iterações",call. = F)
  k=k-1
  }


  sigma=gerar_sigma(theta_val,data)
  media=gerar_mu(theta_val,data)
  Fn=gerar_F(theta_val,data)
  Hn=gerar_H(sigma)
  Sn=gerar_s(media,sigma,theta_val)


  if(simular==F){
  if(bias_correction){
    bias=correcao_vies()
    theta_val[1,]=theta_val_novo-bias[,1]
    theta_val_novo=unlist(theta_val[1,])
    media=gerar_mu(theta_val,data)
    sigma=gerar_sigma(theta_val,data)
    Fn=gerar_F(theta_val,data)
    Hn=gerar_H(sigma)
    Sn=gerar_s(media,sigma,theta_val)
  }}

  corrigido=NA
  if(simular==T){
    bias=correcao_vies()
    corrigido=theta_val_novo-bias[,1]
    corrigido=unlist(corrigido)
  }


  #Information Obs


  G=gerar_G(theta_val,data)
  M2 = gerar_M2(media,sigma)

  obs = Matrix::t(Fn) %*% Hn %*% M2 %*% Hn %*% Fn +
    matrix(t(gerar_s(media,sigma,theta_val,all = F)) %*% Hn %*% G,ncol=length(theta))



  out <- list(
    parameters=theta_val_novo,
    corrigido=corrigido,
    fitted.values = media,
    residuals = resposta-media,
    var = sigma,
    iteracao = k-1,
    score = Matrix::t(Fn)%*%Hn%*%gerar_s(media,sigma,theta_val,all=F),
    fisher = Matrix::t(Fn)%*%Hn%*%Fn,
    fisher_obs = -1*obs,
    loglike = sum(dnorm(resposta,media,sqrt(Matrix::diag(sigma)),log = T)),
    target=resposta,
    funcao_loglike=log_like,
    funcao_mu = gerar_mu,
    funcao_sigma = gerar_sigma,
    funcao_l=gerar_l,
    funcao_F=gerar_F,
    funcao_H=gerar_H,
    funcao_s=gerar_s,
    funcao_D=gerar_D,
    funcao_V=gerar_V,
    funcao_G=gerar_G,
    funcao_C=gerar_C,
    funcao_D2=gerar_D2,
    #matrizes = list(F=Fn,H=Hn,S=Sn),
    data = data,
    call = cl
  )
  class(out) <- "genReg"

  return(out)
}

#' @export
setClass("genReg")


#' @method print genReg
#' @export
print.genReg <- function(x, ...){
  cat("\nCall:\n")
  print(x$call)

  ndec=stringr::str_length(stringr::str_extract(as.character(x$parameters), "\\.[0-9]*")) - 1

  nround=(min(5,max(ndec,na.rm=T)))
  cat("\n\nCoefficients:\n")
  print(round(x$parameters,nround))
  cat("\n")
}

#' @method summary genReg
#' @export
summary.genReg <- function(x, ...){
  cat("\nCall:\n")
  print(x$call)

  ndec=stringr::str_length(stringr::str_extract(as.character(x$parameters), "\\.[0-9]*")) - 1
  nround=(min(5,max(ndec,na.rm=T)))
  cat("\n\nResiduals:\n")
  print(summary(x$res)[-4])
  cat("\n\nCoefficients:\n")
  
  df=length(x$target)-length(x$parameters)+1
  coef=x$parameters
  se=x$fisher_obs %>% solve %>% diag %>% sqrt
  tv=coef/se
  pv=pt(2*abs(tv),df,lower.tail=F)
  #Pr(>|t|)
  print(data.frame(Estimate=round(coef,nround),Std.Error=round(se,nround), t.value=tv, Pv=pv,
             row.names = names(x$parameters)))
  cat("\n")
}



#' @method coef genReg
#' @export
coef.genReg <- function(x, ...){
  x$parameters
}

#' @method predict genReg
#' @export
predict.genReg <- function(x, newdata=NULL, type = "mean", ...){
  if(!type %in% c("mean","var")) stop("type must be 'mean' or 'var'",call. = FALSE)
  parameters=sapply(x$parameters,list)
  if(type == "mean"){
    if(is.null(newdata)) return(x$fitted.values)
    return(x$funcao_mu(parameters,newdata))}
  if(type == "var")
    if(is.null(newdata)) return(Matrix::diag(x$var))
  return(Matrix::diag(x$funcao_sigma(parameters,newdata)))
}

#' @export
likelihood_ratio <- function(x, parameters,correction=FALSE){
  if(class(x)!="genReg")
    stop("x must be a genReg class",call. = F)

  loglike=function(y,media,var){
    sum(dnorm(y,mean=media,sqrt(var),log=T))}
  theta=coef(x) %>% data.frame(nome=names(.),valor=.) %>% tidyr::spread(nome,valor)
  par_teste = dplyr::bind_rows(theta,parameters) %>% dplyr::summarise_all(funs(tail(na.omit(.),1)))
  p=length(theta)
  q=length(parameters)
  data=x$data

  mu=x$fitted.values
  var=Matrix::diag(x$var)
  mu0=x$funcao_mu(par_teste,x$data)
  var0=Matrix::diag(x$funcao_sigma(par_teste,x$data))
  y=x$target

  l1=loglike(y,mu,var)
  l0=loglike(y,mu0,var0)

  LR=2*(l1-l0)
  pv=pchisq(LR,df = q,lower.tail = F)

  if(correction==FALSE){return(list(LR=LR,p_value=pv))}

  ####Skovgaard correction:
  z=x$target-x$fitted.values
  z0=z+mu-mu0
  u0=(z0^2)/var0
  P=sqrt(var)
  P0=Matrix::diag(x$funcao_sigma(par_teste,x$data))^0.5
  Pd=x$funcao_Pd(theta,data)
  Pd0=x$funcao_Pd(par_teste,data)
  a=(y-x$fitted.values)/P
  u00=(a^2)*(P0^2)/var0

  T=x$funcao_D(theta,data)+z*x$funcao_V(theta,data)/var
  T00=x$funcao_D(par_teste,data)+a*P0*x$funcao_V(par_teste,data)/var0

  B=-z*x$funcao_D(theta,data) -0.5*x$funcao_V(theta,data)
  B00=-a*P0*x$funcao_D(par_teste,data) -0.5*x$funcao_V(par_teste,data)
  A=-x$funcao_V(theta,data)/var^2
  A0=-x$funcao_V(par_teste,data)/var0^2

  aux1=rep(1:length(theta),length(theta))
  aux2=rep(1:length(theta),each=length(theta))
  Ad = -2*A[,aux1]*x$funcao_V(theta,data)[,aux2]/var
  Ad = Ad - x$funcao_C(theta,data)/var^2
  Ad0 = -2*A0[,aux1]*x$funcao_V(par_teste,data)[,aux2]/var0
  Ad0 = Ad0 - x$funcao_C(par_teste,data)/var0^2

  E = -0.5*(Ad*(var-z^2)) - x$funcao_D2(theta,data)*z/var
  E00 = -0.5*(Ad0*(var0-(a*P0)^2)) - x$funcao_D2(par_teste,data)*a*P0/var0
  G = apply(B[,aux1]*A[,aux2] + E,2,sum) %>% matrix(nrow=length(theta))
  G00 = apply(B00[,aux1]*A0[,aux2] + E00,2,sum) %>% matrix(nrow=length(theta))

  R = x$funcao_Pd(theta,data)*a + x$funcao_D(theta,data)

  Q = x$funcao_D(theta,data)+z*x$funcao_V(theta,data)/var
  Q00 = x$funcao_D(par_teste,data)+a*P0*x$funcao_V(par_teste,data)/var0
  #Q00 = x$funcao_D(par_teste,data)+z0*x$funcao_V(par_teste,data)/var0


  Ud0 = Matrix::t(Q00)%*%solve(x$funcao_sigma(par_teste,data))%*%R
  U0= Matrix::t(x$funcao_F(par_teste,data)) %*% x$funcao_H(x$funcao_sigma(par_teste,data)) %*% x$funcao_s(mu0,x$funcao_sigma(par_teste,data),par_teste,all = F)
  J00 = Matrix::t(T00)%*%solve(x$funcao_sigma(par_teste,data))%*%x$funcao_D(par_teste,data) + G00
  J = Matrix::t(T)%*%solve(x$var)%*%x$funcao_D(theta,data) + G

  ld = Matrix::t(R)%*%solve(x$var)%*%(-z)
  ld0 = Matrix::t(R)%*%solve(x$funcao_sigma(par_teste,data))%*%(-z0)

  w=which(!names(theta) %in% names(parameters))

  J00=as.matrix(J00)
  p1=det(J)^0.5*det(Ud0)^-1*det(J[w,w])^0.5*det(J00[w,w])^-0.5*det(J00)^0.5
  p2=(t(U0)%*%solve(J00,tol=1e-200) %*% U0)^(p/2)
  p3=LR^(q/2-1)
  p4=t(ld-ld0)%*%Matrix::solve(Ud0,tol=1e-20)%*%U0
  rho=p1*p2/(p3*p4)
  LR2=LR-2*log(rho)[1,1]

  pv2=tryCatch(pchisq(LR2,df = q,lower.tail = F), error=function(e)NA)
  erro=0 + (p1<0) + 2*(p1<0) + 3*(p1<0) + 4*(p4<0)

  out = list(LR=LR,p_value=pv, LR_correction=LR2, p_value_correction=pv2,erro=erro)
  class(out) = "lr_ratio"
  return(out)
}
