# Added by J. Fox 2019-09-06

.onAttach <- function(...){
    packageStartupMessage(gettext(paste(c("This package has some major limitations", 
         "(for example, it does not work reliably when",
         "the number of variables exceeds 30),",
         "and has been superseded by the norm2 package."),
         collapse="\n")))
    return()
}
