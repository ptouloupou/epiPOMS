# package start up message

.onAttach <- function(lib, pkg) {
    info <- packageDescription("epiPOMS")
    packageStartupMessage("Welcome to epiPOMS!\n", paste("Version ", info$Version, ", created on ", 
        info$Date, ".", "\n", sep = ""))
}
