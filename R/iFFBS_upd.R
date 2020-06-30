#' Update the hidden process of an individual using the iFFBS algorithm.
#'
#' @description
#' This function updates the hidden carriage states of an individual in a group
#' conditional on the states of the remaining individuals using the
#' individual-Forward Filtering Backward Sampling (iFFBS) algorithm. This function
#' is only used in conjunction with \code{\link{epiPOMS_mcmc}}.
#'
#' @param
#' ntypes An integer greater than 1 corresponding the number of different strains
#' in the study, \eqn{n_g}.
#' @param
#' nindsgr An integer greater than 1 corresponding to the number of individuals
#' in the group.
#' @param
#' ind A positive integer corresponding to the index of the individual in the
#' group for which the hidden states will be updated.
#' @param
#' tmaxgr An integer greater than 1 corresponding to the time point that the
#' last sample was obtained from all individuals in the group, i.e. the last
#' pre-specified group observation time.
#' @param
#' tmaxind A positive integer corresponding to the time point that the last
#' sample was collected from individual \code{ind}. Note that, this value is
#' less than the time of the last observation in the group, \code{tmaxgr}, only
#' if the individual withdrawn from the study before its completion.
#' @param
#' curhidprocgr A \code{nindsgr} by \code{tmaxgr} matrix containing
#' the current values of the (true) hidden carriage states of individuals in the
#' specific group at each time point over the group observation period. All
#' elements must be non-negative integers, ranging from 0 to \eqn{n_g} (see
#' details, for more information), or \code{NA}. \code{NA} entries corresponds
#' to missing values, for example, individual dropouts. However, we assume that
#' at most one individual from each group can withdraw before the study period
#' ends.
#' @param
#' obsertest1ind A vector of length \code{tmaxgr} containing the
#' observed results of the first test for individual \code{ind} at each time
#' point during the group observation period. All elements must be non-negative
#' integers, ranging from 0 to \eqn{(n_g + 1)}, or \code{NA}. See
#' \code{\link{info_epiPOMSdata}} and \code{\link{as_epiPOMSdata}} for more
#' information.
#' @param
#' obsertest2ind A vector with the same structure as \code{obsertest1ind},
#' but for the second test.
#' @param
#' curalpha A vector of non-negative real numbers, one for each strain,
#' containing the current values of \eqn{\alpha} parameters, i.e. the
#' strain-specific external colonisation rates. This is a vector of \eqn{n_g}
#' values. The first element corresponds to the first strain, the second element
#' to the second strain and so on.
#' @param
#' curbeta A vector of non-negative real numbers, one for each strain,
#' containing the current values of \eqn{\beta} parameters, i.e. the
#' strain-specific within-group colonisation rates. This is a vector of \eqn{n_g}
#' values. The first element corresponds to the first strain, the second element
#' to the second strain and so on.
#' @param
#' curmu A vector of non-negative real numbers, one for each strain,
#' containing the current values of \eqn{\mu} parameters, i.e. the
#' strain-specific clearance rates. This is a vector of \eqn{n_g}
#' values. The first element corresponds to the first strain, the second element
#' to the second strain and so on.
#' @param
#' curdelta A non-negative real number corresponding to the current value of
#' \eqn{\delta} parameter, i.e. the relative colonisation rate in a
#' carrier versus non-carrier individual.
#' @param
#' curgamma A non-negative real number corresponding to the current value of
#' \eqn{\gamma} parameter, i.e. the relative colonisation rate in smaller
#' versus bigger groups in terms of area (in square meters). If there is no
#' difference between groups, set \code{curgamma = 1}.
#' @param
#' curnu A vector of numbers between 0 and 1, containing the current values of
#' \eqn{\nu} parameters, i.e. the probabilities of carriage at the beginning
#' of the study. This is a vector of \eqn{(n_g + 1)} values. The first element
#' corresponds to the non-carriage state. The second and subsequent elements
#' correspond to the carriage of one of the \eqn{n_g} strains. These
#' probabilities should sum up to one.
#' @param
#' curclassmat A list of two matrices, one for each test, containing
#' the current values of the classification probabilities, obtained using function
#' \code{\link{classprob}}.
#'
#' @details
#' Our approach to analyse partially observed epidemic data, involves assuming
#' that the classifications at an observation time are imperfect measures of an
#' underlying true (hidden) epidemic process. Sampling the hidden carriage process
#' is done by using a Gibbs step via the iFFBS algorithm by
#' \insertCite{Touloupou2019;textual}{epiPOMS}, where the hidden states are
#' updated individually per subject conditionally on the states of the remaining
#' subjects.
#'
#' The \code{\link{iFFBS_upd}} is an auxiliary function that uses the iFFBS
#' algorithm to update the hidden states of individual \code{ind} conditional
#' on the current values of the remaining individuals, the model parameters and
#' the observed data. The possible states include being a non-carrier
#' (state 0) or being a carrier of one of the \eqn{n_g} strains (states
#' \eqn{1, 2, \ldots, n_g}), where 1 to \eqn{(n_g - 1)} refers to carriage of
#' one of the main strains and \eqn{n_g} to carriage of the remaining strains
#' (pooled group).
#'
#' A point which is worth emphasising is that individuals withdrawn from the
#' study are removed from the model on their times of drop-out. Therefore,
#' carriage states from these individuals after the drop-out time are not
#' imputed since they do not play further role in the spread of the epidemic.
#'
#' @return
#' A vector of length \code{tmaxgr} containing the updated carriage states
#' of individual \code{ind}.
#'
#' @author
#' Panayiota Touloupou, Simon Spencer
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso
#' \code{\link{epiPOMS_mcmc}} for performing inference in partially
#' observed multi-strain epidemic models.
#'
#' @export
iFFBS_upd <- function(ntypes, nindsgr, ind, tmaxgr, tmaxind, curhidprocgr, obsertest1ind, obsertest2ind, 
    curalpha, curbeta, curmu, curdelta, curgamma, curnu, curclassmat) {
    
    CP1 <- matrix(0, tmaxind, ntypes + 1)
    CP2 <- matrix(0, tmaxind, ntypes + 1)
    CP2rest <- rep(0, ntypes + 1)
    likelihood <- rep(NA, tmaxind)
    curhidprocgr_rest <- curhidprocgr[-ind, ]
    classtest1 <- curclassmat$classprobtest1
    classtest2 <- curclassmat$classprobtest2
    
    if (nindsgr == 2) {
        curhidprocgr_rest <- matrix(curhidprocgr_rest, nrow = 1, ncol = length(curhidprocgr_rest))
    }
    
    # for t=1
    
    CP1[1, ] <- curnu
    prob_star <- rep(1, ntypes + 1)
    LikeContri <- .Call("condlike", as.integer(ntypes + 1), as.integer(rep(obsertest1ind[1], ntypes + 
        1)), as.integer(rep(obsertest2ind[1], ntypes + 1)), as.integer(0:ntypes), as.integer(ntypes + 
        1), as.double(classtest1), as.double(classtest2))
    prob_star <- prob_star * LikeContri
    
    Rest <- .Call("transprobnotNA", as.double(c(1, curalpha)), as.double(c(1, curbeta)), as.double(c(1, 
        curmu)), as.double(curdelta), as.double(curgamma), as.integer(nindsgr - 1), as.integer(as.vector(curhidprocgr_rest[, 
        1:2])), as.integer(ntypes + 1), as.integer(rep(0, ntypes + 1)), as.integer(matrix(0, ntypes + 
        1, ntypes + 1)))
    CP2rest <- Rest[1:(ntypes + 1)]
    
    prob_star <- prob_star * CP1[1, ] * CP2rest
    
    CP2[1, ] <- prob_star/sum(prob_star)
    
    
    if (tmaxind == tmaxgr) {
        
        for (t in 2:(tmaxind - 1)) {
            probTrans <- .Call("transprob", as.double(c(1, curalpha)), as.double(c(1, curbeta)), 
                as.double(c(1, curmu)), as.double(curdelta), as.double(curgamma), as.integer(nindsgr - 
                  1), as.integer(curhidprocgr_rest[, t - 1]), as.integer(ntypes + 1), as.integer(rep(0, 
                  ntypes + 1)))
            probTrans <- matrix(probTrans, ntypes + 1, ntypes + 1, byrow = T)
            
            CP1[t, ] <- CP2[t - 1, ] %*% probTrans
            indsnotNA <- which(!is.na(curhidprocgr_rest[, t + 1]))
            Rest <- .Call("transprobnotNA", as.double(c(1, curalpha)), as.double(c(1, curbeta)), 
                as.double(c(1, curmu)), as.double(curdelta), as.double(curgamma), as.integer(length(indsnotNA)), 
                as.integer(as.vector(curhidprocgr_rest[indsnotNA, t:(t + 1)])), as.integer(ntypes + 
                  1), as.integer(rep(0, ntypes + 1)), as.integer(matrix(0, ntypes + 1, ntypes + 
                  1)))
            CP2rest <- Rest[1:(ntypes + 1)]
            
            if (is.na(obsertest1ind[t]) == TRUE) {
                prob_star <- CP1[t, ] * CP2rest
                CP2[t, ] <- prob_star/sum(prob_star)
                
            } else {
                prob_star <- rep(1, ntypes + 1)
                LikeContri <- .Call("condlike", as.integer(ntypes + 1), as.integer(rep(obsertest1ind[t], 
                  ntypes + 1)), as.integer(rep(obsertest2ind[t], ntypes + 1)), as.integer(0:ntypes), 
                  as.integer(ntypes + 1), as.double(classtest1), as.double(classtest2))
                prob_star <- prob_star * LikeContri
                prob_star <- prob_star * CP1[t, ] * CP2rest
                CP2[t, ] <- prob_star/sum(prob_star)
            }
        }
        
        t <- tmaxind
        probTrans <- .Call("transprob", as.double(c(1, curalpha)), as.double(c(1, curbeta)), as.double(c(1, 
            curmu)), as.double(curdelta), as.double(curgamma), as.integer(nindsgr - 1), as.integer(curhidprocgr_rest[, 
            t - 1]), as.integer(ntypes + 1), as.integer(rep(0, ntypes + 1)))
        probTrans <- matrix(probTrans, ntypes + 1, ntypes + 1, byrow = T)
        
        CP1[t, ] <- CP2[t - 1, ] %*% probTrans
        
        if (is.na(obsertest1ind[t]) == TRUE) {
            CP2[t, ] <- CP1[t, ]
        } else {
            prob_star <- rep(1, ntypes + 1)
            LikeContri <- .Call("condlike", as.integer(ntypes + 1), as.integer(rep(obsertest1ind[t], 
                ntypes + 1)), as.integer(rep(obsertest2ind[t], ntypes + 1)), as.integer(0:ntypes), 
                as.integer(ntypes + 1), as.double(classtest1), as.double(classtest2))
            prob_star <- prob_star * LikeContri
            prob_star <- prob_star * CP1[t, ]
            CP2[t, ] <- prob_star/sum(prob_star)
        }
    } else {
        
        for (t in 2:tmaxind) {
            probTrans <- .Call("transprob", as.double(c(1, curalpha)), as.double(c(1, curbeta)), 
                as.double(c(1, curmu)), as.double(curdelta), as.double(curgamma), as.integer(nindsgr - 
                  1), as.integer(curhidprocgr_rest[, t - 1]), as.integer(ntypes + 1), as.integer(rep(0, 
                  ntypes + 1)))
            probTrans <- matrix(probTrans, ntypes + 1, ntypes + 1, byrow = T)
            
            
            CP1[t, ] <- CP2[t - 1, ] %*% probTrans
            indsnotNA <- which(!is.na(curhidprocgr_rest[, t + 1]))
            Rest <- .Call("transprobnotNA", as.double(c(1, curalpha)), as.double(c(1, curbeta)), 
                as.double(c(1, curmu)), as.double(curdelta), as.double(curgamma), as.integer(length(indsnotNA)), 
                as.integer(as.vector(curhidprocgr_rest[indsnotNA, t:(t + 1)])), as.integer(ntypes + 
                  1), as.integer(rep(0, ntypes + 1)), as.integer(matrix(0, ntypes + 1, ntypes + 
                  1)))
            CP2rest <- Rest[1:(ntypes + 1)]
            
            if (is.na(obsertest1ind[t]) == TRUE) {
                prob_star <- CP1[t, ] * CP2rest
                CP2[t, ] <- prob_star/sum(prob_star)
                
            } else {
                prob_star <- rep(1, ntypes + 1)
                LikeContri <- .Call("condlike", as.integer(ntypes + 1), as.integer(rep(obsertest1ind[t], 
                  ntypes + 1)), as.integer(rep(obsertest2ind[t], ntypes + 1)), as.integer(0:ntypes), 
                  as.integer(ntypes + 1), as.double(classtest1), as.double(classtest2))
                prob_star <- prob_star * LikeContri
                prob_star <- prob_star * CP1[t, ] * CP2rest
                CP2[t, ] <- prob_star/sum(prob_star)
            }
        }
        
        
    }
    
    
    # Backward Sampling
    pbs <- CP2[tmaxind, ]
    thesis <- which(rmultinom(n = 1, size = 1, prob = pbs) == 1)
    curhidprocgr[ind, tmaxind] <- thesis - 1
    likelihood[tmaxind] <- pbs[thesis]
    
    for (t in (tmaxind - 1):1) {
        probTrans <- .Call("transprob", as.double(c(1, curalpha)), as.double(c(1, curbeta)), as.double(c(1, 
            curmu)), as.double(curdelta), as.double(curgamma), as.integer(nindsgr - 1), as.integer(curhidprocgr_rest[, 
            t]), as.integer(ntypes + 1), as.integer(rep(0, ntypes + 1)))
        probTrans <- matrix(probTrans, ntypes + 1, ntypes + 1, byrow = T)
        probTrans <- probTrans[, curhidprocgr[ind, t + 1] + 1]
        pbs <- (probTrans * CP2[t, ])/(CP1[t + 1, curhidprocgr[ind, t + 1] + 1])
        thesis <- which(rmultinom(n = 1, size = 1, prob = pbs) == 1)
        curhidprocgr[ind, t] <- thesis - 1
        likelihood[t] <- pbs[thesis]
    }
    
    return(curhidprocgr[ind, ])
}



