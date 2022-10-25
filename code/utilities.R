

##############################################################################
###
### Some auxiliary functions
###
##############################################################################


hhg <- function(x, y) {
  # to directly provide vectors to hhg.test
  
  hhg.test(as.matrix(dist((x), diag = TRUE, upper = TRUE)),
           as.matrix(dist((y), diag = TRUE, upper = TRUE)),
           nr.perm = 0)
}


show.size <- function(out, out_ind, what = "Circle", level = 0.95){
  cv <- apply(out_ind[[what]], 2, quantile, level) 
  return(list(cv = cv,
              size = rowMeans(t(out[[what]]) > cv) 
  ))
}


#ss <- show.size(out, out_ind, what = exm[1], level = 0.95)$size


print.table.size <- function(out, out_ind, level = 0.95){
  rbind("Circle"     = show.size(out, out_ind, what = exm[1], level)$size,
        "Diamond"    = show.size(out, out_ind, what = exm[2], level)$size,
        "Parabola"   = show.size(out, out_ind, what = exm[3], level)$size,
        "2Parabolas" = show.size(out, out_ind, what = exm[4], level)$size,
        "W"          = show.size(out, out_ind, what = exm[5], level)$size,
        "4indclouds" = show.size(out, out_ind, what = exm[6], level)$size
  )
}


print.table.cv <- function(out, out_ind, level = 0.95){
  rbind("Circle"     = show.size(out, out_ind, what = exm[1], level)$cv,
        "Diamond"    = show.size(out, out_ind, what = exm[2], level)$cv,
        "Parabola"   = show.size(out, out_ind, what = exm[3], level)$cv,
        "2Parabolas" = show.size(out, out_ind, what = exm[4], level)$cv,
        "W"          = show.size(out, out_ind, what = exm[5], level)$cv,
        "4indclouds" = show.size(out, out_ind, what = exm[6], level)$cv
  )
}


print.table.size.noise <- function(out, out_ind, level = 0.95){
  rbind("Circle"     = show.size(out, out_ind, what = exm[1], level)$size,
        "Diamond"    = show.size(out, out_ind, what = exm[2], level)$size,
        "Parabola"   = show.size(out, out_ind, what = exm[3], level)$size #,
        # "2Parabolas" = show.size(out, out_ind, what = exm[4], level)$size,
        # "W"          = show.size(out, out_ind, what = exm[5], level)$size,
        # "4indclouds" = show.size(out, out_ind, what = exm[6], level)$size
  )
}


print.table.simon <- function(out, out_ind, level = 0.95){
  # lin+noise          typ==1
  # parabolic+noise    typ==2
  # cubic+noise        typ==3
  # sin+noise          typ==4
  # x^(1/4) + noise    typ==5
  # circle             typ==6
  # two curves         typ==7
  # X function         typ==8
  # Diamond            typ==9
  rbind("1" = show.size(out, out_ind, what = type[1], level)$size,
        "2" = show.size(out, out_ind, what = type[2], level)$size,
        "3" = show.size(out, out_ind, what = type[3], level)$size,
        "4" = show.size(out, out_ind, what = type[4], level)$size,
        "5" = show.size(out, out_ind, what = type[5], level)$size,
        "6" = show.size(out, out_ind, what = type[6], level)$size,
        "7" = show.size(out, out_ind, what = type[7], level)$size,
        "8" = show.size(out, out_ind, what = type[8], level)$size,
        "9" = show.size(out, out_ind, what = type[9], level)$size 
  )
}







