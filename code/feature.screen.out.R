
out.table.screen <- function(res){
  res.es <- cbind(apply(res, 2, quantile, probs  = 0.05),
                  apply(res, 2, quantile, probs  = 0.25),
                  apply(res, 2, quantile, probs  = 0.5),
                  apply(res, 2, quantile, probs  = 0.75),
                  apply(res, 2, quantile, probs  = 0.95))
  colnames(res.es) <- c("5%","25%", "50%", "75%", "95%")
  return(res.es)
}


library(xtable)

# Example Case 1 --------------------------------------------------------------

load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_1.Rdata")
es1a_1000 <- out.1a
load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_1_ultra.Rdata")
es1a_2500 <- out.1a
load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_1_ultra5.Rdata")
es1a_5000 <- out.1a

a <- rbind(out.table.screen(es1a_1000),
          out.table.screen(es1a_2500),
          out.table.screen(es1a_5000))

Sa <- rbind(apply(es1a_1000,2,function(x)sum(x <= 20 )/1000),
            apply(es1a_2500,2,function(x)sum(x <= 20 )/1000),
            apply(es1a_5000,2,function(x)sum(x <= 20 )/1000))


Sa.t2 <- rbind(apply(es1a_1000,2,function(x) sum(x <= 40 )/1000),
               apply(es1a_2500,2,function(x) sum(x <= 40 )/1000),
               apply(es1a_5000,2,function(x) sum(x <= 40 )/1000))


# Example Case 2 --------------------------------------------------------------

load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_2.Rdata")
es1c_1000 <- out.1a
load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_2_ultra.Rdata")
es1c_2500 <- out.1a
load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_2_ultra5.Rdata")
es1c_5000 <- out.1a


c <- rbind(out.table.screen(es1c_1000),
           out.table.screen(es1c_2500),
           out.table.screen(es1c_5000))

Sc <- rbind(apply(es1c_1000,2,function(x)sum(x <= 20 )/1000),
            apply(es1c_2500,2,function(x)sum(x <= 20 )/1000),
            apply(es1c_5000,2,function(x)sum(x <= 20 )/1000))


Sc.t2 <- rbind(apply(es1c_1000,2,function(x)sum(x <= 40 )/1000),
               apply(es1c_2500,2,function(x)sum(x <= 40 )/1000),
               apply(es1c_5000,2,function(x)sum(x <= 40 )/1000))


# Example Case 3 --------------------------------------------------------------

load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_3.Rdata")
es1e_1000 <- out.1a
load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_3_ultra.Rdata")
es1e_2500 <- out.1a
load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_3_ultra5.Rdata")
es1e_5000 <- out.1a


e <- rbind(out.table.screen(es1e_1000),
           out.table.screen(es1e_2500),
           out.table.screen(es1e_5000))

Se <- rbind(apply(es1e_1000,2,function(x)sum(x <= 20 )/1000),
            apply(es1e_2500,2,function(x)sum(x <= 20 )/1000),
            apply(es1e_5000,2,function(x)sum(x <= 20 )/1000))

Se.t2 <- rbind(apply(es1e_1000,2,function(x)sum(x <= 40 )/1000),
            apply(es1e_2500,2,function(x)sum(x <= 40 )/1000),
            apply(es1e_5000,2,function(x)sum(x <= 40 )/1000))

# Example Case 4 --------------------------------------------------------------

load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_4.Rdata")
es1f_1000 <- out.1a
load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_4_ultra.Rdata")
es1f_2500 <- out.1a
load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_4_ultra5.Rdata")
es1f_5000 <- out.1a


f <- rbind(out.table.screen(es1f_1000),
           out.table.screen(es1f_2500),
           out.table.screen(es1f_5000))

Sf <- rbind(apply(es1f_1000,2,function(x)sum(x <= 20 )/1000),
            apply(es1f_2500,2,function(x)sum(x <= 20 )/1000),
            apply(es1f_5000,2,function(x)sum(x <= 20 )/1000))

Sf.t2 <- rbind(apply(es1f_1000,2,function(x)sum(x <= 40 )/1000),
            apply(es1f_2500,2,function(x)sum(x <= 40 )/1000),
            apply(es1f_5000,2,function(x)sum(x <= 40 )/1000))

print(xtable(cbind(e,f), 
             caption = " "), 
      floating = TRUE, latex.environments = "center")

print(xtable(cbind(Se,Sf), 
             caption = " "), 
      floating = TRUE, latex.environments = "center")



# Example Case 5 --------------------------------------------------------------

load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_5.Rdata")
es2g_1000 <- out.1a
load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_5_ultra.Rdata")
es2g_2500 <- out.1a

load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_5_ultra5.Rdata")
es2g_5000 <- out.1a


g2 <- rbind(out.table.screen(es2g_1000),
            out.table.screen(es2g_2500),
            out.table.screen(es2g_5000))

# g2<- g2[-c(4,8,12),]  ### Bn_pq

g2<- g2[-c(3,7,11),]  ### Bn_pq

Sg2 <- rbind(apply(es2g_1000,2,function(x)sum(x <= 20 )/1000),
             apply(es2g_2500,2,function(x)sum(x <= 20 )/1000),
             apply(es2g_5000,2,function(x)sum(x <= 20 )/1000))

# Sg2 <- Sg2[,-4]
Sg2 <- Sg2[,-3]

## d2 threshold
Sg2.t2 <- rbind(apply(es2g_1000,2,function(x)sum(x <= 40 )/1000),
                apply(es2g_2500,2,function(x)sum(x <= 40 )/1000),
                apply(es2g_5000,2,function(x)sum(x <= 40 )/1000))

# Sg2.t2 <- Sg2.t2[,-4]

Sg2.t2 <- Sg2.t2[,-3]


# Example Case 6 --------------------------------------------------------------

load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_6.Rdata")
es2e_1000 <- out.1a
load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_6_ultra.Rdata")
es2e_2500 <- out.1a
load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_6_ultra5.Rdata")
es2e_5000 <- out.1a


e2 <- rbind(out.table.screen(es2e_1000),
            out.table.screen(es2e_2500),
            out.table.screen(es2e_5000))

e2<- e2[-c(3,7,11),]  ### Bn_pq

Se2 <- rbind(apply(es2e_1000,2,function(x)sum(x <= 20 )/1000),
             apply(es2e_2500,2,function(x)sum(x <= 20 )/1000),
             apply(es2e_5000,2,function(x)sum(x <= 20 )/1000))

Se2 <- Se2[,-3]

## d2 threshold
Se2.t2 <- rbind(apply(es2e_1000,2,function(x)sum(x <= 40 )/1000),
                apply(es2e_2500,2,function(x)sum(x <= 40 )/1000),
                apply(es2e_5000,2,function(x)sum(x <= 40 )/1000))

Se2.t2 <- Se2.t2[,-3]


# Example Case 7 --------------------------------------------------------------

load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_7.Rdata")
es2f_1000 <- out.1a
load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_7_ultra.Rdata")
es2f_2500 <- out.1a

load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_7_ultra5.Rdata")
es2f_5000 <- out.1a


f2 <- rbind(out.table.screen(es2f_1000),
            out.table.screen(es2f_2500),
            out.table.screen(es2f_5000))

f2<- f2[-c(3,7,11),]  ### Bn_pq

Sf2 <- rbind(apply(es2f_1000,2,function(x)sum(x <= 20 )/1000),
             apply(es2f_2500,2,function(x)sum(x <= 20 )/1000),
             apply(es2f_5000,2,function(x)sum(x <= 20 )/1000))

Sf2 <- Sf2[,-3]

## d2 threshold
Sf2.t2 <- rbind(apply(es2f_1000,2,function(x)sum(x <= 40 )/1000),
                apply(es2f_2500,2,function(x)sum(x <= 40 )/1000),
                apply(es2f_5000,2,function(x)sum(x <= 40 )/1000))

Sf2.t2 <- Sf2.t2[,-3]




# Example Case 8 --------------------------------------------------------------
## non è riportato nel testo 

load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_8.Rdata")
es8_1000 <- out.1a
load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_8_ultra.Rdata")
es8_2500 <- out.1a

load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_8_ultra5.Rdata")
es8_5000 <- out.1a


case8 <- rbind(out.table.screen(es8_1000),
            out.table.screen(es8_2500),
            out.table.screen(es8_5000))

case8<- case8[-c(3,7,11),]  ### Bn_pq

Ses8 <- rbind(apply(es8_1000,2,function(x)sum(x <= 20 )/1000),
             apply(es8_2500,2,function(x)sum(x <= 20 )/1000),
             apply(es8_5000,2,function(x)sum(x <= 20 )/1000))

Ses8 <- Ses8[,-3]

## d2 threshold
Ses8.t2 <- rbind(apply(es8_1000,2,function(x)sum(x <= 40 )/1000),
                apply(es8_2500,2,function(x)sum(x <= 40 )/1000),
                apply(es8_5000,2,function(x)sum(x <= 40 )/1000))

Ses8.t2 <- Ses8.t2[,-3]


# Example Case 9 --------------------------------------------------------------
## nel testo compare come case 8

load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_9.Rdata")
es9_1000 <- out.1a

load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_9_ultra.Rdata")
es9_2500 <- out.1a

load("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/case_9_ultra5.Rdata")
es9_5000 <- out.1a


case9 <- rbind(out.table.screen(es9_1000),
               out.table.screen(es9_2500),
               out.table.screen(es9_5000)
               )

case9<- case9[-c(3,7,11),]  ### Bn_pq

Ses9 <- rbind(apply(es9_1000,2,function(x)sum(x <= 20 )/1000),
              apply(es9_2500,2,function(x)sum(x <= 20 )/1000),
              apply(es9_5000,2,function(x)sum(x <= 20 )/1000))

Ses9 <- Ses9[,-3]

## d2 threshold
Ses9.t2 <- rbind(apply(es9_1000,2,function(x)sum(x <= 40 )/1000),
                 apply(es9_2500,2,function(x)sum(x <= 40 )/1000),
                 apply(es9_5000,2,function(x)sum(x <= 40 )/1000)
                 )

Ses9.t2 <- Ses9.t2[,-3]



# ### tabelle finali R ----------------------------------------------------

## cases 1-4
print(xtable(t(cbind(a, c, e, f)), 
             caption = " "), 
      floating = TRUE, latex.environments = "center")

# print(xtable(cbind(e, f), 
#              caption = " "), 
#       floating = TRUE, latex.environments = "center")
# 
# print(xtable(cbind(a2, d2), 
#              caption = " "), 
#       floating = TRUE, latex.environments = "center")


### tabelle finali S

print(xtable(cbind(Sa, Sc, Se, Sf), 
             caption = " "), 
      floating = TRUE, latex.environments = "center")


print(xtable(cbind(Sa.t2, Sc.t2, Se.t2, Sf.t2), 
             caption = " "), 
      floating = TRUE, latex.environments = "center")

# print(xtable(cbind(Sf.t2, Sa2.t2, Sd2.t2), 
#              caption = " "), 
#       floating = TRUE, latex.environments = "center")



## cases 5-7

print(xtable(t(case9),
             caption = " "),
      floating = TRUE, latex.environments = "center")

print(xtable(t(cbind(g2,e2,f2, case9)), 
             caption = " "), 
      floating = TRUE, latex.environments = "center")

# print(xtable(cbind(e, f), 
#              caption = " "), 
#       floating = TRUE, latex.environments = "center")
# 
# print(xtable(cbind(a2, d2), 
#              caption = " "), 
#       floating = TRUE, latex.environments = "center")


### tabelle finali S

print(xtable(cbind(Sg2, Se2, Sf2, Ses9), 
             caption = " "), 
      floating = TRUE, latex.environments = "center")


print(xtable(cbind(Sg2.t2, Se2.t2, Sf2.t2, Ses9.t2), 
             caption = " "), 
      floating = TRUE, latex.environments = "center")



