library("BTYDplus")
library("BTYDplus")
library("BTYD")
library("hypergeo")

#elog2cbs converts the customer Id log and order Transactional dates to customer-by-sufficient-statistic

nop_cbs1  <-elog2cbs(NOP_DATA3_xlsx, unit= "weeks" , T.cal = "2017-01-09" , T.tot = "2017-01-16"  )

#estimate Regularity 

(k.wheat <- estimateRegularity(NOP_DATA2, method = "wheat",
                               plot = TRUE, title = "Wheat & Morrison"))

(k.mle <- estimateRegularity(NOP_DATA2, method = "mle",
                             plot = TRUE, title = "Maximum Likelihood"))

#NBD 

round(params.nbd <- nbd.EstimateParameters(nop_cbs1), 3)
nbd.cbs.LL(params.nbd, nop_cbs)

NBD_star    <- nbd.ConditionalExpectedTransactions(params.nbd,
                                    T.star = 1 , x = nop_cbs1$x , T.cal = nop_cbs1$T.cal )

rbind( Actuals = c(holdout = sum(nop_cbs1$x.star ) ) ,    
        NBD =    c(holdout =  round(sum(NBD_star))) ) 

#Pareto

params.pnbd <- BTYD::pnbd.EstimateParameters(nop_cbs , par.start =  c(1,1,1,1) , max.param.value = 10000  )

names(params.pnbd) <- c("r", "alpha", "s", "beta")
 round(params.pnbd, 3)



#estimate parameters for various models(BGNBD , BGCNBD , MBGNBD , MBGCNBD)

params.bgnbd <- BTYD::bgnbd.EstimateParameters(nop_cbs) # BG/NBD
params.bgcnbd <- bgcnbd.EstimateParameters(nop_cbs) # BG/CNBD-k
params.mbgnbd <- mbgnbd.EstimateParameters(nop_cbs) # MBG/NBD
params.mbgcnbd <- mbgcnbd.EstimateParameters(nop_cbs) # MBG/CNBD-k


row <- function(params, LL) {
  names(params) <- c("k", "r", "alpha", "a", "b")
  c(round(params, 3), LL = round(LL))
}
rbind(`BG/NBD` = row(c(1, params.bgnbd),
                     BTYD::bgnbd.cbs.LL(params.bgnbd, nop_cbs)),
      `BG/CNBD-k` = row(params.bgcnbd,
                        bgcnbd.cbs.LL(params.bgcnbd, nop_cbs)),
      `MBG/NBD` = row(params.mbgnbd,
                      mbgcnbd.cbs.LL(params.mbgnbd, nop_cbs)),
      `MBG/CNBD-k` = row(params.mbgcnbd,
                         mbgcnbd.cbs.LL(params.mbgcnbd, nop_cbs)))


# predict whole customer cohort number of transactions for the test period of week using
# nbd , mbgcnbd ,bgcnbd ) 

nop_cbs$xstar.nbd <- nbd.ConditionalExpectedTransactions(
  params = params.nbd, T.star = 1,
  x = nop_cbs$x,
  T.cal = nop_cbs$T.cal)


nop_cbs$xstar.mbgcnbd <- mbgcnbd.ConditionalExpectedTransactions(
  params = params.mbgcnbd, T.star = 1,
  x = nop_cbs$x, t.x = nop_cbs$t.x,
  T.cal = nop_cbs$T.cal)

nop_cbs$xstar.bgcnbd <- bgcnbd.ConditionalExpectedTransactions(
  params = params.bgcnbd, T.star = 1,
  x = nop_cbs$x, t.x = nop_cbs$t.x,
  T.cal = nop_cbs$T.cal)


# compare predictions with actuals at aggregated level

rbind(`Actuals` = c(`Holdout` = sum(nop_cbs$x.star)),
      `MBG/CNBD-k` = c(`Holdout` = round(sum(nop_cbs$xstar.mbgcnbd))))



rbind(`Actuals` = c(`Holdout` = sum(nop_cbs$x.star)),
      `BG/CNBD-k` = c(`Holdout` = round(sum(nop_cbs$xstar.bgcnbd))))



rbind(`Actuals` = c(`Holdout` = sum(nop_cbs$x.star)),
      `NBD-k` = c(`Holdout` = round(sum(nop_cbs$xstar.nbd))))



nil <- mbgcnbd.PlotTrackingInc(params.mbgcnbd,
                               T.cal = nop_cbs$T.cal,
                               T.tot = max(nop_cbs$T.cal + nop_cbs$T.star),
                               actual.inc.tracking = elog2inc(NOP_DATA3_xlsx))



nil <- nbd.PlotTrackingInc(params.mbgcnbd,
                               T.cal = nop_cbs$T.cal,
                               T.tot = max(nop_cbs$T.cal + nop_cbs$T.star),
                               actual.inc.tracking = elog2inc(NOP_DATA3_xlsx))



# mean absolute error (MAE)
mae <- function(act, est) {
  stopifnot(length(act)==length(est))
  sum(abs(act-est)) / sum(act)
}
mae.nbd <- mae(nop_cbs$x.star, nop_cbs$xstar.nbd)
mae.bgcnbd <- mae(nop_cbs$x.star, nop_cbs$xstar.bgcnbd)
mae.mbgcnbd <- mae(nop_cbs$x.star, nop_cbs$xstar.mbgcnbd)


rbind ( NBD = c(`MAE` = round(mae.nbd, 3)),
      `BG/CNBD-k` = c(`MAE` = round(mae.bgcnbd, 3)),
      `MBG/CNBD-k` = c(`MAE` = round(mae.mbgcnbd, 3))
      )


 
 ###Markov Chains Monte Carlo Simulation Models 
 
 pnbd.draws <- pnbd.mcmc.DrawParameters(nop_cbs)
 #> set param_init: 0.7463, 5.4544, 0.3817, 6.9221
 #> running in parallel on 2 cores
 # generate draws for holdout period
 pnbd.xstar.draws <- mcmc.DrawFutureTransactions(nop_cbs, pnbd.draws)
 # conditional expectations
 nop_cbs$xstar.pnbd.hb <- apply(pnbd.xstar.draws, 2, mean)
 # P(active)
 nop_cbs$pactive.pnbd.hb <- mcmc.PActive(pnbd.xstar.draws)
 # P(alive)
 nop_cbs$palive.pnbd.hb <- mcmc.PAlive(pnbd.draws)
 # show estimates for first few customers
 head( nop_cbs[, c("x", "t.x", "x.star",
                     "xstar.pnbd.hb", "pactive.pnbd.hb",
                     "palive.pnbd.hb")])
 



cohort.draws <- pnbd.draws$level_2
head(as.matrix(cohort.draws), 5)

nop_cbs$xstar.pnbd.hb <- pnbd.hb.ConditionalExpectedTransactions(
  params = pnbd.draws, T.star = 16,
  x = nop_cbs$x, t.x = nop_cbs$t.x,
  T.cal = nop_cbs$T.cal)

      
rbind(`Actuals` = c(`Holdout` = sum(nop_cbs$x.star)),
      `pnhb.hb` = c(`Holdout` = round(sum(nop_cbs$xstar.pnbd.hb))))


mae.pnbd.hb <- mae(nop_cbs$x.star, nop_cbs$xstar.pnbd.hb)


# pggg model 

pggg.draws <- pggg.mcmc.DrawParameters(nop_cbs) # ~2mins on 2015 MacBook Pro
# generate draws for holdout period
pggg.xstar.draws <- mcmc.DrawFutureTransactions(nop_cbs, pggg.draws)
# conditional expectations
nop_cbs$xstar.pggg <- apply(pggg.xstar.draws, 2, mean)
# P(active)
nop_cbs$pactive.pggg <- mcmc.PActive(pggg.xstar.draws)
# P(alive)
nop_cbs$palive.pggg <- mcmc.PAlive(pggg.draws)
# show estimates for first few customers
head(nop_cbs[, c("x", "t.x", "x.star",
                    "xstar.pggg", "pactive.pggg", "palive.pggg")])

round(apply(as.matrix(pggg.draws$level_2), 2, median), 3)


median.est <- sapply(pggg.draws$level_1, function(draw) {
  apply(as.matrix(draw), 2, median)
})
round(apply(median.est, 1, mean), 3)


rbind(`Actuals` = c(`Holdout` = sum(nop_cbs$x.star)),
      `Pareto/GGG` = round(sum(nop_cbs$xstar.pggg)),
      `MBG/CNBD-k` = round(sum(nop_cbs$xstar.mbgcnbd)),
      `Pareto/NBD (HB)` = round(sum(nop_cbs$xstar.pnbd.hb)))
