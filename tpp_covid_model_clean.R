# Draw a large number of parameter combinations
# In an array, put the parameters chosen and TPP, TSP, CI, supply and demand for tests by risk group

t_end = 182
testrateifsymp_l=0.05
ds_inc = 0
testrateifasymp_l=3e-5
da_inc = 0
pili=0.01
pili_inc = 0
sensitivity=0.95
highrisk_prop=0.5
rel_susc_h=1
test_inc_perday=30
reltest_symp=1
reltest_asymp=1
nruns = 1000
results_suffix = 'default'
  

npop = 1000000
gamma=1/6
sigma=1/5.2
initial_I = npop/1000
initial_tests = npop/10000
initial_rt = initial_tmax = npop/1000

t_lockdown = c(21,35)
R0_pre=c(3,5)
R0_post=c(0.8,1)

lowerbounds = c("R0_pre"=head(R0_pre,1),
                "R0_post"=head(R0_post,1),
                "t_lockdown"=head(t_lockdown,1),
                "dal"=head(testrateifasymp_l,1),
                "da_inc" = head(da_inc,1),
                "dsl"=head(testrateifsymp_l,1),
                "ds_inc" = head(ds_inc,1),
                "pili"=head(pili,1),
                "pili_inc"=head(pili_inc,1),
                "sensitivity"=head(sensitivity,1),
                "test_inc_perday"=head(test_inc_perday,1),
                "highrisk_prop"=head(highrisk_prop,1),
                "rel_susc_h"=head(rel_susc_h,1),
                "reltest_symp"=head(reltest_symp,1),
                "reltest_asymp"=head(reltest_asymp,1))

upperbounds = c("R0_pre"=tail(R0_pre,1),
                "R0_post"=tail(R0_post,1),
                "t_lockdown"=tail(t_lockdown,1),
                "dal"=tail(testrateifasymp_l,1),
                "da_inc" = tail(da_inc,1),
                "dsl"=tail(testrateifsymp_l,1),
                "ds_inc" = tail(ds_inc,1),
                "pili"=tail(pili,1),
                "pili_inc"=tail(pili_inc,1),
                "sensitivity"=tail(sensitivity,1),
                "test_inc_perday"=tail(test_inc_perday,1),
                "highrisk_prop"=tail(highrisk_prop,1),
                "rel_susc_h"=tail(rel_susc_h,1),
                "reltest_symp"=tail(reltest_symp,1),
                "reltest_asymp"=tail(reltest_asymp,1))

require(adaptivetau)

test_supply_function = function(demand,supply,max_t_perday) {
  
  supply_maxed = min(supply,max_t_perday)
  
  return(min(demand,supply_maxed))
  
}

testing_t = function(t,initial_tmax,tmax_increase_perday) {
  
  initial_tmax + t*tmax_increase_perday
  
}

transitions_SEITR = list(
  c(SL = -1,EL = +1), # infection, low risk
  c(EL = -1, IL = +1), # exposed to infectious, low risk
  c(EL = -1, R = +1), # quarantine of exposed, low risk individuals
  c(IL = -1,R = +1), # recovery of infectious, low risk
  c(SH = -1,EH = +1), # infection, high risk
  c(EH = -1, IH = +1), # exposed to infectious, high risk
  c(EH = -1, R = +1), # quarantine of exposed, high risk individuals
  c(IH = -1,R = +1), # recovery of infectious, high risk
  c(T = +1), # replenishment of test
  c(T = -1, TPHS=+1), # positive test (high risk, symptomatic)
  c(T = -1, TPLS=+1), # positive test (low risk, symptomatic)
  c(T = -1, TNHS=+1), # negative test (high risk, symptomatic)
  c(T = -1, TNLS=+1), # negative test (low risk, symptomatic)
  c(T = -1, TPHA=+1), # positive test (high risk, asymptomatic)
  c(T = -1, TPLA=+1), # positive test (low risk, asymptomatic)
  c(T = -1, TNHA=+1), # negative test (high risk, asymptomatic)
  c(T = -1, TNLA=+1) # negative test (low risk, asymptomatic)
) 

rates <- function(x, params, t) {
  
  # Time-varying params
  dal = params$dal + min(t,35)*params$dal_inc
  dah = params$dah + min(t,35)*params$dah_inc
  dsl = params$dsl + min(t,35)*params$dsl_inc
  dsh = params$dsh + min(t,35)*params$dsh_inc
  pili = params$pili + min(t,35)*params$pili_inc
  
  test_inc_perday = params$test_inc_perday
  
  avail_test_perday = testing_t(t,params$initial_tmax,test_inc_perday)
  
  # Beta params
  beta_l = (params$beta_pre*(t<params$t_lockdown)+params$beta_post*(t>=params$t_lockdown))
  beta_h = (params$beta_pre_h*(t<params$t_lockdown)+params$beta_post_h*(t>=params$t_lockdown))
  
  # Calculate total demand and supply
  total_demand = dal*x["EL"]+
    (pili*dsl+(1-pili)*dal)*x["SL"] + 
    dsl*x["IL"] + 
    dah*x["EH"]+
    (pili*dsh+(1-pili)*dah)*x["SH"]+
    dsh*x["IH"]
  total_supply = test_supply_function(total_demand,x["T"],avail_test_perday)
  
  # Apportion total supply by relative demand for tests from each group
  supply_els = pili*dsl*x["EL"]/total_demand*total_supply
  supply_sls = pili*dsl*x["SL"]/total_demand*total_supply
  supply_il = dsl*x["IL"]/total_demand*total_supply
  supply_ehs = pili*dsh*x["EH"]/total_demand*total_supply
  supply_shs = pili*dsh*x["SH"]/total_demand*total_supply
  supply_ih = dsh*x["IH"]/total_demand*total_supply
  supply_ela = (1-pili)*dal*x["EL"]/total_demand*total_supply
  supply_sla = (1-pili)*dal*x["SL"]/total_demand*total_supply
  supply_eha = (1-pili)*dah*x["EH"]/total_demand*total_supply
  supply_sha = (1-pili)*dah*x["SH"]/total_demand*total_supply
  
  return(c(beta_l*(x["IL"]+x["IH"])*x["SL"]/(x["SL"]+x["SH"]+x["EL"]+x["EH"]+x["IL"]+x["IH"]+x["R"]),     # rate of infection
           x["EL"]*params$sigma,  # rate of exposed becoming infectious
           (supply_els+supply_ela)*params$sensitivity, # quarantine of exposed individuals who get tested
           x["IL"]*params$gamma + supply_il*params$sensitivity, # rate of recovery plus testing
           beta_h*(x["IL"]+x["IH"])*x["SH"]/(x["SL"]+x["SH"]+x["EL"]+x["EH"]+x["IL"]+x["IH"]+x["R"]),     # rate of infection
           x["EH"]*params$sigma,  # rate of exposed becoming infectious
           (supply_eha+supply_ehs)*params$sensitivity, # quarantine of exposed individuals who get tested
           x["IH"]*params$gamma+supply_ih*params$sensitivity, # rate of recovery plus testing
           testing_t(t,params$initial_rt,test_inc_perday), # rate of replenishment of tests
           (supply_ehs+supply_ih)*params$sensitivity, # rate of tests being positive (high risk, symptomatic)
           (supply_els+supply_il)*params$sensitivity, # rate of tests being positive (low risk, symptomatic)
           supply_shs+(supply_ehs+supply_ih)*(1-params$sensitivity), # rate of tests being negative (high risk, symptomatic)
           supply_sls+(supply_els+supply_il)*(1-params$sensitivity), # rate of tests being negative (low risk, symptomatic)
           supply_eha*params$sensitivity, # rate of tests being positive (high risk, asymptomatic)
           supply_ela*params$sensitivity, # rate of tests being positive (low risk, asymptomatic)
           supply_sha+supply_eha*(1-params$sensitivity), # rate of tests being negative (high risk, asymptomatic)
           supply_sla+supply_ela*(1-params$sensitivity) # rate of tests being negative (low risk, asymptomatic)
  ) 
  )     
}

t_weeks = seq(7,t_end,7)/7
res_array = matrix(NA,nrow=nruns,ncol=length(lowerbounds)+t_end/7*10)
colnames(res_array) = c(names(lowerbounds),
                        paste0("confcaseh_",t_weeks),
                        paste0("confcasel_",t_weeks),
                        paste0("tpr_",t_weeks),
                        paste0("tprh_",t_weeks),
                        paste0("tprl_",t_weeks),
                        paste0("tsr_",t_weeks),
                        paste0("tsrh_",t_weeks),
                        paste0("tsrl_",t_weeks),
                        paste0("supply_",t_weeks),
                        paste0("demand_",t_weeks)
)

for (run in 1:nruns) {
  
  draw_params = runif(length(lowerbounds),min=lowerbounds,max=upperbounds)
  names(draw_params)=names(lowerbounds)
  
  beta_pre = draw_params[['R0_pre']]*gamma/
    (1+draw_params[['highrisk_prop']]*(draw_params[['rel_susc_h']]-1))
  beta_post = draw_params[['R0_post']]*gamma/
    (1+draw_params[['highrisk_prop']]*(draw_params[['rel_susc_h']]-1))
  testrateifsymp_h=draw_params[['dsl']]*draw_params[['reltest_symp']]
  testrateifasymp_h=draw_params[['dal']]*draw_params[['reltest_asymp']]
  
  res_array[run,colnames(res_array) %in% names(draw_params)]=draw_params
  
  params = list(beta_pre=beta_pre,
                beta_post=beta_post,
                t_lockdown=draw_params[['t_lockdown']],
                beta_pre_h=beta_pre*draw_params[['rel_susc_h']],
                beta_post_h=beta_post*draw_params[['rel_susc_h']],
                sigma=sigma,gamma=gamma,
                initial_rt=initial_rt,
                initial_tmax=initial_tmax,
                test_inc_perday=draw_params[['test_inc_perday']],
                dal=draw_params[['dal']],
                dal_inc=draw_params[['da_inc']],
                dah=testrateifasymp_h,
                dah_inc=draw_params[['da_inc']]*draw_params[['reltest_asymp']],
                dsl=draw_params[['dsl']],
                dsl_inc=draw_params[['ds_inc']],
                dsh=testrateifsymp_h,
                dsh_inc=draw_params[['ds_inc']]*draw_params[['reltest_symp']],
                sensitivity=sensitivity,
                pili=draw_params[['pili']],
                pili_inc=draw_params[['pili_inc']]
                )
  
  r = ssa.adaptivetau(c(SL=as.integer((npop-initial_I)*(1-draw_params[['highrisk_prop']])),EL=0,IL=as.integer(initial_I*(1-draw_params[['highrisk_prop']])),
                        SH=as.integer((npop-initial_I)*draw_params[['highrisk_prop']]),EH=0,IH=as.integer(initial_I*draw_params[['highrisk_prop']]),
                        R=0,
                        T=initial_tests,TPHS=0,TPLS=0,TNHS=0,TNLS=0,TPHA=0,TPLA=0,TNHA=0,TNLA=0),
                      transitions_SEITR,rates,params,t_end)
  
  avail_test_perday = testing_t(r[,"time"],params$initial_tmax,params$test_inc_perday)
  
  # Calculate total demand and supply
  total_demand = params$dal*r[,"EL"]+
    (params$pili*params$dsl+(1-params$pili)*params$dal)*r[,"SL"] + 
    params$dsl*r[,"IL"] + 
    params$dah*r[,"EH"]+
    (params$pili*params$dsh+(1-params$pili)*params$dah)*r[,"SH"]+
    params$dsh*r[,"IH"]
  total_supply = sapply(1:nrow(r),function(x) test_supply_function(total_demand[x],r[x,"T"],avail_test_perday[x]))
  
  r = cbind(r,total_demand,total_supply)
  
  # Smooth output and round to nearest day
  tp_hs_loess = loess(r[,"TPHS"]~r[,"time"],span=0.05)
  tp_hs_smooth = round(predict(tp_hs_loess,newdata=1:t_end))
  tp_ls_loess = loess(r[,"TPLS"]~r[,"time"],span=0.05)
  tp_ls_smooth = round(predict(tp_ls_loess,newdata=1:t_end))
  tn_hs_loess = loess(r[,"TNHS"]~r[,"time"],span=0.05)
  tn_hs_smooth = round(predict(tn_hs_loess,newdata=1:t_end))
  tn_ls_loess = loess(r[,"TNLS"]~r[,"time"],span=0.05)
  tn_ls_smooth = round(predict(tn_ls_loess,newdata=1:t_end))
  
  tp_ha_loess = loess(r[,"TPHA"]~r[,"time"],span=0.05)
  tp_ha_smooth = round(predict(tp_ha_loess,newdata=1:t_end))
  tp_la_loess = loess(r[,"TPLA"]~r[,"time"],span=0.05)
  tp_la_smooth = round(predict(tp_la_loess,newdata=1:t_end))
  tn_ha_loess = loess(r[,"TNHA"]~r[,"time"],span=0.05)
  tn_ha_smooth = round(predict(tn_ha_loess,newdata=1:t_end))
  tn_la_loess = loess(r[,"TNLA"]~r[,"time"],span=0.05)
  tn_la_smooth = round(predict(tn_la_loess,newdata=1:t_end))
  
  demand_loess = loess(r[,"total_demand"]~r[,"time"],span=0.05)
  demand_smooth = round(predict(demand_loess,newdata=1:t_end))
  supply_loess = loess(r[,"total_supply"]~r[,"time"],span=0.05)
  supply_smooth = round(predict(supply_loess,newdata=1:t_end))
  demand_weekly = sapply(seq(1,t_end-6,7),function(x) 
    mean(demand_smooth[x:(x+6)]))
  supply_weekly = sapply(seq(1,t_end-6,7),function(x) 
    mean(supply_smooth[x:(x+6)]))
  
  
  new_tp_hs = diff(c(0,tp_hs_smooth))
  new_tn_hs = diff(c(0,tn_hs_smooth))
  new_tp_ls = diff(c(0,tp_ls_smooth))
  new_tn_ls = diff(c(0,tn_ls_smooth))
  
  new_tp_ha = diff(c(0,tp_ha_smooth))
  new_tn_ha = diff(c(0,tn_ha_smooth))
  new_tp_la = diff(c(0,tp_la_smooth))
  new_tn_la = diff(c(0,tn_la_smooth))
  
  new_tp_h = new_tp_hs+new_tp_ha
  new_tn_h = new_tn_hs+new_tn_ha
  new_tp_l = new_tp_ls+new_tp_la
  new_tn_l = new_tn_ls+new_tn_la
  
  new_tp_hs_weekly = sapply(seq(1,t_end-6,7),function(x) 
    mean(new_tp_hs[x:(x+6)]))
  new_tp_ha_weekly = sapply(seq(1,t_end-6,7),function(x) 
    mean(new_tp_ha[x:(x+6)]))
  new_tp_ls_weekly = sapply(seq(1,t_end-6,7),function(x) 
    mean(new_tp_ls[x:(x+6)]))
  new_tp_la_weekly = sapply(seq(1,t_end-6,7),function(x) 
    mean(new_tp_la[x:(x+6)]))
  new_tn_hs_weekly = sapply(seq(1,t_end-6,7),function(x) 
    mean(new_tn_hs[x:(x+6)]))
  new_tn_ha_weekly = sapply(seq(1,t_end-6,7),function(x) 
    mean(new_tn_ha[x:(x+6)]))
  new_tn_ls_weekly = sapply(seq(1,t_end-6,7),function(x) 
    mean(new_tn_ls[x:(x+6)]))
  new_tn_la_weekly = sapply(seq(1,t_end-6,7),function(x) 
    mean(new_tn_la[x:(x+6)]))
  
  new_tests_h_weekly = new_tp_hs_weekly+new_tp_ha_weekly+new_tn_hs_weekly+new_tn_ha_weekly
  new_tests_l_weekly = new_tp_ls_weekly+new_tp_la_weekly+new_tn_ls_weekly+new_tn_la_weekly
  new_tests_weekly = new_tests_h_weekly+new_tests_l_weekly
  
  tpr = (new_tp_hs_weekly+new_tp_ha_weekly+new_tp_ls_weekly+new_tp_la_weekly)/new_tests_weekly
  tpr_h = (new_tp_hs_weekly+new_tp_ha_weekly)/
    new_tests_h_weekly
  tpr_l = (new_tp_ls_weekly+new_tp_la_weekly)/
    new_tests_l_weekly
  
  tsr = (new_tp_hs_weekly+new_tp_ls_weekly+new_tn_hs_weekly+new_tn_ls_weekly)/new_tests_weekly
  tsr_h = (new_tp_hs_weekly+new_tn_hs_weekly)/
    new_tests_h_weekly
  tsr_l = (new_tp_ls_weekly+new_tn_ls_weekly)/
    new_tests_l_weekly
  
  pred_confih_tpr = (params$sensitivity*params$dsh)*
    tpr_h*(params$pili+params$dah/params$dsh*(1-params$pili))/(params$sensitivity-tpr_h*(1-params$dah/params$dsh)*(1-params$pili))
  pred_confih_tsr = (params$sensitivity*params$dsh)*
    (tsr_h*(params$pili+params$dah/params$dsh*(1-params$pili))-params$pili)/((1-params$pili)*(1-tsr_h*(1-params$dah/params$dsh)))
  confih = (new_tp_hs_weekly+new_tp_ha_weekly)/(npop*draw_params[['highrisk_prop']])
  
  pred_confil_tpr = (params$sensitivity*params$dsl)*
    tpr_l*(params$pili+params$dal/params$dsl*(1-params$pili))/(params$sensitivity-tpr_l*(1-params$dal/params$dsl)*(1-params$pili))
  pred_confil_tsr = (params$sensitivity*params$dsl)*
    (tsr_l*(params$pili+params$dal/params$dsl*(1-params$pili))-params$pili)/((1-params$pili)*(1-tsr_l*(1-params$dal/params$dsl)))
  confil = (new_tp_ls_weekly+new_tp_la_weekly)/(npop*(1-draw_params[['highrisk_prop']]))
  
  confi = draw_params[['highrisk_prop']]*confih+(1-draw_params[['highrisk_prop']])*confil
  
  # Output tpr_h, tpr_l, tsr_h, tsr_l, confih, confil, supply and demand
  res_line = c(confih,confil,tpr,tpr_h,tpr_l,tsr,tsr_h,tsr_l,supply_weekly,demand_weekly)
  
  res_array[run,(length(lowerbounds)+1):ncol(res_array)]=
    res_line
  
}

res_array = as.data.frame(res_array)

require(scales)
require(ggplot2)
require(dplyr)
require(grid)
require(gridExtra)
theme_set(theme_bw())

dat_forcurve2 = res_array %>% summarise(sensitivity=mean(sensitivity),
                                        dsl=mean(dsl),
                                        dal=mean(dal),
                                        pili=mean(pili))
datx = data.frame(confcase=c(0,0.002))
dat_forcurve2 = cbind(datx,dat_forcurve2)

g1=ggplot(dat = 
            res_array %>% 
            mutate(confcase_predTPRH = sensitivity*dsl*reltest_symp*tprh_4*(pili+dal*reltest_asymp/(dsl*reltest_symp)*(1-pili))/
                     (sensitivity-tprh_4*(1-dal*reltest_asymp/(dsl*reltest_symp))*(1-pili)))
)+
  geom_point(aes(x=10000*confcaseh_4,y=tprh_4,color=factor(supply_4>=demand_4)),alpha=0.1,size=0.1)+
  xlab('Confirmed incidence per 10,000 at 4 weeks')+
  ylab('TPP at 4 weeks')+
  stat_function(data=dat_forcurve2,aes(x=10000*confcase),
                fun=function(x) x/10000/(dat_forcurve2$dsl[1]*dat_forcurve2$pili[1]+dat_forcurve2$dal[1]*(1-dat_forcurve2$pili[1])+
                                           (1-dat_forcurve2$dal[1]/dat_forcurve2$dsl[1])*(1-dat_forcurve2$pili[1])/
                                           dat_forcurve2$sensitivity[1]*x/10000),color='black')+
  ylim(c(0,1))+
  scale_color_manual(values=c('red','blue'),labels=c('No','Yes'),name='Sufficient supply of tests')+
  theme(
    panel.grid = element_line(color = "#fafafa")
    , strip.background = element_blank()
    , strip.text = element_text(hjust = 0, size = 8, face = "bold")
    , axis.text = element_text(size=8)
    , axis.title = element_text(size=8)
    , legend.text = element_text(size=8)
    , legend.title = element_text(size=8)
    , legend.position = "bottom"
  )+
  guides(
    color = guide_legend(override.aes = list(alpha=1,size=1))
  )

g2=ggplot(dat = 
            res_array %>% 
            mutate(confcase_predTPRH = sensitivity*dsl*reltest_symp*tprh_4*(pili+dal*reltest_asymp/(dsl*reltest_symp)*(1-pili))/
                     (sensitivity-tprh_4*(1-dal*reltest_asymp/(dsl*reltest_symp))*(1-pili)))
)+
  geom_point(aes(x=10000*confcaseh_4,y=tsrh_4,color=factor(supply_4>=demand_4)),alpha=0.1,size=0.1)+
  xlab('Confirmed incidence per 10,000 at 4 weeks')+
  ylab('TSP at 4 weeks')+
  ylim(c(0.9,1))+
  stat_function(data=dat_forcurve2,aes(x=10000*confcase),
                fun=function(x) (dat_forcurve2$dsl[1]*dat_forcurve2$pili[1]+(1-dat_forcurve2$pili[1])/
                                   dat_forcurve2$sensitivity[1]*x/10000)/
                  (dat_forcurve2$dsl[1]*dat_forcurve2$pili[1]+dat_forcurve2$dal[1]*(1-dat_forcurve2$pili[1])+
                     (1-dat_forcurve2$dal[1]/dat_forcurve2$dsl[1])*(1-dat_forcurve2$pili[1])/
                     dat_forcurve2$sensitivity[1]*x/10000),color='black')+
  scale_color_manual(values=c('red','blue'),labels=c('No','Yes'),name='Sufficient supply of tests')+
  theme(
    panel.grid = element_line(color = "#fafafa")
    , strip.background = element_blank()
    , strip.text = element_text(hjust = 0, size = 8, face = "bold")
    , axis.text = element_text(size=8)
    , axis.title = element_text(size=8)
    , legend.text = element_text(size=8)
    , legend.title = element_text(size=8)
    , legend.position = "bottom"
  )+
  guides(
    color = guide_legend(override.aes = list(alpha=1,size=1))
  )
grid.arrange(g1,g2)
