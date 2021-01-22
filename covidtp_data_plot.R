# Hitchings et al - The usefulness of SARS-CoV-2 test positive proportion as a surveillance tool
# Supplementary Code - Processing and plotting data from COVID Tracking Project
require(ggplot2)
require(ggrepel)
require(grid)
require(gridExtra)
require(reshape2)
require(dplyr)
require(httr)
require(jsonlite)

theme_set(theme_bw())

# Get data from COVID Tracking Project
state_tests_historical = GET("https://api.covidtracking.com/v1/states/daily.json")
state_tests_historical = fromJSON(rawToChar(state_tests_historical$content))

# .csv of state populations
state_pops = read.csv('state_pops.csv',header=T)
state_tests_historical$Population = sapply(state_tests_historical$state,function(x) 
  ifelse(length(state_pops$POPESTIMATE2019[as.character(state_pops$StateCode)==as.character(x)])>0,
         state_pops$POPESTIMATE2019[as.character(state_pops$StateCode)==as.character(x)],NA))
state_tests_historical$date = as.Date(as.character(state_tests_historical$date),format='%Y%m%d')

# Crude estimate of prevalence based on detected tests
# To normalize the start time in each state relative to the start of the outbreak
# Mean of 8 days from infection to confirmation
# Use distribution of time from infection to confirmation (log-normal) to estimate how many future detected cases are prevalent on each day
mean_inftoconf=log(8)
sd_inftoconf=log(1.35)

for (s in unique(state_tests_historical$state)) {
  state_tests_historical$prevalent_cases[state_tests_historical$state==s] = 
    sapply(state_tests_historical$date[state_tests_historical$state==s],function(d) 
      round(sum(state_tests_historical$positiveIncrease[state_tests_historical$state==s & state_tests_historical$date>d & 
                                                          state_tests_historical$date<=d+15]*
                  (1-plnorm(1:length(state_tests_historical$positiveIncrease[state_tests_historical$state==s & state_tests_historical$date>d & 
                                                                               state_tests_historical$date<=d+15]),mean_inftoconf,sd_inftoconf)))))
}
state_tests_historical$IP = state_tests_historical$prevalent_cases/state_tests_historical$Population
state_tests_historical$TPR = state_tests_historical$positiveIncrease/
  (state_tests_historical$negativeIncrease+state_tests_historical$positiveIncrease)

# First date on which each state hit a given prevalence threshold
min_prev = 1/10000
firstdates = sapply(unique(state_tests_historical$state),function(s) 
  state_tests_historical$date[max(which(state_tests_historical$state==s & state_tests_historical$IP>=min_prev))])
names(firstdates)=unique(state_tests_historical$state)
firstdates=as.Date(firstdates,origin='1970-01-01')

state_tests_historical$firstdate = sapply(1:nrow(state_tests_historical),function(s) 
  firstdates[names(firstdates)==state_tests_historical$state[s]])
state_tests_historical$firstdate = as.Date(state_tests_historical$firstdate,origin='1970-01-01')

# Aggregate to weekly since first date
for (s in unique(state_tests_historical$state)) {
  state_tests_historical$positivetests_pastweek[state_tests_historical$state==s] = 
    sapply(state_tests_historical$date[state_tests_historical$state==s],function(d) 
      sum(state_tests_historical$positiveIncrease[state_tests_historical$state==s & state_tests_historical$date>=d & 
                                                    state_tests_historical$date<=d+6]))
  state_tests_historical$tests_pastweek[state_tests_historical$state==s] = 
    sapply(state_tests_historical$date[state_tests_historical$state==s],function(d) 
      sum(state_tests_historical$totalTestResultsIncrease[state_tests_historical$state==s & state_tests_historical$date>=d & 
                                                            state_tests_historical$date<=d+6]))
}
state_tests_historical$TPRweekly=state_tests_historical$positivetests_pastweek/state_tests_historical$tests_pastweek

# Function to "fit" TPP to confirmed incidence (Equation (2) in the paper)
max_ir = function(tpr,da,ds,pili,sens) {
  sapply(tpr,function(x) x*(da/ds*(1-pili)+pili)/(sens-x*(1-da/ds)*(1-pili)))
}

expit=function(x){1/(1+exp(-x))}
expit_half = function(x){0.5/(1+exp(-x))}
logit=function(x){log(x/(1-x))}
f = function(par,data) {
  ds=expit(par[['da']])+exp(par[['ds']])
  pred_ci = log(10000*ds*expit(par[['sens']])*
                  max_ir(data$TPR,expit(par[['da']]),ds,expit_half(par[['pili']]),
                         expit(par[['sens']])))
  ci = log(data$CI)
  return(sum((pred_ci-ci)^2))
}
init_params=list("da"=logit(3e-5),"ds"=logit(0.5),"pili"=logit(0.011*2),"sens"=logit(0.95))

# The parameter values from the fit are meaningless as there is not a unique solution
# And note that the "fit" will not work for all states below. Many states have irregularities in the data which means TPP<0 for some dates.
# Some manipulation of the data is required to deal with these cases (e.g. exclusion)
make_states_crosssectional_plot = function(days_since_epidemic_start,states_to_label) {
  # Input the number of days at which to view cross-sectional incidence and TPP (relative to start of epidemic in each state)
  # States to label in the plot
  dat = 
    state_tests_historical %>% 
    filter(!is.na(Population) & state != "PR" & date==firstdate+days_since_epidemic_start & !is.na(firstdate)) %>% 
    mutate(CI=10000*positivetests_pastweek/(7*Population),
           TPR=TPRweekly,
           ) %>% 
    group_by(state) %>%
    summarise(TPR=mean(TPR),
              CI=mean(CI))
  res=optim(init_params,f,data=dat)
  par=expit(res$par[!(names(res$par) %in% c("pili","ds"))])
  par = c(par,'ds'=par[['da']]+exp(res$par[['ds']]))
  par = c(par,'pili'=expit_half(res$par[['pili']]))
  dat$state = as.character(dat$state)
  dat$state[!(dat$state %in% states_to_label)]=""
  
  p=ggplot(dat)+
    scale_y_continuous(breaks=(-3):0,labels=c('0.001','0.010','0.100','1.000'),limits=c(-3,0))+
    scale_x_continuous(breaks=(-3):1,labels=c('0.001','0.010','0.100','1.000','10.000'),limits=c(-2.5,1.5))+
    geom_point(aes(y=log10(TPR),x=log10(CI)))+
    geom_text_repel(aes(y=log10(TPR),x=log10(CI),label=state),hjust=-0.5,size=2,show.legend = F)+
    stat_function(data=data.frame(confcase=1),aes(x=10000*confcase),
                  fun=function(x) pmin(1,log10(10^(x)/10000/(par[['ds']]*par[['pili']]+par[['da']]*(1-par[['pili']])+
                                                               (1-par[['da']]/par[['ds']])*(1-par[['pili']])/
                                                               par[['sens']]*10^(x)/10000))))+
    theme(plot.margin = margin(t=0.1,l=0.5,r=0.2,unit='cm')
          , panel.grid = element_blank()
          , panel.border = element_blank()
          , axis.line = element_line(size = 0.353)
          , strip.background = element_blank()
          , strip.text = element_text(hjust = 0, size = 8, face = "bold")
          , axis.text = element_text(size=8)
          , axis.title = element_text(size=8)
          , plot.title = element_text(size=8)
          , legend.position = ""
    )+
    ylab('Test Positive Proportion')+
    xlab('Confirmed Incidence per 10,000 Persons')+
    ggtitle(paste0(days_since_epidemic_start,' days since detected infectious prevalence reached 0.01%'))
  
  p
}
make_states_crosssectional_plot(250,c("FL","NY","CA","MA","RI","IL","TX","ND","AZ"))

make_state_time_series_plot = function(s,lastdate,sname) {
  dat = 
    state_tests_historical %>% 
    filter(!is.na(Population) & state == s & date<=lastdate & 
             date>=firstdate+21) %>% 
    mutate(CI=10000*positivetests_pastweek/(7*Population),
           TPR=TPRweekly) %>% 
    group_by(date) %>%
    summarise(TPR=mean(TPR),
              CI=mean(CI))
  res=optim(init_params,f,data=dat)
  par=expit(res$par[!(names(res$par) %in% c("pili","ds"))])
  par = c(par,'ds'=par[['da']]+exp(res$par[['ds']]))
  par = c(par,'pili'=expit_half(res$par[['pili']]))
  
  dat_long = melt(dat %>% 
                    mutate(transTPR=10000*par[['ds']]*par[['sens']]*max_ir(TPR,par[['da']],par[['ds']],par[['pili']],par[['sens']])) %>% 
                    select(date,transTPR,CI),id='date')
  dat_long$variable=relevel(dat_long$variable,ref='CI')
  
  
  p=ggplot(dat_long,aes(x=date,y=value,linetype=variable))+
    geom_line(size=0.353)+
    scale_linetype_manual(values=c(1,2),labels=c('Observed','TPP-estimated'),name=NULL)+
    theme(plot.margin = margin(t=0.1,l=0.5,r=0.1,unit='cm')
          , panel.grid = element_blank()
          , panel.border = element_blank()
          , axis.line = element_line(size = 0.353)
          , strip.background = element_blank()
          , strip.text = element_text(hjust = 0, size = 8, face = "bold")
          , axis.text = element_text(size=8)
          , axis.title = element_text(size=8)
          , plot.title = element_text(size=8)
          , legend.position = c(0.8,0.8)
          , legend.margin = margin(r=0.1,l=0.1,unit='cm')
          , legend.key = element_blank()
          , legend.title = element_text(size=8)
          , legend.text = element_text(size=8)
          , legend.background = element_rect(linetype = 1,color='black',size=0.353)
          , legend.title.align = 0.5
          )+
    ylab('Confirmed Incidence per 10,000 Persons')+xlab('Date')+ggtitle(sname)
  
  p
}
make_state_time_series_plot('FL',as.Date('2020-12-31',format='%Y-%m-%d'),'Florida')





