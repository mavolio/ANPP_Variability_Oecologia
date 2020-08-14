library(tidyverse)
library(gridExtra)
library(codyn)
library(rsq)
library(gtable)
library(grid)

setwd("C:\\Users\\mavolio2\\Dropbox\\converge_diverge\\datasets\\LongForm")

theme_set(theme_bw(10)) 

# read in the data ----------------------------------------------------------------

#all ANPP data
anpp<-read.csv("ANPP_Oct2017_2.csv")%>%
  select(-X)%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment, plot_id)%>%
  mutate(numyear=length(treatment_year))%>%
  filter(numyear>5)

#binning the treatments into categories
trtint<-read.csv('treatment interactions_ANPP_datasets_using.csv')%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  select(site_project_comm, treatment, trt_type7, trt_type5, trt_type6, trt_type, trt_type8)

#linking plots to treatments
anpp_expInfo<-read.csv("ExperimentInformation_ANPP_Dec2017.csv")%>%
  select(-X)


#getting site attributes MAP, MAT and rrich
site_info<-read.csv("SiteExperimentDetails_March2019.csv")%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  select(site_project_comm, site_code, project_name, community_type, MAP, MAT, rrich)

#yearly precip data
precip<-read.csv("C:\\Users\\mavolio2\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\real_precip_anppSites.csv")%>%
  mutate(calendar_year=year, precip_mm=precip)%>%
  select(-year, -X, -precip)


# Clean up anpp data and make calculations --------------------------------------------------------


###select the data to use

#for CDR e001/e002 selecting treatments , 6, 8, 9. For BGP dropping mowing treatments 
##all outliers were checked and are actual data.

dat2<-merge(anpp_expInfo, anpp, by=c("site_code","project_name","community_type","treatment"))%>%
  select(-nutrients, -light, -carbon, -water, -other_manipulation, -max_trt, -public, -factorial, -block)%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  mutate(delete=ifelse(site_code=="CDR"&treatment==2|site_code=="CDR"&treatment==3|site_code=="CDR"&treatment==4|site_code=="CDR"&treatment==5|site_code=="CDR"&treatment==7|site_code=="CDR"&treatment=="2_f_u_n"|site_code=="CDR"&treatment=="3_f_u_n"|site_code=="CDR"&treatment=="4_f_u_n"|site_code=="CDR"&treatment=="5_f_u_n"|site_code=="CDR"&treatment=="7_f_u_n"|project_name=="BGP"&treatment=="u_m_c"|project_name=="BGP"&treatment=="u_m_b"|project_name=="BGP"&treatment=="u_m_n"|project_name=="BGP"&treatment=="u_m_p"|project_name=="BGP"&treatment=="b_m_c"|project_name=="BGP"&treatment=="b_m_b"|project_name=="BGP"&treatment=="b_m_n"|project_name=="BGP"&treatment=="b_m_p"|project_name=="RHPs"&calendar_year==2003, 1, 0))%>%
  filter(delete!=1)

##NOTE KBS tilling treatments did not start until 1990, 2 years after the start of the N additions and control data.
# kbs<-dat2%>%
#   filter(site_code=="KBS")%>%
#   group_by(treatment, calendar_year)%>%
#   summarise(anpp=mean(anpp))%>%
#   ungroup%>%
#   group_by(treatment)%>%
#   summarize(n=length(calendar_year))

###fix sev data
nosev<-dat2%>%
  filter(site_project_comm!="SEV_Nfert_0")

sev<-dat2%>%
  filter(site_project_comm=="SEV_Nfert_0")%>%
  select(-treatment_year)%>%
  mutate(treatment_year= ifelse(calendar_year==2004, 10, ifelse(calendar_year==2005, 11, ifelse(calendar_year==2006, 12, ifelse(calendar_year==2007, 13, ifelse(calendar_year==2008, 14, ifelse(calendar_year==2009, 15, ifelse(calendar_year==2010, 16, ifelse(calendar_year==2011, 17, 18)))))))))

all_anpp_dat<-rbind(sev, nosev)

#Calculate temporal CV
anpp_temp_cv<-all_anpp_dat%>%
  group_by(site_code, project_name, community_type, treatment,plot_mani, plot_id)%>%
  summarize(anpp_temp_mean=mean(anpp, na.rm=T),
            anpp_temp_sd=sd(anpp, na.rm=T),
            anpp_temp_cv=(anpp_temp_sd/anpp_temp_mean)*100)

##calculating PD ANPP for each year
meandat<-all_anpp_dat%>%
  group_by(site_project_comm, plot_mani, treatment, treatment_year, calendar_year)%>%
  summarize(manpp=mean(anpp))

mcontrol<-meandat%>%
  filter(plot_mani==0)%>%
  mutate(contanpp=manpp)%>%
  select(-manpp)

mtrt<-meandat%>%
  filter(plot_mani!=0)

PD_anpp_yr<-merge(mtrt, mcontrol, by=c("site_project_comm","treatment_year","calendar_year"))%>%
  mutate(PD=((manpp-contanpp)/contanpp)*100,
         Diff=manpp-contanpp)%>%
  mutate(treatment=treatment.x)%>%
  left_join(site_info)

###Calculating PD of ANPP and CV of ANPP for each treatment plot
Means<-anpp_temp_cv%>%
  group_by(site_code, project_name, community_type, treatment, plot_mani)%>%
  summarize(anpp_temp_cv=mean(anpp_temp_cv), 
            anpp_temp_mean=mean(anpp_temp_mean),
            anpp_temp_sd = mean(anpp_temp_sd))%>%
  ungroup()%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))

cont_temp<-Means%>%
  filter(plot_mani==0)%>%
  rename(cont_temp_cv=anpp_temp_cv,
         cont_temp_sd=anpp_temp_sd,
         cont_temp_mean=anpp_temp_mean)%>%
  select(-plot_mani, -treatment)

CT_comp<-Means%>%
  filter(plot_mani!=0)%>%
  left_join(cont_temp)%>%
  left_join(trtint)%>%
  mutate(PD_CV=((anpp_temp_cv-cont_temp_cv)/cont_temp_cv)*100,
         PD_sd=((anpp_temp_sd-cont_temp_sd)/cont_temp_sd)*100,
         PD_mean=((anpp_temp_mean-cont_temp_mean)/cont_temp_mean)*100)

# Getting site characteristics --------------------------------------------

##getting average production
ave_prod<-all_anpp_dat%>%
  filter(plot_mani==0)%>%
  group_by(site_code, project_name, community_type, calendar_year)%>%
  summarize(anpp = mean(anpp))%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(manpp = mean(anpp))


###calculate community data
#get evenness for each plot
#To calcualte evenness
#getting evenness
anpp_spc<-all_anpp_dat%>%
  select(site_project_comm)%>%
  unique()

#read in community data
community<-read.csv("SpeciesRelativeAbundance_Nov2019.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  right_join(anpp_spc)
spc<-unique(community$site_project_comm)
rich_even<-data.frame()

for (i in 1:length(spc)){
  subset<-community%>%
    filter(site_project_comm==spc[i])
  
  out<-community_structure(subset, time.var = 'calendar_year', abundance.var = 'relcov', replicate.var = 'plot_id')
  out$site_project_comm<-spc[i]
  
  rich_even<-rbind(rich_even, out)
}

trt<-community%>%
  select(site_project_comm, plot_id, treatment)%>%
  unique()
plot_mani<-all_anpp_dat%>%
  select(site_project_comm, treatment, plot_mani)%>%
  unique()

ave_even<-rich_even%>%
  left_join(trt)%>%
  group_by(site_project_comm, plot_id, treatment)%>%
  summarize(Evar=mean(Evar, na.rm=T))%>%
  ungroup()%>%
  group_by(site_project_comm, treatment)%>%
  summarize(Evar = mean(Evar, na.rm=T))%>%
  right_join(plot_mani)

cont_even<-ave_even%>%
  filter(plot_mani==0)%>%
  select(-plot_mani, -treatment)

site_char<-site_info%>%
  right_join(cont_even)%>%
  left_join(ave_prod)

# Analysis 1 PD diff from zero for each treatment? -----------------------

anpp_temp_cv2<-anpp_temp_cv%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

spc<-unique(anpp_temp_cv2$site_project_comm)
ttest_out<-data.frame()

for (i in 1:length(spc)){
  
  subset<-anpp_temp_cv2%>%
    filter(site_project_comm==spc[i])
  
  trts<-unique(subset(subset, plot_mani!=0)$treatment)
  
  sub_controls<-subset(subset, plot_mani==0)%>%
    rename(cont_mean=anpp_temp_mean,
           cont_cv=anpp_temp_cv)%>%
    ungroup()%>%
    select(cont_mean, cont_cv, site_project_comm)
  
  
  for (j in 1:length(trts)){
    
    sub_treat<-subset(subset, treatment==trts[j])%>%
    ungroup()%>%
      select(anpp_temp_mean, anpp_temp_cv, site_project_comm, treatment)
    
    combined<-sub_treat%>%
      bind_cols(sub_controls)
    
    t.test.mean<-t.test(combined$cont_mean, combined$anpp_temp_mean)
    t.test.cv<-t.test(combined$cont_cv, combined$anpp_temp_cv)

    out<-data.frame(site_project_comm=spc[i],
                    treatment=trts[j],
                    t_mean=round(t.test.mean$statistic, digits=3),
                    p_mean=t.test.mean$p.value,
                    c.mean=t.test.mean$estimate[1], 
                    t.mean=t.test.mean$estimate[2],
                    t_cv=round(t.test.cv$statistic, digits=3), 
                    p_cv=t.test.cv$p.value,
                    c.cv=t.test.cv$estimate[1],
                    t.cv=t.test.cv$estimate[2])
    
    ttest_out<-rbind(ttest_out, out)
    }
  
}

ttest_summary<-ttest_out%>%
  mutate(meandiff=t.mean-c.mean,
         cvdiff=t.cv-c.cv)%>%
  mutate(resp_mean=ifelse(p_mean>0.05, "not sig", ifelse(p_mean<0.05&meandiff<0, "dec", ifelse(p_mean<0.05&meandiff>0, "inc", 999))),
         resp_cv=ifelse(p_cv>0.05, "not sig", ifelse(p_cv<0.05&cvdiff<0, "dec", ifelse(p_cv<0.05&cvdiff>0, "inc", 999))))%>%
  left_join(trtint)

mean.overall<-ttest_summary%>%
  group_by(resp_mean)%>%
  summarize(n=length(resp_mean))%>%
  mutate(response="A) Mean ANPP")%>%
  rename(effect=resp_mean)%>%
  mutate(trt_type8="All Trts")

mean.trt<-ttest_summary%>%
  group_by(trt_type8, resp_mean)%>%
  summarize(n=length(resp_mean))%>%
  mutate(response="A) Mean ANPP")%>%
  rename(effect=resp_mean)

cv.overall<-ttest_summary%>%
  group_by(resp_cv)%>%
  summarize(n=length(resp_cv))%>%
  mutate(response="B) CV ANPP")%>%
  rename(effect=resp_cv)%>%
  mutate(trt_type8="All Trts")

cv.trt<-ttest_summary%>%
  group_by(trt_type8, resp_cv)%>%
  summarize(n=length(resp_cv))%>%
  mutate(response="B) CV ANPP")%>%
  rename(effect=resp_cv)

##first overall for PD_CV
t.test(CT_comp$PD_CV, mu=0) # overall No, and not for the difference GCDs
t.test(CT_comp$PD_sd, mu=0) #yes overall sig for SD
t.test(CT_comp$PD_mean, mu=0)#yes overall sig for mean

##num postive or negative PD
sign<-CT_comp%>%
  mutate(pos=ifelse(PD_CV<0, 0,1))

#pos
sum(sign$pos)

45/95 #47% are postive and 53% are negative
50/95
##for not sig for all
irr<-subset(CT_comp, trt_type8=="Water")
t.test(irr$PD_CV, mu=0)
nit<-subset(CT_comp, trt_type8=="Nitrogen")
t.test(nit$PD_CV, mu=0)
nuts<-subset(CT_comp, trt_type8=="Multiple Nutrients")
t.test(nuts$PD_CV, mu=0)
int<-subset(CT_comp, trt_type8=="Interacting Drivers")
t.test(int$PD_CV, mu=0)

## for SD sig for all
# irr<-subset(CT_comp, trt_type6=="Water")
# t.test(irr$PD_sd, mu=0)
# nit<-subset(CT_comp, trt_type6=="Nitrogen")
# t.test(nit$PD_sd, mu=0)
# nuts<-subset(CT_comp, trt_type6=="Multiple Nutrients")
# t.test(nuts$PD_sd, mu=0)
# int<-subset(CT_comp, trt_type8=="Interacting Drivers")
# t.test(int$PD_sd, mu=0)


##for mean sig for all
irr<-subset(CT_comp, trt_type6=="Water")
t.test(irr$PD_mean, mu=0)
nit<-subset(CT_comp, trt_type6=="Nitrogen")
t.test(nit$PD_mean, mu=0)
nuts<-subset(CT_comp, trt_type6=="Multiple Nutrients")
t.test(nuts$PD_mean, mu=0)
int<-subset(CT_comp, trt_type8=="Interacting Drivers")
t.test(int$PD_mean, mu=0)

# Making figure 1 ---------------------------------------------------------
vote.fig<-mean.trt%>%
  bind_rows(cv.trt)%>%
  bind_rows(cv.overall)%>%
  bind_rows(mean.overall)%>%
  mutate(prop=ifelse(trt_type8=="Multiple Nutrients", n/33, ifelse(trt_type8=="Nitrogen", n/11, ifelse(trt_type8=="Water", n/7, ifelse(trt_type8=="Interacting Drivers", n/11, ifelse(trt_type8=="Other GCD", n/33, ifelse(trt_type8=="All Trts", n/95, 999)))))))

vot_mean<-ggplot(data=subset(vote.fig, response=="A) Mean ANPP"), aes(y=prop, x=trt_type8, fill=effect))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("Treatment")+
  ylab("")+
  scale_fill_manual(name="Treatment\n Response", label=c("Not Sig.", "Increase", "Decrease"), limits=c("not sig", "inc", "dec"), values = c("Gray", "lightblue", "darkblue"))+
  scale_x_discrete(limits=c("Other GCD", "Interacting Drivers", "Water", "Nitrogen", "Multiple Nutrients", "All Trts"))+
  geom_vline(xintercept = 5.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  ggtitle("a) Mean ANPP")

vot_cv<-ggplot(data=subset(vote.fig, response=="B) CV ANPP"), aes(y=prop, x=trt_type8, fill=effect))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("Treatment")+
  ylab("Prop. of Trts Diff. from Control")+
  scale_fill_manual(name="Treatment Response", label=c("Not Sig.", "Increase", "Decrease"), limits=c("not sig", "inc", "dec"), values = c("Gray", "skyblue", "darkblue"))+
  scale_x_discrete(limits=c("Other GCD", "Interacting Drivers", "Water", "Nitrogen", "Multiple Nutrients", "All Trts"))+
  geom_vline(xintercept = 5.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  ggtitle("b) CV ANPP")

#draw this and then remove legend position above.
legend=gtable_filter(ggplot_gtable(ggplot_build(vot_cv)), "guide-box") 
grid.draw(legend)

fig1a<-
  grid.arrange(vot_mean, legend, vot_cv, ncol=1, heights=c(1,0.5,1))



##making a bar graph of this
PD_bargraph_trt<-CT_comp%>%
  group_by(trt_type8)%>%
  summarize(cv=mean(PD_CV),
            sd_cv=sd(PD_CV),
            sd=mean(PD_sd),
            sd_sd=sd(PD_sd),
            mn=mean(PD_mean),
            sd_mn=sd(PD_mean),
            num=length(PD_CV))%>%
  mutate(se_cv=sd_cv/sqrt(num),
         se_sd=sd_sd/sqrt(num),
         se_mn=sd_mn/sqrt(num))%>%
  filter(trt_type8=="Nitrogen"|trt_type8=="Multiple Nutrients"|trt_type8=="Water"| trt_type8=="Interacting Drivers")

PD_bargraph_all<-CT_comp%>%
  summarize(cv=mean(PD_CV),
            sd_cv=sd(PD_CV),
            sd=mean(PD_sd),
            sd_sd=sd(PD_sd),
            mn=mean(PD_mean),
            sd_mn=sd(PD_mean),
            num=length(PD_CV))%>%
  mutate(se_cv=sd_cv/sqrt(num),
         se_sd=sd_sd/sqrt(num),
         se_mn=sd_mn/sqrt(num))%>%
  mutate(trt_type8="All Treatments")

PD_bargraph<-rbind(PD_bargraph_trt, PD_bargraph_all)


cv_fig<-ggplot(data=PD_bargraph, aes(x=trt_type8, y=cv, fill=trt_type8))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=cv-se_cv, ymax=cv+se_cv),position= position_dodge(0.9), width=0.2)+
  ylab("")+
  ylab("Difference (%) CV ANPP")+
  scale_x_discrete(limits = c("All Treatments",'Multiple Nutrients','Nitrogen','Water', "Interacting Drivers"),labels = c("All Trts", "Multiple\n Nutrients", "Nitrogen","Water", "Interacting"))+
  xlab("")+
  scale_fill_manual(values=c("darkred","orange","green3","blue", "black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  ggtitle("d) CV ANPP")

##this is figure S2
sd_fig<-ggplot(data=PD_bargraph, aes(x=trt_type8, y=sd, fill=trt_type8))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=sd-se_sd, ymax=sd+se_sd),position= position_dodge(0.9), width=0.2)+
  ylab("")+
  ylab("Percent Difference\nSD of ANPP")+
  scale_x_discrete(limits = c("All Treatments",'Multiple Nutrients','Nitrogen','Water', "Interacting Drivers"),labels = c("All Trts", "Multiple\n Nutrients", "Nitrogen","Water", "Interacting"))+
  xlab("")+
  scale_fill_manual(values=c("darkred","orange","green3","blue", "black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  geom_text(x=1, y=35, label="*", size=8)+
  geom_text(x=2, y=75, label="*", size=8)+
  geom_text(x=5, y=50, label="*", size=8)+
  scale_y_continuous(limits=c(0, 80))

mn_fig<-ggplot(data=PD_bargraph, aes(x=trt_type8, y=mn, fill=trt_type8))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=mn-se_mn, ymax=mn+se_mn),position= position_dodge(0.9), width=0.2)+
  ylab("")+
  ylab("Difference (%) Mean ANPP")+
  scale_x_discrete(limits = c("All Treatments",'Multiple Nutrients','Nitrogen','Water', "Interacting Drivers"),labels = c("All Trts", "Multiple\n Nutrients", "Nitrogen","Water", "Interacting"))+
  xlab("")+
  scale_fill_manual(values=c("darkred","orange","green3","blue", "black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  geom_text(x=1, y=35, label="*", size=8)+
  geom_text(x=2, y=55, label="*", size=8)+
  geom_text(x=3, y=35, label="*", size=8)+
  geom_text(x=4, y=48, label="*", size=8)+
  geom_text(x=5, y=35, label="*", size=8)+
  scale_y_continuous(limits=c(0, 60))+
  ggtitle("c) Mean ANPP")


fig1b<-
  grid.arrange(arrangeGrob(mn_fig+theme(legend.position="none"),
                           cv_fig+theme(legend.position="none"),
                           ncol=1))

figure1<-grid.arrange(fig1a,fig1b, ncol=2)

ggsave("C:\\Users\\mavolio2\\Dropbox\\Manuscripts\\Corre - ANPP\\Submit oecologia\\Revision\\Figure1.jpeg", figure1, dpi=800, units="mm", width=84*2, height=84*2)

# Analysis 2 abiotic and biotic drivers of PD ANPP and CV of ANPP --------------------


###what correlates with PD_CV?

PD_cor<-CT_comp%>%
  left_join(site_char)

#shouldn't use site ANPP as predictive, because it alone is correlated with CV of anpp for the treated and control plots
with(PD_cor, plot(manpp, cont_temp_cv))
with(PD_cor, plot(manpp, anpp_temp_cv))
with(PD_cor, cor.test(manpp, cont_temp_cv))
with(PD_cor, cor.test(manpp, anpp_temp_cv))
with(PD_cor, cor.test(MAT, MAP))

#library(MASS) # MASS masks select in tidyverse, so only load this when doing mutliple regressions

##how correlated are the predictor variables?
pairs(PD_cor[,c(21:24)])

stepAIC(lm(PD_CV~MAP+MAT+rrich+Evar, data=PD_cor))
summary(model.cv<-lm(PD_CV~MAP+Evar+MAT, data=PD_cor))
rsq.partial(model.cv, adj = T)

# stepAIC(lm(PC_sd~MAT+MAP+anpp+sdppt+cont_rich+Evar, data=PC_cor))
# summary(model.sd<-lm(PC_CV~MAP+anpp+sdppt, data=PC_cor))
# rsq.partial(model.sd, adj =T)

# stepAIC(lm(PD_mean~MAT+MAP+rrich+Evar, data=PD_cor))
# summary(model.mn<-lm(PD_CV~MAT+rrich+Evar, data=PD_cor))
# rsq.partial(model.mn)


# Making figure 2 ---------------------------------------------------------

tograph_cor<-PD_cor%>%
  select(site_project_comm, treatment,PD_CV, PD_sd, PD_mean, MAP, MAT, rrich, Evar)%>%
  gather(parm, value, MAP:Evar)%>%
  gather(vari_metric, vari_value, PD_CV:PD_mean)%>%
  mutate(parm_group=factor(parm, levels = c("rrich", "Evar","MAP","MAT")),
         vari_group=factor(vari_metric, levels=c("PD_mean","PD_sd","PD_CV")))

rvalues <- tograph_cor %>% 
  group_by(vari_group, parm_group) %>%
  summarize(r.value = round((cor.test(vari_value, value)$estimate), digits=3),
            p.value = (cor.test(vari_value, value)$p.value))

ttest_sig<-ttest_out%>%
  mutate(PD_mean=ifelse(p_mean>0.05,0,1),
         PD_CV=ifelse(p_cv>0.05, 0, 1))%>%
  gather(vari_metric, sig, PD_mean:PD_CV)%>%
  select(site_project_comm, treatment, vari_metric, sig)

tograph_cor2<-tograph_cor%>%
  filter(vari_metric!="PD_sd")%>%
  left_join(ttest_sig)
rvalues2<-rvalues %>% 
  filter(vari_group!="PD_sd")

srgraph<-
ggplot(data=subset(tograph_cor2, vari_metric=="PD_CV"&parm_group=="rrich" ), aes(x = value, y = vari_value))+
  geom_point(aes(color=as.factor(sig)))+
  scale_color_manual(name="", values=c("darkgray", 'blue'))+
  xlab("Sp. Richness")+
  ylab("Percent Difference in Temporal Variability")+
  annotate("text", x=Inf, y = Inf, hjust=1.05, vjust=1.5, label="-0.151")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  ggtitle("a) Sp. Richness")

evengraph<-
  ggplot(data=subset(tograph_cor2, vari_metric=="PD_CV"&parm_group=="Evar" ), aes(x = value, y = vari_value))+
  geom_point(aes(color=as.factor(sig)))+
  scale_color_manual(name="", values=c("darkgray", 'blue'))+
  xlab("Evenness")+
  ylab("")+
  geom_smooth(method="lm", se=F, color = "black")+
  annotate("text", x=Inf, y = Inf, hjust=1.05, vjust=1.5, label="0.259")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  ggtitle("b) Evenness")

MAPgraph<-
  ggplot(data=subset(tograph_cor2, vari_metric=="PD_CV"&parm_group=="MAP" ), aes(x = value, y = vari_value))+
  geom_point(aes(color=as.factor(sig)))+
  scale_color_manual(name="", values=c("darkgray", 'blue'))+
  xlab("Annual Precip.")+
  ylab("")+
  geom_smooth(method="lm", se=F, color = "black")+
  annotate("text", x=Inf, y = Inf, hjust=1.05, vjust=1.5, label="-0.497")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  ggtitle("c) MAP")

MATgraph<-
  ggplot(data=subset(tograph_cor2, vari_metric=="PD_CV"&parm_group=="MAT" ), aes(x = value, y = vari_value))+
  geom_point(aes(color=as.factor(sig)))+
  scale_color_manual(name="", values=c("darkgray", 'blue'))+
  xlab("Annual Temp.")+
  ylab("")+
  geom_smooth(method="lm", se=F, color = "black")+
  annotate("text", x=Inf, y = Inf, hjust=1.05, vjust=1.5, label="-0.612")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  ggtitle("d) MAT")

figure2<-grid.arrange(srgraph, evengraph, MAPgraph, MATgraph, ncol=4)

ggsave("C:\\Users\\mavolio2\\Dropbox\\Manuscripts\\Corre - ANPP\\Submit oecologia\\Revision\\Figure2.jpeg", figure2, dpi=800, units="mm", width=84*2, height=84)


# Analysis 3, sensitivity of anpp to precip --------------------------------------------------

#precip analysis #1986 in CDR has no precip data, this one year is being dropped.

#drop irrigation treatments becuase it is confusing to add water and then test against precip OR maybe just KNZ because not the same amount of water each year.

#look at average change in sensitivity by treatments.

anpp_precip<-merge(all_anpp_dat, precip, by=c("site_code","calendar_year"))%>%
  mutate(trt=ifelse(plot_mani==0,"C","T"))%>%
  group_by(site_project_comm, trt, calendar_year, precip_mm, treatment, plot_mani)%>%
  summarize(anpp=mean(anpp))%>%
  mutate(spc_trt=paste(site_project_comm, treatment, sep="::"))

# #visually inspecting the data
# ggplot(data=anpp_precip, aes(x=precip_mm, y=anpp, group=treatment, color=trt))+
#   geom_point()+
#   geom_smooth(method="lm", se=F)+
#   facet_wrap(~site_project_comm, ncol=8, scales = "free")


##get slopes for each treatment including controls
spc<-unique(anpp_precip$spc_trt)
lm.slopes<-data.frame()
for (i in 1:length(spc)){
  subset<-anpp_precip%>%
    filter(spc_trt==spc[i])
  test.lm<-lm(anpp~precip_mm, data=subset)
  output.lm<-data.frame(site_project_comm=unique(subset$site_project_comm), 
                        treatment=unique(subset$treatment), 
                        plot_mani=unique(subset$plot_mani), 
                        est=summary(test.lm)$coef["precip_mm", c("Estimate")], 
                        st.er=summary(test.lm)$coef["precip_mm", c("Std. Error")], 
                        p.val=summary(test.lm)$coef["precip_mm","Pr(>|t|)"])
  lm.slopes<-rbind(lm.slopes, output.lm)
}


#graphing this
c.slope<-lm.slopes%>%
  filter(plot_mani==0)%>%
  rename(c_est=est, c_se=st.er)%>%
  select(site_project_comm, c_est, c_se)

t.slope<-lm.slopes%>%
  filter(plot_mani!=0)%>%
  select(-p.val)

slopes_tograph<-c.slope%>%
  left_join(t.slope)%>%
  left_join(site_char)%>%
  mutate(diff=est-c_est)%>%
  separate(site_project_comm, into=c("site_code","project_name","community_type"), sep="_", remove=F)%>%
  left_join(trtint)

#overall ttest
t.test(slopes_tograph$diff, mu=0)

#do t-test do the slopes differ from zero?
irr<-subset(slopes_tograph, trt_type8=="Water")
t.test(irr$diff, mu=0)

nit<-subset(slopes_tograph, trt_type8=="Nitrogen")
t.test(nit$diff, mu=0)

nuts<-subset(slopes_tograph, trt_type8=="Multiple Nutrients")
t.test(nuts$diff, mu=0)

int<-subset(slopes_tograph, trt_type8=="Interacting Drivers")
t.test(int$diff, mu=0)


# Making figure 3 ---------------------------------------------------------
slopes_bar_trt<-slopes_tograph%>%
  left_join(trtint)%>%
  group_by(trt_type8)%>%
  summarise(mdiff=mean(diff),
            ndiff=length(diff),
            sddiff=sd(diff))%>%
  mutate(sediff=sddiff/sqrt(ndiff))%>%
  filter(trt_type8=="Nitrogen"|trt_type8=="Multiple Nutrients"|trt_type8=="Water"|trt_type8=="Interacting Drivers")

slopes_bar_overall<-slopes_tograph%>%
  summarise(mdiff=mean(diff),
            ndiff=length(diff),
            sddiff=sd(diff))%>%
  mutate(sediff=sddiff/sqrt(ndiff))%>%
  mutate(trt_type8="All Treatments")


slopes_bar<-rbind(slopes_bar_overall, slopes_bar_trt)


figure3<-ggplot(data=slopes_bar, aes(x=trt_type8, y=mdiff, fill=trt_type8))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=mdiff-sediff, ymax=mdiff+sediff),position= position_dodge(0.9), width=0.2)+
  xlab("")+
  ylab("Difference in Slopes")+
  #theme(axis.text.x=element_text(angle=45, hjust=1))+
  scale_x_discrete(limits = c("All Treatments",'Multiple Nutrients','Nitrogen','Water', "Interacting Drivers"),labels = c("All Trts", "Multiple\n Nutrients", "Nitrogen","Water", "Interacting"))+
  scale_fill_manual(values=c("black","darkred","orange","green3", "blue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  geom_text(x=1, y=0.255, label="*", size=8)+
  geom_text(x=2, y=0.24, label="*", size=8)


ggsave("C:\\Users\\mavolio2\\Dropbox\\Manuscripts\\Corre - ANPP\\Submit oecologia\\Revision\\Figure3.jpeg", figure3, dpi=800, units="mm", width=84, height=85)

# appendix Figure S1 analysis PD anpp over time ---------------------------------------------------------
fig<-PD_anpp_yr%>%
  left_join(trtint)%>%
  mutate(id=paste(site_project_comm, treatment, sep="::"))

ggplot(data=fig, aes(x=treatment_year, y=PD))+
  geom_smooth(method="lm", formula = y ~  poly(x, 2), size=0.5, color="gray", aes(group=id), se=F)+
  geom_point(size=0.05, color="gray")+
  geom_smooth(method="lm", formula = y ~  poly(x, 2), size=1, color="black", se=F)+
  geom_hline(yintercept=0)+
  xlab("Treatment Year")+
  ylab("Percent Difference in ANPP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~trt_type5, ncol=5, scales="free")

# appendix Figure S3 correlating trt amount with reponses ------------------------------------
trtdetails<-all_anpp_dat%>%
  select(-anpp, -plot_id, -calendar_year, -numyear, -treatment_year)%>%
  filter(plot_mani!=0)%>%
  unique()

CT_comp2<-CT_comp%>%
  left_join(trtdetails)%>%
  group_by(site_project_comm, treatment)%>%
  mutate(multnuts=sum(n, p, k))

#nitrogen - not sig for mean
with(subset(CT_comp2, trt_type6=="Nitrogen"), cor.test(n, PD_mean))# p 0.66
with(subset(CT_comp2, trt_type6=="Nitrogen"), plot(n, PD_mean))
#not sig for CV
with(subset(CT_comp2, trt_type6=="Nitrogen"), cor.test(n, PD_CV))# p 0.32
with(subset(CT_comp2, trt_type6=="Nitrogen"), plot(n, PD_CV))

#precip - not sig for mean
with(subset(CT_comp2, trt_type6=="Water"), cor.test(precip, PD_mean))# p 0.09
with(subset(CT_comp2, trt_type6=="Water"), plot(precip, PD_mean))
#not sig for CV
with(subset(CT_comp2, trt_type6=="Water"), cor.test(precip, PD_CV))# p 0.60
with(subset(CT_comp2, trt_type6=="Water"), plot(precip, PD_CV))

#temp - all only 1 C increase, but lots of vari
with(subset(CT_comp2, trt_type6=="Heat"), cor.test(temp, PD_mean))
with(subset(CT_comp2, trt_type6=="Heat"), plot(temp, PD_mean))

#mult nuts - very sig for mean
with(subset(CT_comp2, trt_type6=="Multiple Nutrients"), cor.test(multnuts, PD_mean))#SIG
with(subset(CT_comp2, trt_type6=="Multiple Nutrients"), plot(multnuts, PD_mean))
## not sig for CV
with(subset(CT_comp2, trt_type6=="Multiple Nutrients"), cor.test(multnuts, PD_CV))#0.15
with(subset(CT_comp2, trt_type6=="Multiple Nutrients"), plot(multnuts, PD_CV))

Nit_mean<-
  ggplot(data=subset(CT_comp2, trt_type6=='Nitrogen'), aes(x=n, y=PD_mean))+
  geom_point(size = 2)+
  xlab("N added (g m-2)")+
  ylab("Difference (%) Mean ANPP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Nit_cv<-
  ggplot(data=subset(CT_comp2, trt_type6=='Nitrogen'), aes(x=n, y=PD_CV))+
  geom_point(size = 2)+
  xlab("N added (g m-2)")+
  ylab("Difference (%) CV ANPP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Wat_mean<-
  ggplot(data=subset(CT_comp2, trt_type6=='Water'), aes(x=precip, y=PD_mean))+
  geom_point(size = 2)+
  xlab("% Increase in Precipitation")+
  ylab("Difference (%) Mean ANPP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Wat_cv<-
  ggplot(data=subset(CT_comp2, trt_type6=='Water'), aes(x=precip, y=PD_CV))+
  geom_point(size = 2)+
  xlab("% Increase in Precipitation")+
  ylab("Difference (%) CV ANPP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
mult_mean<-
  ggplot(data=subset(CT_comp2, trt_type6=='Multiple Nutrients'), aes(x=multnuts, y=PD_mean))+
  geom_point(size = 2)+
  xlab("Summed NPK added (g m-2)")+
  ylab("Difference (%) Mean ANPP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(method="lm", se=F, color="black")
mult_cv<-
  ggplot(data=subset(CT_comp2, trt_type6=='Multiple Nutrients'), aes(x=multnuts, y=PD_CV))+
  geom_point(size = 2)+
  xlab("Summed NPK added (g m-2)")+
  ylab("Difference (%) CV ANPP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

grid.arrange(Nit_mean, Nit_cv,Wat_mean, Wat_cv,mult_mean,   mult_cv)




# Making graphical abstract -----------------------------------------------

boxfig<-CT_comp%>%
  select(site_project_comm, treatment, PD_CV, PD_mean)%>%
  gather(measure, value, PD_CV:PD_mean)

ggplot(data=boxfig, aes(x=measure, y=value))+
  geom_boxplot(outlier.color="black", outlier.size=2)+
  geom_jitter(aes(color=measure))+
  ylab("Difference (% GCD Treatment Effect)")+
  xlab("")+
  scale_x_discrete(limits=c("PD_mean", "PD_CV"), labels = c("Mean ANPP", "ANPP CV"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_hline(yintercept = 0, size = 1)