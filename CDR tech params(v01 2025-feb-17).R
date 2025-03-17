#' Plot assumptions and parameters for CDR tech
#' 

#### Libraries ####
#' Technically gdxtools, but since I first run some part of plotgdx_witch. R it loaded them
#' 
#' 

#------------------------------------------------------------------------------------------------

gdx_path<- '/Users/cindyazuero/Documents/WITCH model/witch-master/all_data_temp_ssp2_bau.gdx'
gdx_in_R<-gdx(gdx_path)

graphsPath<-'/Users/cindyazuero/Dropbox (Personal)/Otros CGA/RFF-CMCC EIEE/Projects/UPTAKE/WP3/Technologies info/CCS'

#' Get CCS info from gdx
ccs_stor_cap_max<-gdx_in_R["ccs_stor_cap_max"]
ccs_floor_cost<-gdx_in_R["ccs_floor_cost"]
mcost_inv0<-gdx_in_R["mcost_inv0"]
ccs_wcum0 <-gdx_in_R["ccs_wcum0"]
ccs_learn<- gdx_in_R["ccs_learn"]

graphMax_CCS_storaage<- function(){
  dataGraph<-ccs_stor_cap_max
  
  png(filename = paste(graphsPath, "maxCCSstorage_master.png", sep="/"), width=8, height=4, units = "in", res=200 )
  p<-ggplot()+
    geom_bar(data=dataGraph[dataGraph$ccs_stor_estim=="best",], aes(x=ccs_stor, y=value, fill= ccs_stor),stat="identity")+
    geom_point(data=dataGraph[dataGraph$ccs_stor_estim=="low",], aes(x=ccs_stor, y=value), color="red")+
    geom_point(data=dataGraph[dataGraph$ccs_stor_estim=="high",], aes(x=ccs_stor, y=value), color="blue")+
    labs(x="",y="GtCO2", fill= "Storage type")+
    facet_wrap(~n, scales="free_y")+
    theme( axis.text.y = element_text(size = 14),
           strip.text = element_text(size=12))
  print(p)
  dev.off()  
}



graph_MCOST_INV_curve<-function(jlccs, t=30, n)
{
  browser()
  #' Filter data by chosen options
  ccs_floor_cost_g<-ccs_floor_cost[ccs_floor_cost$V1==jlccs & ccs_floor_cost$n==n,"value"]
  mcost_inv0_g<-mcost_inv0[mcost_inv0$j==jlccs & mcost_inv0$n==n & mcost_inv0$t==t,"value"]
  ccs_wcum0_g<-ccs_wcum0[ccs_wcum0$jlccs==jlccs,"value"]
  ccs_learn_g<-ccs_learn[ccs_learn$jlccs==jlccs & ccs_learn$t==t,"value"]
  
  my_function<-function(x){
    max(ccs_floor_cost_g, mcost_inv0_g*(x/ccs_wcum0_g)**(-ccs_learn_g))
  }
  
  #' Graph function
  ggplot(data.frame(x = c(0, 2000)), aes(x)) + 
    stat_function(fun =my_function, color = "blue") +
    labs(title = "MCOST_INV(ccs_wcum_spill))", x = "TW", y = "T$/TW")
}

graph_MCOST_INV_curve("elbigcc", 30, "brazil")
