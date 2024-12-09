#' ---
#' title: "Constructing CDR graphs"
#' author: "Cindy Azuero"
#' date: "Oct 4, 2024"
#' 
#' Meant to construct the graphs for the CDR Curves
#' 
#' This file must be called inside plotgdx_witch.R or at least ran after plotgdx_witch.R
#' ---

library(roxygen2)
#### q_emi_co2_beccs ####

#-> Calculate the parameter outside gams

#--- Obtain data---
stop_nash<-get_witch("stop_nash") # If this parameter is =1 it means the nash loop converged!
Q_IN<- get_witch("Q_IN")
emi_sys<- get_witch("emi_sys")
ccs_emi_capt<- get_witch("ccs_emi_capt")
ccs_capture_eff<-get_witch("ccs_capture_eff")
model_ctax<-get_witch("ctax")

#--- Filter to the ones I need---
Q_IN<-Q_IN[Q_IN$fuel=="wbio" & Q_IN$jfed=="elbigcc",]
emi_sys<-emi_sys[emi_sys$sys=="co2_ffi",]
ccs_emi_capt<-ccs_emi_capt[ccs_emi_capt$fuel=="wbio",]
ccs_capture_eff<-ccs_capture_eff[ccs_capture_eff$jfed=="elbigcc"]
model_ctax<-model_ctax[model_ctax$e=="co2_ffi"]

#--- Change name of the variables---
colnames(Q_IN)[colnames(Q_IN)=="value"]="Q_IN"
colnames(emi_sys)[colnames(emi_sys)=="value"]="emi_sys"
#colnames(ccs_emi_capt)[colnames(ccs_emi_capt)=="value"]="ccs_emi_capt"
#colnames(ccs_capture_eff)[colnames(ccs_capture_eff)=="value"]="ccs_capture_eff"
colnames(model_ctax)[colnames(model_ctax)=="value"]="ctax"

q_emi_co2_beccs<-merge(Q_IN,emi_sys,all.x=TRUE)
q_emi_co2_beccs<-merge(q_emi_co2_beccs,model_ctax,all.x=TRUE) # Those years and scenarios without ctax are with NA

#Since ccs_emi_capt and ccs_capture_eff are a constant across models, t, n I'll just take the value
# for BAU
ccs_emi_capt_value<-ccs_emi_capt[1,"value"]
ccs_capture_eff_value<- ccs_capture_eff[1,"value"]

q_emi_co2_beccs$ccs_emi_capt<-ccs_emi_capt_value
q_emi_co2_beccs$ccs_capture_eff<-ccs_capture_eff_value

#--- Calculate var ---
q_emi_co2_beccs$value<-q_emi_co2_beccs[,"Q_IN"]*q_emi_co2_beccs[,"emi_sys"]*q_emi_co2_beccs[,"ccs_emi_capt"]*q_emi_co2_beccs[,"ccs_capture_eff"]

# Add the model ctax value

# Attempt to use an existing function to create the global aggregation
#test<-make_global_tr(q_emi_co2_beccs,cols = c("t", "file")) # didn't work

# Obtain world as the aggregation of all regions n for each pair of (t, file)
world_q_emi_co2_beccs <- q_emi_co2_beccs %>% group_by_at(c("t", "file")) %>% summarize(value=sum(value)) %>% mutate(n="World")
# Since in this ran, all regions have the same tax, the tax used for world will be the one of Brazil
world_q_emi_co2_beccs<-merge(world_q_emi_co2_beccs, model_ctax[model_ctax$n=="brazil",  c("t", "file", "ctax")],all.x=TRUE)
world_q_emi_co2_beccs<-world_q_emi_co2_beccs[,c(1,4,2,5,3)]

# Create simplified dataset
q_emi_co2_beccs2<-q_emi_co2_beccs[,c("t", "n", "file", "ctax", "value")]

q_emi_co2_beccs2<-rbind(q_emi_co2_beccs2, world_q_emi_co2_beccs)

# Replace ctax of NA with 0
q_emi_co2_beccs2<- q_emi_co2_beccs2 %>% mutate(across(ctax, ~ replace(., is.na(.), 0)))

# Add a column with ctax_initial level

# Correct the BAU file rows
q_emi_co2_beccs2<- q_emi_co2_beccs2 %>% mutate(file=case_when(
  file == "results_results_ssp2_bau" ~ "results_ssp2_bau",  # Replace the specific value
  TRUE ~ file  # Keep the original value for all other cases
))


q_emi_co2_beccs2<- q_emi_co2_beccs2 %>% mutate(ctax_final= case_when(
  file == "results_ssp2_bau" ~ 0,
  TRUE ~ as.integer(substring(file,19, length(file)))
  
))

#' @description
#' Removes from the dataset the rows that correspond to those runs (i.e file) for which the nash loop
#' did not converge
#' 
filterRuns_NashNotConverged<-function(dataset){
  dataset<-dataset[apply(as.matrix(dataset$file),1,function(x) all (x %in% as.matrix(stop_nash[stop_nash$value==1,"file"]))),]
  dataset
}


#--- Make the graph---

graph_q_emi_co2_beccs_pathway<-function(){
  
  #Filter those runs that didn't converge
  graphData<-q_emi_co2_beccs2[apply(as.matrix(q_emi_co2_beccs2$file),1,function(x) all (x %in% as.matrix(stop_nash[stop_nash$value==1,"file"]))),]
  ggplot(data=graphData, aes(x=ttoyear(t), y=value, group=ctax_final))+
    geom_line(aes(color=ctax_final))+
    labs(y= "[GtC/yr]", title = "CO2 emissions from bioenergy with CCS" )+
    facet_wrap(~n,scales="free_y") +
     theme(plot.title=element_text(size=14),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
   #       panel.background = element_blank(),
           axis.title.x = element_blank(), #element_text(size=12),
           axis.title.y = element_text(size = 14),
           axis.text = element_text(size=14),
           legend.text = element_text(size=12),
           legend.title = element_text(size=12),
           legend.position = "bottom",
    #       #legend.spacing = unit(3, "mm"),
           legend.key.width = unit(3, "cm"),
           strip.text = element_text(size=14)
     )
    
}

graph_q_emi_co2_beccs_pathway()
saveplot("q_emi_co2_beccs pathway", width = 15, height=10)

graph_ctax_vs_q_emi_beccs_co2<-function(){
  graphData<- q_emi_co2_beccs2
  graphData<-graphData[apply(as.matrix(graphData$file),1,function(x) all (x %in% as.matrix(stop_nash[stop_nash$value==1,"file"]))),]

  ggplot(data=graphData, aes(x=ctax, y=value))+
    geom_point(aes(color=t))+
    labs(y= "[GtC/yr]", x= "ctax [T$/GTonC]", title = "CO2 emissions from bioenergy with CCS" )+
    geom_smooth(method = "loess", 
                formula = y ~ x)+ 
    facet_wrap(~n,scales="free_y") +
    theme(plot.title=element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #       panel.background = element_blank(),
          axis.title.x = element_text(size=12),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size=14),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          legend.position = "bottom",
          #       #legend.spacing = unit(3, "mm"),
          legend.key.width = unit(3, "cm"),
          strip.text = element_text(size=14)
    )
}

graph_ctax_vs_q_emi_beccs_co2()
saveplot("q_emi_co2_beccs vs ctax", width = 15, height=10)

#' @param dataVar the dataframe of the var I want to create the world value
#' @return dataframe with columns t, n, file, ctax, value where n="World"
getWorld<-function(dataVar, idIndeces){
  idIndeces2<-idIndeces[idIndeces!="n"]
  indecesFilter<-c(idIndeces2,"file")
  # Obtain world as the aggregation of all regions n for each pair of (t, file)
  world_var <- dataVar %>% group_by_at(indecesFilter) %>% summarize(value=sum(value)) %>% mutate(n="World")
  # Since in this ran, all regions have the same tax, the tax used for world will be the one of Brazil
  world_var<-merge(world_var, model_ctax[model_ctax$n=="brazil",  c("t", "file", "ctax")],all.x=TRUE)
  indecesOrder<-c(idIndeces, "file", "ctax", "value")
  world_var<-world_var[,indecesOrder]
  world_var
}

#' @param indices text with the indices. e.g.: ""t","n""
#' @param varName the name of the variable that I want to graph
#' @return a dataframe with the data required to create a graph using ctax and the value of the var, including the "World" region
getVarGraphData<-function(varName, indicesM){
  var<-get_witch(varName)
  var<-merge(var, model_ctax[,c("t", "n","file", "ctax")], all.x=TRUE) # Those years and scenarios without ctax are with NA
  #browser()
  cols<-c(indicesM,"file", "ctax", "value")
  varSimplified<-var[,..cols]
  #varSimplified<-var[,c("t", "n", "file", "ctax", "value")]
  
  world_var<-getWorld(var, indicesM)
  
  varSimplified<-rbind(varSimplified,world_var)
  
  # Replace ctax of NA with 0
  varSimplified<- varSimplified %>% mutate(across(ctax, ~ replace(., is.na(.), 0)))
  
  # Add a column with ctax_initial level
  
  # Correct the BAU file rows
  varSimplified<- varSimplified %>% mutate(file=case_when(
    file == "results_results_ssp2_bau" ~ "results_ssp2_bau",  # Replace the specific value
    TRUE ~ file  # Keep the original value for all other cases
  ))
  
  
  varSimplified<- varSimplified %>% mutate(ctax_final= case_when(
    file == "results_ssp2_bau" ~ 0,
    TRUE ~ as.integer(substring(file,19, length(file)))
    
  ))
  
  varSimplified
}

#' @description
#' Use this function to construct a graph with ctax in the x-axis and the var of interest in the y-axis.
#' Use this one if the variable does not have additional indices to (t,n) that need a filter, if need a filter use
#' graph_ctax_vs_var_data()
#' 
#'@param varName the name of the variable that I want to graph
#'@return creates a graph with ctax on the x-axis, the value of the var in the y-axis,
#' all years t as data points with different color
graph_ctax_vs_var<-function(varName, graphTitle, varUnit){
  dataGraph<- getVarGraphData(varName, c("t", "n"))
  dataGraph<-dataGraph[apply(as.matrix(dataGraph$file),1,function(x) all (x %in% as.matrix(stop_nash[stop_nash$value==1,"file"]))),]
  dataGraph<-dataGraph[dataGraph$t==30]
  
  
  ggplot(data=dataGraph, aes(x=ctax, y=value))+
    geom_point(aes(color=t))+ 
    labs(y= varUnit, , x= "ctax [T$/GTonC]", title = graphTitle )+
    facet_wrap(~n,scales="free_y") +
    geom_smooth(method = "loess", 
                formula = y ~ x)+ 
    theme(plot.title=element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #       panel.background = element_blank(),
          axis.title.x = element_text(size=12),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size=14),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          legend.position = "bottom",
          #       #legend.spacing = unit(3, "mm"),
          legend.key.width = unit(3, "cm"),
          strip.text = element_text(size=14)
    )
}
# https://www.geeksforgeeks.org/add-regression-line-to-ggplot2-plot-in-r/

#' @description
#' Use this function to construct a graph with ctax in the x-axis and the var of interest in the y-axis.
#' Use this one if the variable DO have additional indices to (t,n) that need a filter, if no need for additional filter use
#' graph_ctax_vs_var()
#' 
#'@param dataGraph the dataset that needs to be graphed, already with necessary filters
#'@return creates a graph with ctax on the x-axis, the value of the var in the y-axis,
#' all years t as data points with different color
graph_ctax_vs_var_data<-function(dataGraph, graphTitle, varUnit){
  ggplot(data=dataGraph, aes(x=ctax, y=diff(value)/5))+
    geom_point(aes(color=t))+ 
    labs(y= varUnit, , x= "ctax [T$/GTonC]", title = graphTitle )+
    facet_wrap(~n,scales="free_y") +
    geom_smooth(method = "loess", 
                formula = y ~ x)+ 
    theme(plot.title=element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #       panel.background = element_blank(),
          axis.title.x = element_text(size=12),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size=14),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          legend.position = "bottom",
          #       #legend.spacing = unit(3, "mm"),
          legend.key.width = unit(3, "cm"),
          strip.text = element_text(size=14)
    )
}

graph_ctax_vs_var("K_DAC", "DAC installed capacity", "[GtC]")
saveplot("K_DAC vs ctax", width = 15, height=10)



graph_var_pathway<- function(varName, graphTitle, varUnit){
  dataGraph<- getVarGraphData(varName, c("t", "n"))
  dataGraph<-dataGraph[apply(as.matrix(dataGraph$file),1,function(x) all (x %in% as.matrix(stop_nash[stop_nash$value==1,"file"]))),]
  
  ggplot(data=dataGraph, aes(x=ttoyear(t), y=value, group=ctax_final))+
    geom_line(aes(color=ctax_final))+
    labs(y= varUnit, title = graphTitle)+
    facet_wrap(~n,scales="free_y") +
    theme(plot.title=element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #       panel.background = element_blank(),
          axis.title.x = element_blank(), #element_text(size=12),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size=14),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          legend.position = "bottom",
          #       #legend.spacing = unit(3, "mm"),
          legend.key.width = unit(3, "cm"),
          strip.text = element_text(size=14)
    )
}

graph_var_pathway_data<- function(dataGraph, graphTitle, varUnit){
  ggplot(data=dataGraph, aes(x=ttoyear(t), y=value, group=ctax_final))+
    geom_line(aes(color=ctax_final))+
    labs(y= varUnit, title = graphTitle)+
    facet_wrap(~n,scales="free_y") +
    theme(plot.title=element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #       panel.background = element_blank(),
          axis.title.x = element_blank(), #element_text(size=12),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size=14),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          legend.position = "bottom",
          #       #legend.spacing = unit(3, "mm"),
          legend.key.width = unit(3, "cm"),
          strip.text = element_text(size=14)
    )
}

graph_var_pathway("K_DAC", "DAC installed capacity pathway", "[GtC]")
saveplot("K_DAC pathway", width = 15, height=10)


#'--------------------------------
#' For K_En(elbigcc)

K_EN_elbigcc<- getVarGraphData("K_EN", c("jreal", "t", "n"))
K_EN_elbigcc<-filterRuns_NashNotConverged(K_EN_elbigcc) #slow, takes time
K_EN_elbigcc<-K_EN_elbigcc[K_EN_elbigcc$jreal=="elbigcc",]

graph_ctax_vs_var_data(K_EN_elbigcc, "Capital in electricity with biomass and CCS","[TW]")
saveplot("K_EN(elbigcc) vs ctax", width = 15, height=10)

graph_ctax_vs_var_data(K_EN_elbigcc, "Capital in electricity with biomass and CCS-diff","[TW]")
saveplot("K_EN(elbigcc) vs ctax -diff", width = 15, height=10)

graph_var_pathway_data(K_EN_elbigcc,"Capital in electricity with biomass and CCS","[TW]")
saveplot("K_EN(elbigcc) pathway", width = 15, height=10)

#'--------------------------------
#' For Q_EMI(ccs, co2_daccs, co2_plant_ccs)

Q_EMI<-getVarGraphData("Q_EMI", c("e", "t", "n"))
Q_EMI<-filterRuns_NashNotConverged(Q_EMI)

Q_EMI_ccs<-Q_EMI[Q_EMI$e=="ccs",]
graph_ctax_vs_var_data(Q_EMI_ccs, "Emissions from ccs",  "[GtCe/year]")
saveplot("Q_EMI(ccs) vs ctax", width=15, height=10)

graph_var_pathway_data(Q_EMI_ccs, "Emissions from ccs",  "[GtCe/year]")
saveplot("Q_EMI(ccs) pathway", width=15, height=10)

Q_EMI_co2_daccs<-Q_EMI[Q_EMI$e=="co2_daccs",]
graph_ctax_vs_var_data(Q_EMI_co2_daccs, "Emissions from co2_daccs",  "[GtCe/year]")
saveplot("Q_EMI(co2_daccs) vs ctax", width=15, height=10)

graph_var_pathway_data(Q_EMI_co2_daccs, "Emissions from co2_daccs",  "[GtCe/year]")
saveplot("Q_EMI(co2_daccs) pathway", width=15, height=10)

Q_EMI_co2_plant_ccs<-Q_EMI[Q_EMI$e=="co2_plant_ccs",]
graph_ctax_vs_var_data(Q_EMI_co2_plant_ccs, "Emissions from co2_plant_ccs",  "[GtCe/year]")
saveplot("Q_EMI(co2_plant_ccs) vs ctax", width=15, height=10)

graph_var_pathway_data(Q_EMI_co2_plant_ccs, "Emissions from co2_plant_ccs",  "[GtCe/year]")
saveplot("Q_EMI(co2_plant_ccs) pathway", width=15, height=10)
