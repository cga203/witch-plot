#Auxiliary Functions




ttoyear <- function(t){year=(as.numeric(t) * 5 + 2000); return(year);}
yeartot <- function(year){t=(as.numeric(year) - 2000) / 5; return(t);}



saveplot <- function(plotname, width=7, height=5, text_size=10){
  if(!exists("legend_position")){legend_position = "bottom"}
  if(legend_position=="bottom"){legend_direction="horizontal"}else{legend_direction="vertical"}
  print(last_plot() + labs(title=plotname)); 
  ggsave(filename=paste(graphdir,as.character(gsub(" ", "_", plotname)),".png", sep=""), plot = last_plot() + labs(title=plotname) + theme(text = element_text(size=text_size), legend.position=legend_position, legend.direction = legend_direction), width=width, height=height)
}



ssptriple <- function(df) #Function converts a single "file" columns to three with SSP, RCP, SPA
{
  scenario <- df$file
  triple <- as.data.frame(matrix(0, ncol = 0, nrow = length(scenario)))
  triple$ssp=substr(scenario, 1, 4)
  triple$rcp=substr(scenario, 6, 9)
  triple$spa=substr(scenario, 11, 14)  
  triple$spa <- str_replace(triple$spa, "spa[1-5]", "spaX")
  #special cases for BAU
  if(length(triple[str_detect(triple$rcp, "bau"),1])>0){triple[str_detect(triple$rcp, "bau"),]$rcp <- "bau"}
  if(length(triple[str_detect(triple$rcp, "bau"),1])>0){triple[str_detect(triple$rcp, "bau"),]$spa <- "spa0"}
  df_new <- cbind(df, triple)
  df_new$file <- NULL
  return(df_new)
}

readkey <- function()
{
  cat ("Press [enter] to continue")
  line <- readline()
}