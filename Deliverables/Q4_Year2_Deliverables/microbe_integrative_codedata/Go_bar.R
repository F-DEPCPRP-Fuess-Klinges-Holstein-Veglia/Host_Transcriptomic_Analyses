##this code generates a plot of associations between microbial families of interest and immune GO terms
##last updated LEF 6/12/25

##load packages
library(gridExtra)
library(grid)
library(ggplot2)
##add data##
data<-read.csv("go_bar.csv")
##sort families and roles as desired##
data$family = factor(data$family, levels=c( "Kiloniellaceae", "Lentisphaeraceae", "Oligoflexaceae", "Polixibacteraceae", 
                                            "Desulfobacteraceae","Rhodocylcaceae","Terasakiellaceae"))
data$role = factor(data$role, levels=c("positive", "neutral", "negative"))

##make the plot

ggplot(data, aes(x=reorder(name, pval),y=pval, fill=role)) + geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#bd004a", "#dee6f6", "#5a7bac"), name = "Immune Role") +
  facet_wrap(~family, ncol=4, scales="free") +
  coord_flip() +
  geom_col(color = "black", size = .2) +
  ylab("-log(padj)")+xlab(NULL) + theme_bw() + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(color="black", face="bold", size = 15),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        strip.background = element_rect(fill = "#F0F0F0", size=1, linetype="solid", 
                                        colour ="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        panel.background = element_rect(fill = "#F0F0F0"),
        legend.position = c(.93, .15),
        legend.background = element_rect(fill="#F0F0F0",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        legend.key = element_rect(fill="#F0F0F0"),
        legend.title = element_text(colour="black", size=16, 
                                    face="bold"),
        legend.text=element_text(size=14),
        legend.key.size = unit(.8,"cm"))
