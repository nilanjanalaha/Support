plot.create <- function(pow, sd, tau, what)
{
  if(missing(sd)) sd <- matrix(0, nrow(pow), ncol(pow))
#Reshape the data
longpow <- melt(pow, id.vars = 1:2, variable.name = "p", value.name="estimate")
longsd <- melt(sd, id.vars = 1:2, variable.name = "p", value.name="sd")
mydata <- data.frame(longpow, sd=longsd$sd)

#errorbar
lb <- mydata$estimate - 2*mydata$sd
lb <- ifelse(lb>0, lb, 0)
ub <- mydata$estimate + 2*mydata$sd
ub <- ifelse(ub<1, ub, 1)
mydata <- cbind(mydata, lb, ub)
# factoring Method
mydata$Method <- factor(mydata$Method, levels=c("SCCA","Cleaned_SCCA","CT"))


g <- ggplot(data=mydata, aes(x=frac, y=estimate, color=p)) 
g <- g + geom_line() + facet_grid(cols = vars(Method))
g <- g+ labs(x=expression(s/sqrt(n)), y="Proportion of recoverd support")
#error-bar
g <- g + geom_errorbar(data=mydata, aes(ymin=lb, ymax=ub), width=0.02)
# change theme and text size 
tts <- 15
ts <- 15
tl <- 15
g <- g + theme_bw() + theme(legend.position = "top", axis.title= element_text(size= tts),
                            axis.text.y=element_text(size= ts),
                            #axis.text.x=element_text(size= ts),
                            legend.text=element_text(size= tl),
                            #strip.text.x=element_text(size= ts),
                            legend.title=element_text(size= tts))
g <- g + ggtitle(paste("Simulation for ", what, " when tau","=",tau))
print(g)
#dest3 <- paste(fp, "/Ramu/support/main_", where,"/",where,"plot.pdf", sep="")
#ggsave(dest3, device="pdf", width = 9, height=3.5, units = "in")
}
