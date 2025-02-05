

#function to colour background
colour_background <- function() { 
  # Paint background where ... in green
  #rect(xleft, ybottom, xright, ytop, density = NULL,
  rect(0, 1, 1, 2, col = rgb(0, 0.5, 0, alpha = 0.2), border = NA)
  # Paint background where ... in green
  rect(1, 0, 2, 1, col = rgb(0, 0.5, 0, alpha = 0.2), border = NA)
  
  # Paint background for non rectangular zones
  #Define vertices for the blue non-rectangular zones
  vertices_blue_x <- c(1, 2, 2) # Example vertices for a blue zone
  vertices_blue_y <- c(1, 1, 2)
  # Paint background where y < x in transparent blue
  polygon(vertices_blue_x, vertices_blue_y, col = rgb(0, 0, 1, alpha = 0.2), border = NA)
  
  vertices2_blue_x <- c(1, 0, 0) 
  vertices2_blue_y <- c(1, 1, 0)
  # Paint background where y < x in transparent blue
  polygon(vertices2_blue_x, vertices2_blue_y, col = rgb(0, 0, 1, alpha = 0.2), border = NA)
  
  # Define vertices for the pink non-rectangular zones
  vertices_pink_x <- c(1, 2, 1) # Example vertices for a pink zone
  vertices_pink_y <- c(1, 2, 2) # Example vertices for a pink zone
  # Paint background where y > x in transparent pink
  polygon(vertices_pink_x,vertices_pink_y, col = rgb(1, 0, 1, alpha = 0.2), border = NA)
  
  # Define vertices for the pink non-rectangular zones
  vertices2_pink_x <- c(1, 0, 1) # Example vertices for a pink zone
  vertices2_pink_y <- c(1, 0, 0) # Example vertices for a pink zone
  # Paint background where y > x in transparent pink
  polygon(vertices2_pink_x,vertices2_pink_y, col = rgb(1, 0, 1, alpha = 0.2), border = NA)
}


###########################
#agora y=meanWPT23/meanPT23 scatter plot (used in the paper)
###########################

n <- 6
par(oma = c(4,1,1,2), mfrow = c(2, 3), mar = c(5, 4, 4, 1))


# 6 figures arranged in 2 rows and 3 columns
#par(mfrow=c(2,3))


#get data from csv file
data1 <- read.table('343genes_lowLatitude_v1.csv', header = TRUE, sep = ";", dec = ".")
#or
#if keeping the simbol "/" in header
#data1 <- read.csv('343genes_lowLatitude.csv', header = TRUE, sep = ";", dec = ".", check.names = FALSE)

#select columns
data2 <- subset(data1, select = c("Gene", "meanWPT9DIVmeanPT9", "meanWPT23DIVmeanPT23"))

common_genes_G9_G23 <- c("gene-LOC117903562","gene-LOC117891490","gene-LOC117899978","gene-LOC117902328","gene-LOC117896701","gene-LOC117895424","gene-LOC117892300","gene-LOC117899984","gene-LOC117897558","gene-LOC117897194","gene-LOC117903409","gene-LOC117895433","gene-LOC117894537","gene-LOC117896003","gene-LOC117897553","gene-LOC117903196","gene-LOC117901354","gene-LOC117890238","gene-LOC117896755","gene-LOC117898322","gene-LOC117903798","gene-LOC117896651","gene-LOC117896465","gene-LOC117897087","gene-LOC117902204","gene-LOC117891988","gene-LOC117894984","gene-LOC117893817","gene-LOC117902716","gene-LOC117901234","gene-LOC117896958","gene-LOC117889634","gene-LOC117897353","gene-LOC117896697","gene-LOC117902565","gene-LOC117895229","gene-LOC117898119","gene-LOC117899573","gene-LOC117898615","gene-LOC117896543","gene-LOC117895186","gene-LOC117891341","gene-LOC117893991")

# Simple Scatterplot
plot(data2$meanWPT9DIVmeanPT9, data2$meanWPT23DIVmeanPT23, main="candidate genes at generation 9 (343 genes)",
     xlab="meanWPT9/meanPT9", ylab="meanWPT23/meanPT23", pch=19, col = ifelse(data2$Gene %in% common_genes_G9_G23, "red", "black"),xlim = c(0, 2), ylim = c(0, 2), cex.main=0.9)


# Add a vertical line at x = 5
abline(v = 1, col = "red", lty = 2)

# Add a horizontal line at y = 25
abline(h = 1, col = "red", lty = 2)

# Add a diagonal line intercept=0, slope=1
abline(a=0,b=1,col = "black", lwd=2)

colour_background()



#scatter plots 1173 genes 

#get data from csv file
data1 <- read.table('1173genes_lowLatitude.csv', header = TRUE, sep = ";", dec = ".")
#or
#if keeping the simbol "/" in header
#data1 <- read.csv('343genes_lowLatitude.csv', header = TRUE, sep = ";", dec = ".", check.names = FALSE)

#select columns
data2 <- subset(data1, select = c("Gene", "meanWPT9DIVmeanPT9", "meanWPT23DIVmeanPT23"))

common_genes_G9_G23 <- c("gene-LOC117903562","gene-LOC117891490","gene-LOC117899978","gene-LOC117902328","gene-LOC117896701","gene-LOC117895424","gene-LOC117892300","gene-LOC117899984","gene-LOC117897558","gene-LOC117897194","gene-LOC117903409","gene-LOC117895433","gene-LOC117894537","gene-LOC117896003","gene-LOC117897553","gene-LOC117903196","gene-LOC117901354","gene-LOC117890238","gene-LOC117896755","gene-LOC117898322","gene-LOC117903798","gene-LOC117896651","gene-LOC117896465","gene-LOC117897087","gene-LOC117902204","gene-LOC117891988","gene-LOC117894984","gene-LOC117893817","gene-LOC117902716","gene-LOC117901234","gene-LOC117896958","gene-LOC117889634","gene-LOC117897353","gene-LOC117896697","gene-LOC117902565","gene-LOC117895229","gene-LOC117898119","gene-LOC117899573","gene-LOC117898615","gene-LOC117896543","gene-LOC117895186","gene-LOC117891341","gene-LOC117893991")

# Simple Scatterplot
plot(data2$meanWPT9DIVmeanPT9, data2$meanWPT23DIVmeanPT23, main="candidate genes at generation G23 (1173 genes)",
     xlab="meanWPT9/meanPT9 ", ylab="meanWPT23/meanPT23", pch=19, col = ifelse(data2$Gene %in% common_genes_G9_G23, "red", "black"),xlim = c(0, 2), ylim = c(0, 2), cex.main=0.9)

# Add a vertical line at x = 5
abline(v = 1, col = "red", lty = 2)

# Add a horizontal line at y = 25
abline(h = 1, col = "red", lty = 2)

# Add a diagonal line intercept=0, slope=1
abline(a=0,b=1, col = "black", lwd=2)


colour_background()

#scatter plots 127 genes 

#get data from csv file
data1 <- read.table('127genes_lowLatitude.csv', header = TRUE, sep = ";", dec = ".")
#or
#if keeping the simbol "/" in header
#data1 <- read.csv('343genes_lowLatitude.csv', header = TRUE, sep = ";", dec = ".", check.names = FALSE)

#select columns
data2 <- subset(data1, select = c("Gene", "meanWPT9DIVmeanPT9", "meanWPT23DIVmeanPT23"))

common_genes_G9_G23 <- c("gene-LOC117903562","gene-LOC117891490","gene-LOC117899978","gene-LOC117902328","gene-LOC117896701","gene-LOC117895424","gene-LOC117892300","gene-LOC117899984","gene-LOC117897558","gene-LOC117897194","gene-LOC117903409","gene-LOC117895433","gene-LOC117894537","gene-LOC117896003","gene-LOC117897553","gene-LOC117903196","gene-LOC117901354","gene-LOC117890238","gene-LOC117896755","gene-LOC117898322","gene-LOC117903798","gene-LOC117896651","gene-LOC117896465","gene-LOC117897087","gene-LOC117902204","gene-LOC117891988","gene-LOC117894984","gene-LOC117893817","gene-LOC117902716","gene-LOC117901234","gene-LOC117896958","gene-LOC117889634","gene-LOC117897353","gene-LOC117896697","gene-LOC117902565","gene-LOC117895229","gene-LOC117898119","gene-LOC117899573","gene-LOC117898615","gene-LOC117896543","gene-LOC117895186","gene-LOC117891341","gene-LOC117893991")
sig_G9<- c("gene-LOC117893766","gene-LOC117895462","gene-LOC117896231","gene-LOC117897170","gene-LOC117892067","gene-LOC117899010","gene-LOC117892918","gene-LOC117902005","gene-LOC117903020","gene-LOC117899796","gene-LOC117896252","gene-LOC117901317","gene-LOC117893840","gene-LOC117898493","gene-LOC117900137","gene-LOC117889790","gene-LOC117903562","gene-LOC117898135","gene-LOC117897827","gene-LOC117899750","gene-LOC117899528","gene-LOC117900208","gene-LOC117891490","gene-LOC117899978","gene-LOC117897033","gene-LOC117896145","gene-LOC117902328","gene-LOC117895887","gene-LOC117891203","gene-LOC117902444","gene-LOC117896701","gene-LOC117889935","gene-LOC117903536","gene-LOC117890579","gene-LOC117897864","gene-LOC117897942","gene-LOC117895424","gene-LOC117897565","gene-LOC117899557","gene-LOC117893129","gene-LOC117891888","gene-LOC117899307","gene-LOC117896079","gene-LOC117898021","gene-LOC117892070","gene-LOC117892584","gene-LOC117892300","gene-LOC117894448","gene-LOC117903176","gene-LOC117897656","gene-LOC117903056","gene-LOC117899984","gene-LOC117903862","gene-LOC117891118","gene-LOC117900175","gene-LOC117897558","gene-LOC117894775",
           "gene-LOC117897194","gene-LOC117891037","gene-LOC117891288","gene-LOC117898607","gene-LOC117900326","gene-LOC117901666","gene-LOC117896396","gene-LOC117890903","gene-LOC117895329","gene-LOC117894839","gene-LOC117903409","gene-LOC117900834","gene-LOC117901613","gene-LOC117894151","gene-LOC117900458","gene-LOC117891944","gene-LOC117901678","gene-LOC117892505","gene-LOC117902551","gene-LOC117903346","gene-LOC117890653","gene-LOC117902338","gene-LOC117893621","gene-LOC117903621","gene-LOC117894684","gene-LOC117903827","gene-LOC117900926","gene-LOC117891771","gene-LOC117895433","gene-LOC117894537","gene-LOC117894314","gene-LOC117894853","gene-LOC117890750","gene-LOC117898469","gene-LOC117892110","gene-LOC117896505","gene-LOC117894353","gene-LOC117896003","gene-LOC117894568","gene-LOC117892143","gene-LOC117900170","gene-LOC117897319","gene-LOC117897798","gene-LOC117897767","gene-LOC117895461","gene-LOC117890818","gene-LOC117895091","gene-LOC117903750","gene-LOC117897350","gene-LOC117898790","gene-LOC117898029","gene-LOC117899043","gene-LOC117900322","gene-LOC117892069","gene-LOC117902462","gene-LOC117895370","gene-LOC117892971","gene-LOC117903161","gene-LOC117897255","gene-LOC117900080","gene-LOC117891850","gene-LOC117898767","gene-LOC117889641","gene-LOC117891538","gene-LOC117891066","gene-LOC117891276","gene-LOC117898014","gene-LOC117895787","gene-LOC117895294","gene-LOC117895964","gene-LOC117898443","gene-LOC117901194","gene-LOC117891569","gene-LOC117897553","gene-LOC117894362","gene-LOC117900309","gene-LOC117899065","gene-LOC117900083","gene-LOC117894766","gene-LOC117895917","gene-LOC117898203","gene-LOC117903196","gene-LOC117902839","gene-LOC117902358","gene-LOC117891219","gene-LOC117895844","gene-LOC117899358","gene-LOC117892608","gene-LOC117899034","gene-LOC117901760","gene-LOC117900846","gene-LOC117892707","gene-LOC117891694","gene-LOC117902257","gene-LOC117902001","gene-LOC117901354","gene-LOC117896559","gene-LOC117890238","gene-LOC117897386","gene-LOC117896755","gene-LOC117894347","gene-LOC117892190","gene-LOC117890883","gene-LOC117891308","gene-LOC117898322","gene-LOC117903798","gene-LOC117902309","gene-LOC117890069","gene-LOC117890631","gene-LOC117903011","gene-LOC117899301","gene-LOC117896395","gene-LOC117896651","gene-LOC117890384","gene-LOC117901237","gene-LOC117898656","gene-LOC117901579","gene-LOC117896465","gene-LOC117895947","gene-LOC117892064","gene-LOC117903038","gene-LOC117903108","gene-LOC117893485","gene-LOC117891945","gene-LOC117900034","gene-LOC117903019","gene-LOC117894326","gene-LOC117895538","gene-LOC117903144","gene-LOC117892282","gene-LOC117899079","gene-LOC117903164","gene-LOC117898181","gene-LOC117897087","gene-LOC117892679","gene-LOC117903110","gene-LOC117899162","gene-LOC117902204","gene-LOC117896616","gene-LOC117899113","gene-LOC117897854","gene-LOC117891988","gene-LOC117903767","gene-LOC117897867","gene-LOC117896607","gene-LOC117896874","gene-LOC117892527","gene-LOC117897844","gene-LOC117897871","gene-LOC117898936","gene-LOC117898437","gene-LOC117895604","gene-LOC117890134","gene-LOC117902661","gene-LOC117893563","gene-LOC117896707","gene-LOC117899585","gene-LOC117890503","gene-LOC117894984","gene-LOC117899109","gene-LOC117893667","gene-LOC117892357","gene-LOC117892961","gene-LOC117894381","gene-LOC117893817","gene-LOC117902716","gene-LOC117898663","gene-LOC117898222","gene-LOC117902714","gene-LOC117901234","gene-LOC117900043",
           "gene-LOC117894644","gene-LOC117890234","gene-LOC117891956","gene-LOC117895285","gene-LOC117898072","gene-LOC117894213","gene-LOC117891884","gene-LOC117898627","gene-LOC117894793","gene-LOC117901048","gene-LOC117890120","gene-LOC117903228","gene-LOC117890679","gene-LOC117896958","gene-LOC117899342","gene-LOC117890998","gene-LOC117892315","gene-LOC117889634","gene-LOC117893519","gene-LOC117897353","gene-LOC117900017","gene-LOC117899823","gene-LOC117899552","gene-LOC117902561","gene-LOC117890742","gene-LOC117895199","gene-LOC117896697","gene-LOC117892552","gene-LOC117892058","gene-LOC117896764","gene-LOC117894663","gene-LOC117898201","gene-LOC117898013","gene-LOC117892329","gene-LOC117902565","gene-LOC117899052","gene-LOC117898906","gene-LOC117903515","gene-LOC117896127","gene-LOC117895012","gene-LOC117899758","gene-LOC117889791","gene-LOC117897263","gene-LOC117902028","gene-LOC117892325","gene-LOC117901691","gene-LOC117895229","gene-LOC117895248","gene-LOC117896278","gene-LOC117900213","gene-LOC117889731","gene-LOC117889932","gene-LOC117899376","gene-LOC117896308","gene-LOC117900402","gene-LOC117896065","gene-LOC117898119","gene-LOC117895841","gene-LOC117894799","gene-LOC117903338","gene-LOC117893989","gene-LOC117903120","gene-LOC117890320","gene-LOC117901311","gene-LOC117889623","gene-LOC117893609","gene-LOC117895524","gene-LOC117896762","gene-LOC117897379","gene-LOC117899573","gene-LOC117890708","gene-LOC117896460","gene-LOC117892358","gene-LOC117889993","gene-LOC117893522","gene-LOC117898615","gene-LOC117901825","gene-LOC117901859","gene-LOC117893711","gene-LOC117897307","gene-LOC117901033","gene-LOC117896402","gene-LOC117902057","gene-LOC117896543","gene-LOC117897127","gene-LOC117901687","gene-LOC117896571","gene-LOC117902046","gene-LOC117893819","gene-LOC117897213","gene-LOC117895186","gene-LOC117895562","gene-LOC117901762","gene-LOC117895429","gene-LOC117898089","gene-LOC117891341","gene-LOC117889987","gene-LOC117900074","gene-LOC117894054","gene-LOC117900463","gene-LOC117903535","gene-LOC117903152","gene-LOC117902061","gene-LOC117891835","gene-LOC117901261","gene-LOC117900432","gene-LOC117892319","gene-LOC117889610","gene-LOC117896185","gene-LOC117890412","gene-LOC117894007","gene-LOC117895731","gene-LOC117899263","gene-LOC117893991","gene-LOC117899394")

sig_G23<- readLines("book2.csv")

# Simple Scatterplot
plot(data2$meanWPT9DIVmeanPT9, data2$meanWPT23DIVmeanPT23, main="genes with sig. Generation x Selection interaction (127 genes)",
     xlab="meanWPT9/meanPT9 ", ylab="meanWPT23/meanPT23", pch=19,col = ifelse(data2$Gene %in% common_genes_G9_G23, "red", 
                                                                              ifelse(data2$Gene %in% sig_G9, "orange",
                                                                                     ifelse(data2$Gene %in% sig_G23, "purple", "black"))),xlim = c(0, 2), ylim = c(0, 2), cex.main=0.9)



# Add a vertical line at x = 5
abline(v = 1, col = "red", lty = 2)

# Add a horizontal line at y = 25
abline(h = 1, col = "red", lty = 2)

# Add a diagonal line intercept=0, slope=1
abline(a=0,b=1, col = "black", lwd=2)


colour_background()

########
#analysing high latitude
########

# 3 figures arranged in 1 rows and 3 columns
#par(mfrow=c(1,3))

#scatter plots 1145 genes 

#get data from csv file
data1 <- read.table('1145genes_highLatitude.csv', header = TRUE, sep = ";", dec = ".")
data1[data1 == '#DIV/0!'] <- NA
#or
#if keeping the simbol "/" in header
#data1 <- read.csv('343genes_lowLatitude.csv', header = TRUE, sep = ";", dec = ".", check.names = FALSE)

#select columns
data2 <- subset(data1, select = c("Gene", "meanWNL9DIVmeanNL9", "meanWNL23DIVmeanNL23"))

common_genes_G9_G23 <- c("gene-LOC117893872","gene-LOC117891585","gene-LOC117895118","gene-LOC117893414","gene-LOC117901706","gene-LOC117902136","gene-LOC117900128","gene-LOC117890945","gene-LOC117896319","gene-LOC117901040","gene-LOC117902041","gene-LOC117890167","gene-LOC117900060","gene-LOC117902929","gene-LOC117897413","gene-LOC117899736","gene-LOC117902321","gene-LOC117890305","gene-LOC117893939","gene-LOC117902115","gene-LOC117890853","gene-LOC117900068","gene-LOC117902348","gene-LOC117901442","gene-LOC117901954","gene-LOC117893922","gene-LOC117903476","gene-LOC117898554","gene-LOC117897757","gene-LOC117901941","gene-LOC117890977","gene-LOC117898207","gene-LOC117902159","gene-LOC117893835","gene-LOC117891589","gene-LOC117902413","gene-LOC117896621","gene-LOC117896579","gene-LOC117900249","gene-LOC117894626","gene-LOC117897609","gene-LOC117892331","gene-LOC117895223","gene-LOC117901292","gene-LOC117903743","gene-LOC117895393","gene-LOC117899719","gene-LOC117890900","gene-LOC117902339","gene-LOC117898587","gene-LOC117902237","gene-LOC117897765","gene-LOC117903339","gene-LOC117889954","gene-LOC117899025","gene-LOC117902480","gene-LOC117900528","gene-LOC117896332","gene-LOC117894044","gene-LOC117890974","gene-LOC117894098","gene-LOC117895067","gene-LOC117893465","gene-LOC117899981","gene-LOC117897549","gene-LOC117895502","gene-LOC117894596","gene-LOC117891916","gene-LOC117892254","gene-LOC117896450","gene-LOC117891283","gene-LOC117897352","gene-LOC117901888","gene-LOC117892904","gene-LOC117890735","gene-LOC117893113","gene-LOC117896706","gene-LOC117890543","gene-LOC117897798","gene-LOC117903764","gene-LOC117894236","gene-LOC117898339","gene-LOC117903395","gene-LOC117891166","gene-LOC117891201","gene-LOC117890829","gene-LOC117898465","gene-LOC117902067","gene-LOC117894588","gene-LOC117895213","gene-LOC117897061","gene-LOC117897433","gene-LOC117897570","gene-LOC117890283","gene-LOC117897879","gene-LOC117899655","gene-LOC117895449","gene-LOC117891031","gene-LOC117894663","gene-LOC117894602","gene-LOC117897764","gene-LOC117900499","gene-LOC117898764","gene-LOC117892381","gene-LOC117895783","gene-LOC117896814","gene-LOC117901926","gene-LOC117898955","gene-LOC117893363","gene-LOC117891219","gene-LOC117896768","gene-LOC117896739","gene-LOC117898425","gene-LOC117889619","gene-LOC117893460","gene-LOC117892592","gene-LOC117903342","gene-LOC117894412","gene-LOC117903890","gene-LOC117898032","gene-LOC117900644","gene-LOC117902319","gene-LOC117890610","gene-LOC117897946","gene-LOC117899366","gene-LOC117901081","gene-LOC117894075","gene-LOC117902233","gene-LOC117903344","gene-LOC117898446","gene-LOC117899506","gene-LOC117893241","gene-LOC117896561","gene-LOC117897389","gene-LOC117893792","gene-LOC117893527","gene-LOC117895507","gene-LOC117901790","gene-LOC117900721","gene-LOC117892970","gene-LOC117899965","gene-LOC117903436","gene-LOC117901390","gene-LOC117898371","gene-LOC117898744","gene-LOC117901082","gene-LOC117891717","gene-LOC117896126","gene-LOC117894065","gene-LOC117896075","gene-LOC117890107","gene-LOC117892819","gene-LOC117892991","gene-LOC117891878","gene-LOC117903769","gene-LOC117897084","gene-LOC117897559","gene-LOC117890307","gene-LOC117893165","gene-LOC117895341","gene-LOC117901549","gene-LOC117898086","gene-LOC117890179","gene-LOC117897524","gene-LOC117898223","gene-LOC117892357","gene-LOC117895455","gene-LOC117890411","gene-LOC117893475","gene-LOC117898915","gene-LOC117899406","gene-LOC117903121","gene-LOC117895569","gene-LOC117902597","gene-LOC117893837","gene-LOC117895201","gene-LOC117903217","gene-LOC117892244","gene-LOC117891060","gene-LOC117900805","gene-LOC117902756","gene-LOC117896562","gene-LOC117890987","gene-LOC117890113","gene-LOC117898390","gene-LOC117898952","gene-LOC117891511","gene-LOC117892053","gene-LOC117892042","gene-LOC117896837","gene-LOC117891246","gene-LOC117894924","gene-LOC117893814","gene-LOC117897872","gene-LOC117898393","gene-LOC117893102","gene-LOC117901152","gene-LOC117901026","gene-LOC117901689") 

# Simple Scatterplot
plot(data2$meanWNL9DIVmeanNL9, data2$meanWNL23DIVmeanNL23, main="candidate genes at generation 9 (1145 genes)",
     xlab="meanWNL9/meanNL9 ", ylab="meanWNL23/meanNL23", pch=19, col = ifelse(data2$Gene %in% common_genes_G9_G23, "red", "black"),xlim = c(0, 2), ylim = c(0, 2), cex.main=0.9)

# Add a vertical line at x = 5
abline(v = 1, col = "red", lty = 2)

# Add a horizontal line at y = 25
abline(h = 1, col = "red", lty = 2)

# Add a diagonal line intercept=0, slope=1
abline(a=0,b=1, col = "black", lwd=2)


colour_background()

#scatter plots 1939 genes 

#get data from csv file
data1 <- read.table('1939genes_highLatitude.csv', header = TRUE, sep = ";", dec = ".")
data1[data1 == '#DIV/0!'] <- NA
#or
#if keeping the simbol "/" in header
#data1 <- read.csv('343genes_lowLatitude.csv', header = TRUE, sep = ";", dec = ".", check.names = FALSE)

#select columns
data2 <- subset(data1, select = c("Gene", "meanWNL9DIVmeanNL9", "meanWNL23DIVmeanNL23"))

common_genes_G9_G23 <- c("gene-LOC117893872","gene-LOC117891585","gene-LOC117895118","gene-LOC117893414","gene-LOC117901706","gene-LOC117902136","gene-LOC117900128","gene-LOC117890945","gene-LOC117896319","gene-LOC117901040","gene-LOC117902041","gene-LOC117890167","gene-LOC117900060","gene-LOC117902929","gene-LOC117897413","gene-LOC117899736","gene-LOC117902321","gene-LOC117890305","gene-LOC117893939","gene-LOC117902115","gene-LOC117890853","gene-LOC117900068","gene-LOC117902348","gene-LOC117901442","gene-LOC117901954","gene-LOC117893922","gene-LOC117903476","gene-LOC117898554","gene-LOC117897757","gene-LOC117901941","gene-LOC117890977","gene-LOC117898207","gene-LOC117902159","gene-LOC117893835","gene-LOC117891589","gene-LOC117902413","gene-LOC117896621","gene-LOC117896579","gene-LOC117900249","gene-LOC117894626","gene-LOC117897609","gene-LOC117892331","gene-LOC117895223","gene-LOC117901292","gene-LOC117903743","gene-LOC117895393","gene-LOC117899719","gene-LOC117890900","gene-LOC117902339","gene-LOC117898587","gene-LOC117902237","gene-LOC117897765","gene-LOC117903339","gene-LOC117889954","gene-LOC117899025","gene-LOC117902480","gene-LOC117900528","gene-LOC117896332","gene-LOC117894044","gene-LOC117890974","gene-LOC117894098","gene-LOC117895067","gene-LOC117893465","gene-LOC117899981","gene-LOC117897549","gene-LOC117895502","gene-LOC117894596","gene-LOC117891916","gene-LOC117892254","gene-LOC117896450","gene-LOC117891283","gene-LOC117897352","gene-LOC117901888","gene-LOC117892904","gene-LOC117890735","gene-LOC117893113","gene-LOC117896706","gene-LOC117890543","gene-LOC117897798","gene-LOC117903764","gene-LOC117894236","gene-LOC117898339","gene-LOC117903395","gene-LOC117891166","gene-LOC117891201","gene-LOC117890829","gene-LOC117898465","gene-LOC117902067","gene-LOC117894588","gene-LOC117895213","gene-LOC117897061","gene-LOC117897433","gene-LOC117897570","gene-LOC117890283","gene-LOC117897879","gene-LOC117899655","gene-LOC117895449","gene-LOC117891031","gene-LOC117894663","gene-LOC117894602","gene-LOC117897764","gene-LOC117900499","gene-LOC117898764","gene-LOC117892381","gene-LOC117895783","gene-LOC117896814","gene-LOC117901926","gene-LOC117898955","gene-LOC117893363","gene-LOC117891219","gene-LOC117896768","gene-LOC117896739","gene-LOC117898425","gene-LOC117889619","gene-LOC117893460","gene-LOC117892592","gene-LOC117903342","gene-LOC117894412","gene-LOC117903890","gene-LOC117898032","gene-LOC117900644","gene-LOC117902319","gene-LOC117890610","gene-LOC117897946","gene-LOC117899366","gene-LOC117901081","gene-LOC117894075","gene-LOC117902233","gene-LOC117903344","gene-LOC117898446","gene-LOC117899506","gene-LOC117893241","gene-LOC117896561","gene-LOC117897389","gene-LOC117893792","gene-LOC117893527","gene-LOC117895507","gene-LOC117901790","gene-LOC117900721","gene-LOC117892970","gene-LOC117899965","gene-LOC117903436","gene-LOC117901390","gene-LOC117898371","gene-LOC117898744","gene-LOC117901082","gene-LOC117891717","gene-LOC117896126","gene-LOC117894065","gene-LOC117896075","gene-LOC117890107","gene-LOC117892819","gene-LOC117892991","gene-LOC117891878","gene-LOC117903769","gene-LOC117897084","gene-LOC117897559","gene-LOC117890307","gene-LOC117893165","gene-LOC117895341","gene-LOC117901549","gene-LOC117898086","gene-LOC117890179","gene-LOC117897524","gene-LOC117898223","gene-LOC117892357","gene-LOC117895455","gene-LOC117890411","gene-LOC117893475","gene-LOC117898915","gene-LOC117899406","gene-LOC117903121","gene-LOC117895569","gene-LOC117902597","gene-LOC117893837","gene-LOC117895201","gene-LOC117903217","gene-LOC117892244","gene-LOC117891060","gene-LOC117900805","gene-LOC117902756","gene-LOC117896562","gene-LOC117890987","gene-LOC117890113","gene-LOC117898390","gene-LOC117898952","gene-LOC117891511","gene-LOC117892053","gene-LOC117892042","gene-LOC117896837","gene-LOC117891246","gene-LOC117894924","gene-LOC117893814","gene-LOC117897872","gene-LOC117898393","gene-LOC117893102","gene-LOC117901152","gene-LOC117901026","gene-LOC117901689") 

# Simple Scatterplot
plot(data2$meanWNL9DIVmeanNL9, data2$meanWNL23DIVmeanNL23, main="candidate genes at generation 23 (1939 genes)",
     xlab="meanWNL9/meanNL9 ", ylab="meanWNL23/meanNL23", pch=19, col = ifelse(data2$Gene %in% common_genes_G9_G23, "red", "black"),xlim = c(0, 2), ylim = c(0, 2), cex.main=0.9)

# Add a vertical line at x = 5
abline(v = 1, col = "red", lty = 2)

# Add a horizontal line at y = 25
abline(h = 1, col = "red", lty = 2)

# Add a diagonal line intercept=0, slope=1
abline(a=0,b=1, col = "black", lwd=2)

colour_background()

#scatter plots 409 genes 

#get data from csv file
data1 <- read.table('409genes_highLatitude.csv', header = TRUE, sep = ";", dec = ".")
data1[data1 == '#DIV/0!'] <- NA
#or
#if keeping the simbol "/" in header
#data1 <- read.csv('343genes_lowLatitude.csv', header = TRUE, sep = ";", dec = ".", check.names = FALSE)

#select columns
data2 <- subset(data1, select = c("Gene", "meanWNL9DIVmeanNL9", "meanWNL23DIVmeanNL23"))

common_genes_G9_G23 <- c("gene-LOC117893872","gene-LOC117891585","gene-LOC117895118","gene-LOC117893414","gene-LOC117901706","gene-LOC117902136","gene-LOC117900128","gene-LOC117890945","gene-LOC117896319","gene-LOC117901040","gene-LOC117902041","gene-LOC117890167","gene-LOC117900060","gene-LOC117902929","gene-LOC117897413","gene-LOC117899736","gene-LOC117902321","gene-LOC117890305","gene-LOC117893939","gene-LOC117902115","gene-LOC117890853","gene-LOC117900068","gene-LOC117902348","gene-LOC117901442","gene-LOC117901954","gene-LOC117893922","gene-LOC117903476","gene-LOC117898554","gene-LOC117897757","gene-LOC117901941","gene-LOC117890977","gene-LOC117898207","gene-LOC117902159","gene-LOC117893835","gene-LOC117891589","gene-LOC117902413","gene-LOC117896621","gene-LOC117896579","gene-LOC117900249","gene-LOC117894626","gene-LOC117897609","gene-LOC117892331","gene-LOC117895223","gene-LOC117901292","gene-LOC117903743","gene-LOC117895393","gene-LOC117899719","gene-LOC117890900","gene-LOC117902339","gene-LOC117898587","gene-LOC117902237","gene-LOC117897765","gene-LOC117903339","gene-LOC117889954","gene-LOC117899025","gene-LOC117902480","gene-LOC117900528","gene-LOC117896332","gene-LOC117894044","gene-LOC117890974","gene-LOC117894098","gene-LOC117895067","gene-LOC117893465","gene-LOC117899981","gene-LOC117897549","gene-LOC117895502","gene-LOC117894596","gene-LOC117891916","gene-LOC117892254","gene-LOC117896450","gene-LOC117891283","gene-LOC117897352","gene-LOC117901888","gene-LOC117892904","gene-LOC117890735","gene-LOC117893113","gene-LOC117896706","gene-LOC117890543","gene-LOC117897798","gene-LOC117903764","gene-LOC117894236","gene-LOC117898339","gene-LOC117903395","gene-LOC117891166","gene-LOC117891201","gene-LOC117890829","gene-LOC117898465","gene-LOC117902067","gene-LOC117894588","gene-LOC117895213","gene-LOC117897061","gene-LOC117897433","gene-LOC117897570","gene-LOC117890283","gene-LOC117897879","gene-LOC117899655","gene-LOC117895449","gene-LOC117891031","gene-LOC117894663","gene-LOC117894602","gene-LOC117897764","gene-LOC117900499","gene-LOC117898764","gene-LOC117892381","gene-LOC117895783","gene-LOC117896814","gene-LOC117901926","gene-LOC117898955","gene-LOC117893363","gene-LOC117891219","gene-LOC117896768","gene-LOC117896739","gene-LOC117898425","gene-LOC117889619","gene-LOC117893460","gene-LOC117892592","gene-LOC117903342","gene-LOC117894412","gene-LOC117903890","gene-LOC117898032","gene-LOC117900644","gene-LOC117902319","gene-LOC117890610","gene-LOC117897946","gene-LOC117899366","gene-LOC117901081","gene-LOC117894075","gene-LOC117902233","gene-LOC117903344","gene-LOC117898446","gene-LOC117899506","gene-LOC117893241","gene-LOC117896561","gene-LOC117897389","gene-LOC117893792","gene-LOC117893527","gene-LOC117895507","gene-LOC117901790","gene-LOC117900721","gene-LOC117892970","gene-LOC117899965","gene-LOC117903436","gene-LOC117901390","gene-LOC117898371","gene-LOC117898744","gene-LOC117901082","gene-LOC117891717","gene-LOC117896126","gene-LOC117894065","gene-LOC117896075","gene-LOC117890107","gene-LOC117892819","gene-LOC117892991","gene-LOC117891878","gene-LOC117903769","gene-LOC117897084","gene-LOC117897559","gene-LOC117890307","gene-LOC117893165","gene-LOC117895341","gene-LOC117901549","gene-LOC117898086","gene-LOC117890179","gene-LOC117897524","gene-LOC117898223","gene-LOC117892357","gene-LOC117895455","gene-LOC117890411","gene-LOC117893475","gene-LOC117898915","gene-LOC117899406","gene-LOC117903121","gene-LOC117895569","gene-LOC117902597","gene-LOC117893837","gene-LOC117895201","gene-LOC117903217","gene-LOC117892244","gene-LOC117891060","gene-LOC117900805","gene-LOC117902756","gene-LOC117896562","gene-LOC117890987","gene-LOC117890113","gene-LOC117898390","gene-LOC117898952","gene-LOC117891511","gene-LOC117892053","gene-LOC117892042","gene-LOC117896837","gene-LOC117891246","gene-LOC117894924","gene-LOC117893814","gene-LOC117897872","gene-LOC117898393","gene-LOC117893102","gene-LOC117901152","gene-LOC117901026","gene-LOC117901689") 

sig_G9<- read.table("1145genes_highLatitude.csv", sep=";",header = TRUE)
sig_G9<- sig_G9$Gene
sig_G23<- read.table("1939genes_highLatitude.csv", sep=";",header = TRUE)
sig_G23<- sig_G23$Gene

# Simple Scatterplot
plot(data2$meanWNL9DIVmeanNL9, data2$meanWNL23DIVmeanNL23, main="genes with sig. Generation x Selection interaction (409 genes)",
     xlab="meanWNL9/meanNL9 ", ylab="meanWNL23/meanNL23", pch=19, col = ifelse(data2$Gene %in% common_genes_G9_G23, "red", 
                                                                               ifelse(data2$Gene %in% sig_G9, "orange",
                                                                                      ifelse(data2$Gene %in% sig_G23, "purple", "black"))),xlim = c(0, 2), ylim = c(0, 2), cex.main=0.9)

# Add a vertical line at x = 5
abline(v = 1, col = "red", lty = 2)

# Add a horizontal line at y = 25
abline(h = 1, col = "red", lty = 2)

# Add a diagonal line intercept=0, slope=1
abline(a=0,b=1, col = "black", lwd=2)

colour_background()


par(fig = c(0, 1, 0, 1), oma = c(1, 1, 1, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("candidates in both generations", "candidates at G9", "candidates at G23"), col = c("red","orange","purple"), pch=19, xpd = TRUE, horiz = TRUE, cex = 1, bty = 'n')
# xpd = TRUE makes the legend plot to the figure



############
#Empty plots
############

# 6 figures arranged in 2 rows and 3 columns
par(mfrow=c(2,3))

ylim<- c(0,2)
xlim<- c(0,2)

# Adjust figure margins
par(mar = c(2, 2, 2, 2) + 0.1)  # c(bottom, left, top, right)

# Create an empty plot without points
plot(NA, xlim=xlim, ylim=ylim, xlab="X-axis", ylab="Y-axis")
  
# Add a vertical line at x = 5
abline(v=1,col="red",lty=2)

# Add a horizontal line at y = 25
abline(h = 1, col = "red", lty = 2)

# Draw the diagonal line with thicness 2
abline(a = 0, b = 1, col = "black", lwd=2)

colour_background()








