
## prepare the input data ---------------------
   
   Climat <- read.table("grapevine_climat.txt",dec='.',header=T, sep='\t')

   Climat$Doy <- as.numeric(format(as.Date(substr(Climat$TIMESTAMP,1,10),"%d/%m/%Y"),'%j')) 
   Climat$Date <- as.Date(substr(Climat$TIMESTAMP,1,10),"%d/%m/%Y")
   Climat$Hour <- as.numeric(substr(Climat$TIMESTAMP,12,13))
   Climat$min <- as.numeric(substr(Climat$TIMESTAMP,15,16))
   Climat$Temp<- Climat$AirTC_Avg
   Climat$Rh <- (Climat$RH_Max + Climat$RH_Min)/2
  
   DATA <- aggregate(cbind(Temp,Rh,PPFD) ~ Hour+Doy,data=Climat,mean)
   DATA$dap <- DATA$Doy - 126                            ## assume that the pollination date is 05/05/2016
   DATA$Tps <- DATA$dap *24 + DATA$Hour
   DATA$RH <- DATA$Rh
   DATA$vpd <- 0.611*exp(17.27*DATA$Temp/(DATA$Temp+237.3))*(1-DATA$RH/100)

   RGmax <- aggregate(PPFD ~ Doy,data=DATA,max);   RGmin <- aggregate(PPFD ~ Doy,data=DATA,min)
   RGmax$dap <- RGmax$Doy - 126;   RGmin$dap <- RGmin$Doy - 126

   DATA$RGmin <- RGmin$PPFD[match(DATA$dap,RGmin$dap)]
   DATA$RGmax <- RGmax$PPFD[match(DATA$dap,RGmax$dap)]
   DATA <- subset(DATA,Tps >=0,select=c(dap,Tps,Temp,RH,PPFD,vpd,RGmin,RGmax))    
 
   DATA$Cpmin <- 75                                      ## mmol/L
   DATA$Cpmax <- 275                                     ## mmol/L
   DATA$Wpmax <- -1.5                                    ##  bar 
   DATA$Wpmin <- -7.5                                    ##  bar  

   DATA$Cp <- DATA$Cpmin+ (DATA$PPFD- DATA$RGmin)/(DATA$RGmax - DATA$RGmin) *(DATA$Cpmax- DATA$Cpmin)
   DATA$PTLx <- DATA$Wpmax+ (DATA$PPFD- DATA$RGmin)/(DATA$RGmax - DATA$RGmin) * (DATA$Wpmin- DATA$Wpmax)

   if(F){
   DATA$daily <- ifelse(sin((DATA$Tps%%24-5)*pi/15)>0,sin((DATA$Tps%%24-5)*pi/15),0)
   DATA$Cp <- DATA$Cpmin + DATA$daily *(DATA$Cpmax - DATA$Cpmin )
   DATA$PTLx <- DATA$Wpmax + DATA$daily * (DATA$Wpmin - DATA$Wpmax)
         }
   
   DATA_input <- subset(DATA,dap >=4 & dap < 106,select=c(Tps,Temp,RH,Cp,PTLx))

   if(F)
   {
     write.table(DATA_input,"DATA_input.txt",sep='\t',row.names=F)
   }


