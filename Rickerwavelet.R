rickerwavelet<-function(dt,nt,f){
  
  t<<-seq((-(nt/2.0)*dt),(((nt-1)/2.0)*dt),dt) #creates time array 
  velocityinmeter_nanosecond<-(1.13)*1e8
  t<-seq((-(nt/2.0)*dt),(((nt-1)/2.0)*dt),dt) #creates time array 
  depthinmeters<-(velocityinmeter_nanosecond*t)/2
  b<-(pi*f)^2
  bt2=b*(t^2)
  w<- (1-2*bt2)*exp(-bt2)
}

conversiontodepthinmeters<-function(timearray,velocityinmeterperns){
  depthinmeters<-(velocityinmeterperns*1e8*timearray)/2
  return(depthinmeters)
  }
