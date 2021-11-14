# R script to compare UC and gaussian anamorphosis from rgeostat and pygslib

library(RGeostats)
setwd("C:/OG_python/pygslib/pygslib/Ipython_templates")

print('working directory:')    
print (getwd())

# get the data

#import drillhole data
cluster <- read.csv('../data/cluster.csv' , header=TRUE, na.strings="")

# create Rgeostat database
data = db.create(cluster)
data = db.locate(data,c("Xlocation","Ylocation"),"x")         
data = db.locate(data,c("Primary"),"z")
data = db.locate(data,c("Declustering.Weight"),"w")


#plot the data
plot(data,name.post="Primary",title="Information")

# get anamorphosis
anam = anam.fit(data, name="Primary", 
                type="gaus", nbpoly = 40,
                ndisc=5000)
print(anam)

# create a frame and save
df =  data.frame(anam@hermite)
write.csv(df, "hpoly_fromrgeostat.csv", row.names = FALSE)
