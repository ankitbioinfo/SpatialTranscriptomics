
Sys.setenv(LANG="en")
require(foreign)
require(ggplot2)
require(MASS)

dat <- read.dta("https://stats.idre.ucla.edu/stat/stata/dae/nb_data.dta")

write.csv(dat,"sampleNBdata2.dat",quote=FALSE,row.names=TRUE)


dat <- within(dat, {
    prog <- factor(prog, levels = 1:3, labels = c("General", "Academic", "Vocational"))
    id <- factor(id)
})

summary(dat)

write.csv(dat,"sampleNBdata.dat",quote=FALSE,row.names=TRUE)

summary(m1 <- glm.nb(daysabs ~ math + prog, data = dat))
