# This is a script to make the famous Badass Banners first conceptualized and produced by Darshi Banan a graduate student in Andrew Leakey's Lab
library(ggplot2)

# Read in arguments
args <- commandArgs(TRUE)

# Directory program launched from
directory<-getwd()
# Base directory of trait
job <- args[1]
comp<-args[2]


dir.job<-paste(directory, job, sep="/")

st.name<-paste(job, 'concatenated_summary_table.csv', sep="_")
st.path<-paste(dir.job, st.name, sep="/")
st<-read.csv(st.path, header=F)
colnames(st)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')

# Need to add grepl statement here
if(comp == 'y'){
    st.diff<-subset(st, grepl("comp", type))
}
st<-st[st$type == 'raw',]

# Genetic map is in base directory (Lets hard code in v0.96 for the moment) 
# Will need to be flexibly wired in the future for other RIL populations

map<-read.table(paste(directory, "GBS_map_A10xB100_v0.96.csv", sep="/"), header=T, sep=",")
chrs<-t(map[1,2:ncol(map)])
pos<-t(map[2,2:ncol(map)])

genome<-cbind(chrs, pos)
colnames(genome)<-c('chrs', 'pos')
genome<-as.data.frame(genome)
genome$chrs<-as.numeric(as.character(genome$chrs))
genome$pos<-as.numeric(as.character(genome$pos))

c1.max<-max(genome[genome$chrs == 1, 'pos'])
c2.max<-max(genome[genome$chrs == 2, 'pos'])
c3.max<-max(genome[genome$chrs == 3, 'pos'])
c4.max<-max(genome[genome$chrs == 4, 'pos'])
c5.max<-max(genome[genome$chrs == 5, 'pos'])
c6.max<-max(genome[genome$chrs == 6, 'pos'])
c7.max<-max(genome[genome$chrs == 7, 'pos'])
c8.max<-max(genome[genome$chrs == 8, 'pos'])
c9.max<-max(genome[genome$chrs == 9, 'pos'])

# Here you are making an empty map. You'll plot the QTL on this scaffold plot
blank_data<-data.frame(chr=c('1','1','2','2','3','3','4','4','5','5','6','6','7','7','8','8','9','9'), x=c(0,c1.max,0,c2.max,0,c3.max,0,c4.max,0,c5.max,0,c6.max,0,c7.max,0,c8.max,0,c9.max), y=0)
st$chr<-factor(st$chr, levels=c(1,2,3,4,5,6,7,8,9))
if (comp == 'y') {
  st.diff$chr<-factor(st.diff$chr, levels=c(1,2,3,4,5,6,7,8,9))
}

# If effect size > 0 make plotting character an arrow up or if < 0 make arrow down
fx.size<-st$additive.fx
fx.size<-as.numeric(as.character(fx.size))

plot.char<-c()
for(i in 1:length(fx.size)){
  if (fx.size[i] > 0) {plot.char<-c(plot.char, '24')}
  if (fx.size[i] < 0) {plot.char<-c(plot.char, '25')}
}

if (nrow(st) == 0) {
  write("No QTL detected in the traits.", file=paste(job, 'badass_banner.pdf', sep="_"))
  quit()
} 

st$plot.char<-plot.char
st$plot.char<-as.factor(st$plot.char)

if (comp == 'y') {
# Color plotting character by experiment
treatments<-st$treatment
treatment.name<-unique(treatments)
plot.col<-c()
for(i in 1:length(treatments)){
  logical<-treatments[i] == treatment.name
  col<-which(logical, arr.ind=TRUE)
  plot.col<-c(plot.col, col)
}

st$plot.col<-plot.col
st$plot.col<-as.factor(st$plot.col)

# make the plot
setwd(dir.job)
pdf(paste(job, 'badass_banner.pdf', sep="_"))
p<-ggplot() + geom_point(data = st, aes(x = pos, y = prop.var, shape=plot.char, colour=plot.col, fill=plot.col),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(name=c("B100 allelic effect"), values=c(24,25), labels=c("Up", "Down"))
print(p + geom_errorbarh(data = st, aes(y=prop.var, xmax = R.CI_pos, xmin = L.CI_pos, height=5, x=pos)) + scale_color_manual(name=c("Treatment"), values=c(unique(plot.col)),labels=c(unique(as.character(treatments))[1], unique(as.character(treatments))[2])) + scale_fill_manual(values=c(unique(plot.col)), guide=FALSE) + ylab("% Variance")  + xlab("Genome Position"))
dev.off()

# Now lets make a badass banner for the difference
if (nrow(st.diff) > 0) {
  fx.size<-st.diff$additive.fx
  fx.size<-as.numeric(as.character(fx.size))

  plot.char<-c()
  for(i in 1:length(fx.size)){
    if (fx.size[i] > 0) {plot.char<-c(plot.char, '24')}
    if (fx.size[i] < 0) {plot.char<-c(plot.char, '25')}
  }

  st.diff$plot.char<-plot.char
  st.diff$plot.char<-as.factor(st.diff$plot.char)

  pdf(paste(job, 'diff_badass_banner.pdf', sep="_"))
  p<-ggplot() + geom_point(data = st.diff, aes(x = pos, y = prop.var, colour = type, shape=plot.char),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(name=c("B100 allelic effect"), values=c(24,25), labels=c("Up", "Down")) + scale_color_manual(values = c("comp_diff" = "red", "comp_rel_diff" = "blue", "comp_ratio" = "green"))
  print(p + geom_errorbarh(data = st.diff, aes(y=prop.var, xmax = R.CI_pos, xmin = L.CI_pos, height=5, x=pos)) + ylab("% Variance")  + xlab("Genome Position"))
  dev.off()
  }
}
if (nrow(st.diff) == 0) {
  write("No QTL detected in the difference.", file=paste(job, 'badass_banner.txt', sep="_"))
} 

if (comp == 'n') {

  # Color plotting character by experiment
  traits<-st$trait
  trait.name<-unique(traits)
  plot.col<-c()
  for(i in 1:length(traits)){
    logical<-traits[i] == trait.name
    col<-which(logical, arr.ind=TRUE)
    plot.col<-c(plot.col, col)
  }
  
st$plot.col<-plot.col
st$plot.col<-as.factor(st$plot.col)
  
# make the plot
setwd(dir.job)
pdf(paste(job, 'no_treatment_badass_banner.pdf', sep="_"))
p<-ggplot() + geom_point(data = st, aes(x = pos, y = prop.var, shape=plot.char, colour=plot.col, fill=plot.col),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(name=c("B100 allelic effect"), values=c(24,25), labels=c("Up", "Down"))
print(p + geom_errorbarh(data = st, aes(y=prop.var, xmax = R.CI_pos, xmin = L.CI_pos, height=5, x=pos)) + scale_color_manual(values=c(unique(plot.col))) + scale_fill_manual(values=c(unique(plot.col))) + ylab("% Variance")  + xlab("Genome Position"))
dev.off()

}
