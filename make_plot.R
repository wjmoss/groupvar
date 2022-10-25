.libPaths("/home/wuju/linux/R/x86_64-pc-linux-gnu-library/4.1")
library(ggplot2)
library(tidyr)
library(dplyr)
library(plyr)
source('theme_for_plots.R')


# initialize empty data frame
df = data.frame(
  p=integer(),
  n=integer(),
  sp=character(),
  gev=double(),
  pc=double(),
  ges_nv=double()
)

# load data and fill in the data frame
folder <- 'empty_stdbic'
blk <- 'npart=2'
#blk <- ''
for (p in c(5,10,20,40)){
  for (n in c(100,500,1000)){
    for (sp in c('s','d')){
      fnm = paste0('./', folder, '/p=', p, 'n=', n, 'sp=', sp, 'seed=1', blk, '.RData')
      load(fnm)
      tmp = data.frame(t(err2))
      names(tmp) = c("gev", "pc", "ges")
      tmp$p = p
      tmp$n = n
      tmp$sp = sp
      df = rbind(df, tmp)
    }
  }
}

#mean(sapply(1:100, function (i) sum(res[[i]]$ges.fit)))

# change columns into var and val
df0 = df %>% pivot_longer(., cols = c('gev','pc','ges'), names_to = "Var", values_to = "Val")

# set axis, corresponding to values in grid columns
df0$p <- mapvalues(df0$p, from = c(5,10,20,40),
                   to = c("p=5","p=10", "p=20", "p=40"))
df0$n <- mapvalues(df0$n, from = c(100,500,1000),
                   to = c("n=100","n=500", "n=1000"))

# create factor columns for fixing orders
df0$p_f = factor(df0$p, levels=c("p=5","p=10", "p=20", "p=40"))
df0$n_f = factor(df0$n, levels=c("n=100","n=500", "n=1000"))
df0$var_f = factor(df0$Var, levels=c('gev','pc','ges'))


# make plot
#theme_update(plot.title = element_text(hjust = 0.5))
ggplot(df0, aes(x=var_f, y=Val, fill=sp)) +
  geom_boxplot() +
  facet_grid(p_f~n_f, scales = 'free_y') +
  ylab("SHD to true CPDAG") + xlab("Methods") 


# sparse
df1 = df0 %>% filter(sp == 's')
ggplot(df1, aes(x=var_f, y=Val)) +
  geom_boxplot() +
  facet_grid(p_f~n_f, scales = 'free_y') +
  ylab("SHD to true CPDAG") + xlab("Methods") + 
  #ggtitle("Sparse graphs, ") + 
  theme(strip.text = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        plot.title = element_text(size = 20, face = "bold"))

# dense
df1 = df0 %>% filter(sp == 'd')
ggplot(df1, aes(x=var_f, y=Val)) +
  geom_boxplot() +
  facet_grid(p_f~n_f, scales = 'free_y') +
  ylab("SHD to true CPDAG") + xlab("Methods") + 
  #ggtitle("Plot of length \n by dose") + 
  theme(strip.text = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        plot.title = element_text(size = 20, face = "bold"))
  

# save plot?
