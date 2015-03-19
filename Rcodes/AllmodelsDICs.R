# Supplementary Material
# Applying Bayesian spatio-temporal models to fisheries bycatch in the Canadian Arctic 
# R codes for all models tested with final mesh

# Notes: For some model tested, there will be R-INLA warnings, this is normal, some models tested were inappropriate; moreover you are working with a similar but simulated dataset due to Fisheries and Oceans Canada data privacy policy.

rm(list = ls(all = TRUE))
require(xtable)
library(INLA)
### INLA website: www.r-inla.org
### the next comand line is used if the 
### inla computations is made on another machine
### if you use your computer, just comment the next line
###inla.setOption(inla.call='remote') 
inla.setOption(num.threads=8) 

#load Greenland shark simulated data
gr1 <- read.table('fackdata.txt', header=TRUE)
gr1$month <- factor(gr1$month, 7:11)
head(gr1)

# \section{The models}
# 
# We need a model that take into account 
# the spatial and the temporal domains. 
# 
# For spatial domain we consider a Gaussian Field. 
# 
# For time we consider that we have 20 months 
# and these are nested within four years. 
# We consider different temporal evolutions 
# of the spatial random effect.
# 
# \begin{enumerate} 
#   \item $m1$: the spatial random effect is constant 
#     over time
#   \item $m2$: the spatial random effect is replicated 
#     each of the four years, constant within year 
#   \item $m3a$: the spatial random effect is correlated 
#     on consecutive years, constant within year 
#   \item $m3b$: the spatial random effect is correlated 
#     on years not considering the order of year, i.e., 
#     the correlation between year 1 and 2 is considered 
#     same that correlation between year 1 and 4. 
#     Also, the spatial effect is constant within year. 
#   \item $m4$: the spatial random effect is replicated 
#     each of the 20 months 
#   \item $m5a:$ the spatial random effect is replicated 
#     each year and correlated in consecutive months, 
#     within each year 
#   \item $m5b:$ the spatial random effect is replicated 
#     each year and correlated in mounths within each year, 
#     for example, 2009-07 has same correlation with 2009-08, 
#     2009-10 and 2009-11 but correlation zero between 
#     months in another years.
# \end{enumerate}

# Testing different families 
fams <- c(paste(c('', rep('zeroinflated', 3)), 
                'poisson', c('', 0:2), sep=''), 
          paste(c('', rep('zeroinflated', 3)), 
                'nbinomial', c('', 0:2), sep=''))
names(fams) <- gsub('ial', '', gsub('eroinflated', '', fams))
fams

# covariates
names(gr1)
covars <- gr1[,c(4:7)]
head(covars)

#Formulae
form0 <-  y ~ 0 + b0 + I(duration) + Ngillnets + I(TC.tspp) + 
    f(inla.group(bathymetry), model='rw1') #inla.group reduce the number to unique values

# model 0: no spatial effect
dics.0 <- sapply(fams, function(fam) {
    cat(fam, '\n') ### just to see the progress...
    r <- inla(form0, family=fam, 
               data=data.frame(y=gr1$y, covars,b0=1), 
               control.predictor=list(compute=TRUE),
               control.inla=list(strategy='laplace'), #note that we are here using laplace, default in R-INLA is the simplified laplace approximation (run faster)
               control.compute=list(dic=TRUE, cpo=TRUE))
     c(r$dic$dic, r$cpu[4])
})
dics.0
dicsfile <- 'outputs/dic-cpu-table.txt' ##need to create folder outputs or remove the extension
cat('', file=dicsfile)
cat(dics.0[1,], dics.0[2,], '\n', file=dicsfile, append=TRUE)

##############################
### Model with spatial effects
##############################
# STEP 1: MESH CONSTRUCTION
(mesh.a <- inla.mesh.2d(cbind(gr1$X, gr1$Y), max.edge=c(40, 80), cut=5))$n
##
plot(mesh.a, asp=1,main=NULL,sub=NULL)

# STEP 2: CREATE SPDE OBJECT
spde.a <- inla.spde2.matern(mesh.a) 

# STEP 3: PREPARE DATA {STACK}
### m1: 1 spatial effect constant over time
A.1 <- inla.spde.make.A(mesh.a, loc=cbind(gr1$X, gr1$Y))
ind.1 <- inla.spde.make.index('s', mesh.a$n)
stk.1 <- inla.stack(data=list(y=gr1$y), A=list(A.1,1),
                    effects=list(ind.1, list(data.frame(b0=1,covars))))

# STEP 4: Formulae
form1 <-  y ~ 0 + b0 + I(duration) + Ngillnets + I(TC.tspp) + f(inla.group(bathymetry), model='rw1')  +  f(s, model=spde.a)

# STEP 5: FIT THE MODEL; HERE FIT TO SEVERAL FAMILIES 
dics.1 <- sapply(fams, function(fam) {
    cat(fam, '\n') ### just to see the progress...
    r <- inla(form1, family=fam, data=inla.stack.data(stk.1),
               control.predictor=list(A=inla.stack.A(stk.1), compute=TRUE),
               control.inla=list(strategy='laplace'), 
               control.compute=list(dic=TRUE, cpo=TRUE))
    c(r$dic$dic, r$cpu[4])
})
dics.1
cat(dics.1[1,], dics.1[2,], '\n', file=dicsfile, append=TRUE)

### m2: Models with random effect is replicated each year 
# prepare the data {stack}
table(repl.2 <- gr1$year-2007)
dim(A.2 <- inla.spde.make.A(mesh.a, repl=repl.2, loc=cbind(gr1$X, gr1$Y)))
ind.2 <- inla.spde.make.index(name='s', n.spde=mesh.a$n, n.repl=4)
stk.2 <- inla.stack(data=list(y=gr1$y), A=list(A.2,1),
                     effects=list(ind.2, list(data.frame(b0=1,covars))))

# formulae
form.2 <-  y ~ 0 + b0 + I(duration) + Ngillnets + I(TC.tspp) + f(inla.group(bathymetry), model='rw1')  +  f(s, model=spde.a, replicate=s.repl)

# fit the models to all families
dics.2 <- sapply(fams, function(fam) {
    cat(fam, '\n') ### just to see the progress...
    r <- inla(form.2, family=fam, data=inla.stack.data(stk.2),
               control.predictor=list(A=inla.stack.A(stk.2), compute=TRUE),
               control.inla=list(strategy='laplace'), 
               control.compute=list(dic=TRUE, cpo=TRUE))
     c(r$dic$dic, r$cpu[4])
})
dics.2
cat(dics.2[1,], dics.2[2,], '\n', file=dicsfile, append=TRUE)

#model 3: Models with the random effect is correlated between years. 
# prepare the data {stack}
table(gr.3 <- gr1$year-2007)
dim(A.3 <- inla.spde.make.A(mesh.a, group=gr.3, loc=cbind(gr1$X, gr1$Y)))
ind.3 <- inla.spde.make.index(name='s', n.spde=mesh.a$n, n.group=4)
stk.3 <- inla.stack(data=list(y=gr1$y), A=list(A.3,1),
                     effects=list(ind.3, list(data.frame(b0=1,covars))))

# formula 
# autoregressive correlatin ($m3a$) 
form.3a <- y ~ 0 + b0 + I(duration) + Ngillnets + I(TC.tspp) + f(inla.group(bathymetry), model='rw1')  + f(s, model=spde.a, group=s.group, control.group=list(model='ar1'))

# exchangeable ($m3b$)
form.3b <- y ~ 0 + b0 + I(duration) + Ngillnets + I(TC.tspp) + f(inla.group(bathymetry), model='rw1') + f(s, model=spde.a, group=s.group, control.group=list(model='exchangeable'))

# fit models with all families
# m3a
dics.3a <- sapply(fams, function(fam) {
    cat(fam, '\n') ### just to see the progress...
    r <- inla(form.3a, family=fam, data=inla.stack.data(stk.3),
               control.predictor=list(A=inla.stack.A(stk.3), compute=T),
               control.inla=list(strategy='laplace'), 
               control.compute=list(dic=TRUE, cpo=TRUE))
     c(r$dic$dic, r$cpu[4])
})
dics.3a
cat(dics.3a[1,], dics.3a[2,], '\n', file=dicsfile, append=TRUE)

# m3b
dics.3b <- sapply(fams, function(fam) {
    cat(fam, '\n') ### just to see the progress...
    r <- inla(form.3b, family=fam, data=inla.stack.data(stk.3),
               control.predictor=list(A=inla.stack.A(stk.3), compute=T),
               control.inla=list(strategy='laplace'), 
               control.compute=list(dic=TRUE, cpo=TRUE))
     c(r$dic$dic, r$cpu[4])
})
dics.3b
cat(dics.3b[1,], dics.3b[2,], '\n', file=dicsfile, append=TRUE)

#m4:Replication each of the 20 months
# prepare the data {stack}
table(repl.4 <- (gr1$year-2008)*5 + 
      as.integer(as.character(gr1$month))-6)
dim(A.4 <- inla.spde.make.A(mesh.a, repl=repl.4, loc=cbind(gr1$X, gr1$Y)))
ind.4 <- inla.spde.make.index(name='s', n.spde=mesh.a$n, n.repl=20)
stk.4 <- inla.stack(data=list(y=gr1$y), A=list(A.4,1), tag='e',
                     effects=list(ind.4, list(data.frame(b0=1,covars))))
# Formula
form.4 <-y ~ 0 + b0 + I(duration) + Ngillnets + I(TC.tspp) + f(inla.group(bathymetry), model='rw1') + f(s, model=spde.a, replicate=s.repl) 

# fit models with all families
dics.4 <- sapply(fams, function(fam) {
    cat(fam, '\n') ### just to see the progress...
    r <- inla(form.4, family=fam, data=inla.stack.data(stk.4),
               control.predictor=list(A=inla.stack.A(stk.4), compute=TRUE),
               control.inla=list(strategy='laplace'), 
               control.compute=list(dic=TRUE, cpo=TRUE))
     c(r$dic$dic, r$cpu[4])
})
dics.4
cat(dics.4[1,], dics.4[2,], '\n', file=dicsfile, append=TRUE)

#m5:Replication between years. Correlation between months within year. 
# prepare the data {stack}
table(repl.5 <- gr1$year-2007)
table(gr.5 <- as.integer(as.character(gr1$month))-6)
dim(A.5 <- inla.spde.make.A(mesh.a, loc=cbind(gr1$X, gr1$Y), 
                            repl=repl.5, group=gr.5))
ind.5 <- inla.spde.make.index(name='s', n.spde=mesh.a$n, n.repl=4, n.group=5)
stk.5 <- inla.stack(data=list(y=gr1$y), A=list(A.5,1),
                     effects=list(ind.5, list(data.frame(b0=1,covars))))

#Formula 
# autoregressive correlatin ($m5a$) 
form.5a <- y ~ 0 + b0 + I(duration) + Ngillnets + I(TC.tspp) + f(inla.group(bathymetry), model='rw1') + f(s, model=spde.a, group=s.group, control.group=list(model='ar1'), replicate=s.repl)

# exchangeable ($m5b$)
form.5b <- y ~ 0 + b0 + I(duration) + Ngillnets + I(TC.tspp) + f(inla.group(bathymetry), model='rw1') + f(s, model=spde.a, group=s.group, control.group=list(model='exchangeable'), replicate=s.repl) 

# fit models with all families
dics.5a <- sapply(fams, function(fam) {
    cat(fam, '\n') ### just to see the progress...
    r <- inla(form.5a, family=fam, data=inla.stack.data(stk.5),
               control.predictor=list(A=inla.stack.A(stk.5), compute=T),
               control.inla=list(strategy='laplace'), 
               control.compute=list(dic=TRUE, cpo=TRUE))
     c(r$dic$dic, r$cpu[4])
})
dics.5a
cat(dics.5a[1,], dics.5a[2,], '\n', file=dicsfile, append=TRUE)

dics.5b <- sapply(fams, function(fam) {
    cat(fam, '\n') ### just to see the progress...
    r <- inla(form.5b, family=fam, data=inla.stack.data(stk.5),
               control.predictor=list(A=inla.stack.A(stk.5), compute=T),
               control.inla=list(strategy='laplace'), 
               control.compute=list(dic=TRUE, cpo=TRUE))
    c(r$dic$dic, r$cpu[4])
})
dics.5b
cat(dics.5b[1,], dics.5b[2,], '\n', file=dicsfile, append=TRUE)

dic.t <- rbind(dics.0, dics.1, dics.2, dics.3a, dics.3b, dics.4, dics.5a, dics.5b)[seq(1, 16, 2),]
rownames(dic.t) <- paste('m', c(0:2, '3a', '3b', 4, '5a', '5b'), sep='')
dic.t

require(xtable)
xtable(dic.t, caption='DIC values for all model tested. Poisson and negative binomial with 0, 1, 2 are referring to the zero-inflated types.',label="tab:DIC") 

cpu.t <- rbind(dics.0, dics.1, dics.2, dics.3a, dics.3b, dics.4, dics.5a, dics.5b)[seq(2, 16, 2),]
rownames(cpu.t) <- paste('m', c(0:2, '3a', '3b', 4, '5a', '5b'), sep='')
cpu.t

xtable(cpu.t, caption='CPU values for all model tested. Poisson and negative binomial with 0, 1, 2 are referring to the zero-inflated types.',label="tab:CPU")


