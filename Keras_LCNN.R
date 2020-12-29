require(keras)

# Toy data set
require(SoyNAM)
dataset = BLUP(rm.rep=F)
y = dataset$Phen
gen = dataset$Gen
fam = dataset$Fam

# Data
set.seed(0)
x = sample(nrow(gen),0.50*nrow(gen))
train_data = gen[x,]
train_labels = y[x]
test_data = gen[-x,]
test_labels = y[-x]
test_fam = fam[x]

dim(train_data) = c(dim(train_data),1)
dim(test_data) = c(dim(test_data),1)

# build model
nSnp = ncol(train_data)
inputs <- layer_input(shape = c(nSnp,1))
outputs <- inputs %>% 
  layer_locally_connected_1d(filters=1,kernel_size=10,strides=10) %>%
  layer_flatten() %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 1, activation = 'linear')
model <- keras_model(inputs, outputs) %>%
  compile( loss = "mse", optimizer = 'adam')
summary(model)

# Display training progress by printing a single dot for each completed epoch
print_dot_callback <- callback_lambda(
  on_epoch_end = function(epoch, logs) {
    if (epoch %% 5 == 0) cat("\n"); cat(".")} )    

# Fit the model and store training stats
epochs <- 25
history <- model %>% fit(
  train_data,
  train_labels,
  epochs = epochs,
  validation_split = 0.1,
  verbose = 0,
  callbacks = list(print_dot_callback)
)

# Plot outcome
library(ggplot2)
plot(history, metrics = "loss", smooth = FALSE) 

# Check predictions via LCNN
z = predict(model,test_data)
af_lcnn = cor(z,test_labels)
wf_lcnn = c(by(cbind(z,test_labels),test_fam,function(x)cor(x)[1,2] ))

# Check predictions via GBLUP
require(bWGR)
dim(train_data) = dim(train_data)[1:2]
dim(test_data) = dim(test_data)[1:2]
gblup = emML(train_labels,train_data)
g = gblup$mu + c(test_data%*%gblup$b)
af_gblup = cor(g,test_labels)
wf_gblup = c(by(cbind(g,test_labels),test_fam,function(x)cor(x)[1,2] ))

# Plot
barplot(c(LCNN=af_lcnn,GBLUP=af_gblup),ylab='PA across family',col=c(2,3),args.legend = list(x='bottomright'))
barplot(rbind(`LCNN    (avg = 0.63)`=wf_lcnn,`GBLUP (avg = 0.71)`=wf_gblup),ylab='PA across family',las=2,beside = T,legend=T,col=c(2,3),args.legend = list(x='bottomright'),xlab='Family')
