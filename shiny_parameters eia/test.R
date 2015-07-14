rm(list = ls())
load('test.params.ws.rdata')
test.params <- runSim(test.params)
save(test.params, file = 'testOutput.rdata')
